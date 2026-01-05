#!/usr/bin/env python3

##############################################################################################
# Copyright 2022 The Johns Hopkins University Applied Physics Laboratory LLC
# All rights reserved.
# Permission is hereby granted, free of charge, to any person obtaining a copy of this
# software and associated documentation files (the "Software"), to deal in the Software
# without restriction, including without limitation the rights to use, copy, modify,
# merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
# PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
# OR OTHER DEALINGS IN THE SOFTWARE.
#

"""Provide a command line tool to fetch a list of refseq genome ids to a single file, useful for kraken2 database building or alignment purposes"""
import time
from Bio import SeqIO, Entrez
from xmlrpc.client import Boolean
from functools import partial
from mimetypes import guess_type
from typing import List
from tokenize import String
from tabnanny import filename_only
from contextlib import closing
import gzip

import argparse
import csv
import logging
import sys
from collections import Counter, defaultdict
from pathlib import Path
import re
import os
import shutil
import urllib.request as request
import ssl
from determine_priority_assembly import determine_priority_assembly, format_description

ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE


# import requests
logger = logging.getLogger()


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and transform a tabular samplesheet.",
        epilog="Example: python check_samplesheet.py samplesheet.csv samplesheet.valid.csv",
    )
    parser.add_argument(
        "-i",
        "--input",
        metavar="INPUT",
        nargs="+",
        help="List of taxIDs or names",
    )
    parser.add_argument(
        "-b",
        "--db",
        metavar="DB",
        default="nuccore",
        help="Database of choice to pull IDs from",
    )
    parser.add_argument(
        "-s",
        "--ftp_path",
        type=int,
        default=19,
        help="ftp path column, to be used instead of esummary when using the assembly file as reference",
    )
    parser.add_argument(
        "-r",
        "--refresh",
        action='store_true',
        help="Dont index already pulled fasta files",
    )
    parser.add_argument(
        "-f",
        "--type",
        default='file',
        help="Input type, can be a list of taxids or names or from a file",
    )
    parser.add_argument(
        "-m",
        "--missingfile",
        default=None,
        required = False,
        help="Missing taxids or names to a file",
    )
    parser.add_argument(
        "-g",
        "--gcf_map",
        default='file',
        required = False,
        help="Output file that is a mapping of gcf & chr accession to taxid or name. 3 columns with tab separator is made",
    )
    parser.add_argument(
        "-c",
        "--colnumber_file_hits",
        type=int,
        default=5, # 5 is for taxid, 6 for name
        help="Column number to get taxids or names from if using a file. Starts at 1st index",
    )
    parser.add_argument(
        "-a",
        "--assembly_names",
        type=str,
        default="5,6",
        help="Assembly refseq accession and taxid or name to be matched to the imported taxids",
    )
    parser.add_argument(
        "-y",
        "--name_col_assembly",
        default=7,
        type=int,
        help="Name column in assembly file you'd like to make in the header",
    )
    parser.add_argument(
        "-k",
        "--kraken2output",
        action='store_true',
        default = False,
        help="reformat header for each fasta to a kraken:taxid/name|id parsing. Requires the setup of the file to be kraken:taxid|taxid|refseqId (whatever text can come after this, separated by space(s)",
    )
    parser.add_argument(
        "-p",
        "--assembly_map_idx",
        default=0,
        help="Assembly refseq accession and taxid or name",
    )
    parser.add_argument(
        "-t",
        "--assembly_metadata",
        type=Path,
        nargs="+",
        help="Assembly refseq to pull taxids or names from instead of pulling straight from a set of IDs. Map the taxid to accession",
    )
    parser.add_argument(
        '--taxids',
        action='store_true',
        default=False,
        help='Add taxids as a 5th column to the output gcfmapping file',
    )
    parser.add_argument(
        "-e",
        "--email",
        metavar="EMAIL",
        type=str,
        help="Email for entrez querying, optional",
    )
    parser.add_argument(
        "-d",
        "--dry_run",
        action='store_true',
        default=False,
        help="Only Dry run the total filesizes of download",
    )
    parser.add_argument(
        "--skip_incomplete",
        action='store_true',
        default=False,
        help="Skip any genome that doesnt have a complete genome",
    )
    parser.add_argument(
        "-o",
        "--file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Name of the output FASTA file to put all fasta references into",
    )

    return parser.parse_args(argv)


def import_genome_file(filename, kraken2output):
    refs = dict()
    with open(filename, "r") as f:
        line = f.readline()
        for line in f:
            line = line.strip()
            try:
                if (kraken2output):
                    firstline = line.split(" ")[0]
                    header = firstline.split(" | ")[1]
                    refs[header] = line
                else:
                    header = line.split(" ")[0]
                    refs[header] = line
            except Exception as er:
                pass
    print("Done")
    return refs

#

def import_assembly_file(
    input,
    filename,
    matchcol,
    idx,
    nameidx,
    index_ftp,
    missingfile=None,
    skip_incomplete=False,
):
    assemblies = {}
    first = {}

    # Normalize input -> list
    if not isinstance(input, list):
        input = input.split(" ")

    # Normalize filename -> list (ordered)
    if isinstance(filename, str):
        filenames = [filename]
    else:
        filenames = list(filename)

    matchcol = str(matchcol)
    assembly_cols = [int(x) for x in matchcol.split(",")]

    def get_url(utl, id_):
        bb = os.path.basename(utl)
        return utl + "/" + bb + "_genomic.fna.gz"

    # priorities[matchval][prio] = (file_index, obj)
    priorities = {}
    seencols = set()

    def maybe_set(priorities_dict, matchval, prio_key, file_i, obj):
        """
        Keep the earliest-file entry for a given (matchval, prio_key).
        """
        cur = priorities_dict[matchval].get(prio_key)
        if cur is None or file_i < cur[0]:
            priorities_dict[matchval][prio_key] = (file_i, obj)

    for file_i, fn in enumerate(filenames):
        with open(fn, "r") as f:
            # Preserve your header handling style (you had readline + next)
            _ = f.readline()
            try:
                next(f)
            except StopIteration:
                continue

            for line in f:
                line = line.strip()
                if not line:
                    continue

                linesplit = line.split("\t")
                # sanity: need indexes
                if len(linesplit) <= max(idx, nameidx, index_ftp, 11, 4):
                    continue

                gcfidx = linesplit[idx]
                matchidces = [linesplit[x] for x in assembly_cols if x < len(linesplit)]
                namecol = linesplit[nameidx]
                taxidcol = linesplit[5]
                species_taxidcol = linesplit[6]
                urlcol = linesplit[index_ftp]

                formatted_header = str(namecol).replace(" ", "_")

                # find first match in input order
                match_index = -1
                for match in matchidces:
                    if match in input:
                        match_index = input.index(match)
                        break
                if match_index < 0:
                    continue

                # optional filter
                if skip_incomplete and linesplit[11] != "Complete Genome":
                    continue

                matchval = input[match_index]
                seencols.add(matchval)

                if matchval not in priorities:
                    priorities[matchval] = {}

                obj = dict(
                    id="{}|{}".format(gcfidx, formatted_header),
                    accession=gcfidx,
                    characteristic=None,
                    chrs=[],
                    organism=namecol,
                    reference=get_url(urlcol, gcfidx),
                    name=formatted_header,
                    taxidcol=taxidcol,
                    species_taxidcol=species_taxidcol,
                )

                # assign priority bucket
                # NOTE: your numbering comments say 1-4 but your code uses '0','1','2','3'
                # We'll keep your existing meaning:
                # 0 = representative genome
                # 1 = reference genome
                # 2 = complete genome
                # 3 = other

                if linesplit[4] == "representative genome":
                    obj["characteristic"] = "representative"
                    maybe_set(priorities, matchval, "0", file_i, obj)

                elif linesplit[4] == "reference genome":
                    obj["characteristic"] = "reference"
                    maybe_set(priorities, matchval, "1", file_i, obj)

                elif linesplit[11] == "Complete Genome":
                    obj["characteristic"] = "complete genome"
                    maybe_set(priorities, matchval, "2", file_i, obj)

                else:
                    # preserve your "first time seen" behavior
                    if formatted_header not in first:
                        first[formatted_header] = True
                        obj["characteristic"] = "other"
                        maybe_set(priorities, matchval, "3", file_i, obj)

    # pick best per matchval: lowest prio, then earliest file index
    for matchval, prmap in priorities.items():
        for prio in ("0", "1", "2", "3"):
            if prmap.get(prio):
                assemblies[matchval] = prmap[prio][1]  # obj
                break
        else:
            print("No reference genome found for", matchval)

    # format as you did: keyed by organism
    assembliesformat = {}
    for _, item in assemblies.items():
        assembliesformat[item["organism"]] = dict(
            accession=item["accession"],
            id=item["id"],
            name=item["name"],
            reference=item["reference"],
            chrs=item["chrs"],
            taxid=item["taxidcol"],   # taxid
            species_taxid=item["species_taxidcol"],   # species taxid
        )
    missing = list(set(input) - set(seencols))
    if missing:
        print("Missing", missing)
        if missingfile:
            with open(os.path.join(missingfile), "w") as w:
                for m in missing:
                    w.write(f"{m}\n")

    return assembliesformat


def get_assembly_summary(id):
    """Get esummary for an entrez id"""
    esummary_handle = Entrez.esummary(db="assembly", id=id, report="full")
    esummary_record = Entrez.read(esummary_handle)
    return esummary_record


def download_fasta(ftp_site):
    encoding = guess_type(ftp_site)[1]   # uses file extension
    _open = partial(
        gzip.open, mode='rt') if encoding == 'gzip' else open
    try:
        with closing(request.urlopen(ftp_site, context=ctx)) as r:
            # file.gz with basename of the inputfilename is opened
            with open('file.gz', 'wb') as f:
                shutil.copyfileobj(r, f)
            f.close()
        r.close()
        return _open
    except Exception as ex:
        print("Could not download", ftp_site, ex)
        raise ex

def get_assemblies(assemblies, outfile, GCFs_to_skip, refresh):
    """Download genbank assemblies for a given search term.
    Args:
        term: search term, usually organism name
        download: whether to download the results
        path: folder to save to
    """
    # get the value.accession from the refs dict asa list
    ids = list(assemblies.keys())
    accessions = [assemblies[id]['accession'] for id in ids]
    # filter accessions on those present in GCFS_to_skip
    accessions = [x for x in accessions if x not in GCFs_to_skip]
    caught_ncs = [ ]
    new_mapping = dict()
    typee = "a"
    if refresh:
        typee = "w"


    with open(outfile, typee) as w:
        for id in ids:
            try:
                accession = assemblies[id]['accession']
                if accession in GCFs_to_skip:
                    assemblies[id]['chrs'] = GCFs_to_skip[accession]
                    print("key already seen:", id, "; skipping")
                else:
                    ftp_site = assemblies[id]['reference']
                    obj = assemblies[id]['id']
                    name = assemblies[id]['name']
                    _open = download_fasta(ftp_site)
                    print(f"Downloading {ftp_site} to for {accession}: {id}")
                    with _open('file.gz') as uncompressed:
                        for record in SeqIO.parse(uncompressed, "fasta"):
                            # pattern = f"{record.id}\s*"
                            # record.description = re.sub(pattern, "", record.description)
                            # record.description = record.description.replace(" ", "_")
                            # newobj = f"{record.id} {accession} {record.description}"
                            if record.id not in caught_ncs:
                                caught_ncs.append(record.id)
                            if (not refresh and record.id not in new_mapping) or refresh:
                                print("writing", record.id, "to file")

                                SeqIO.write(record, w, "fasta")
                            elif not refresh and record.id in new_mapping:
                                print("already seen", record.id, "skipping")
                            new_mapping[record.id] = dict(accession=accession, name=obj)
                            if record.id not in assemblies[id]['chrs']:
                                fid, description = format_description(record.id, record.description)
                                assemblies[id]['chrs'].append(dict(acc=record.id, name=description) )
                    uncompressed.close()
                    try:
                        os.remove('file.gz')
                    except Exception as ex:
                        print(f"Could not remove file.gz {ex}")
                        pass
            except Exception as ex:
                print(ex)
                pass
    w.close()
    return new_mapping

def download(refs, db, outfile, seen):
    # if refs is not empty
    if len(refs.items()) == 0:
        return
    with open(outfile, "a") as w:
        i = 0
        maxt = 30
        next_ = []
        retry_max = 3
        for key, value in refs.items():
            try:

                if seen and key in seen:
                    print("key already seen:", key, "; skipping")
                else:
                    next_.append(key)
                if i % maxt == 0 and len(next_) > 0:
                    print(str(i), " th iteration of ids to submit..", next_, db)
                    for retry_count in range(retry_max):
                        try:
                            handle = Entrez.efetch(
                                db=db, rettype="fasta", retmode="fasta", id=",".join(next_), idtype="acc")
                        except Exception as e:
                            if retry_count < retry_max - 1:
                                print(f"Failed to fetch records (attempt {retry_count + 1})), retrying...")
                                time.sleep(10) # Wait for a few seconds before retrying
                            else:
                                print(f"Failed to fetch records after {retry_max} attempts. Error: {str(e)}")
                    seq_records = SeqIO.parse(handle, 'fasta')
                    for seq_record in seq_records:
                        if (seq_record):
                            seq_record.id = str(value.replace(">", ""))
                            SeqIO.write(seq_record, w, "fasta")
                    handle.close()
                    next_ = []

            except Exception as err:
                print("No seq record found", next_, err)
                next_ = []
                pass
            i = i+1


def get_hits_from_file(filenames, colnumber):
    taxids = []
    for filename in filenames:
        with open(filename, "r") as f:
            lines = f.readlines()
            for line in lines:
                line = line.rstrip()
                name = line.split("\t")[colnumber-1]
                if name not in taxids:
                    taxids.append(name)
    return taxids
def return_format_size(size_in_bytes):
    # convert bytes to gb, mb or tb
    format = "bytes"
    if size_in_bytes > 1000000000000:
        format = "tb"
        size_in_bytes = size_in_bytes/1000000000000
    elif size_in_bytes > 1000000000:
        format = "gb"
        size_in_bytes = size_in_bytes/1000000000
    elif size_in_bytes > 1000000:
        format = "mb"
        size_in_bytes = size_in_bytes/1000000
    elif size_in_bytes > 1000:
        format = "kb"
        size_in_bytes = size_in_bytes/1000
    return size_in_bytes, format
def check_size(assemblies, GCFs_to_skip):
    final_size = 0
    i = 0
    for key, value in assemblies.items():
        # retrieve the expected size of the file at value.reference
        ftp_site = value['reference']
        try:
            with closing(request.urlopen(ftp_site, context=ctx)) as r:
                final_size = final_size + int(r.info()['Content-Length'])
                r.close()
        except Exception as ex:
            print("Could not get info for file", ftp_site, ex)
            raise ex
        finally:
            i+=1
            if i % 100== 0:
                print("Checked", i, "files")
                format_size, format = return_format_size(final_size)
                print("Current size is", format_size, format)
    format_size, format = return_format_size(final_size)

    print("Total size of files to download is", format_size, format)

    return final_size
def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)


    gcf_mapping = dict()

    if args.gcf_map:
        if os.path.exists(args.gcf_map):
            #import the file, and save them to a dict of column 2 as key, column 1 as value
            with open(args.gcf_map, "r") as f:
                lines = f.readlines()
                for line in lines:
                    line = line.strip()
                    linesplit = line.split("\t")
                    if len(linesplit) > 1:
                        key = linesplit[0]
                        value = linesplit[1]
                        gcf_mapping[key] = value
            f.close()
    # logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if args.type == 'file':
        #colnumber file hits is the column from the input top hits file you want to match to the args.assembly_names
        seen_in_tops = get_hits_from_file(args.input, args.colnumber_file_hits)
    else:
        seen_in_tops = args.input

    assemblies = import_assembly_file(
        seen_in_tops, args.assembly_metadata, args.assembly_names, args.assembly_map_idx, args.name_col_assembly, args.ftp_path, args.missingfile, args.skip_incomplete
    )



    i = 0
    seen_in_fasta = dict()
    if os.path.exists(args.file_out):
        for seq_record in SeqIO.parse(args.file_out, "fasta"):
            line = str(seq_record.id)
            desc = str(seq_record.description)
            if i % 1000 == 0:
                print("grabbed the " + str(i+1) +
                "the reference from existing fasta")
            i = i+1
            try:
                # Splitting on space or pipe
                delimiters = "[ ]"  # Split on underscore or comma
                linesplit = re.split(delimiters, desc)
                acc, desc = format_description(line, desc)
                if not args.refresh:
                    seen_in_fasta[acc] = desc
            except Exception as ex:
                print(ex)
                pass
    GCFs_to_skip = dict()
    if not args.refresh:
        for key, value in seen_in_fasta.items():
            if key in gcf_mapping:
                if gcf_mapping[key] not in GCFs_to_skip:
                    GCFs_to_skip[gcf_mapping[key]] = [dict(acc=key, name=value)]
                    print("GCF found in mapping file, skipping", key, "from", gcf_mapping[key])
                else:
                    GCFs_to_skip[gcf_mapping[key]].append(dict(acc=key, name=value))

            else:
                print("No mapping for chr/contig to GCF. Will need to redownload to get gcf mapping for: ", key)
    ## Now, check what contigs/chrs in seen_in_fasta are present in the gcf_mapping
    if args.email:
        Entrez.email = args.email
    if not (args.refresh):
        print(len(seen_in_fasta.keys()), "already seen reference ids")
    # Now use the assembly refseq file to get the ftp path and download the fasta files
    print("Get assemblies now")

    if args.dry_run:
        check_size(assemblies, GCFs_to_skip)
    else:
        get_assemblies(
            assemblies, #this is the top hits from the input file you want to retrieve
            args.file_out, # this is the file you want to write to
            GCFs_to_skip, # this is the list of GCFs you want to skip
            args.refresh # this is the boolean to check if you want to refresh the file
        )
        if args.gcf_map:
            with open(args.gcf_map, "w") as w:
                for organism_key, value in assemblies.items():
                    # organism_key is the dict key you used (organism name)
                    for chr in value["chrs"]:
                        if args.taxids:
                            taxid = value.get("taxid", "")
                            species_taxid = value.get("species_taxid", "")
                            outstring = f"{chr['acc']}\t{value['accession']}\t{organism_key}\t{chr['name']}\t{taxid}\t{species_taxid}"
                        else:
                            outstring = f"{chr['acc']}\t{value['accession']}\t{organism_key}\t{chr['name']}"
                        w.write(outstring + "\n")
            w.close()

if __name__ == "__main__":
    sys.exit(main())
