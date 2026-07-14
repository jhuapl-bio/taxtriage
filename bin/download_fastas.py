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
    parser.add_argument(
        "--accession_col",
        type=int,
        default=None,
        help="1-based column in the input file holding a specific assembly accession "
             "(e.g. GCF_000002415.2) to use FIRST for each row. When set, the row is "
             "resolved by matching this accession against column --assembly_map_idx of "
             "the local assembly summary; if empty or not found it falls back to taxid "
             "matching. Starts at 1st index.",
    )
    parser.add_argument(
        "--accession_map",
        type=Path,
        default=None,
        help="Optional CSV (e.g. assets/pathogen_sheet.csv) with 'taxid' and "
             "'assembly_accession' columns (and optionally 'name'). Provides a "
             "curated taxid/name -> assembly accession lookup used FIRST for each "
             "input entry, before falling back to local taxid matching. Complements "
             "--accession_col.",
    )
    parser.add_argument(
        "--ncbi_fallback",
        action='store_true',
        default=False,
        help="For entries not resolved by the local assembly summary, query NCBI Entrez "
             "(by accession first, then by taxid) to find and download an assembly. "
             "Recommended to also pass -e/--email.",
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
    requested_accessions=None,
    taxid_to_accession=None,
):
    """Resolve each input taxid to an assembly using the local assembly summary.

    Resolution order per input entry (accession-first):
      1. If the entry has a specific assembly accession (via taxid_to_accession)
         and that accession is present in the local assembly summary
         (column `idx`), use that exact assembly.
      2. Otherwise fall back to the best-priority assembly whose taxid/name
         column (`matchcol`) matches the entry (original behavior).
    Entries resolved by neither are returned as `unresolved` so the caller can
    optionally attempt an NCBI lookup.
    """
    if requested_accessions is None:
        requested_accessions = set()
    if taxid_to_accession is None:
        taxid_to_accession = {}

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
        # remove any trailing / that is empty from utl
        utl_formatted = utl.rstrip("/")
        bb = os.path.basename(utl_formatted)
        return utl_formatted + "/" + bb + "_genomic.fna.gz"

    def build_obj(gcfidx, formatted_header, namecol, urlcol, taxidcol,
                  species_taxidcol, characteristic):
        return dict(
            id="{}|{}".format(gcfidx, formatted_header),
            accession=gcfidx,
            characteristic=characteristic,
            chrs=[],
            organism=namecol,
            reference=get_url(urlcol, gcfidx),
            name=formatted_header,
            taxidcol=taxidcol,
            species_taxidcol=species_taxidcol,
        )

    # priorities[matchval][prio] = (file_index, obj)   -> taxid fallback buckets
    priorities = {}
    first = {}
    accession_hits = {}   # accession -> obj (direct accession-first match)
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
                # Strip only the line ending (not spaces) so that multi-space
                # fields such as asm_submitter are preserved, then split
                # strictly on tab.
                line = line.rstrip("\r\n")
                if not line.strip():
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

                # optional filter
                if skip_incomplete and linesplit[11] != "Complete Genome":
                    continue

                # --- Accession-first capture ---
                # If this row's accession was explicitly requested, remember it
                # regardless of whether its taxid is in the input list.
                if gcfidx in requested_accessions and gcfidx not in accession_hits:
                    accession_hits[gcfidx] = build_obj(
                        gcfidx, formatted_header, namecol, urlcol,
                        taxidcol, species_taxidcol, "accession",
                    )

                # --- Taxid/name priority buckets (fallback) ---
                # find first match in input order
                match_index = -1
                for match in matchidces:
                    if match in input:
                        match_index = input.index(match)
                        break
                if match_index < 0:
                    continue

                matchval = input[match_index]
                seencols.add(matchval)

                if matchval not in priorities:
                    priorities[matchval] = {}
                obj = build_obj(
                    gcfidx, formatted_header, namecol, urlcol,
                    taxidcol, species_taxidcol, None,
                )
                # assign priority bucket
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

    # pick best per taxid matchval: lowest prio, then earliest file index
    taxid_best = {}
    for matchval, prmap in priorities.items():
        for prio in ("0", "1", "2", "3"):
            if prmap.get(prio):
                taxid_best[matchval] = prmap[prio][1]  # obj
                break

    # Per-input selection: accession-first, then local taxid fallback.
    chosen = []
    unresolved = []
    for taxid in input:
        acc = taxid_to_accession.get(taxid)
        if acc and acc in accession_hits:
            chosen.append(accession_hits[acc])
        elif taxid in taxid_best:
            chosen.append(taxid_best[taxid])
        else:
            print("No local assembly found for", taxid,
                  ("(accession " + acc + " not in summary)") if acc else "")
            unresolved.append(dict(taxid=taxid, accession=acc or ""))

    # format as you did: keyed by organism
    assembliesformat = {}
    for item in chosen:
        assembliesformat[item["organism"]] = dict(
            accession=item["accession"],
            id=item["id"],
            name=item["name"],
            reference=item["reference"],
            chrs=item["chrs"],
            taxid=item["taxidcol"],   # taxid
            species_taxid=item["species_taxidcol"],   # species taxid
        )
    if unresolved:
        print("Missing locally:", [u["taxid"] for u in unresolved])
        if missingfile:
            with open(os.path.join(missingfile), "w") as w:
                for u in unresolved:
                    w.write(f"{u['taxid']}\n")

    return assembliesformat, unresolved


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


def get_accession_map_from_file(filenames, taxid_colnumber, accession_colnumber):
    """Build {taxid: accession} and the set of requested accessions from the
    input file. Rows with a blank accession are skipped (taxid-only fallback)."""
    mapping = {}
    accessions = set()
    for filename in filenames:
        with open(filename, "r") as f:
            for line in f:
                line = line.rstrip("\r\n")
                if not line.strip():
                    continue
                cols = line.split("\t")
                if len(cols) < taxid_colnumber:
                    continue
                taxid = cols[taxid_colnumber - 1].strip()
                acc = ""
                if accession_colnumber and len(cols) >= accession_colnumber:
                    acc = cols[accession_colnumber - 1].strip()
                if taxid and acc and taxid not in mapping:
                    mapping[taxid] = acc
                    accessions.add(acc)
    return mapping, accessions


def load_accession_map_csv(path):
    """Build {taxid|name: accession} and the set of accessions from a curated
    CSV (name-based) such as assets/pathogen_sheet.csv. Keyed by BOTH taxid and
    name so it works whether the download input matches on taxid or on name.
    Rows with a blank accession are skipped."""
    mapping = {}
    accessions = set()
    try:
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                acc = (row.get("assembly_accession") or "").strip()
                if not acc:
                    continue
                taxid = (row.get("taxid") or "").strip()
                name = (row.get("name") or "").strip()
                if taxid:
                    mapping.setdefault(taxid, acc)
                if name:
                    mapping.setdefault(name, acc)
                accessions.add(acc)
    except Exception as ex:
        print(f"Could not read accession map {path}: {ex}")
    return mapping, accessions


def _assembly_obj_from_docsum(doc, fallback_taxid=None):
    """Convert an NCBI assembly esummary DocumentSummary into the internal
    assembly object shape used by get_assemblies / the gcf mapping writer."""
    ftp = doc.get("FtpPath_RefSeq") or doc.get("FtpPath_GenBank") or ""
    if not ftp:
        return None
    acc = doc.get("AssemblyAccession", "")
    org = doc.get("SpeciesName") or doc.get("Organism") or acc
    taxid = str(doc.get("Taxid") or doc.get("SpeciesTaxid") or fallback_taxid or "")
    species_taxid = str(doc.get("SpeciesTaxid") or taxid)
    formatted = str(org).replace(" ", "_")
    base = os.path.basename(ftp.rstrip("/"))
    reference = ftp.rstrip("/") + "/" + base + "_genomic.fna.gz"
    return dict(
        organism=org,
        accession=acc,
        id="{}|{}".format(acc, formatted),
        name=formatted,
        reference=reference,
        chrs=[],
        taxidcol=taxid,
        species_taxidcol=species_taxid,
    )


def ncbi_lookup(accession=None, taxid=None):
    """Query NCBI Entrez assembly db by accession (preferred) or taxid and
    return an internal assembly object, or None if nothing usable is found."""
    try:
        if accession:
            term = accession
        elif taxid:
            term = "txid{}[Organism:exp]".format(taxid)
        else:
            return None
        handle = Entrez.esearch(db="assembly", term=term, retmax=1)
        rec = Entrez.read(handle)
        handle.close()
        ids = rec.get("IdList", [])
        if not ids:
            return None
        summary = get_assembly_summary(ids[0])
        docs = summary["DocumentSummarySet"]["DocumentSummary"]
        if not docs:
            return None
        return _assembly_obj_from_docsum(docs[0], fallback_taxid=taxid)
    except Exception as ex:
        print("NCBI lookup failed for", accession or taxid, ex)
        return None


def resolve_via_ncbi(unresolved):
    """For each unresolved entry, try NCBI by accession first, then by taxid.
    Returns a dict keyed by organism (same shape as import_assembly_file)."""
    resolved = {}
    for entry in unresolved:
        acc = entry.get("accession")
        taxid = entry.get("taxid")
        obj = None
        if acc:
            obj = ncbi_lookup(accession=acc)
        if obj is None and taxid:
            obj = ncbi_lookup(taxid=taxid)
        if obj is None:
            print("Could not resolve via NCBI:", entry)
            continue
        resolved[obj["organism"]] = dict(
            accession=obj["accession"],
            id=obj["id"],
            name=obj["name"],
            reference=obj["reference"],
            chrs=obj["chrs"],
            taxid=obj["taxidcol"],
            species_taxid=obj["species_taxidcol"],
        )
        time.sleep(0.34)  # be polite to NCBI (~3 req/s)
    return resolved
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
    requested_accessions = set()
    taxid_to_accession = {}
    if args.type == 'file':
        #colnumber file hits is the column from the input top hits file you want to match to the args.assembly_names
        seen_in_tops = get_hits_from_file(args.input, args.colnumber_file_hits)
        if args.accession_col:
            taxid_to_accession, requested_accessions = get_accession_map_from_file(
                args.input, args.colnumber_file_hits, args.accession_col
            )
    else:
        seen_in_tops = args.input

    # Optional curated taxid/name -> accession lookup (e.g. the pathogen sheet).
    # Merged with anything from --accession_col; the per-row accession column
    # takes precedence when both provide a value for the same key.
    if args.accession_map:
        map_from_csv, acc_from_csv = load_accession_map_csv(args.accession_map)
        for key, acc in map_from_csv.items():
            taxid_to_accession.setdefault(key, acc)
        requested_accessions |= acc_from_csv

    if taxid_to_accession:
        print(f"Accession-first enabled: {len(requested_accessions)} curated accessions available")

    assemblies, unresolved = import_assembly_file(
        seen_in_tops, args.assembly_metadata, args.assembly_names, args.assembly_map_idx,
        args.name_col_assembly, args.ftp_path, args.missingfile, args.skip_incomplete,
        requested_accessions=requested_accessions, taxid_to_accession=taxid_to_accession,
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

    # NCBI fallback: for entries the local assembly summary could not resolve,
    # query NCBI Entrez (by accession first, then taxid) and merge results in.
    if args.ncbi_fallback and unresolved:
        if not args.email:
            print("Warning: --ncbi_fallback set without -e/--email; NCBI may throttle requests")
        print(f"Attempting NCBI fallback for {len(unresolved)} unresolved entries...")
        ncbi_assemblies = resolve_via_ncbi(unresolved)
        for org, obj in ncbi_assemblies.items():
            if org not in assemblies:
                assemblies[org] = obj
        print(f"Resolved {len(ncbi_assemblies)} of {len(unresolved)} via NCBI")

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
