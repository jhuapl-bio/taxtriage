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
from collections import Counter
from pathlib import Path
import re
import os
import shutil
import urllib.request as request
import ssl
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
        help="List of taxIDs",
    )
    parser.add_argument(
        "-d",
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
        "-f",
        "--type",
        default='file',
        help="Input type, can be a list of taxids or from a file",
    )
    parser.add_argument(
        "-c",
        "--colnumber_file_taxids",
        type=int,
        default=5,
        help="Column number to get taxids from if using a file. Starts at 1st index",
    )
    parser.add_argument(
        "-k",
        "--kraken2output",
        action='store_true',
        help="reformat header for each fasta to a kraken:taxid|id parsing. Requires the setup of the file to be kraken:taxid|taxid|refseqId (whatever text can come after this, separated by space(s)",
    )
    parser.add_argument(
        "-p",
        "--assembly_map_idx",
        default=[0, 5],
        help="Assembly refseq accession and taxid",
    )
    parser.add_argument(
        "-t",
        "--assembly_refseq_file",
        type=Path,
        help="Assembly refseq to pull taxids from instead of pulling straight from a set of IDs. Map the taxid to accession",
    )
    parser.add_argument(
        "-e",
        "--email",
        metavar="EMAIL",
        type=str,
        help="Email for entrez querying, optional",
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
                    header = firstline.split("|")[2]
                    refs[header] = line
                else:
                    header = line.split(" ")[0]
                    refs[header] = line
            except Exception as er:
                pass
    print("Done")
    return refs

#
def import_assembly_file(input, filename, idx):
    refs = dict()
    seen = dict()
    first = dict()


    print("--------------")
    if (not isinstance(input, list)):
        input = input.split(" ")
    with open(filename, "r") as f:
        line = f.readline()
        for line in f:
            line = line.strip()
            linesplit = line.split("\t")

            if len(linesplit) >= 12 and (linesplit[idx[1]] in input) and linesplit[11] == "Complete Genome" and linesplit[idx[1]] not in seen:

                #If the refseq_category column in the assembly.txt is reference genome
                if linesplit[4] == "reference genome":
                    #Set taxid as seen
                    seen[linesplit[idx[1]]] = True
                    #Save reference to dict
                    refs[linesplit[idx[0]]] = dict(id="kraken:taxid|{}|{}".format(
                        linesplit[idx[1]], linesplit[idx[0]]), fulline=linesplit)
                #If there is no reference genome
                else:
                    #If this is the first time the taxa without reference genome is seen
                    if linesplit[idx[1]] not in first:
                        #Save reference to dict
                        refs[linesplit[idx[0]]] = dict(id="kraken:taxid|{}|{}".format(
                            linesplit[idx[1]], linesplit[idx[0]]), fulline=linesplit)
                        #Save taxa as the first to be seen
                        first[linesplit[idx[1]]] = True
                    #If the taxa without reference genome has already been seen previously, pass (save first seen only)
                    elif linesplit[idx[1]] in first:
                        pass
            #If no complete genome found, pass
            else:
                pass
    return refs


def get_assembly_summary(id):
    """Get esummary for an entrez id"""
    esummary_handle = Entrez.esummary(db="assembly", id=id, report="full")
    esummary_record = Entrez.read(esummary_handle)
    return esummary_record


def get_assemblies(refs, outfile, seen, index_ftp):
    """Download genbank assemblies for a given search term.
    Args:
        term: search term, usually organism name
        download: whether to download the results
        path: folder to save to
    """
    ids = refs.keys()
    if index_ftp and len(ids) > 0:
        with open(outfile, "a") as w:
            for id in ids:
                try:
                    if seen and id in seen  :
                        print(id)
                        print("key already seen:", id, "; skipping")
                    else:
                        ftp_site = refs[id]['fulline'][index_ftp]
                        obj = refs[id]['id']
                        fullid = os.path.basename(ftp_site) + '_genomic.fna.gz'
                        ftp_site = ftp_site+'/'+fullid
                        print(ftp_site, id)
                        encoding = guess_type(fullid)[1]   # uses file extension
                        _open = partial(
                            gzip.open, mode='rt') if encoding == 'gzip' else open
                        with closing(request.urlopen(ftp_site, context=ctx)) as r:
                            with open('file.gz', 'wb') as f:
                                shutil.copyfileobj(r, f)
                            f.close()
                        r.close()
                        with _open('file.gz') as uncompressed:
                            for record in SeqIO.parse(uncompressed, "fasta"):
                                if (len(record.seq) > 1000):
                                    newobj = obj+"|"+record.id
                                    record.id = newobj
                                    SeqIO.write(record, w, "fasta")
                        uncompressed.close()
                except Exception as ex:
                    print(ex)
                    pass
        w.close()
    return

    # provide your own mail here
    ids = refs.keys()
    handle = Entrez.efetch(db="assembly", id=ids, retmax='200')
    record = Entrez.read(handle)
    ids = [record[i] for i in range(len(record)) if i % 2 == 0]
    print(f'found {len(ids)} ids')
    links = []
    for id in ids:
        # get summary
        print(id)
        summary = get_assembly_summary(id)
        print("_")
        # get ftp link
        url = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
        if url == '':
            continue
        # label = os.path.basename(url)
        # #get the fasta link - change this to get other formats
        # link = os.path.join(url,label+'_genomic.fna.gz')
        # print (link)
        # links.append(link)
        # if download == True:
        #     #download link
        #     urllib.request.urlretrieve(link, f'{label}.fna.gz')
    return links


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


def get_taxids_from_file(filename, colnumber):
    taxids = []
    with open(filename, "r") as f:
        lines = f.readlines()
        for line in lines:
            line = line.rstrip()
            taxids.append(line.split("\t")[colnumber-1])
    return taxids


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    # logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if args.type == 'file':
        taxids = get_taxids_from_file(args.input, args.colnumber_file_taxids)
    else:
        taxids = args.input
    if (args.assembly_refseq_file):
        refs = import_assembly_file(
            taxids, args.assembly_refseq_file, args.assembly_map_idx)
    else:
        refs = import_genome_file(taxids, args.kraken2output)
    seen = dict()

    i = 0
    if os.path.exists(args.file_out):
        for seq_record in SeqIO.parse(args.file_out, "fasta"):
            line = str(seq_record.id)
            if i % 1000 == 0:
                print("grabbed the " + str(i+1) +
                "the reference from existing fasta")
            i = i+1
            try:
                if (args.kraken2output):
                    firstline = line.split(" ")[0]
                    header = firstline.split("|")[2].replace(">", "")
                    refs[header] = line.replace(">", "")
                else:
                    header = line.split(" ")[0].replace(">", "")
                    firstline = line.replace(">", "")
                    refs[header] = firstline
                seen[header] = True
            except Exception as ex:
                print(ex)
                pass
    if args.email:
        Entrez.email = args.email
    print(len(seen.keys()), "already seen reference ids")
    if (not args.assembly_refseq_file):
        print("downloading refseq file")
        download(refs, args.db, args.file_out, seen)
    else:
        print("get assemblies")
        get_assemblies(refs, args.file_out, seen, args.ftp_path)


if __name__ == "__main__":
    sys.exit(main())
