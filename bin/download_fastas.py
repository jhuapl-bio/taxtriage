#!/usr/bin/env python3

"""Provide a command line tool to fetch a list of refseq genome ids to a single file, useful for kraken2 database building or alignment purposes"""

import gzip

import argparse
import csv
import logging
import sys
from collections import Counter
from pathlib import Path
import re
import os
from tabnanny import filename_only
from tokenize import String
from xmlrpc.client import Boolean 
logger = logging.getLogger()
from Bio import SeqIO
from Bio import Entrez
from mimetypes import guess_type
from functools import partial

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
        help="List of refseq IDs",
    )
    parser.add_argument(
        "-d",
        "--db",
        metavar="DB",
        default="nuccore",
        help="Database of choice to pull IDs from",
    )
    parser.add_argument(
        "-k",
        "--kraken2output",
        action='store_true',
        help="reformat header for each fasta to a kraken:taxid|id parsing. Requires the setup of the file to be kraken:taxid|taxid|refseqId (whatever text can come after this, separated by space(s)",
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
    with open(filename,"r") as f:
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
def download(refs, db, outfile, seen):
    with open(outfile, "a") as w:
        i = 0
        maxt = 30
        next_ = []
        for key, value in refs.items():
            try:
                if i % maxt == 0 and len(next_)>0:
                    print( str(i), " th iteration of ids to submit..", next_)
                    handle = Entrez.efetch(db=db, rettype="fasta", retmode="fasta", id=",".join(next_), idtype="acc")
                    seq_records = SeqIO.parse(handle, 'fasta') 
                    for seq_record in seq_records:
                        if (seq_record):
                            seq_record.id = str(value.replace(">", ""))
                            SeqIO.write(seq_record, w, "fasta")
                    next_ = []
                    handle.close()
                if seen and key  in seen:
                    print("key already seen:", key, "; skipping")
                    # exit()
                else:
                    next_.append(key)
            except Exception as err:
                print("No seq record found", next_, err)
                next_ = []
                pass
            i = i+1
                
def main(argv=None):
    """Coordinate argument parsing and program execut      ion."""
    args = parse_args(argv)   
    # logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    refs = import_genome_file(args.input, args.kraken2output)
    seen = dict()
    
    if os.path.exists(args.file_out):
        for seq_record in SeqIO.parse(args.file_out, "fasta"):
            line = str(seq_record.id)
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
                pass;
    if args.email:
        Entrez.email = args.email
    download(refs, args.db, args.file_out, seen)
    
if __name__ == "__main__":
    sys.exit(main())
