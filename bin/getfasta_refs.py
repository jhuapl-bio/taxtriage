#!/usr/bin/env python3

"""Provide a command line tool to pull one or more references from either a list of genome ids (kraken) or a file, separated by newline per id."""

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
        help="List of taxids to pull from reference",
    )
    parser.add_argument(
        "-r",
        "--reference",
        metavar="REFERENCE",
        type=Path,
        nargs="+",
        help="Reference fasta files",
    )
    parser.add_argument(
        "-s",
        "--taxid_header_sep",
        metavar="SEP",
        default="|",
        help="Separator to pull the taxids from",
    )
    parser.add_argument(
        "-d",
        "--reads",
        metavar="reads",
        default=None,
        type=Path,
        nargs="+",
        help="Reads file, can be compressed",
    )
    parser.add_argument(
        "-a",
        "--assignment_reads",
        metavar="assignment",
        default=None,
        type=Path,
        help="Reads assignment file, pulled out and sent to a map list",
    )
    parser.add_argument(
        "-x",
        "--pos_taxid_header",
        metavar="pos_taxid_header",
        default=None,
        type=int,
        help="Position header is in, default is first index AFTER kraken:taxid",
    )
    parser.add_argument(
        "-q",
        "--file_reads_out",
        action="store_true",
        help="Output classified reads out filtered",
    )
    parser.add_argument(
        "-o",
        "--file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Name of the output tsv file containing a mapping of top n organisms at individual taxa levels",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    parser.add_argument(
        "-t",
        "--type",
        help="The type of input, can be a file or a list, contained with a newline and separated by spaces.",
        choices=("file", "list"),
        type=str,
        default="file",
    )
    return parser.parse_args(argv)

def import_taxids(filename):
    taxids = []
    with open(filename,"r") as f:
        taxids = [line.strip() for line in f]
    return taxids
def import_filter_fasta(taxids, fastafile, sep, pos):
    
    mapping  = dict()
    if isinstance(taxids, str):
        taxids = [taxids]
    for file in fastafile:
        with open(file,"r") as f:
            for seq_record in SeqIO.parse(file, "fasta"):
                grabbed = seq_record.id.split(sep)
                if pos:
                    idx = grabbed[pos]
                else:
                    idx = grabbed[grabbed.index('kraken:taxid') + 1]
                idx = str(idx)
                if idx in taxids and idx not in mapping:
                    mapping[idx] = seq_record.seq
        f.close()
    return mapping
def write_filtered(outfile, record_dict):
    try:
        # if not os.path.isdir(outfile):
        #     os.mkdir(outfile)
        with open( os.path.join(outfile),"w") as f:
            for taxid, seq in record_dict.items():
                f.write(">" + str(taxid) + "\n")
                f.write(str(seq) + "\n")
        f.close()
        # for taxid, seq in record_dict.items():
        #     with open( os.path.join(outfile, str(taxid)+".fa"),"w") as f:
        #         f.write(">" + str(taxid) + "\n")
        #         f.write(str(seq) + "\n")
        #     f.close()
        
    except OSError as error:
        print(error)
def get_fastq_filtered(filtered_taxids, reads, assignment_reads, file_reads_out):
    
    classified_reads  = dict()
    with open(assignment_reads,"r") as f:
        for line in f.readlines():
            splitline = line.rstrip().split("\t")
            classified_reads[(splitline[1])]= splitline[2]
    f.close()
    i=0
    for read in reads:
        encoding = guess_type(read)[1]  # uses file extension
        _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
        filename = os.path.join(os.path.dirname(read), str(i)+"_filtered.fastq")
        with open(filename, "w") as w:
            with _open(read) as f:
                seen = dict()
                g = 0
                for seq_record in SeqIO.parse(f, "fastq"):
                    if seq_record.id in classified_reads and str(classified_reads[seq_record.id]) in filtered_taxids:
                        SeqIO.write(seq_record,w,"fastq")
                        taxid = str(classified_reads[seq_record.id])
                        if taxid not in seen:
                            seen[taxid] = []
                        seen[taxid].append(g)
                        g = g + 1
            f.close()
            # for taxid,value in seen.items():
            #     with open(os.path.join(file_reads_out[i], taxid+".fastq"), "w") as f:
            #         o = 0
            #         with _open(read) as t:
            #             for seq_record in SeqIO.parse(t, "fastq"):
            #                 if o in value:
            #                     SeqIO.write(seq_record,f,"fastq")
            #                 o = 1+o
            #         t.close()
                    
            #     f.close()

        w.close()
        i+=1    
def main(argv=None):
    """Coordinate argument parsing and program execut      ion."""
    args = parse_args(argv)   
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if  args.type == 'file' and not args.input.is_file() :
        logger.error(f"The given input file {args.input} was not found!")
        sys.exit(2)
    elif args.type == 'file' and args.input.is_file():
        logger.info("File exists, importing and filtering")
        taxids = import_taxids(args.input)
        filtered_taxids = import_filter_fasta(taxids, 
            args.reference, 
            args.taxid_header_sep, 
            args.pos_taxid_header
        )
        write_filtered(args.file_out, filtered_taxids)
        if args.reads and args.assignment_reads:
            get_fastq_filtered(filtered_taxids, args.reads, args.assignment_reads, args.file_reads_out)
    elif args.type == 'list' and len(args.input) > 0:
        filtered_taxids = import_filter_fasta(
            args.input.split(" "), 
            args.reference, 
            args.taxid_header_sep, 
            args.pos_taxid_header
        )
        
        write_filtered(args.file_out, filtered_taxids)
        if args.reads and args.assignment_reads:
            get_fastq_filtered(filtered_taxids, args.reads, args.assignment_reads, args.file_reads_out)
    elif not args.input and args.type == 'list':
        logger.error(f"The given input list of taxids: {args.input} was not found!")
        sys.exit(2)
    
if __name__ == "__main__":
    sys.exit(main())
