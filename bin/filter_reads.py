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

"""Provide a command line tool to pull one or more references from either a list of genome ids (kraken) or a file, separated by newline per id."""

from functools import partial
from mimetypes import guess_type
from Bio import SeqIO
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
        "--dir_out",
        metavar="DIR_OUT",
        type=Path,
        help="Name of the output directory to place the filtered fastq file(s) into",
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
    with open(filename, "r") as f:
        taxids = [line.strip().split("\t")[4] for line in f]
    return taxids


def write_filtered(outfile, record_dict):
    try:
        # if not os.path.isdir(outfile):
        #     os.mkdir(outfile)
        with open(os.path.join(outfile), "w") as f:
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


def get_fastq_filtered(filtered_taxids, reads, assignment_reads, outputdir):

    classified_reads = dict()
    with open(assignment_reads, "r") as f:
        for line in f.readlines():
            splitline = line.rstrip().split("\t")
            classified_reads[(splitline[1])] = splitline[2]
    f.close()
    i = 0
    for read in reads:
        encoding = guess_type(read)[1]  # uses file extension
        sample_base = Path(read).stem.split(".")[0]
        filehandles = dict()
        matched = dict()
        g = 0
        for taxid in filtered_taxids:
            matched[taxid] = []
        _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
        os.makedirs(outputdir, exist_ok=True)
        filename = os.path.join(
            outputdir, sample_base+"_"+str(i+1)+"_filtered.fastq")
        with open(filename, "w") as w:
            try:
                with _open(read) as f:
                    for seq_record in SeqIO.parse(f, "fastq"):
                        if seq_record.id in classified_reads and str(classified_reads[seq_record.id]) in filtered_taxids:
                            taxid = str(classified_reads[seq_record.id])
                            SeqIO.write(seq_record, w, "fastq")
                        g = g + 4
                f.close()
            except Exception as ex:
                print(ex, "failed with file")
                pass
        i += 1


def main(argv=None):
    """Coordinate argument parsing and program execut      ion."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level,
                        format="[%(levelname)s] %(message)s")
    # if  args.type == 'file' and not args.input.is_file() :
    #     logger.error(f"The given input file {args.input} was not found!")
    #     sys.exit(2)
    taxids = import_taxids(args.input)
    if args.reads and args.assignment_reads:
        get_fastq_filtered(taxids, args.reads,
            args.assignment_reads, args.dir_out)
    return


if __name__ == "__main__":
    sys.exit(main())
