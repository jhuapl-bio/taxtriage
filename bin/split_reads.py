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
        description="Split a tsv-given fastq file into individual fastq files",
        epilog="",
    )
    parser.add_argument(
        "-m",
        "--input_metadata",
        metavar="INPUT",
        help="2 column tsv file that is Read ID AND second column as the category it belongs to. ",
    )
    parser.add_argument(
        "-f",
        "--genome_file",
        metavar="INPUT",
        help="Split the FASTA file",
    )
    parser.add_argument(
        "-q",
        "--reads",
        nargs="+",
        metavar="READS",
        help="Fastq File(s)",
    )
    parser.add_argument(
        "-k",
        "--kraken2_genome_format",
        action="store_true",
        help="If called, the second column from the metadata is in kraken2 format so will pull the 2nd column of that format e.g. Kraken:taxid|1366|GCF|NZ",
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

    return parser.parse_args(argv)


def split_fastqs(metadata, reads, outputdir):

    classified_reads = dict()
    taxids = []
    with open(metadata, "r") as f:
        for line in f.readlines():
            splitline = line.rstrip().split("\t")
            if splitline[1] not in taxids:
                taxids.append(splitline[1])
            classified_reads[(splitline[0])] = splitline[1]
    f.close()
    outputfilehandles = dict()

    i = 0
    for read in reads:
        encoding = guess_type(read)[1]  # uses file extension
        sample_base = Path(read).stem.split(".")[0]
        g = 0
        _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
        os.makedirs(outputdir, exist_ok=True)
        for f in taxids:
            outpath = os.path.join(
                outputdir, f+"_"+sample_base+"_"+str(i+1)+"_filtered.fastq")
            if not os.path.exists(os.path.dirname(outpath)):
                os.makedirs(os.path.dirname(outpath))
            outputfilehandles[f] = open(outpath, "w")
        try:
            with _open(read) as f:
                for seq_record in SeqIO.parse(f, "fastq"):
                    if seq_record.id in classified_reads:
                        taxid = str(classified_reads[seq_record.id])
                        if (taxid in taxids):
                            SeqIO.write(
                                seq_record, outputfilehandles[taxid], "fastq")
                    g = g + 4
            f.close()
        except Exception as ex:
            print(ex, "failed with file")
            pass
        finally:
            for f in outputfilehandles.values():
                f.close()
        i += 1


def main(argv=None):
    """Coordinate argument parsing and program execut      ion."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level,
                        format="[%(levelname)s] %(message)s")
    # if  args.type == 'file' and not args.input.is_file() :
    #     logger.error(f"The given input file {args.input} was not found!")
    #     sys.exit(2)
    split_fastqs(args.input_metadata, args.reads, args.dir_out)
    return


if __name__ == "__main__":
    sys.exit(main())
