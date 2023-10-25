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
        "-r",
        "--reference",
        metavar="REFERENCE",
        type=Path,
        nargs="+",
        help="Reference fasta files",
    )
    parser.add_argument(
        "-p",
        "--samplename",
        metavar="SAMPLENAME",
        default=None,
        help="Append Samplename to output filenames of taxids",
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
    parentsdict = dict()
    with open(filename, "r") as f:
        for line in f:
            splitt = line.strip().split("\t")
            taxid = splitt[4]
            if len(splitt) >= 7:
                parents = splitt[6]
            else:
                parents = ""
            parentsdict[taxid] = parents.split(";")
            taxids.append(taxid)
    return taxids, parentsdict


def import_filter_fasta(taxids, fastafile, sep, pos, parents):
    mapping = dict()
    if isinstance(taxids, str):
        taxids = [taxids]
    for file in fastafile:
        with open(file, "r") as f:
            for seq_record in SeqIO.parse(file, "fasta"):
                grabbed = seq_record.id.split(sep)
                if pos:
                    idx = grabbed[pos]
                else:
                    idx = grabbed[grabbed.index('kraken:taxid') + 1]
                idx = str(idx)
                if idx in taxids:
                    if not idx in mapping:
                        mapping[idx] = []
                    mapping[idx].append([seq_record.id, seq_record.seq])
        f.close()
    return mapping


def write_filtered(outdir, record_dict, sample_base):
    try:
        filehandles = dict()
        if not sample_base:
            sample_base = ""
        elif len(record_dict.keys()) > 0:
            sample_base = "."+sample_base
        for unique_taxids in record_dict.keys():
            outfile = os.path.join(outdir, unique_taxids+sample_base+".fasta")

            filehandles[unique_taxids] = open(outfile, "w")
        if len(record_dict.keys()) <= 0:
            with open(os.path.join(outdir, sample_base+".fasta"), 'w') as fp:
                pass
        for taxid, entries in record_dict.items():
            for seq in entries:
                filehandles[taxid].write(">" + str(seq[0]) + "\n")
                filehandles[taxid].write(str(seq[1]) + "\n")
        for taxid in filehandles.values():
            taxid.close()

    except OSError as error:
        print(error)


def get_fastq_filtered(filtered_taxids, reads, assignment_reads, output, parents):
    classified_reads = dict()
    lowers = dict()
    with open(assignment_reads, "r") as f:
        for line in f.readlines():
            splitline = line.rstrip().split("\t")
            classified_reads[(splitline[1])] = splitline[2]
    f.close()
    # print(classified_reads.keys())
    i = 0
    classified_read_ids = list(classified_reads.keys())
    for read in reads:
        encoding = guess_type(read)[1]  # uses file extension

        sample_base = Path(read).stem.split(".")[0]
        _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
        filename = os.path.join(output, sample_base +
                                "_"+str(i+1)+"_filtered.fastq")
        key_list = list(filtered_taxids.keys())
        with open(filename, "w") as w:
            try:
                with _open(read) as f:
                    seen = dict()
                    g = 0

                    i = 0
                    for seq_record in SeqIO.parse(f, "fastq"):
                        # regex replace everything after the first "." for the read id
                        newid = ""
                        try:
                            newid = re.sub(r'\/.*$', '', seq_record.id)
                        except:
                            newid = seq_record.id
                        if newid in classified_reads or seq_record.id in classified_reads:
                            seen = False

                            # if (classified_reads[newid] in parents):

                            i += 1

                            SeqIO.write(seq_record, w, "fastq")

                            g = g + 1
                f.close()
            except Exception as ex:
                print(ex, "failed with file", filename, read)
                pass
        w.close()
        i += 1


def main(argv=None):
    """Coordinate argument parsing and program execut      ion."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level,
                        format="[%(levelname)s] %(message)s")
    # if  args.type == 'file' and not args.input.is_file() :
    #     logger.error(f"The given input file {args.input} was not found!")
    #     sys.exit(2)
    os.makedirs(args.dir_out, exist_ok=True)
    if args.type == 'file':
        logger.info("File exists, importing and filtering")
        taxids, parents = import_taxids(args.input)
        filtered_taxids = import_filter_fasta(taxids,
            args.reference,
            args.taxid_header_sep,
            args.pos_taxid_header,
            parents
        )
        # print(taxids, "\n\n\n", parents)
        write_filtered(args.dir_out, filtered_taxids, args.samplename)
        # if args.reads and args.assignment_reads:
        #     get_fastq_filtered(filtered_taxids, args.reads, args.assignment_reads, args.dir_out, taxids)
    elif args.type == 'list' and len(args.input) > 0:
        filtered_taxids = import_filter_fasta(
            args.input.split(" "),
            args.reference,
            args.taxid_header_sep,
            args.pos_taxid_header
        )
        write_filtered(args.dir_out, filtered_taxids, args.samplename)
        # if args.reads and args.assignment_reads:
        #     get_fastq_filtered(filtered_taxids, args.reads, args.assignment_reads, args.file_reads_out,[])
    elif not args.input and args.type == 'list':
        logger.error(
            f"The given input list of taxids: {args.input} was not found!")
        sys.exit(2)
    return


if __name__ == "__main__":
    sys.exit(main())
