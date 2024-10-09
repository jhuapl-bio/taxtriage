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

import argparse
import sys
import gzip

def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Import a list of taxids, find all classified k2 reads from outfile with at least one occurrence of that taxid",
    )
    parser.add_argument(
        "-i",
        "--taxids",
        required=True,
        metavar="TAXIDS",
        help="Taxid list, one taxid per line",
    )
    parser.add_argument(
        "-k",
        "--k2",
        required=True,
        metavar="K2OUT",
        help="Kraken2 output per-read file",
    )
    parser.add_argument(
        "-1",
        "--in1",
        metavar="READ1",
        help="First pair of paired-end fastq.gz file or single-end fastq.gz file",
    )
    parser.add_argument(
        "-2",
        "--in2",
        metavar="READ2",
        help="Second pair of paired-end fastq.gz file (optional)",
    )
    parser.add_argument(
        "-o1",
        "--out1",
        required=True,
        metavar="OUTPUT1",
        help="Output file for single-end or first pair of paired-end",
    )
    parser.add_argument(
        "-o2",
        "--out2",
        metavar="OUTPUT2",
        help="Output file for second pair of paired-end (optional)",
    )
    return parser.parse_args(argv)
def open_file(file):
    """Open file, with support for gzip if the file ends with .gz."""
    if file.endswith(".gz"):
        return gzip.open(file, "rt")  # Open as text mode for .gz files
    else:
        return open(file, "r")

def parse_kraken2(k2_file, taxid_list):
    """Parse Kraken2 output and return a set of read IDs that match the taxid list."""
    read_ids = set()
    with open(k2_file, 'r') as infile:
        for line in infile:
            columns = line.strip().split()
            taxid_info = columns[-1]
            read_id = columns[1]  # Assuming second column is the read ID
            if any(taxid in taxid_info for taxid in taxid_list):
                read_ids.add(read_id)
    return read_ids

def extract_reads(fastq_file, output_file, read_ids):
    """Extract matching reads from fastq.gz file and write to new file."""
    with gzip.open(output_file, "wt") as outgz:
        with gzip.open(fastq_file, "rt") as infq:
            write_read = False
            for line in infq:
                if line.startswith('@'):
                    # Extract the read ID and remove /1 or /2 suffix for paired-end reads
                    raw_read_id = line.split()[0][1:]  # Remove '@'
                    read_id = raw_read_id.split('/')[0]  # Remove /1 or /2 suffix
                    write_read = read_id in read_ids
                if write_read:
                    outgz.write(line)

def extract_paired_reads(fastq1, fastq2, out1, out2, read_ids):
    """Extract paired reads from fastq.gz files and write to new files."""
    with gzip.open(out1, "wt") as outgz1, gzip.open(out2, "wt") as outgz2:
        with gzip.open(fastq1, "rt") as infq1, gzip.open(fastq2, "rt") as infq2:
            write_read1 = False
            write_read2 = False
            for line1, line2 in zip(infq1, infq2):
                if line1.startswith('@') and line2.startswith('@'):
                    # Extract the read IDs and remove /1 or /2 suffix for paired-end reads
                    raw_read_id1 = line1.split()[0][1:]  # Remove '@'
                    raw_read_id2 = line2.split()[0][1:]  # Remove '@'
                    read_id1 = raw_read_id1.split('/')[0]  # Remove /1 or /2 suffix
                    read_id2 = raw_read_id2.split('/')[0]  # Remove /1 or /2 suffix
                    write_read1 = write_read2 = read_id1 in read_ids and read_id2 in read_ids
                if write_read1:
                    outgz1.write(line1)
                if write_read2:
                    outgz2.write(line2)

def main(argv=None):
    args = parse_args(argv)

    # Read the taxid list
    taxidfile = args.taxids
    taxid_list = []
    with open(taxidfile, 'r') as f:
        taxid_list = [line.strip() for line in f]

    # Parse Kraken2 classified reads file and get matching read IDs
    read_ids = parse_kraken2(args.k2, taxid_list)

    # If paired-end files provided, extract reads for both R1 and R2
    if args.in2 and args.out2:
        print("Extracting paired reads...")
        extract_paired_reads(args.in1, args.in2, args.out1, args.out2, read_ids)
    else:
        # Otherwise, handle single-end or R1-only reads
        extract_reads(args.in1, args.out1, read_ids)


if __name__ == "__main__":
    sys.exit(main())
