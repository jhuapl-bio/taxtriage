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

import pysam
import argparse

import sys
def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Split a vcf file on the individual taxa, and call the index",
        epilog="Example: ",
    )
    parser.add_argument(
        "-i",
        "--file_in",
        metavar="FILE_IN",
        type=str,
        help="BAMFILE input",
    )
    parser.add_argument(
        "-o",
        "--output_prefix",
        metavar="OUTPREFIX",
        type=str,
        help="Output prefix",
    )

    return parser.parse_args(argv)

import gzip
import os

# Open the input BAM file
def split_vcf(input_vcf_path, output_prefix):
    seen_files = {}
    is_compressed = input_vcf_path.endswith('.gz')
    uncompressed_file_paths = []
    openfile = None
    lasttaxid = None
    try:
        vcf = pysam.VariantFile(input_vcf_path)
        # with (gzip.open(input_vcf_path, 'rt') if is_compressed else open(input_vcf_path, 'r')) as vcf_file:
            # Open the VCF with pysam
        header = vcf.header
        for record in vcf:
            # Extract the header
            header = vcf.header
            for record in vcf:
                taxid = record.chrom.split('|')[0]
                if lasttaxid != taxid:
                    if taxid in seen_files:
                        print("already seen:", taxid)
                    # print(lasttaxid, ">",taxid)
                    lasttaxid = taxid
                # Open new file if taxid not encountered before
                if taxid not in seen_files:
                    openfile = pysam.VariantFile(f"{output_prefix}{taxid}.vcf.gz", "w", header=header)
                    seen_files[taxid] = 1
                # else:
                    # openfile.close()
                    # openfile = pysam.VariantFile(f"{output_prefix}{taxid}.vcf.gz", "w", header=header)
                # Write record to appropriate VCF file
                # if taxid == '83333':
                openfile.write(record)

    finally:
        if openfile:
            openfile.close()





def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    split_vcf(args.file_in, args.output_prefix)

if __name__ == "__main__":
    sys.exit(main())
