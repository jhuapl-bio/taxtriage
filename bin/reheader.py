

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
import argparse
import sys
import pysam  # Ensure pysam is installed


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
        help="INPUT fasta file",
    )
    parser.add_argument(
        "-m",
        "--mapping",
        metavar="MAP",
        help="MAPPING FILE",
    )
    return parser.parse_args(argv)

def main():
    args = parse_args()
    input_bam = args.input
    mapping_tsv = args.mapping

    # Load the mapping from the TSV file
    mapping = {}
    with open(mapping_tsv, 'r') as f:
        for line in f:
            old, new, chr = line.strip().split('\t')
            mapping[old] = chr

    # Read the current header from the BAM file
    bam = pysam.AlignmentFile(input_bam, "rb")
    header_dict = bam.header.to_dict()

    # Update the reference names in the header
    for ref in header_dict['SQ']:
        old_name = ref['SN']
        if old_name in mapping:
            ref['SN'] = mapping[old_name]
    bam.close()
    header = pysam.AlignmentHeader.from_dict(header_dict)
    print(header)
    # # Write the new header to a SAM file
    # with open(new_header_sam, 'w') as f:
    #     f.write(pysam.AlignmentHeader.from_dict(header_dict).to_string())

if __name__ == "__main__":
    sys.exit(main())
