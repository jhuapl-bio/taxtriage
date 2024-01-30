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

import pandas as pd
import argparse


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
        help="List of features with positions and descriptions/text associated",
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="OUTPUT",
        help="BED File output",
    )
    return parser.parse_args(argv)


args = parse_args()

# Load the gene feature table
# Adjust the file path and column names according to your file's format
gene_feature_table = pd.read_csv(args.input, sep='\t')
gene_feature_table = gene_feature_table[gene_feature_table['# feature'] == 'CDS']

gene_feature_table=gene_feature_table[['genomic_accession','start','end', 'name' ]]
# Select the relevant columns (chromosome, start, end, strand, gene name)
# Adjust the column indices (0, 1, 2, ...) to match the layout of your file
# bed_data = gene_feature_table[[6, 7, 8, 9, 13]]  # Example: chromosome, start, end, strand, gene_name
# filter on CDS in # feature col
print(gene_feature_table)
# Convert to 0-based start position for BED format
# bed_data[1] = bed_data[1] - 1

# print(bed_data.head(20))
# Save to BED file
gene_feature_table.to_csv(args.output, sep='\t', header=False, index=False)
