#!/usr/bin/env python3

##############################################################################################
# Copyright 2024 The Johns Hopkins University Applied Physics Laboratory LLC
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
import pandas as pd
from pathlib import Path

def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Map GCF accession to taxids using input and reference files.",
        epilog="Example: python3 map_gcf_to_taxid.py -i input.txt -r reference.txt -o output.tsv -c taxid",
    )
    parser.add_argument(
        "-i",
        "--file_in",
        metavar="FILE_IN",
        type=Path,
        required=True,
        help="Tabular input txt file",
    )
    parser.add_argument(
        "-r",
        "--ref_file",
        metavar="REF_FILE",
        type=Path,
        required=True,
        help="Reference file with GCF accession and taxid",
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="OUTPUT",
        type=Path,
        required=False,
        help="Output file to write the mapping",
    )
    parser.add_argument(
        "-c",
        "--column",
        metavar="COLUMN",
        type=str,
        required=False,
        default="taxid",
        help="Column name in the reference file to map to (e.g., taxid, species_taxid)",
    )
    return parser.parse_args(argv)

def read_input_file(input_file):
    return pd.read_csv(input_file, sep='\t', header=None, names=["Acc", "Assembly", "Organism_Name", "Description"])

def read_reference_file(ref_file):
    # skip the first line for ref_file
    return pd.read_csv(ref_file, sep='\t', skiprows=1)

def map_gcf_to_taxid(input_df, ref_df, column):
    # Create a dictionary from the reference file for mapping
    gcf_to_taxid = ref_df.set_index('#assembly_accession')[column].to_dict()

    # Map the GCF to the desired column (taxid or specified column)
    input_df['Mapped_Value'] = input_df['Assembly'].map(gcf_to_taxid)

    return input_df

def main(argv=None):
    args = parse_args(argv)

    input_df = read_input_file(args.file_in)
    ref_df = read_reference_file(args.ref_file)

    mapped_df = map_gcf_to_taxid(input_df, ref_df, args.column)
    # Write the output to the specified file
    if args.output:
        mapped_df[['Acc', 'Assembly',  'Organism_Name', 'Description', 'Mapped_Value']].to_csv(args.output, sep='\t', index=False)
        print(f"Output written to {args.output}")
    else:
        print("Output sample, no output file specified:")
        print(mapped_df[['Acc', 'Assembly', 'Organism_Name', 'Description', 'Mapped_Value']].head())
if __name__ == "__main__":
    main()
