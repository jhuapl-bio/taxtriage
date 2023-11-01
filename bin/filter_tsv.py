#!/usr/bin/env python3

import sys
import argparse
import os
import csv

parser = argparse.ArgumentParser(
    description='Filter merged tsv file by taxonomic rank')
parser.add_argument('--input', '-i', type=str,  nargs='+',
                    required=True, help='Input merged tsv file')
parser.add_argument('--output', '-o', type=str,
                    required=True, help='Output TSV files')

# Parse the arguments
args = parser.parse_args()


def filter_tsv(filepath):
    # Define the classes and ranks
    omit = {'U', 'D', 'S1', 'S2', 'K'}
    rank_to_dict = {'P': 'phylum', 'C': 'class', 'O': 'order',
                    'F': 'family', 'G': 'genus', 'S': 'species'}

    # Read the TSV file line by line and create filtered files
    with open(filepath, 'r') as input_file:
        reader = csv.DictReader(input_file, delimiter='\t')

        for row in reader:
            rank = row['rank']

            if rank not in omit:
                class_name = rank_to_dict.get(rank)
                if class_name:
                    output_file_path = f"{class_name}_filtered_mqc.tsv"

                    # Create a new dictionary without the 'rank' column
                    filtered_row = {key: value for key,
                                    value in row.items() if key != 'rank'}

                    with open(output_file_path, 'a') as output_file:
                        writer = csv.DictWriter(
                            output_file, fieldnames=filtered_row.keys(), delimiter='\t')

                        if output_file.tell() == 0:
                            writer.writeheader()

                        writer.writerow(filtered_row)


if __name__ == "__main__":
    print(args.input)
    filter_tsv(args.input[0])
