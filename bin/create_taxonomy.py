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
import re
import argparse

def parse_kraken2_report(report_file, rank_mapping):
    tax_info = {}
    k2_regex = re.compile(r"^(\s+)(.+)")
    depth = dict()
    lastparents = dict()

    with open(report_file, 'r') as infile:
        for line in infile:
            fields = line.strip().split('\t')
            tax_id = fields[4].strip()  # Taxonomic ID
            name = fields[5]    # Scientific name
            rank_code = fields[3].strip()    # Taxonomic rank code
            rank = get_rank_name(rank_code, rank_mapping)  # Get the rank name
            # Extract hierarchy information
            match = k2_regex.search(name)

            if match:
                current_depth = int(len(match.group(1)) / 2)
                name = match.group(2)
            else:
                current_depth = 0

            # Determine parent ID
            if current_depth > 0:
                parent_id = lastparents.get(current_depth - 1, "0")
            else:
                parent_id = "0"
            # Update depth and parent tracking
            depth[tax_id] = current_depth
            lastparents[current_depth] = tax_id

            # Store taxonomic information
            tax_info[tax_id] = {
                "name": name,
                "rank": rank,
                "parent": parent_id,
                "depth": current_depth
            }

    return tax_info

def get_rank_name(rank_code, rank_mapping):
    # Check if the rank_code matches exactly
    if rank_code in rank_mapping:
        return rank_mapping[rank_code]
    # If the rank code has a trailing number, extract the prefix and the number
    else:
        return None
        # base_rank = rank_code[0]
        # suffix = rank_code[1:] if len(rank_code) > 1 else ""
        # if base_rank in rank_mapping:
        #     return f"{rank_mapping[base_rank]}{suffix}"
        # else:
        #     return "no rank"

def write_names_dmp(tax_info, output_file):
    with open(output_file, 'w') as out:
        for tax_id, info in tax_info.items():
            # Write to names.dmp format
            out.write(f"{tax_id}\t|\t{info['name']}\t|\t\t|\tscientific name\t|\n")

def write_nodes_dmp(tax_info, output_file):
    with open(output_file, 'w') as out:
        for tax_id, info in tax_info.items():
            # Write to nodes.dmp format with static fields
            if info.get('rank', None):
                out.write(f"{tax_id}\t|\t{info['parent']}\t|\t{info['rank']}\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\n")

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description='Generate names.dmp and nodes.dmp from a Kraken2 report.')
    parser.add_argument('-k', '--report', required=True, help='Input Kraken2 report file')
    parser.add_argument('--names', default='names.dmp', help='Output names.dmp file (default: names.dmp)')
    parser.add_argument('--nodes', default='nodes.dmp', help='Output nodes.dmp file (default: nodes.dmp)')
    args = parser.parse_args()

    # Rank mapping based on your specific nodes.dmp file codes
    rank_mapping = {
        'S': 'species',
        'S1': 'subspecies',
        'G': 'genus',
        'G1': 'subgenus',
        'F': 'family',
        'F1': 'subfamily',
        'O': 'order',
        'O1': 'suborder',
        'C': 'class',
        'C1': 'subclass',
        'P': 'phylum',
        'P1': 'subphylum',
        'K': 'kingdom',
        'U': 'no rank',
        'R': 'no rank'
    }

    # Parse the Kraken2 report
    tax_info = parse_kraken2_report(args.report, rank_mapping)

    # Write names.dmp and nodes.dmp
    write_names_dmp(tax_info, args.names)
    write_nodes_dmp(tax_info, args.nodes)

if __name__ == "__main__":
    main()
