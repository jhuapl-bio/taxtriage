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

"""Provide a command line tool to validate and transform tabular samplesheets."""


import argparse
import csv
import logging
import sys
from collections import Counter
from pathlib import Path
import re
import os
logger = logging.getLogger()


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Species merge/count of all subspecies and strains",
        epilog="Example: python3 merge_subspecies.py -i input.krakenreport.txt -o output.csv",
    )
    parser.add_argument(
        "-i",
        "--file_in",
        metavar="FILE_IN",
        type=Path, nargs="+", required=True,
        help="Tabular input Kraken2 report",
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="OUTPUT",
        type=str, default=['S2', 'S1', 'S', 'C', 'O', 'F', 'G', 'P', 'K', 'D', 'U'],
        help="Filter on only showing specific ranks ",
    )
    parser.add_argument(
        "-e",
        "--level",
        metavar="OUTPUT",
        type=str, default="S", choices=['S2', 'S1', 'S', 'C', 'O', 'F', 'G', 'P', 'K', 'D', 'U'],
        help="Filter merge on this rank. Default is Species (S)",
    )
    parser.add_argument(
        "-p",
        "--matchpath",
        metavar="MATCH",
        type=str, default=None, required=False,
        help="Sheet to match a taxid set or list to taxids in the kraken2 report. Annotates if there is a match",
    )
    parser.add_argument(
        "-c",
        "--colmatch",
        metavar="COL",
        type=int, default=2,
        help="Column to match the matchpath to for kraken2 taxid match, Default is 2",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)
from collections import defaultdict

def count_double_spaces(line):
    double_space_count = 0
    index = 0
    while index < len(line) - 1:
        if line[index] == ' ' and line[index + 1] == ' ':
            double_space_count += 1
            index += 2
        else:
            break
    return double_space_count

def process_kraken2_report(report_file, target_rank):
    with open(report_file, 'r') as file:
        lines = file.readlines()

    name_children_count = []
    current_parent_name = None
    current_parent_indent_level = -1

    for line in lines:
        columns = line.split('\t')

        if len(columns) < 6:
            continue

        current_rank = columns[3].strip()
        current_taxid = columns[4].strip()
        current_name = columns[5].strip()
        current_indent_level = count_double_spaces(columns[5])

        if current_rank == target_rank:
            # Start tracking a new "target_rank" and its children
            current_parent_name = current_name
            current_parent_indent_level = current_indent_level
            name_children_count.append(dict(
                taxid=current_taxid,
                name=current_name,
                match=False,
                ranks=defaultdict(list)
            ))
        elif current_parent_name and current_indent_level > current_parent_indent_level:
            # We are within the children of the current "target_rank" rank
            matchdict = dict(
                taxid=current_taxid,
                name=current_name,
                match=False
            )
            name_children_count[-1]['ranks'][current_rank].append(matchdict)
        elif current_indent_level <= current_parent_indent_level:
            # We moved out of the current "target_rank" rank's children
            current_parent_name = None
            current_parent_indent_level = -1

    return name_children_count

def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    # read in the kraken report
    report_files = args.file_in
    target_rank = args.level
    colmatch = args.colmatch
    matchpath = args.matchpath

    all_name_children_counts = []
    if len(report_files) == 0:
        print("Error: No report files provided, exiting...")
        return 1
    else:
        for report_file in report_files:
            name_children_count = process_kraken2_report(report_file, target_rank)

            if matchpath:
                # read in the matchpath
                if not os.path.exists(matchpath):
                    print(f"Error: File {matchpath} does not exist")
                    continue

                with open(matchpath, 'r', encoding='utf-8', errors='replace') as file:
                    lines = file.readlines()
                    for line in lines:
                        columns = line.split(',')
                        if len(columns) < colmatch:
                            continue
                        taxid = columns[colmatch-1].strip()
                        for entry in name_children_count:
                            if entry['taxid'] == taxid:
                                entry['match'] = True
                            for rank, organisms in entry['ranks'].items():
                                for organism in organisms:
                                    if organism['taxid'] == taxid:
                                        organism['match'] = True

            # Get the base name of the file
            sample_name = os.path.basename(report_file).replace('.filtered.report', '')
            for entry in name_children_count:
                entry['sample'] = sample_name

            all_name_children_counts.extend(name_children_count)

        generate_output_csv(all_name_children_counts, args.output, target_rank)

def generate_output_csv(name_children_count, combined_output_file, target_rank):
    output_data = []

    # Collect data for combined output and individual samples
    sample_groups = defaultdict(list)
    for entry in name_children_count:
        rank_counts = defaultdict(int)
        matchthis = 1 if entry['match'] else 0
        rank_details = []

        for rank, organisms in entry['ranks'].items():
            rank_counts[rank] = len(organisms)
            rank_details.append(f"{rank}({len([org for org in organisms if org['match']])})")
        record = {
            'Sample': entry['sample'],
            'Taxid': entry['taxid'],
            'Name': entry['name'],
            'Match': matchthis + sum(1 for org in entry['ranks'].values() for item in org if item['match']),
            'Lower Ranks': ', '.join(rank_details)
        }
        output_data.append(record)
        sample_groups[entry['sample']].append(record)

    # Write combined CSV output
    with open(combined_output_file, 'w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=['Sample', 'Taxid', 'Name', 'Match', 'Lower Ranks'])
        writer.writeheader()
        writer.writerows(output_data)

    # Write individual CSV files for each sample
    dirname = os.path.dirname(combined_output_file)
    for sample, records in sample_groups.items():
        individual_output_file = os.path.join(dirname, f"{sample}.single.csv")
        with open(individual_output_file, 'w', newline='') as file:
            writer = csv.DictWriter(file, fieldnames=["Sample", 'Taxid', 'Name', 'Match', 'Lower Ranks'])
            writer.writeheader()
            writer.writerows(records)

if __name__ == "__main__":
    sys.exit(main())

