#!/usr/bin/env python3

import re
import csv
import argparse


# Create the parser
parser = argparse.ArgumentParser(description='Process kraken report.')
parser.add_argument('--input', '-i', type=str, required=True,
                    help='Input kraken report file')
parser.add_argument('--output', '-o', type=str,
                    required=True, help='Output TSV file')
parser.add_argument('--id', '-d', type=str, default="N/A",
                    required=False, help='Sample id of the sample')
# Parse the arguments
args = parser.parse_args()


k2_regex = re.compile(
    r"^\s{0,2}(\d{1,3}\.\d{1,2})\t(\d+)\t(\d+)\t([\dUDKRPCOFGS-]{1,3})\t(\d+)(\s+)(.+)")

t_ranks = {
    # 'D': 'Domain',
    # 'K': 'Kingdom',
    # 'P': 'Phylum',
    # 'C': 'Class',
    # 'O': 'Order',
    # 'F': 'Family',
    # 'G': 'Genus',
    'S': 'Species',
    'S1': 'Species 1',
    'S2': 'Species 2',
    'S3': 'Species 3',
    'S4': 'Species 4',
    'S5': 'Species 5',
}
limited = {
    # 'D': 1,
    # 'K': 1,
    # 'P': 1,
    # 'C': 1,
    # 'O': 1,
    # 'F': 2,
    'G': 1,
    'S': 5
}

top_n_counts = {}
top_n_baseAll = 10  # adjust this as needed

with open(args.input, 'r') as file:
    data = []
    for line in file:
        match = k2_regex.search(line)
        if match:
            row = {
                "percent": float(match.group(1)),
                "counts_rooted": int(match.group(2)),
                "counts_direct": int(match.group(3)),
                "rank_code": match.group(4),
                "tax_id": int(match.group(5)),
                "num_spaces": len(match.group(6)),
                "classif": match.group(7).strip('\t'),
            }
            data.append(row)

    table_data = []
    i = 0
    keys = ['taxid_sample', 'rank', 'classification', 'percent', 'depth']
    for index, row in enumerate(data, start=1):
        if row["rank_code"] in limited:
            top_n_base = limited[row['rank_code']]
        else:
            top_n_base = top_n_baseAll

        if row["rank_code"] not in t_ranks:
            continue
        if row["rank_code"] not in top_n_counts:
            top_n_counts[row["rank_code"]] = 0
        if top_n_counts[row["rank_code"]] >= top_n_base:
            continue
        top_n_counts[row["rank_code"]] = top_n_counts.get(
            row["rank_code"], 0) + 1
        table_row = {
            "taxid_sample": str(row['tax_id'])+"_"+args.id,
            "rank": t_ranks.get(row["rank_code"], row["rank_code"]),
            "classification": row["classif"],
            "percent": row["percent"],
            "depth": row["num_spaces"],
        }
        i += 1
        table_data.append(table_row)

# Write table data to a tsv file
with open(args.output, 'w', newline='') as file:
    print("Writing to file: {}".format(args.output))
    fieldnames = keys
    writer = csv.DictWriter(file, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()
    for index, row in enumerate(table_data, start=1):
        writer.writerow(row)
