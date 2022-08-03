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
        description="Validate and transform a tabular samplesheet.",
        epilog="Example: python check_samplesheet.py samplesheet.csv samplesheet.valid.csv",
    )
    parser.add_argument(
        "-i",
        "--file_in",
        metavar="FILE_IN",
        type=Path,
        help="Tabular input Kraken2 report",
    )
    parser.add_argument(
        "-o",
        "--file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Name of the output tsv file containing a mapping of top n organisms at individual taxa levels",
    )
    parser.add_argument(
        "-t",
        "--top_per_rank",
        default=15,
        metavar="TOP_PER_RANK",
        type=int,
        help="Max Top per rank code",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def import_file(input):
    tsv_file = open(input)
    read_tsv = csv.reader(tsv_file, delimiter="\t")
    mapping=dict()
    header = ['abundance', 'clade_fragments_covered', 'number_fragments_assigned', 'rank', 'taxid','name']
    total = []
    for row in read_tsv:
        
        entry = dict()
        if row[3] not in mapping:
            mapping[row[3]] = []
        for x in range(0, len(header)):
            if (header[x] != 'name' and header[x] != 'rank' and header[x] != 'taxid'):
                entry[header[x]] = float(row[x])
            else:
                if header[x] == 'taxid':
                    entry[header[x]] = int(row[x])
                else:
                    entry[header[x]] = row[x]
                
        mapping[row[3]].append(entry)
    # k2_regex = re.compile(r"^\s{0,2}(\d{1,3}\.\d{1,2})\t(\d+)\t(\d+)\t([\dUDKRPCOFGS-]{1,3})\t(\d+)(\s+)(.+)")
    k2_regex = re.compile(r"^(\s+)(.+)")
    data = []
    depth = dict()
    for key, value in mapping.items():
        for  l in value:
            match = k2_regex.search(l['name'])
            if match:
                depth[l['taxid']] = int(len(match.group(1))/2)
                l['name'] = match.group(2)
            else:
                depth[l['taxid']] = 0
            l['depth'] = depth[l['taxid']]
            
    return mapping
def top_hit(mapping):
    uniq_ranks = list(mapping.keys())
    sorted_mapping = dict()
    i = 0
    for rank in uniq_ranks:
        sorted_specific_rank  = sorted(mapping[rank], key=lambda d: d['abundance'], reverse =True) 
        mapping[rank] = sorted_specific_rank
    header = ['abundance', 'clade_fragments_covered', 'number_fragments_assigned', 'rank', 'taxid','name']
    return mapping
def make_files(mapping,outdir,top):
    import json
    # try:
    #     os.mkdir(os.path.dirname(outdir))
    # except OSError as e:
    #     print(e)
    header = ['abundance', 'clade_fragments_covered', 'number_fragments_assigned', 'rank', 'taxid','name']
    path = str(outdir)+"_top_report.tsv"
    path  = open(path, "w")
    writer = csv.writer(path, delimiter='\t')
    writer.writerow(header )
    for rank_code, value in mapping.items():
        i = 0
        for row in value:
            out = []
            if i >= top:
                break
            else:
                for head_item in header:
                    out.append(row[head_item])
                writer.writerow(out)
            i = i + 1
    path.close()
    print("done")
    print("Done exporting files")
def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    print("yes")
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    print("yes")
    mapping = import_file(args.file_in)
    print("yes")
    mapping = top_hit(mapping)
    make_files(mapping, args.file_out, args.top_per_rank)
if __name__ == "__main__":
    sys.exit(main())
