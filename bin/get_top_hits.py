#!/usr/bin/env python


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
def make_files(mapping,outdir):
    import json
    try:
        os.mkdir(outdir)
    except OSError as e:
        print(e)
    header = ['abundance', 'clade_fragments_covered', 'number_fragments_assigned', 'rank', 'taxid','name']
    
    for rank_code, value in mapping.items():
        path = os.path.join(outdir, rank_code+".tsv")
        path  = open(path, "w")
        writer = csv.writer(path, delimiter='\t')
        writer.writerow(header)
        for row in value:
            out = []
            for head_item in header:
                out.append(row[head_item])
            writer.writerow(out)
        path.close()
    print("Done exporting files")
def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    mapping = import_file(args.file_in)
    mapping = top_hit(mapping)
    make_files(mapping, args.file_out)
if __name__ == "__main__":
    sys.exit(main())
