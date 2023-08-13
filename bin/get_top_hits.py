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
        "-f",
        "--filter_ranks",
        metavar="TOP_HITS_STRING",
        type=str, nargs="+", default=['S2', 'S1', 'S', 'C', 'O', 'F', 'G', 'P', 'K', 'D', 'U' ],
        help="Filter on only showing specific ranks ",
    )
    parser.add_argument(
        "-s",
        "--top_hits_string",
        metavar="TOP_HITS_STRING",
        type=str,
        help="Top Hits ",
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
        default=5,
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


def import_file(input, filter_ranks):
    tsv_file = open(input,newline='' )
    read_tsv = csv.reader(tsv_file, delimiter="\t")
    mapping=[]
    
    taxids=dict()
    header = ['abundance', 'clade_fragments_covered', 'number_fragments_assigned', 'rank', 'taxid','name', 'parents']
    total = []
    for row in read_tsv:
        entry = dict()
        for x in range(0, len(header)):
            if (header[x] != 'name' and header[x] != 'rank' and header[x] != 'taxid'and x < len(row)) :
                entry[header[x]] = float(row[x])
            elif x < len(row):
                if header[x] == 'taxid':
                    entry[header[x]] = int(row[x])
                else:
                    entry[header[x]] = row[x]
            else:
                entry[header[x]]=""
        entry['name'] = entry['name'].strip()
        mapping.append(entry)
        taxids[entry['taxid']]=entry['name']
        
    # k2_regex = re.compile(r"^\s{0,2}(\d{1,3}\.\d{1,2})\t(\d+)\t(\d+)\t([\dUDKRPCOFGS-]{1,3})\t(\d+)(\s+)(.+)")
    k2_regex = re.compile(r"^(\s+)(.+)")
    depth = dict()
    lastparents = dict()
    for  l in mapping:
        match = k2_regex.search(l['name'])
        if match:
            depth[l['taxid']] = int(len(match.group(1))/2)
            l['name'] = match.group(2)
        else:
            depth[l['taxid']] = 0
        parents=[]
        
        for i in range (depth[l['taxid']]-1,0,-1):
            parents.append(int(lastparents[i]))
        lastdepth=depth[l['taxid']]
        lastparents[depth[l['taxid']]] = l['taxid']
        l['depth'] = depth[l['taxid']]
        l['parents'] =  parents
    if len(filter_ranks) > 0: 
        mapping = [m  for m in mapping if m['rank'] in filter_ranks]
        return mapping
    else: 
        return mapping
def top_hit(mapping, specific_limits, top_per_rank):
    uniq_ranks = list([x['rank'] for x in mapping])
    sorted_mapping = dict()
    for x in mapping:
        if not x['rank'] in sorted_mapping:
            sorted_mapping[x['rank']]=[]
        sorted_mapping[x['rank']].append(x)
    i = 0
    for rank in uniq_ranks:
        sorted_specific_rank  = sorted(sorted_mapping[rank], key=lambda d: d['abundance'], reverse =True) 
        sorted_mapping[rank] = sorted_specific_rank
    # header = ['abundance', 'clade_fragments_covered', 'number_fragments_assigned', 'rank', 'taxid','name','parents']
    newdata_seen = dict()
    countranks = dict()
    newdata  = dict()
    for specifics, value in specific_limits.items():
        rank = value['rank']
        limit = value['limit']
    
        if rank in sorted_mapping:
            i=0
            for row in sorted_mapping[rank]:
                if specifics in row['parents']:
                    if i >= limit:
                        break
                    else:
                        newdata_seen[row['taxid']] =  True
                        if not rank in countranks:
                            countranks[rank] = 1
                        else:
                            countranks[rank] += 1
                        i+=1
                        newdata[row['taxid']] = row
    for key, value in sorted_mapping.items():
        for row in value:
            if not key in countranks:
                countranks[key] = 0
            if countranks[key] >= top_per_rank:
                break
            else:
                if not row['taxid'] in newdata:
                    newdata[row['taxid']] = row
                    countranks[key] +=1
    for key, value in newdata.items():
        rank = value['rank']
    return newdata


def make_files(mapping, outpath):
    header = ['abundance', 'clade_fragments_covered', 'number_fragments_assigned', 'rank', 'taxid','name']
    path = str(outpath)
    path  = open(path, "w")
    writer = csv.writer(path, delimiter='\t')
    writer.writerow(header)
    for taxid, row in mapping.items():
        out = []
        for head_item in header:
            if head_item == 'parents':
                out.append(";".join([str(x) for x in row[head_item]]))
            else:
                out.append(row[head_item])
        writer.writerow(out)
    path.close()
    print("Done exporting the top hits report ")
def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    specific_limits = dict()
    if args.top_hits_string:
        for x in args.top_hits_string.split(";"):
            fulllist = x.split(":")
            if len(fulllist) >= 3:
                specific_limits[int(fulllist[0])] = dict( limit=int(fulllist[1]), rank=fulllist[2] )
    
    mapping = import_file(args.file_in, args.filter_ranks)
    mapping = top_hit(mapping, specific_limits, args.top_per_rank)
    make_files(mapping, args.file_out)
if __name__ == "__main__":
    sys.exit(main())
