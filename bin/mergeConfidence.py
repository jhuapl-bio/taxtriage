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
        description="Merge a paf to confidence output mqc.tsv with the kraken report it corresponds to include tax info",
        epilog="Example: ",
    )
    parser.add_argument(
        "-i",
        "--file_in",
        metavar="FILE_IN",
        type=Path,
        help="Tabular mqc.tsv output of paf_to_confidence.sh",
    )
    parser.add_argument(
        "-s",
        "--sample",
        metavar="SAMPLE",
        type=str,
        default=None,
        help="OPTIONAL, samplename analyzed",
    )
    parser.add_argument(
        "-k",
        "--kraken_report",
        metavar="FILE_OUT",
        type=Path,
        help="Kraken report file to merge information with",
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


def import_file(tsv_filename, kraken_report, samplename):
    tsv_file = open(tsv_filename)
    rows = dict()
    confidence_tsv = csv.reader(tsv_file, delimiter="\t")
    for row in confidence_tsv:
        rows[row[0]] = row
        if samplename:
            row[0] = row[0] +"_" + samplename
    tsv_file.close()
    report_file = open(kraken_report)
    
    report = csv.reader(report_file, delimiter="\t")
    for row in report:
        taxid = row[4]
        rank=row[3]  
        
        name = row[5].strip()
        if taxid in rows:
            rows[taxid].insert(1,rank)      
            rows[taxid].insert(1,name)      
    report_file.close()
    
    return rows

def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    rows = import_file(args.file_in, args.kraken_report, args.sample)
    path  = open(args.file_out, "w")
    writer = csv.writer(path, delimiter='\t')
    writer.writerow(['Acc', 'Name', 'Rank', 'Length', 'Reads Aligned', 'Abu Total Aligned', 'Abu Total Bases', 'Total Bases Adjusted', 'Breadth Genome Coverage', 'Depth Coverage Mean', 'Depth Coverage Stdev', 'Depth Coef. Variation'])
    if len(rows.keys()) == 0:
        empty = ['None']
        for x in range(0, 11):
            empty.append("N/A")
        writer.writerow(empty)

    else:
        for key, value in rows.items():
            writer.writerow(value)
    path.close()    
if __name__ == "__main__":
    sys.exit(main())
