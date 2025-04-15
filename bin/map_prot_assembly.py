#!/usr/bin/env python


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

"""Provide a command line tool to fetch a list of refseq genome ids to a single file, useful for kraken2 database building or alignment purposes"""
from collections import defaultdict
import sys
import argparse
import re
import csv
import pysam
from math import log2



def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and transform a tabular samplesheet.",
        epilog="Example: python check_samplesheet.py samplesheet.csv samplesheet.valid.csv",
    )
    parser.add_argument(
        "-c",
        "--column",
        required=False,
        default="accession.version",
        metavar="FEATURE_COLUMN",
        help="Name of column to map feature. Default is accession.version.",
    )
    parser.add_argument(
        "-f",
        "--features",
        required=True,
        metavar="FEATURES",
        help="Features map file of accession.version of protein to taxid. Can be format from prot.accession2taxid.gz",
    )
    parser.add_argument(
        "-p",
        "--mapnames",
        required=False,
        metavar="MAPNAMES",
        help="Optional Name mapping of name of organism to a given taxid. Matches the Mapped_Vaue column to the taxid column in the output ",
    )
    parser.add_argument(
        "-n",
        "--names",
        required=False,
        metavar="NAMES",
        help="Optional Name mapping of a accession (2nd colun of --features input). Also, match patterns based on toxin, virulence, AMR using regex",
    )
    parser.add_argument(
        "-d",
        "--doutput",
        required=True,
        metavar="DIAMOND",
        help="DIAMOUND outfmt:6  File to process",
    )
    parser.add_argument(
        "-a",
        "--assembly",
        required=False,
        metavar="ASSEMBLY",
        help="Optional assembly file to retrieve assembly files from taxid column",
    )
    parser.add_argument(
        "-m",
        "--minidentity",
        required=False,
        default=60,
        metavar="MINIDENTITY",
        help="Filter on Minimum Identity for all queries to a contig/accessions/sequence - Default is 60",
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="OUTPUT",
        type=str,
        help="Name of output file that contains diamond outputs mapped to taxid/gcf",
    )


    return parser.parse_args(argv)
# make main function


def main(argv=None):
    """Run the main script."""
    args = parse_args(argv)
    name_mapping = dict()
    taxid_map = dict()
    if args.mapnames:
        with open(args.mapnames, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                taxid_map[row['Mapped_Value']] = row['Organism_Name']
        f.close()
    if args.names:
        header = ['Type', 'with', 'Assembly', 'Assembly Type', 'subtype',
          'empty', 'NUC Accession', 'length', 'spec', 'strand', 'Prot Accession',
          'Prot', 'Name', 'Description', 'empty2', 'empty3', 'empty4', 'empty5', 'empty6', 'empty7']
        # set regex patterns to only fine toxin, not anti-toxin, virulence, or Antimicrobial resistance or AMR
        # for toxin, make sure to ignore antitoxin if it is present, only if it is toxin, like ((?!anti)toxin
        patterns = [
            re.compile(r'(?i)toxin'),
            re.compile(r'(?i)virulence'),
            re.compile(r'(?i)antimicrobial resistance'),
            re.compile(r'(?i)amr')
        ]
        filterpattern = [
            # dont match antitoxin
            re.compile(r'(?i)antitoxin'),
        ]
        with open(args.names, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t', fieldnames=header)
            for row in reader:
                # if description in row matches any pattern, print it
                if any(pattern.search(row['Description']) for pattern in patterns):
                    # chekc if the filterpattern any doesnt match any
                    if not any(filter.search(row['Description']) for filter in filterpattern):
                        name_mapping[row['Prot Accession']] = row['Description']
                # name_mapping[row['Prot Accession']] = row['Description']
        f.close()
    # Load the features map
    features = {}
    with open(args.features, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            features[row[args.column]] = row['taxid']
    # Load the assembly map
    assembly = defaultdict(set)
    if args.assembly:
        # Open the file and read the data
        with open(args.assembly, mode='r') as file:
            reader = csv.reader(file, delimiter='\t')

            # Skip the first row
            next(reader)

            # Read the second row as the header and remove the starting '#'
            header = next(reader)
            header = [col.lstrip('#') for col in header]

            # Read the rest of the rows
            for row in reader:
                assembly[row[5]].add(row[0])
                assembly[row[6]].add(row[0])
    # Load the diamond output
    diamond = defaultdict(list)
    with open(args.doutput, 'r') as f:
        # Define the header manually since the file doesn't contain one
        header = ['Query ID', 'Subject ID', '% Identity', 'Alignment Length', 'Mismatches',
          'Gap Openings', 'Query Start', 'Query End', 'Subject Start',
          'Subject End', 'E-value', 'Bit Score']
        reader = csv.DictReader(f, delimiter='\t', fieldnames=header)
        for row in reader:
            diamond[row.get('Subject ID')].append({
                "Bit Score": row.get('Bit Score', 0),
                "E-value": row.get('E-value', 0),
                "Identity": row.get('% Identity', 0),
                'Mismatches': row.get('Mismatches', 0),
                'Gap Openings': row.get('Gap Openings', 0),
                'Alignment Length': row.get('Alignment Length', 0),
                "Query ID": row.get('Query ID', None),
            })

    # for each key, value, sort based on the Identity, and filter on Identity > args.minIdentity
    for key in diamond:
        if name_mapping.get(key):
            print('>', name_mapping.get(key), diamond[key])
        diamond[key] = [x for x in diamond[key] if float(x['Identity']) > args.minidentity ]
        diamond[key] = sorted(diamond[key], key=lambda x: x['Identity'], reverse=True)
        for entry in diamond[key]:
            entry['Description'] = name_mapping.get(key, "")
    # Map the diamond output to taxid
    taxid = defaultdict(list)
    query_lis = defaultdict(int)
    for key in diamond:
        if key in features:
            # Assign the result of features.get(key) to a different variable
            feature_taxid = features.get(key, None)
            if feature_taxid is not None:
                # Merge all elements of diamond[key] into taxid[key]
                for ele in diamond.get(key, []):
                    ele['Accession'] = key
                    taxid[feature_taxid].append(ele)  # Now taxid refers to the defaultdict



    # Write the output
    with open(args.output, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['taxid', "Name", '# contigs', '# CDS', 'medianIdentity', 'medianlength', 'medianmismatched', 'medianeval', 'desc'])
        for key, val in taxid.items():
            uniqueproteins = set([x['Accession'] for x in val])
            unique_contigs = set([x['Query ID'] for x in val])
            if key in taxid_map:
                name = taxid_map[key]
            else:
                name = "Unknown Name"
            desc = set([x['Description'] for x in val])
            # remove any ''
            desc = ", ".join(list(filter(None, desc)))
            # get the median of the Identity Score
            medianIdentity = 0
            medianlength = 0
            medianmismatched = 0
            medianeval = 0
            if len(val) > 0:
                medianIdentity = val[len(val)//2]['Identity']
                medianlength = val[len(val)//2]['Alignment Length']
                medianmismatched = val[len(val)//2]['Mismatches']
                medianeval = val[len(val)//2]['E-value']
            writer.writerow([key, name, len(unique_contigs),  len(uniqueproteins), medianIdentity, medianlength, medianmismatched, medianeval, desc])

if __name__ == "__main__":
    main()
