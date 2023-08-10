#!/usr/bin/env python3

import sys
import argparse 
parser = argparse.ArgumentParser(description='Process and merge multiple tsv files')
parser.add_argument('--input', '-i', type=str,  nargs="+", required=True, help='Input  tsv files')
# add argument boolean for if header is true or not 
parser.add_argument('--header', '-d', type=bool, default=True, required=False, help='If header is present, extract first row first file, ignore for rest')
parser.add_argument('--output', '-o', type=str, required=True, help='Output TSV file that is merged')
# Parse the arguments
args = parser.parse_args()

def merge_tsvs(file_paths):
    with open(args.output, 'w') as outfile:
        for i, file_path in enumerate(file_paths):
            with open(file_path, 'r') as infile:
                if args.header:
                    if i == 0:
                        # write header for the first file
                        outfile.write(infile.readline())
                    else:
                        # skip header for subsequent files
                        infile.readline()
                else:
                    infile.readline()
                # copy the rest of the file content
                for line in infile:
                    outfile.write(line)

if __name__ == "__main__":
    merge_tsvs(args.input)