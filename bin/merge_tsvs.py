#!/usr/bin/env python3

import sys
import argparse
import os
parser = argparse.ArgumentParser(
    description='Process and merge multiple tsv files')
parser.add_argument('--input', '-i', type=str,  nargs="+",
                    required=True, help='Input  tsv files')
# add argument boolean for if header is true or not
parser.add_argument('--header', '-d',  action='store_true', default=True,
                    help='If header is present, extract first row first file, ignore for rest')
parser.add_argument('--output', '-o', type=str, required=True,
                    help='Output TSV file that is merged')
parser.add_argument('--append_name', '-a', action='store_true', default=False,
                    help='Append the basename of the file as the first column')
# Parse the arguments
args = parser.parse_args()


def merge_tsvs(file_paths):
    with open(args.output, 'w') as outfile:
        idx = 0
        for i, file_path in enumerate(file_paths):
            with open(file_path, 'r') as infile:
                if args.header:
                    if i == 0:
                        # write header for the first file
                        if args.append_name:
                            outfile.write("Sample\t"+infile.readline())
                        else:
                            outfile.write(infile.readline())
                    else:
                        # skip header for subsequent files
                        infile.readline()
                else:
                    infile.readline()
                # copy the rest of the file content
                p = 0
                for line in infile:
                    # cange the first index (on tab separated) to be index +1
                    linesplit = line.split("\t")
                    linesplit[0] = str(idx)
                    line = "\t".join(linesplit)
                    if args.append_name:
                        parsed = os.path.basename(file_path.split(".")[0])
                        line = parsed+"_"+str(p) + "\t" + line
                        p += 1
                    idx +=1
                    outfile.write(line)


if __name__ == "__main__":
    merge_tsvs(args.input)
