#!/bin/bash

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


usage() {
    cat << EOF

Usage: $0 [options]

This script runs samtools coverage on a BAM file and replaces reference names
in the output according to a provided mapping file.

OPTIONS:
    -h          Show this help message
    -b FILE     Input BAM file
    -m FILE     Reference name mapping file (original and new names, tab-separated)
    -o FILE     Output file for the modified samtools coverage results
    -c ARGS     Which column to use for mapping. Default is 2
    -o FILE     Output file for the modified samtools coverage results

The mapping file should have two columns:
1. Original reference name (as in the BAM file)
2. New reference name

Each line in the mapping file corresponds to one reference name mapping.

Example:
bash $0 -b input.bam -m mapping.txt -o coverage_output.txt



EOF
}

# Default values for arguments
bam_file=""
mapping_file=""
output_file=""
args=""
mapcol=2
format="tsv"
# Parse command-line arguments
while getopts "hb:m:o:a:c:f:" opt; do
    case ${opt} in
        h )
            usage
            exit 0
            ;;
        b )
            bam_file=${OPTARG}
            ;;
        f )
            format=${OPTARG}
            ;;
        c )
            mapcol=${OPTARG}
            ;;
        m )
            mapping_file=${OPTARG}
            ;;
        o )
            output_file=${OPTARG}
            ;;
        a )
            args=${OPTARG}
            ;;
        \? )
            usage
            exit 1
            ;;
    esac
done

# Check for mandatory options
if [ -z "${bam_file}" ] || [ -z "${mapping_file}" ] || [ -z "${output_file}" ]; then
    echo "Error: Missing required arguments."
    usage
    exit 1
fi

# Check if samtools is installed
if ! command -v samtools &> /dev/null; then
    echo "Error: samtools could not be found. Please install samtools and try again."
    exit 1
fi

if [[ $format != "fasta" ]]; then
    # Run samtools coverage and process the output
    samtools coverage $args "${bam_file}"  | awk -F '\t' -v mapcol=$mapcol -v mapFile="${mapping_file}" '
        BEGIN {
            FS=OFS="\t"
            # Load the mapping from the mapping file
            while ((getline < mapFile) > 0) {
                map[$1] = $mapcol;
                # print $1,$2
            }
        }
        /^[^>]/ {
            # Replace the reference name if it exists in the mapping
            # split on space, get first index
            split($1, a, " ");

            if (a[1] in map) {
                $1 =a[1]" "map[a[1]];
            }
        } {print}
    ' > "${output_file}"
else
    # Run samtools coverage and process the output
    # get the headers from fasta file and

    echo "header coverage"

    samtools coverage $args "${bam_file}" | awk -F '\t' -v mapFile="${mapping_file}" '
        BEGIN {
            FS = OFS = "\t"
            while ((getline < mapFile) > 0) {
                if ($1 ~ /^>/) {
                    split($0, a, " ")
                    mapped = ""
                    # loop through all elemnts start at 2 to end for a and append to mapped

                    n = split($0, a, " ")
                    for (i = 2; i <= n; i++) {
                        mapped = mapped""a[i]" "
                    }
                    # remove > from a[1]
                    a[1] = substr(a[1], 2)
                    map[a[1]] = mapped
                }
            }
            close(mapFile)
        }
        /^[^>]/ {
            split($1, a, " ")
            if (a[1] in map) {
                # split on " " and get first index
                $1 = $1" "map[a[1]]
            }
        }
        { print }'  > "${output_file}"
fi
echo "Modified samtools coverage output saved to ${output_file}"
