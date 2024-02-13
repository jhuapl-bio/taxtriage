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

Help message for fasta_reheader.sh:

DESCRIPTION:
    Processes a FASTA file to modify its headers. In each header:
    - Replaces the first space with a "|"
    - Replaces all subsequent spaces with "_"
    - Ignores everything after the first comma

NOTES:
    - This script is designed to work with FASTA files that have descriptive headers.

OPTIONS:
    -h          Show this help message
    -i FILE     Input FASTA file
    -o FILE     Output FASTA file where the processed headers will be written

USAGE:
    bash fasta_reheader.sh -i <input_fasta_file> -o <output_fasta_file>
    Example: bash fasta_reheader.sh -i input.fasta -o output.fasta

EOF
}

while getopts "hi:o:" OPTION
do
    case $OPTION in
        h) usage; exit 0 ;;
        i) INPUT_FASTA="$OPTARG" ;;
        o) OUTPUT_FASTA="$OPTARG" ;;
        ?) usage; exit 1 ;;
    esac
done

if [ -z "$INPUT_FASTA" ] || [ -z "$OUTPUT_FASTA" ]; then
    echo "Both input and output files must be specified."
    usage
    exit 1
fi

awk -F, 'BEGIN{OFS=""} /^>/{split($1, a, " "); gsub(" ", "_", $1); print a[1]"|"substr($1, length(a[1])+2)} !/^>/{print}' "$INPUT_FASTA" > "$OUTPUT_FASTA"

echo "Processed FASTA file saved to $OUTPUT_FASTA"
