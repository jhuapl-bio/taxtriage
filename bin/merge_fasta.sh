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
#    FUNCTIONS
usage() {
    cat <<EOF
Help message for merge_fasta.sh

DESCRIPTION:
    Returns a merged fasta file from multiple

NOTES:
    - All echo prints are from STDERR (will not be redirected with '>', can be with '2>')
    - alignments are filtered for:
        1) single best alignment per read
        2) reads with a [block length]:[read length] ratio > 0.80


OPTIONS:
    -h    help    show this message
    -i    FILE    input fasta filename (should include full path)
    -s    STRING    Text for samplename to add to each header
    -o    FILE    output fasta filename (should include full path)

USAGE:
bash merge_fasta.sh -i <fasta file(s)> -o <path_to_final_fasta_file>
i="example1.fasta example2.fasta"
o="merged_fasta.fasta"
bash merge_fasta.sh -i "$i" -o "$o"

EOF
}

#    ARGUMENTS
# parse args
while getopts "ho:i:s:" OPTION; do
    case $OPTION in
    h)
        usage
        exit 1
        ;;
    i) fasta=$OPTARG ;;
    s) samplename=$OPTARG ;;
    o) output=$OPTARG ;;
    ?)
        usage
        exit
        ;;
    esac
done
# check args
if [[ -z "$fasta" ]]; then
    printf "%s\n" "Please specify an input filename (-i) [include full path]."
    exit
fi
if [[ -z "$samplename" ]]; then
    printf "%s\n" "The input samplename (-s) does not exist. Exiting."
    exit
fi
if [[ -z "$output" ]]; then
    printf "%s\n" "Please specify an output filename (-o) [include full path]."
    exit
fi

# setup other variables
absolute_path_x="$(
    readlink -fn -- "$0"
    echo x
)"
absolute_path_of_script="${absolute_path_x%x}"
scriptdir=$(dirname "$absolute_path_of_script")
runtime=$(date +"%Y%m%d%H%M%S%N")
tmp="$scriptdir/tmp-$runtime"
#    MAIN

for line in ${fasta}; do
    gawk -v sample=$samplename '{
        if ($0 ~ /^>/){
            gsub(/>/,"",$1); print ">"sample"|"$_
        } else {
            print $0
        }
    }' $line
done >$output
