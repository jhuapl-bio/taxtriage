#!/usr/bin/env bash
##############################################################################################
# Copyright 2023 The Johns Hopkins University Applied Physics Laboratory LLC
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
#	FUNCTIONS
usage() {
    cat <<EOF
Help message for split_vcf.sh:

DESCRIPTION:
    Splits the VCF file on multiple reference accessions into individual by the kraken taxid. All references are split into a single file based on kraken:taxid column

OUTPUT FORMAT:
    Multiple .vcf files AND their respective .vcf.idx files
NOTES:
    Requires AWK 1.32 or higher
    Requires bcftools 1.12 or higher

OPTIONS:
    -h	help	show this message
    -i	FILE	Input VCF file
    -o	FILE	output VCF prefix. Files will be <prefix>.<taxid>.vcf and <prefix>.<taxid>.vcf.idx

USAGE:
bash split_vcf.sh -i </full/path/to/variants.vcf> -o <prefixString>
i="/data/inputs/variants.vcf"
o="/data/outputs/prefixString"
bash split_vcf.sh -i "$i" -o "$o"

EOF
}

#	ARGUMENTS
# parse args
while getopts "hi:o:" OPTION; do
    case $OPTION in
    h)
        usage
        exit 1
        ;;
    i) VCF=$OPTARG ;;
    o) OUTPUTPREFIX=$OPTARG ;;
    ?)
        usage
        exit
        ;;
    esac
done
# check args
if [[ -z "$VCF" ]]; then
    printf "%s\n" "Please specify an input filename (-i) [include full path to vcf]."
    exit
fi
if [[ ! -f "$VCF" ]]; then
    printf "%s\n" "The input filename (-i) does not exist. Exiting."
    exit
fi
if [[ -z "$OUTPUTPREFIX" ]]; then
    printf "%s\n" "Please specify an output prefix name (-o) [include full string name]."
    exit
fi

echo $VCF "into prefix of files > " $OUTPUTPREFIX

bcftools view $VCF | awk -v prefix="$OUTPUTPREFIX" -F'\t' '{
    if($0 ~ /^#/){
        print $0
    } else {
        split($1, parts, "|");
        taxid = parts[2];
        rowname = $0;

        if (!seen[taxid]++) {
            # Overwrite the file for the first encounter of the taxid
            print taxid "\t" rowname > prefix""taxid ".vcf";
        } else {
            # Append to the file for subsequent encounters of the taxid
            print taxid "\t" rowname >> prefix""taxid ".vcf";
        }
    }
} END {
    # for (k in hash) {
    #     print "Key:", k;
    #     print "Rows:";
    #     print hash[k];
    #     print "----------";
    # }

}'
