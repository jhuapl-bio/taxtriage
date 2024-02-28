#!/usr/bin/env bash
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
usage()
{
cat << EOF

Help message for make_gcf_mapping.sh:

DESCRIPTION:
    Optionally merge assembly refseq on a gcf assembly pull FASTA file. FASTA file must be in format >chr GCF description (3 columns) where first 2 are a space and the third is anything after

NOTES:


OPTIONS:
    -h    help    show this message
    -i    FILE    a fasta file of gcf references pulled. should be in chr \s gcf \s description format
    -a    FILE    a file of the ncbi assembly refseq. Optional. If optional, the 3rd column of the mapping file will be pasted to col 4
    -o    DIR        full path to output  file

USAGE:

bash make_gcf_mapping.sh -i ~/Downloads/combined_genomes.fasta -o test_output/test.tsv -a test_output/get/assembly_summary_refseq.txt

EOF
}




#    ARGUMENTS
# parse args
while getopts "hi:a:o:" OPTION
do
    case $OPTION in
        h) usage; exit 1 ;;
        i) fasta=$OPTARG ;;
        a) assembly=$OPTARG ;;
        o) output=$OPTARG ;;
        ?) usage; exit ;;
    esac
done


# check args
if [[ -z $fasta ]] || [[ -z $output ]]
then
    usage
    exit 1
fi
# print something if assembly is empty
if [[ -z $assembly ]]
then
    echo "No assembly file provided. Will only output gcf \t chr \t description \t description (repeat)"
fi


grep '^>' $fasta   | awk -F ' ' '{
    # Initialize the description variable
    description = "";

    # Concatenate fields to form the description, ignoring fields after the first comma
    for (i = 2; i <= NF; i++) {
        if (index($i, ",") != 0) {
            # If theres a comma in the current field, split it and take the part before the comma
            split($i, parts, ",");
            description = description (description == "" ? "" : " ") parts[1];
            break;  # Stop processing further fields after encountering the first comma
        } else {
            description = description (description == "" ? "" : " ") $i;
        }
    }
    # Remove the ">" from the first column
    gsub("^>", "", $1);
    print $1 "\t" description;

}' > $output.tmp

echo "Assembly references of names completed from FASTA file. Moving on to adding the organism name"

if [[  -z $assembly ]]
then
    paste $output.tmp <(awk -F '\t' '{print $3}' $output.tmp) > $output
else
    # use awk, import both output.tmp and assembly, and print the 1st, 2nd, 3rd, and map
    awk -F '\t' 'BEGIN {FS=OFS="\t"}
    NR==FNR {map[$1]=$8; next}
    {
        if ($2 == "") {
            print $0, $3  # If assembly is empty, duplicate the third column
        } else {
            print $0, (map[$2] != "" ? map[$2] : $3)  # Otherwise, use mapped value if present, else default to third column
        }
    }' $assembly $output.tmp > $output
fi





