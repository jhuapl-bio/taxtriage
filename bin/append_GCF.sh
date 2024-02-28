
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
    -i    FILE    a fasta file/dir of gcf references pulled. Appends the filename's GCF to the middle of the accession (between chr accession and description)
    -d    FILE    optional. If called, it is a directory of GCF files (*fna, *fasta, *faa, *fa)
    -o    FILE        full path to output  file

USAGE:

bash append_GCF.sh -i GCF_020735385.1.fasta -o test.fasta -t (looks at all fasta files in dir)
bash append_GCF.sh -i GCF_020735385.1.fasta -o test.fasta  (appends base assembly accession to center of contig/chromosome header)

EOF
}

DIR="FALSE"

#    ARGUMENTS
# parse args
while getopts "hi:do:" OPTION
do
    case $OPTION in
        h) usage; exit 1 ;;
        i) fasta=$OPTARG ;;
        d) DIR="TRUE" ;;
        o) output=$OPTARG ;;
        ?) usage; exit ;;
    esac
done


# make a function to append the GCF to the middle of the accession
# if the file is GCF_102310239.1_aSDGM I want only GCF. Split on the _ and get first 2 elements
# now get all lines that start with ">" in the filename, and append the gcf to the middle of the accession between NC_12123 and "djaskdjasdklj" (random description). It should be $1 then $gcf then everything AFTER $1
# use sed after first space
append_f(){
    filename=$1
    # echo "Name BEGINS with GCF: ${basen}"
    # get the basename of the file removing extension
    basen=$(basename $filename | sed 's/\.[^.]*$//')
    # if the file is GCF_102310239.1_aSDGM I want only GCF. Split on the _ and get first 2 elements
    gcf=$(echo $basen | awk -F '_' '{print $1"_"$2}')
    # now get all lines that start with ">" in the filename, and append the gcf to the middle of the accession between NC_12123 and "djaskdjasdklj" (random description). It should be $1 then $gcf then everything AFTER $1
    # use sed after first space
    sed -e "s/ / $gcf /" $filename
}


if [[ $DIR == "TRUE" ]]; then
    echo "Directory detected"
    for filename in $(find $fasta  \( -name "*.fasta" -o -name "*fna" -o -name "*faa" -o -name "*fa" \) -type f ); do
        basen=$(basename $filename)
        if [[ $basen =~ "GCF" ]]; then
            append_f $filename
        fi

    done

else
    append_f $fasta
fi

