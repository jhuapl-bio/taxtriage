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
usage() {
    cat <<EOF
Help message for remove_taxids.sh:

DESCRIPTION:
    Removes 1 or more taxids from kraken2 report output and adjusts relative abundances

NOTES:
    - All echo prints are from STDERR (will not be redirected with '>', can be with '2>')
    - alignments are filtered for:
        1) single best alignment per read
        2) reads with a [block length]:[read length] ratio > 0.80
OUTPUT FORMAT:
    Kraken2 Report File

OPTIONS:
    -h    help    show this message
    -i    FILE    input sam filename (should include full path)
    -o    FILE    output tsv filename (should include full path)
    -t    STRING    1 or more taxids, separated by spaces and in double quotes

USAGE:
remove_taxids.sh -i </full/path/to/kraken2/report> -o </full/path/to/kraken2/filtered/report> -t "<taxids_separated_by_space>"
i="BC01.kraken.report"
o="BC01.filtered.kraken.report"
t="1 2 9609"
remove_taxids.sh -i "$i" -o "$o" -t "$t"

EOF
}

OUTPUT_READS="false"

#    ARGUMENTS
# parse args
while getopts "hi:o:t:" OPTION; do
    case $OPTION in
    h)
        usage
        exit 1
        ;;
    i) REPORT=$OPTARG ;;
    t) TAXIDS=$OPTARG ;;
    o) OUTPUT=$OPTARG ;;
    ?)
        usage
        exit
        ;;
    esac
done
# check args
if [[ -z "$REPORT" ]]; then
    printf "%s\n" "Please specify an input kraken2 report filename (-i) [include full path]."
    exit
fi
if [[ ! -f "$REPORT" ]]; then
    printf "%s\n" "The input filename (-i) does not exist. Exiting."
    exit
fi
if [[ -z "$TAXIDS" ]]; then
    printf "%s\n" "Please specify taxids count of 1 or more (-t) [include full path]."
    exit
fi
if [[ -z "$OUTPUT" ]]; then
    printf "%s\n" "No output"
    exit
fi
#    MAIN

# find single best alignment per read
echo >&2 "filtering taxids"
#        based on MAPQ

gawk -v filtertaxids="${TAXIDS[@]}" -F '\t' '
    function join(array, start, end, sep,    result, i)
    {
        if (sep == "")
        sep = " "
        else if (sep == SUBSEP) # magic value
        sep = ""
        result = array[start]
        for (i = start + 1; i <= end; i++)
            result = result sep array[i]
        return result
    }
    BEGIN {
        OFS="\t";
        split(filtertaxids,arr," ")
        for (f in arr){
            seen[arr[f]]=1
        }
        start=0
        idx=1+0
        globalidx=0+0
        markedpos=0+0
        lastdepth=0+0
        subtract = 0
    }
    {
        taxidDepth=$6
        n=gsub(/\s{2}|\t/, "_", taxidDepth )
        if (seen[$5]==1){
            start=1
            subtract=$3
            first=1
            markedpos=n
        } else if (start == 1 && markedpos < n){
            start=1
            first=0
            dontremove[$5]=1
            subtract=$3
        } else {
            start=0
            n=n+0
            first=0
            latest[n]=$5
            saved[globalidx]=n
            parents[n][$5]=$2
            taxids[globalidx]=$5
            fulline[$5]=$0
            idx+=1
            total_remain_raw+=$3
        }

        globalidx+=1
        if (start == 1 && markedpos <= n  ){
            for (depth in parents){
                if (depth < markedpos){
                    latesttaxidatdepth=latest[depth]
                    parents[depth][latesttaxidatdepth]-=subtract

                }
            }
        }
        total_raw+=$3
        lastdepth=n
        lasttaxid=$5

    } END{
        divided=total_raw
        PROCINFO["sorted_in"] = "@ind_num_asc"
        for (position in taxids){
            taxid=taxids[position]
            # print "taxid: ", taxid, "line number: ", position, "Value new: ", parents[saved[position]][taxid]
            full=fulline[taxid]
            split(full, linee, "\t")
            linee[2]=parents[saved[position]][taxid]
            linee[1]=sprintf("%.2f",linee[2]*100/divided)
            g=join(linee,1,length(linee),"\t")
            print g
        }
        # print "@TOTAL_RAW:"total_raw"\t@TOTAL_REMAIN_RAW:"total_remain_raw

    }

    ' $REPORT 2>&1 | tee $OUTPUT
#  |     awk  -F '\t' '
#     BEGIN {
#         OFS="\t";
#     }
#     {

#     }
#     '
