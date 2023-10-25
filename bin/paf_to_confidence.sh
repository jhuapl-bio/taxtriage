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
Help message for taxid_to_reflist.sh:

DESCRIPTION:
    Returns a list of reference accessions with associated metadata (importantly, the ftp link for download), for the input tax ID (taxid, taxonomic identifier). Running this script will download the following file at the script's location (if it does not already exist there):
        ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt

NOTES:
    - All echo prints are from STDERR (will not be redirected with '>', can be with '2>')
    - alignments are filtered for:
        1) single best alignment per read
        2) reads with a [block length]:[read length] ratio > 0.80

OUTPUT FORMAT:
1    reference accession number
2    reference accession length (bp)
3    total reads aligned to accession sequence
4    abundance based on total aligned reads (does NOT include unaligned reads in denominator)
5    abundance based on total bases of aligned reads
6    abundance based on total bases of aligned reads (adjusted for genome size)
7    breadth of genome coverage (proportion of genome where position depth >0)
8    depth of coverage mean (uses total ref length, i.e. includes positions where depth=0)
9    depth of coverage standard deviation (stdev)
10    depth coefficient of variation (stdev/mean, smaller is better)

OPTIONS:
    -h    help    show this message
    -i    FILE    input paf filename (should include full path)
    -o    FILE    output tsv filename (should include full path)

USAGE:
bash paf_to_confidence.sh -i </full/path/to/alignment.paf> -o </full/path/to/output.tsv>
i="/data/projects/aphl_basestack/module-paf_to_confidence/examples/b_thuringensis_test1.paf"
o="/data/sandbox/paf_to_conf-b_thuringensis_test1.tsv"
bash paf_to_confidence.sh -i "$i" -o "$o"

EOF
}

#    ARGUMENTS
# parse args
while getopts "hi:o:" OPTION; do
    case $OPTION in
    h)
        usage
        exit 1
        ;;
    i) PAF=$OPTARG ;;
    o) OUTPUT=$OPTARG ;;
    ?)
        usage
        exit
        ;;
    esac
done
# check args
if [[ -z "$PAF" ]]; then
    printf "%s\n" "Please specify an input filename (-i) [include full path]."
    exit
fi
if [[ ! -f "$PAF" ]]; then
    printf "%s\n" "The input filename (-i) does not exist. Exiting."
    exit
fi
if [[ -z "$OUTPUT" ]]; then
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
mkdir -p "$tmp"
#    MAIN

# find single best alignment per read
echo >&2 "finding best single alignment per read"
#        based on MAPQ
awk -F'\t' '{
    if($1 in b){
        if($12>best[$1]){
            b[$1]=$12;
            best[$1]=$0;
        }
    }else{
        b[$1]=$12;
        best[$1]=$0;
    }
}END{
    for(r in best){
        print(best[r]);
    }
}' "$PAF" >"$tmp/tmp1.paf"

# include all alignments

#cat "$outdir/newtarget_filtered.paf" > "$outdir/seeker.paf"
bsa_seqcount=$(awk 'END{print(NR)}' "$tmp/tmp1.paf")
bsa_acccount=$(cut -f6 "$tmp/tmp1.paf" | sort | uniq | awk 'END{print(NR)}')
echo >&2 "  $bsa_seqcount single best alignments"
echo >&2 "  to $bsa_acccount unique accessions"

# get reads with a "block length"/"read length" greater than 80%
echo >&2 "finding reads with a [block length]:[read length] ratio > 0.80"
awk -F'\t' '{if(($10/$2)>0.0){print($0)}}' "$tmp/tmp1.paf" >"$tmp/tmp2.paf"

awk -F'\t' '{if(($10/$2)>0.0){count[$1]++; a[$1]=$6}}END{for(r in count){if(count[r]==1){print(a[r])}}}' "$tmp/tmp1.paf" | sort | uniq -c | sed -e 's/^ \+//' -e 's/ /\t/' >"$tmp/uniq.refs"

umseqs_count=$(awk -F'\t' '{x+=$1}END{print(x)}' "$tmp/uniq.refs")
umrefs_count=$(awk -F'\t' 'END{print(NR)}' "$tmp/uniq.refs")
echo >&2 "  $umseqs_count uniquely mapping reads"
echo >&2 "  to $umrefs_count accessions"
#>&2 echo "gathering alignment stats"
awk -F'\t' '{
    ilen[$6]=$7;
    ireads[$6]+=1;
    for(i=$8+1;i<=$9+1;i++){
        dep[$6][i]+=1
    };
}END{
    for(i in ireads){
        aligned+=ireads[i];
        for(j in dep[i]){
            sum[i]+=dep[i][j];
            sumsq[i]+=dep[i][j]^2;
            total_bases+=dep[i][j];
            cov[i]+=1;
        }
        tlen+=ilen[i];
    }


    for(i in dep){
        sfactor=(ilen[i]/tlen);
        adjsum[i]=(sfactor*sum[i]);
        adjtotal+=(sfactor*sum[i]);
    }

    # print output row per $6
    for(i in dep){
        o1=ireads[i];
        o3=o1/aligned;
        o4=(sum[i]/total_bases);
        o5=(adjsum[i]/adjtotal);
        o6=(cov[i]/ilen[i]);

        mean=sum[i]/ilen[i];
        stdev=sqrt((sumsq[i]-sum[i]^2/ilen[i])/ilen[i]);
        o9=stdev/mean;
        printf("%s\t%.0f\t%.0f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\n", i, ilen[i], o1, o3, o4, o5, o6, mean, stdev, o9);
    }
}' "$tmp/tmp2.paf" >$OUTPUT
# > $OUTPUT
# head $OUTPUT
#    i        ref accession
#    ilen[i]    ref accession length (bp) [assembly length for .assemblies output]
#    o1        total reads aligned to accession
#    o3        abundance based on total aligned reads (does NOT include unaligned reads in denominator)
#    o4        abundance based on total bases of aligned reads
#    o5        abundance based on total bases of aligned reads (adjusted for genome size)
#    o6        breadth of genome coverage (proportion of genome where position depth >0)
#    o7        depth of coverage mean (using total ref length, i.e. positions where depth=0)
#    o8        depth of coverage stdev
#    o9        stdev/mean depth ratio (coefficient of variation, smaller is better)

# clean up tmp dir
rm -r "$tmp"
