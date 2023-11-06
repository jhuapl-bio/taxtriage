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
    -i    FILE    input sam filename (should include full path)
    -o    FILE    output tsv filename (should include full path)

USAGE:
bash sam_to_confidence.sh -i </full/path/to/alignment.sam> -o </full/path/to/output.tsv>
i="/data/projects/aphl_basestack/module-sam_to_confidence/examples/b_thuringensis_test1.sam"
o="/data/sandbox/sam_to_conf-b_thuringensis_test1.tsv"
bash sam_to_confidence.sh -i "$i" -o "$o"

EOF
}




OUTPUT_READS="false"

#    ARGUMENTS
# parse args
while getopts "hi:o:r:m:" OPTION
do
    case $OPTION in
        h) usage; exit 1 ;;
        i) SAM=$OPTARG ;;
        m) PILEUP=$OPTARG ;;
        o) OUTPUT=$OPTARG ;;
        r) OUTPUT_READS=$OPTARG;;
        ?) usage; exit ;;
    esac
done
# check args
if [[ -z "$SAM" ]]; then printf "%s\n" "Please specify an input filename (-i) [include full path]."; exit; fi
if [[ ! -f "$SAM" ]]; then printf "%s\n" "The input filename (-i) does not exist. Exiting."; exit; fi
if [[ -z "$PILEUP" ]]; then printf "%s\n" "Please specify an mpileup filename (-m) [include full path]."; exit; fi
if [[ -z "$OUTPUT" ]]; then printf "%s\n" "Please specify an output filename (-o) [include full path]."; exit; fi
if [[  ! "$OUTPUT_READS" == "false" ]]; then printf "%s\n" "READS to be output in a 2 column file to $OUTPUT_READS"; fi
# setup other variables
absolute_path_x="$(readlink -fn -- "$0"; echo x)"
absolute_path_of_script="${absolute_path_x%x}"
scriptdir=$(dirname "$absolute_path_of_script")
runtime=$(date +"%Y%m%d%H%M%S%N")
tmp="$scriptdir/tmp-$runtime-$RANDOM"
mkdir -p "$tmp"
#    MAIN

# find single best alignment per read
>&2 echo "finding best single alignment per read"
#        based on MAPQ
regex="NM:i:([0-9]+)"


gawk -v regex=$regex -F'\t' 'BEGIN{OFS="\t"}{
    if ($1 ~ /^@/){
        print $0
    }
    match($0, /LN:([0-9]+)/, ary )
    if (ary[1]){
        match($0, /SN:([^\t]+)/, ars )
        t[ars[1]]=ary[1]
    }
    match($0, /NM:i:([0-9])+/, arx )
    if (length(arx) >0){
        mm=arx[1]

        match($6, /([0-9]+)M/, s )

        size=s[1]
        if($1 in b){
            if($5>best[$1]){
                b[$1]=$5;
                best[$1]=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"size"\t"$10"\t"t[$3]"\t"mm;
            }
        }else{
            b[$1]=$5;
                best[$1]=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"size"\t"$10"\t"t[$3]"\t"mm;
        }
    }

}END{
    for(r in best){
        print(best[r]);
    }
}' "$SAM"  > "$tmp/tmp1_included.sam"

# include all alignments
bsa_seqcount=$(awk 'END{print(NR)}' "$tmp/tmp1_included.sam")
bsa_acccount=$(cut -f3 "$tmp/tmp1_included.sam" | sort | uniq | awk 'END{print(NR)}')
>&2 echo "  $bsa_seqcount single best alignments"
>&2 echo "  to $bsa_acccount unique accessions"

# get reads with a "block length"/"read length" greater than 80%
>&2 echo "finding reads with a [block length]:[read length] ratio > 0.80"

awk -F'\t' 'function abs(x){return +sqrt(x*x)}{if($1 ~ /^@/ || abs(($9-$13)/$9)>=0.01){print($0)}}' "$tmp/tmp1_included.sam"   > "$tmp/tmp2.sam"
awk -F'\t' 'function abs(x){return +sqrt(x*x)}{if( $1 ~ /^[^@]/ && abs(($9-$13)/$9)>0.01){count[$1]++; a[$1]=$3}}END{for(r in count){if(count[r]==1){print(a[r])}}}' "$tmp/tmp1_included.sam" \
    | sort \
    | uniq -c \
    | sed -e 's/^ \+//' -e 's/ /\t/' > "$tmp/uniq.refs"

if [[ ! $OUTPUT_READS == 'false' ]]; then
    echo "adding reads to $OUTPUT_READS as 2 column file"
    cut -f 1,3 -d $'\t' "$tmp/tmp2.sam" | grep -v '^@' > $OUTPUT_READS
    # awk -v output=$OUTPUT_READS -F "\t" '{ print>output"_"$3".reads" }' "$tmp/tmp2.sam"
    echo "added reads to necessary file"
fi


umseqs_count=$(awk -F'\t' '{x+=$3}END{print(x)}' "$tmp/uniq.refs")
umrefs_count=$(awk -F'\t' 'END{print(NR)}' "$tmp/uniq.refs")
>&2 echo "  $umseqs_count uniquely mapping reads"
>&2 echo "  to $umrefs_count accessions"
>&2 echo "gathering alignment stats"
echo "______"

awk 'BEGIN {
    printf "%s\t%s\t%s\t%s\t", "Accession", "Covered", "Mean Depth", "Median Depth\n"
}
{
    if ($4 > 0) {
        covered[$1]++;
        sum_depth[$1] += $4;
        depths[$1] = depths[$1] " " $4;  # store depths in an array for median computation
    }
    total[$1]++;
}
END {
    for (genome in total) {
        # compute mean
        mean = sum_depth[genome]/total[genome];
        printf "%s\tCovered: %d\tMean Depth: %.2f\t", genome, 100*covered[genome] / total[genome ], mean;

        # compute median
        split(depths[genome], arr, " ");
        n = asort(arr);
        if (n % 2)
            print "Median Depth:", arr[(n+1)/2];
        else
            print "Median Depth:", (arr[n/2] + arr[n/2+1]) / 2.0;
    }
}' $PILEUP
# > $OUTPUT



# gawk -F'\t'  '
#     function round2(num) {
#         return int(num * 100 + (num < 0 ? -0.5 : 0.5)) / 100
#     }
#     function abs(x){
#         return ((x < 0.0) ? -x : x)
#     };
# BEGIN {
#     # Define the coverage thresholds
#     thresholds[1]=1
#     thresholds[10]=1
#     thresholds[50]=1
#     thresholds[100]=1
# }
# {
#     if(NR==FNR){
#         a[$1][$2]=$4;
#         next
#     } else {
#         if ($1 ~ /^[^@]/ ){
#             ilen[$3]=$11;
#             ireads[$3]+=1;
#             if ($9 < 0){
#                 start = $4-abs($9)
#                 end=$4+1
#             } else {
#                 end=$4+abs($9)
#                 start=$4
#             }
#             end=$4+abs($9)
#             start=$4
#             for (t = start; t <= end; t++) {
#                 countpos[t] = countpos[t] + 1
#             }

#             dep[$3][start][end]+=1
#             for(i=start;i<end;i++){
#                 dep2[$3][i]+=1
#             };
#         }

#     };

# }
# END{
#     for (i in ireads){
#         for (pos in a[i]){
#             sumsq[i]+=a[i][pos]^2
#         }

#         tlen+=ilen[i];

#         aligned+=ireads[i];
#         s=0
#         delete f
#         min=0
#         delete count
#         delete deleteMarks
#         delete seen
#          for(j in dep[i]){
#             max=0

#             j=j+0
#             for (y in dep[i][j]){
#                 y=y+0

#                 plus=dep[i][j][y] * (y  - j )
#                 sum[i]+=plus;
#                 total_bases+=dep[i][j][y];

#                 if (max < y){
#                     max = y
#                 }

#             }
#             seen[j]=max
#         }


#         lastt=0
#         nextt=0
#         delete marks

#         PROCINFO["sorted_in"] = "@ind_num_asc"
#         for (t in seen)  {
#                 left=t+0
#                 right=seen[t]+0
#                 result=0
#                 if (nextt >= left && right > nextt){
#                     result=1
#                     nextt=right
#                     if (left < lastt){
#                         lastt=left

#                     }
#                 } else if (nextt < left){
#                     result=2
#                     if (left < lastt || lastt == 0){
#                         lastt=left
#                     } else if (left > nextt){
#                         lastt=left
#                     }
#                     nextt=right
#                 }
#                 if(result>0){
#                     main=right-left+1
#                     new=nextt-lastt+1
#                     marks[lastt]=nextt
#                 }



#         }
#         for(k in marks){
#             cov[i]+=(marks[k]-k)
#         }
#     }


#     for(i in ireads){
#         aligned+=ireads[i];
#         for(j in dep2[i]){

#             sum2[i]+=dep2[i][j];
#             sumsq2[i]+=dep2[i][j]^2;
#             total_bases+=dep2[i][j];
#             cov2[i]+=1;
#         }

#         tlen+=ilen[i];
#     }

#     for(i in ireads){
#         sfactor=(ilen[i]/tlen);
#         adjsum[i]=(sfactor*sum[i]);
#         adjtotal+=(sfactor*sum[i]);
#     }

#      for(i in ireads){
#         o1=ireads[i];
#         o3=o1/aligned;
#         num_keys = 0
#         for (p in countpos) {
#             num_keys++
#          }
#         if (num_keys == 0) {
#             o6 = 0
#         } else {
#             o6 = (num_keys/ilen[i])
#         }
#         # result_str = ""
#         # for (t in thresholds) {
#         #     threshold = t+0
#         #     threshold_count = 0
#         #     for (p in countpos) {
#         #         amount = countpos[p]+0
#         #         if ( amount  >= threshold) {
#         #             threshold_count++
#         #         }
#         #     }
#         #     result_str = result_str "," threshold "x:" round2(threshold_count/ilen[i])
#         # }
#         o4=(sum[i]/total_bases);
#         o5=(adjsum[i]/adjtotal);
#         # o8=(substr(result_str, 2))
#         # o6=(cov[i]/ilen[i]);
#         mean=sum[i]/ilen[i];
#         stdev=sqrt(abs(sumsq[i]-sum[i]^2/ilen[i])/ilen[i]);
#         o9=stdev/mean;
#         printf("%s\t%.0f\t%.0f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\n", i, ilen[i], o1, o3, o4, o5, o6, mean, stdev, o9);
#     }
# }'  $PILEUP "$tmp/tmp2.sam"    > $OUTPUT
# cat $OUTPUT
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
