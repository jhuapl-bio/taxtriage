#!/bin/bash
#  **********************************************************************
#  Copyright (C) 2023 Johns Hopkins University Applied Physics Laboratory
#
#  All Rights Reserved.
#  For any other permission, please contact the Legal Office at JHU/APL.
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#  **********************************************************************

#	FUNCTIONS
#===============================================================================
usage() {
  cat <<EOF
NOTES:
	- Comparisons are only valid down to the *highest* (not lowest) taxonomic rank in the ground truth.
		For example, if your ground truth consists of 9 strans and 1 species, the comparison
		metrics calculated only for species level and above are reliable.
	- This will be up to the user to track and be aware of when they are designing their input tsv.

DEPENDENCIES:
	GNU Parallel
		O. Tange (2011): GNU Parallel - The Command-Line Power Tool,
		;login: The USENIX Magazine, February 2011:42-47.
	NCBI Taxonomy data
		will be downloaded and formatted on first run
		ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz


REQUIRED:
	-h	help	show this message
	-t	INT		number of threads to GNU parallel over
	-b	FILE	truth abundance profile tsv (the input tsv for META evaluation)
				 *mode 1 (simulation), this is the same file used for simulating reads
				 *there will be a different 'evaluation script' for modes 2 and 3
	-i	DIR		directory containing formatted classifier output
	-o	DIR		directory for tmp files and final comparison output files for visualization

OPTIONAL:
	N/A





EXAMPLE:	using kraken and meseeks classifiers (meseeks is one I've been
			working on for another project, don't worry about including this as part of META...
			yet, or probably ever since we'll want to maintain our unbiased-ness

# format classifier output
w="/data/apps/src/meta_system_system_metrics_evaluation_parsers"
o="$w/classifier_output_formatted/"
i="$w/classifier_output/kraken.report"
bash $w/parse_kraken.sh -i "$i" -o "$o"
i="$w/classifier_output/kraken2.report"
bash $w/parse_kraken2.sh -i "$i" -o "$o"
i="$w/classifier_output/krakenuniq.report"
bash $w/parse_krakenuniq.sh -i "$i" -o "$o"
i="$w/classifier_output/diamond.report"
bash $w/parse_diamond.sh -i "$i" -o "$o"
i="$w/classifier_output/mash.report"
bash $w/parse_mash.sh -i "$i" -o "$o"
i="$w/classifier_output/centrifuge.report"
bash $w/parse_centrifuge.sh -i "$i" -o "$o"
i="$w/classifier_output/kraken2_bracken.report"
bash $w/parse_bracken.sh -i "$i" -o "$o"


# run comparison
t="10"; w="/data/apps/src/meta_system_system_metrics_evaluation_parsers"
b="$w/sim_meta_95gg9031_05vv10245.tsv"
i="$w/classifier_output_formatted/"
o="$w/metacompare_test_output"
bash $w/metacompare.sh -t "$t" -b "$b" -i "$i" -o "$o"




EOF
}

#	DEFAULTS & INPUTS & CHECKS
#===============================================================================
#	notes:
#		echo $? (0 = successful execution)
# absolute path to script dir
absolute_path_x="$(
  readlink -fn -- "$0"
  echo x
)"


absolute_path_of_script="${absolute_path_x%x}"
scriptdir=$(dirname "$absolute_path_of_script")
bin="$scriptdir"

# parse args
while getopts "ht:b:i:o:d:f" OPTION; do
  case $OPTION in
  h)
    usage
    exit 1
    ;;
  t) THREADS=$OPTARG ;;
  b) BASELINE=$OPTARG ;;
  f) FORCE="1";;
  i) INDIR=$OPTARG ;;
  d) TAXDUMP=$OPTARG;;
  o) OUTPUT=$OPTARG ;;
  ?)
    usage
    exit
    ;;
  esac
done
# check args
if [[ -z "$THREADS" ]]; then
  printf "%s\n" "Please specify number of threads (-t)."
  exit
fi
if [[ -z "$BASELINE" ]]; then
  printf "%s\n" "Please specify META abundance profile tsv (-b)."
  exit
fi
if [[ ! -f "$BASELINE" ]]; then
  printf "%s\n" "The input (-b) $BASELINE file does not exist."
  exit
fi
if [[ -z "$INDIR" ]]; then
  printf "%s\n" "Please specify input directory (-i)."
  exit
fi
if [[ ! -d "$INDIR" ]]; then
  printf "%s\n" "The input (-i) $INDIR directory does not exist."
  exit
fi
if [[ -z $OUTPUT ]]; then
  printf "%s\n" "Please specify a final output directory (-o)."
  exit
fi
if [[ -z $TAXDUMP ]]; then
  printf "%s\n" "Taxdump not specified, downloading automatically."
fi


if [[ ! -d "$OUTPUT" ]]; then mkdir -p "$OUTPUT"; fi

# setup other variables
outdir="$OUTPUT"
tmp="$outdir/tmp"
if [[ ! -d "$tmp" ]]; then
  mkdir -p "$tmp"
fi
echo $TAXDUMP


# setup for taxid2taxstring
if [[ ! -f "$TAXDUMP/nodes.dmp" ]] || [[ ! -f $TAXDUMP'/names.dmp' ]]; then
  echo >&2 "getting ncbi taxonomy files"
  mkdir -p "$scriptdir/taxdump"
  wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -O "$scriptdir/taxdump/taxdump.tar.gz"
  gunzip "$scriptdir/taxdump/taxdump.tar.gz"
  tar -C "$scriptdir/taxdump" -xf "$scriptdir/taxdump/taxdump.tar"
  TAXDUMP="$scriptdir/taxdump"
fi

echo $TAXDUMP


# grep "scientific name" "$TAXDUMP/names.dmp" | cut -f1,3 | sort -k1,1 > $TAXDUMP/1.txt
# cut -f1,3,5 "$TAXDUMP/nodes.dmp" | sort -k1,1 > $TAXDUMP/2.txt


# join -j 1 -t $'\t'  $TAXDUMP/1.txt $TAXDUMP/2.txt
echo $FORCE
if [[ ! -f "$TAXDUMP/joined.dmp" ]] || ! [[  -z $FORCE ]]; then
  join -t $'\t' -1 1 -2 1 \
    <(grep "scientific name" "$TAXDUMP/names.dmp" | cut -f1,3 | sort -k1,1) \
    <(cut -f1,3,5 "$TAXDUMP/nodes.dmp" | sort -k1,1) | sort -n > "$TAXDUMP/joined.dmp"
fi


profile_prep() {
  bn=$(basename "$1")
  echo >&2 "profile prep for: $bn"
  tmpp="$tmp/${bn}_dir"
  mkdir -p "$tmpp"
  cut -f1 "$1" >"$tmpp/taxid.list"
  # find taxstrings of each reported taxid
  $bin/taxid2taxstring.sh -i "$tmpp/taxid.list" -t "$scriptdir/taxdump" -o "$tmpp/taxid.list.ts"
  # make 3 col file of taxid, abundance, taxstring
  exit 1
  awk -F'\t' '{
		if(NR==FNR){
			abu[$1]=$2;
		}else{
			printf("%s\t%s\t%s\n",$1,abu[$1],$2);
		}
	}' "$1" "$tmpp/taxid.list.ts" >"$tmpp/taxid.abu.ts"

  # clean up taxid2taxstring of P,C,O,F,G,S per taxid to separate files
  #	first line is abundance of tip taxid
  while read tat; do
    printf "$tat" | awk -F'\t' -v tmp="$tmpp" '{
			split($3,ts,"|");
			print($2) > tmp"/"$1".taxlist";
			for(i in ts){
				print(ts[i]) > tmp"/"$1".taxlist";
			}
		}'
  done <"$tmpp/taxid.abu.ts"

  #	20200528 - adding in additional output format for sunburst where missing ranks are padded in (taxid.abu.ts.padded)
  find "$tmpp" -type f -name "*.taxlist" | sort | while read list; do
    taxid=$(basename "$list" | sed 's/\..*//')
    abu=$(head -1 "$list")
    printf "%s\t%s\t%s|" "$taxid" "$abu" "1;root(no rank)"
    for rank in "superkingdom" "phylum" "class" "order" "family" "genus" "species"; do
      # 20200730 fix - some ranks may appear on multiple lines, however the first will be the highest level (i.e. everything else is a subrank that may be ignored)
      check=$(grep "($rank)" "$list" | head -1)
      if [[ "$check" != "" ]]; then
        printf "$check|"
      else
        printf "x;unassigned($rank)|"
      fi
    done
    if [[ $(grep -A10 "(species)" "$list" | wc -l) -gt "1" ]]; then
      last=$(tail -1 "$list" | sed 's/(no rank)/(strain)/')
      printf "%s\n" "$last"
    else
      printf "x;unassigned(strain)\n"
    fi
  done | sed 's/\(|x;unassigned([A-Za-z]\+)\)\+$//' >"$tmpp/taxid.abu.ts.padded"

  # make abundance files per rank (also roll up abundance for same rank lines)
  for rank in "superkingdom" "phylum" "class" "order" "family" "genus" "species"; do
    find "$tmpp" -type f -name "*.taxlist" | sort | while read list; do
      # ensure rank exists before printing anything
      # 20200608, bug fix... some 'unclassified' intra-rank ranks... (only return first match in *.taxlist)
      row=$(grep -m1 "($rank)" "$list")
      if [[ "$row" != "" ]]; then
        abu=$(head -1 "$list")
        printf "$abu\t$row\n"
      fi
    done | awk -F'\t' '{
			abu[$2]+=$1;
		}END{
			for(ts in abu){
				printf("%.9f\t%s\n",abu[ts],ts);
			}
		}' >"$tmpp/byrank-$rank.tsv"
    # 20200527, might also need to adjust abundances here to output 100% at each rank
    # 20200608, no, this would throw off abudance prediction for correct taxid assignments

    # make file with abundance and org name per taxid
    sed -e 's/;/\t/' -e 's/(.*//' "$tmpp/byrank-$rank.tsv" >"$tmpp/taxid_abu_org-$rank.tsv"

  done
  # check if rank below species
  find "$tmpp" -type f -name "*.taxlist" | sort | while read list; do
    # first check if species row even exists
    if [[ $(grep "(species)" "$list") != "" ]]; then
      ranktip=$(grep -A1000 "(species)" "$list" | tail -1 | sed -e 's/.*(//' -e 's/)//')
      if [[ "$ranktip" != "species" ]]; then
        printf "%s\t" $(head -1 "$list")
        grep -A1000 "(species)" "$list" | tail -1
      fi
    fi
  done | awk -F'\t' '{
		abu[$2]+=$1;
	}END{
		for(ts in abu){
			printf("%.9f\t%s\n",abu[ts],ts);
		}
	}' >"$tmpp/byrank-strain.tsv"
  # make file with abundance and org name per taxid
  sed -e 's/;/\t/' -e 's/(.*//' "$tmpp/byrank-strain.tsv" >"$tmpp/taxid_abu_org-strain.tsv"

}
export tmp bin scriptdir
export -f profile_prep


# get tax info of input and each parsed classifier output
cp "$BASELINE" "$tmp/BASELINE.tsv"
# re-establish the abundance profile based on 3rd column of input tsv
#	i.e. only among rows where $3==1
awk -F'\t' '{
	if($3==1){
		t+=$2;
		abu[$1]=$2;
	}
}END{
	for(tax in abu){
		printf("%s\t%.9f\n",tax,abu[tax]/t);
	}
}' "$tmp/BASELINE.tsv" > "$tmp/BASELINE1.tsv"


find "$INDIR" -maxdepth 1 -type f | sort >>"$tmp/parallel.taxinfo"
parallel --arg-file "$tmp/parallel.taxinfo" --jobs="$THREADS" profile_prep

exit 1
#	new file added from 'metacompare_realdata.sh' script
# make intermediate consensus matrix
for rank in "superkingdom" "phylum" "class" "order" "family" "genus" "species" "strain"; do
  find "$tmp" -type f -name "taxid_abu_org-$rank.tsv" | sort -V | while read profile; do
    name=$(basename $(dirname "$profile") | sed -e 's/parsed_//' -e 's/_dir//')
    awk -F'\t' -v name="$name" -v rank="$rank" '{printf("%s\t%s\t%s\n",name,rank,$0)}' "$profile"
  done
done >"$tmp/classifier_rank_abu_taxid_org.tsv"

# add comma separated list of classifiers that identiy the taxid per line
awk -F'\t' '{
	class[$1]+=1;
	taxid[$4]+=1;
	ct[$1][$4]=$0;
}END{
	for(t in taxid){
		for(c in class){
			if(ct[c][t]!=""){
				printf("%s\t",ct[c][t]);
				for(cc in class){
					if(ct[cc][t]!=""){
						printf("%s,",cc);
					}
				}
			}
			printf("\n");
		}
	}
}' "$tmp/classifier_rank_abu_taxid_org.tsv" | awk '{if($0!=""){print($0)}}' | sed 's/,$//' >"$outdir/classifier_rank_abu_taxid_org_inclusion.tsv"

calc() {
  bn=$(printf "$1" | cut -f1)
  rank=$(printf "$1" | cut -f2)
  tmpp="$tmp/${bn}_dir"
  echo >&2 "processing $bn data - for rank $rank"
  # among all rank taxids in the pair (BASELINE v CLASSIFIER OUTPUT)
  # get list of unique taxids
  cut -f2 "$bltmp/byrank-$rank.tsv" >"$tmpp/ts.list-$rank"
  cut -f2 "$tmpp/byrank-$rank.tsv" >>"$tmpp/ts.list-$rank"
  sort "$tmpp/ts.list-$rank" | uniq >"$tmpp/ts.list-$rank.uniq"
  # and make ordered (sorted) vectors of abundances (for L2 calc)
  # col1 for BASELINE ($b_abu), col2 for CLASSIFIER OUTPUT ($c_abu)
  while read tax; do
    b_abu=$(grep "$tax" "$bltmp/byrank-$rank.tsv" | cut -f1)
    if [[ "$b_abu" == "" ]]; then b_abu="0.000000000"; fi
    c_abu=$(grep "$tax" "$tmpp/byrank-$rank.tsv" | cut -f1)
    if [[ "$c_abu" == "" ]]; then c_abu="0.000000000"; fi
    printf "$tax\t$b_abu\t$c_abu\n"
  done <"$tmpp/ts.list-$rank.uniq" | sort -t$'\t' -rnk2 >"$tmpp/abundance_vectors-$rank.tsv"
  # calculate the L2 distance between BASELINE abu abundance (p) and
  # CLASSIFIER OUTPUT abundance vector (q) at each rank
  #	vectors are of equal length
  # d(p,q) = sqrt( (q1-p1)^2 + ... + (qn-pn)^2 )
  awk -F'\t' -v rank="$rank" '{
		sum_diff_sq+=($3-$2)^2;
	}END{
		printf("%s\t%0.9f\n",rank,sqrt(sum_diff_sq));
	}' "$tmpp/abundance_vectors-$rank.tsv" >>"$tmpp/L2_table.tsv"

  # There are two AUPRC methods below (MIT [presence/absence] and APL [takes abundance into account])
  #		WILL ONLY OUTPUT MIT METHOD
  #		as it is a better complement to the L2 method
  # MIT method
  #	output format:
  #	1	rank
  #	2	>cutoff(t)
  #	3	recall
  #	4	precision
  #	5	TP
  #	6	FN
  #	7	FP
  #	8	TN
  # threshold from 0.000 to 1.000 in increments of 0.001
  seq 0.000 0.001 1.000 | while read cutoff; do
    awk -F'\t' -v rank="$rank" -v cutoff="$cutoff" '{
			# check if truth is positive
			if($2>cutoff){
				# check if prediction is positive wrt truth
				if($3>0){TP+=1;PP+=1;CP+=1}else{FN+=1;PN+=1;CP+=1};
			}else{
				if($3>0){FP+=1;PP+=1;CN+=1}else{TN+=1;PN+=1;CN+=1};
			}
		}END{
			# calc precision/recall
			if(CP==0){recall=0}else{recall=TP/CP};
			if(PP==0){precision=0}else{precision=TP/PP};
			printf("%s\t%.3f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\n",rank,cutoff,recall,precision,TP,FN,FP,TN);
		}' "$tmpp/abundance_vectors-$rank.tsv"
  done >"$tmpp/AUPR-mit_table-$rank.tsv"
  # get AUC (trapezoid method)
  #	from highest to lowest cutoff
  seq 1001 -1 2 | while read from; do
    head -$from "$tmpp/AUPR-mit_table-$rank.tsv" | tail -2 | awk -F'\t' '
		    function abs(v) {return v < 0 ? -v : v};
		{
			x[NR]=$3;
			y[NR]=$4;
		}END{
			h=abs(x[1]-x[2]);
			s=y[1]+y[2];
			aut=(h*s)/2;
			print(aut);
		}'
  done | awk -v rank="$rank" '{
		auc+=$0;
	}END{
		if(auc<0){auc=0};
		printf("%s\t%0.9f\n",rank,auc)}' >>"$tmpp/AUPRC-mit_table.tsv"

  # APL method
  #	output format:
  #	1	rank
  #	2	>cutoff(t)
  #	3	recall
  #	4	precision
  #	5	TP
  #	6	FN
  #	7	FP
  #	8	TN
  # threshold from 0.000 to 1.000 in increments of 0.001
  #	seq 0.001 0.001 1.000 | while read cutoff; do
  #		awk -F'\t' -v rank="$rank" -v cutoff="$cutoff" '{
  # check if truth is positive
  #			if($2>cutoff){
  # check if prediction is positive wrt truth
  #	&& if it is also above the cutoff (primary diff between MIT and APL methods)
  #				if($3>cutoff){TP+=1;PP+=1;CP+=1}else{FN+=1;PN+=1;CP+=1};
  #			}else{
  #				if($3>cutoff){FP+=1;PP+=1;CN+=1}else{TN+=1;PN+=1;CN+=1};
  #			}
  #		}END{
  # calc precision/recall
  #			if(CP==0){recall=0}else{recall=TP/CP};
  #			if(PP==0){precision=0}else{precision=TP/PP};
  #			printf("%s\t%.3f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\n",rank,cutoff,recall,precision,TP,FN,FP,TN);
  #		}' "$tmpp/abundance_vectors-$rank.tsv"
  #	done > "$tmpp/AUPR-apl_table-$rank.tsv"
  # get AUC (trapezoid method)
  #	from highest to lowest cutoff
  #	seq 1000 -1 2 | while read from; do
  #		head -$from "$tmpp/AUPR-apl_table-$rank.tsv" | tail -2 | awk -F'\t' '
  #           function abs(v) {return v < 0 ? -v : v};
  #       {
  #			x[NR]=$3;
  #			y[NR]=$4
  #		}END{
  #			h=abs(x[1]-x[2]);
  #			s=y[1]+y[2];
  #			aut=(h*s)/2;
  #			print(aut)
  #		}'; done | awk -v rank="$rank" '{auc+=$0}END{printf("%s\t%0.9f\n",rank,auc)}' >> "$tmpp/AUPRC-apl_table.tsv"

  # Kurskal-Wallis, H test (one-way ANOVA on ranks)
  #	https://statistics.laerd.com/spss-tutorials/kruskal-wallis-h-test-using-spss-statistics.php
  #	dependent variable:		taxa rank (ordinal)
  #				which has n groups, where n is the number of taxa identified
  #	independent variable:	classification tool (category)
  # distribution of taxa ranks per category are not the same, so only mean ranks may be compared.
  #
  #	20200323 -	leaving this out for now, I'm not sure this is necessary or useful to have
  #				the L2 and AUPRC are comprehensive for comparison

}

export tmp bltmp
export -f calc

# for each prepped classifier output, pair with prepped BASELINE at each tax level
blbn=$(basename $(head -1 "$tmp/parallel.taxinfo"))
bltmp="$tmp/${blbn}_dir"

# clear L2 and AUPRC tables before calc
find "$tmp" -name "L2_table.tsv" -exec rm {} +
find "$tmp" -name "AUPRC-mit_table.tsv" -exec rm {} +
# run calc on each rank for each classifier output
find "$INDIR" -maxdepth 1 -type f | sort | while read output; do
  bn=$(basename "$output")
  for rank in "superkingdom" "phylum" "class" "order" "family" "genus" "species" "strain"; do
    printf "$bn\t$rank\n"
  done
done >"$tmp/parallel.calc"
parallel --arg-file "$tmp/parallel.calc" --jobs="$THREADS" calc

# make final matrix
#printf "classifier_name\trank\tL2\tAUPRC-mit\tAUPRC-apl\n"
printf "classifier_name\trank\tL2\tAUPRC\n" >$outdir/eval.tsv
find "$tmp" -maxdepth 1 -type d -name "parsed_*" | sort | while read d; do
  name=$(basename "$d" | sed -e 's/parsed_//' -e 's/_dir//')
  for rank in "superkingdom" "phylum" "class" "order" "family" "genus" "species" "strain"; do
    L2=$(grep "^$rank" "$d/L2_table.tsv" | cut -f2)
    AUPRC_mit=$(grep "^$rank" "$d/AUPRC-mit_table.tsv" | cut -f2)
    #		AUPRC_apl=$(grep "^$rank" "$d/AUPRC-apl_table.tsv" | cut -f2)
    printf "$name\t$rank\t$L2\t$AUPRC_mit\n" >>$outdir/eval.tsv
  done
done

# cleanup
#rm -rf "$tmp"
