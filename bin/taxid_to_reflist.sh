#!/bin/bash

#	FUNCTIONS
usage()
{
cat << EOF
Help message for taxid_to_reflist.sh:

DESCRIPTION:
	Returns a list of reference accessions with associated metadata (importantly, the ftp link for download), for the input tax ID (taxid, taxonomic identifier). Running this script will download the following file at the script's location (if it does not already exist there):
		ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt

NOTES:
	- refseq summary file will be downloaded if does not already exist at \$scriptdir
	- All echo prints are from STDERR (will not be redirected with '>')
	- a file containing a single row from this script's output is to be used as input to 'refseq_download_single.sh' to fetch the reference fasta.gz file
	- test taxid 1428 is for the species "Bacillus thuringiensis"
	- for 'refseq_download_single.sh' script, input will be for 1257079 "Bacillus thuringiensis DAR 81934"

OPTIONS:
	-h	help	show this message
	-i	INT		taxid to find refseq assembly lines for
	-a  FILE	assembly_refseq.txt file (optional)
	-o	FILE	output filename (should include full path)

USAGE:
bash taxid_to_reflist.sh -i <taxid integer> -i </full/path/to/output.txt>
bash taxid_to_reflist.sh -i 1428 -o "/data/sandbox/test_output-1428.txt"
bash taxid_to_reflist.sh -i 1257079 -o "/data/sandbox/test_output-1257079.txt"

EOF
}


assembly_provided="false"

#	ARGUMENTS
# parse args
while getopts "hi:o:a:" OPTION
do
	case $OPTION in
		h) usage; exit 1 ;;
		i) TAXID=$OPTARG ;;
		o) OUTPUT=$OPTARG ;;
		a) ASSEMBLY=$OPTARG ;;
		?) usage; exit ;;
	esac
done
# check args
if [[ -z "$TAXID" ]]; then printf "%s\n" "Please specify a taxid INT (-i)."; exit; fi
if [[ -z "$OUTPUT" ]]; then printf "%s\n" "Please specify a filename (include full path) (-o)."; exit; fi



# setup other variables
absolute_path_x="$(readlink -fn -- "$0"; echo x)"
absolute_path_of_script="${absolute_path_x%x}"
scriptdir=$(dirname "$absolute_path_of_script")
runtime=$(date +"%Y%m%d%H%M%S%N")

#	MAIN

# get updated refseq assembly summary (ras)
#	NOTE: delete the 'assembly_summary_refseq.txt' file at `$scriptdir/` to force download
if [[ ! -f $ASSEMBLY ]]; then
	>&2 echo "downloading `assembly_summary_refseq.txt`"
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt -O $scriptdir/"assembly_summary_refseq.txt"
	ASSEMBLY="$scriptdir/assembly_summary_refseq.txt"
else
	mod_datetime=$(date -r "$scriptdir/assembly_summary_refseq.txt")
	>&2 echo "skipping download of 'assembly_summary_refseq.txt'"
	>&2 echo "    $mod_datetime"
fi

echo "><<<<<"
echo  $TAXID $ASSEMBLY $OUTPUT
echo "><<<<<"
awk -F'\t' -v taxid="$TAXID" '{if($6==taxid){printf("%s\n",$0)}}' $ASSEMBLY >> "$OUTPUT"
# echo "done" > $OUTPUT




exit 0