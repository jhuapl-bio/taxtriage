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
    Downloads the genome (*.fasta.gz) from a single row pulled from the "assembly_refseq_summary.txt" file:
        ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt

NOTES:
    - Output fasta.gz at (-o) will have the format:
        "\$OUTDIR/\${acc}-\${taxid}-\${version_status}-\${assembly_level}-\${release_type}-\${name}.fasta.gz"
    - All echo prints are from STDERR (will not be redirected with '>')

OPTIONS:
    -h    help    show this message
    -i    FILE    a file containing a single row from assembly_summary_refseq.txt
    -o    DIR        full path to output directory

USAGE:
bash refseq_download_single.sh -i <file with single row from 'assembly_summary_refseq.txt'> -o </full/path/to/output.txt>
bash refseq_download_single.sh -i "/data/sandbox/test_output-1257079.txt" -o "/data/sandbox"


EOF
}




#    ARGUMENTS
# parse args
while getopts "hi:o:" OPTION
do
    case $OPTION in
        h) usage; exit 1 ;;
        i) REFROW=$OPTARG ;;
        o) OUTDIR=$OPTARG ;;
        ?) usage; exit ;;
    esac
done
# check args
if [[ -z "$REFROW" ]]; then printf "%s\n" "Please specify a file containing a single row from assembly_summary_refseq.txt (-i). Exiting."; exit; fi
if [[ ! -f "$REFROW" ]]; then printf "%s\n" "File (-i) does not exist. Exiting."; exit; fi
if [[ -z "$OUTDIR" ]]; then printf "%s\n" "Please specify an output directory (include full path) (-o)."; exit; fi
if [[ ! -d "$OUTDIR" ]]; then printf "%s\n" "Directory (-o) does not exist. Exiting."; exit; fi

# setup other variables
absolute_path_x="$(readlink -fn -- "$0"; echo x)"
absolute_path_of_script="${absolute_path_x%x}"
scriptdir=$(dirname "$absolute_path_of_script")
runtime=$(date +"%Y%m%d%H%M%S%N")
#    MAIN

# pull ftp paths and download reference genomes for accessions in input tsv
#    1    assembly_accession
#    6    taxid                <- strain if available, otherwise is species taxid
#    7    species_taxid
#    8    organism_name
#    9    infraspecific_name
#    20    ftp_path


acc=$(cut -f1 "$REFROW")
taxid=$(cut -f6 "$REFROW")
version_status=$(cut -f11 "$REFROW")
assembly_level=$(cut -f12 "$REFROW")
release_type=$(cut -f13 "$REFROW")
name=$(cut -f8 "$REFROW" | sed 's/ /_/g')
path=$(cut -f20 "$REFROW")
bn=$(basename "$path")
>&2 echo "wgetting..."
>&2 echo "    accession:      $acc"
>&2 echo "    taxid:          $taxid"
>&2 echo "    version_status: $version_status"
>&2 echo "    assembly_level: $assembly_level"
>&2 echo "    release_type:   $release_type"
>&2 echo "    org name:       $name"


echo "$path/${bn}_genomic.fna.gz"
wget "$path/${bn}_genomic.fna.gz" --output-document "$OUTDIR/${acc}-${taxid}-${version_status}-${assembly_level}-${release_type}-${name}.fasta.gz"
# addend `2> /dev/null` if you want to suppress the progress stats











#    STOPPING HERE FOR NOW.
#    below is some cleanup utilized in the META platform:
#    1. removing linebreaks... super annoying ncbi!
#    2. remove smaller sequence fragments that could cause issues when simulating reads
exit
# some cleanup if desired later
gunzip -f "$outdir/$taxid.fasta.gz"
# put all sequence strings under each header into a single line
awk '{if(NR == 1){printf("%s\n", $0)}else{if(substr($0,1,1) == ">"){printf("\n%s\n", $0)} else {printf("%s", $0)}}}END{printf("\n")}' "$outdir/$taxid.fasta" > "$outdir/$taxid.fasta.tmp"
# only retain contigs/assemblies >1 Kbp
sed $'$!N;s/\\\n/\t/' "$outdir/$taxid.fasta.tmp" | awk -F'\t' '{if(length($2)>=1000){printf("%s\n%s\n",$1,$2)}}' > "$outdir/$taxid.fasta"
rm "$outdir/$taxid.fasta.tmp"









