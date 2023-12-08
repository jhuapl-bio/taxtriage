#!/bin/bash
# **********************************************************************
# Copyright (C) 2019 Johns Hopkins University Applied Physics Laboratory
#
# All Rights Reserved.
# For any other permission, please contact the Legal Office at JHU/APL.
# **********************************************************************

# This script takes a list of taxids, and a taxonomy path with names.dmp and nodes.dmp
# and creates a two-column file with taxid, and full taxonomic string
#      by: Thomas Mehoke
# written: July 14, 2015
# edited:	20190327 - RP (fix for binftoolkit-specific pathing)

runtime=$(date +"%Y%m%d%H%M%S%N")

#-------------------------------------------------
usage() {
  cat <<EOF
usage: $0 -i <input> -t <taxonomy path> -o <output>


OPTIONS:
   -h      show this message
   -i      path to file containing a list of taxids, one per line
   -t      path to taxonomy folder containing names.dmp and nodes.dmp
   -o      file to place text output file (default: STDOUT)
   -w      working directory (default: /tmp)

EOF
}

#-------------------------------------------------
create() {
  cat <<EOF

Create joined.dmp in your taxonomy folder
from names.dmp and nodes.dmp as follows:

     join -j 1 -t $'\t' -o 0 2.2 1.2 2.3 \\
       <(grep "scientific name" names.dmp | cut -f1,3 | sort -k1,1) \\
       <(cut -f1,3,5 nodes.dmp | sort -k1,1) | sort -n > joined.dmp

EOF
}

#-------------------------------------------------
# set default values
tempdir="/tmp"
outputfile=""
# RP add 20200214
# absolute path to script dir
absolute_path_x="$(
  readlink -fn -- "$0"
  echo x
)"
absolute_path_of_script="${absolute_path_x%x}"
scriptdir=$(dirname "$absolute_path_of_script")

# parse input arguments
while getopts "hi:t:o:d:w:" OPTION; do
  case $OPTION in
  h)
    usage
    exit 1
    ;;
  i) input=$OPTARG ;;
  t) taxpath=$OPTARG ;;
  o) outputfile=$OPTARG ;;
  d) THREADS=$OPTARG ;;
  w) tempdir=$OPTARG ;;
  ?)
    usage
    exit
    ;;
  esac
done

# if necessary arguments are not present, display usage info and exit
if [[ -z "$input" ]]; then
  echo "Specify a taxid list with -i" >&2
  usage
  exit 2
fi
if [[ -z "$taxpath" ]]; then
  echo "Select a path to the taxonomy folder with -t" >&2
  usage
  exit 2
elif ! [[ -d "$taxpath" ]]; then
  echo "Error: taxonomy path \"$taxpath\" does not exist." >&2
  usage
  exit 2
elif ! [[ -s "$taxpath/joined.dmp" ]]; then
  echo "Error: Taxonomy path \"$taxpath\" does not contain joined.dmp" >&2
  create
  exit 3
fi

#===================================================================================================

workdir="$tempdir/krakenreport_fullstring.sh-$runtime"
mkdir -p "$workdir"
OUTPUT="$workdir/output"

# get full string for all taxids
gawk -F $'\t' -f $scriptdir/get_fullstring.awk <(cut -f1 "$input") "$taxpath/joined.dmp" >"$OUTPUT"

# output to STDOUT
if [[ -z "$outputfile" ]]; then
  paste "$input" "$OUTPUT"
# or output to output file specified by -o argument
else
  paste "$input" "$OUTPUT" >"$outputfile"
fi

rm -rf "$workdir"

#~~eof~~#
