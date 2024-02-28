#!/usr/bin/env python3

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

"""Provide a command line tool to fetch a list of refseq genome ids to a single file, useful for kraken2 database building or alignment purposes"""

import sys
import os
import gzip
import argparse
import pysam



def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and transform a tabular samplesheet.",
        epilog="Example: python check_samplesheet.py samplesheet.csv samplesheet.valid.csv",
    )
    parser.add_argument(
        "-i",
        "--input",
        metavar="INPUT",
        help="BAM File to process",
    )
    parser.add_argument(
        "-m",
        "--match",
        metavar="MATCH",
        help="Reference accessions match, 2 cols min: accession (NZ/NC), name",
    )
    parser.add_argument(
        "-a",
        "--accessioncol",
        metavar="ACCCOL",
        default=0,
        help="Index of the column in mapfile (if specified) to match to the reference accession. 0 index start",
    )
    parser.add_argument(
        "-n",
        "--namecol",
        default=2,
        metavar="NAMECOL",
        help="Index of the column in mapfile (if specified) to match to the name. 0 index start",
    )
    parser.add_argument(
        "-t",
        "--sampletype",
        metavar="SAMPLETYPE",
        default='Unknown',
        help="Sample Type to process. If Empty, defaults to null",
    )
    parser.add_argument(
        "-s",
        "--samplename",
        metavar="SAMPLENAME",
        default="No_Name",
        help="Name of the sample to process. If Empty, defaults to 'No_Name'",
    )
    parser.add_argument(
        "-p",
        "--pathogens",
        metavar="PATHOGENS",
        help="TXT File to process. Must be in header format: ",
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="OUTPUT",
        type=str,
        help="Name of output directory for feature table files",
    )


    return parser.parse_args(argv)


def import_pathogens(pathogens):
    """Import the pathogens from the input file
    """
    pathogensdct = dict()
    with open(pathogens, 'r') as f:
        for line in f:
            splitline = line.split('\t')
            if len(splitline) > 0:
                pathogenname = splitline[0]
            else:
                pathogenname = None
            if len(splitline) > 1:
                taxid = splitline[1]
            else:
                taxid = None
            if len(splitline) > 2:
                callclass = splitline[2]
            else:
                callclass = None
            if len(splitline) > 3:
                sites = splitline[3]
            else:
                sites = None
            if len(splitline) > 4:
                commensal = splitline[4]
            else:
                commensal = None
            if len(splitline) > 5:
                status = splitline[5]
            else:
                status = None
            if len(splitline) > 6:
                pathology = splitline[6]
            else:
                pathology = None
            # assign these values into a dict where key is the pathogenname
            pathogensdct[pathogenname] = {
                'taxid': taxid,
                'callclass': callclass,
                'sites': sites,
                'commensal': commensal,
                'status': status,
                'pathology': pathology
            }
    f.close()
    return pathogensdct

def identify_pathogens(inputfile, pathogens):
    """Identify the pathogens in the input file"""
    with open(inputfile, 'r') as f:
        for line in f:
            if line.startswith('pathogen'):
                pathogens = line.split('\t')
                return pathogens
    f.close()
    return None

def count_reference_hits(bam_file_path, matchdct):
    """
    Count the number of reads aligned to each reference in a BAM file.

    Args:
    bam_file_path (str): Path to the BAM file.

    Returns:
    dict: A dictionary with reference names as keys and counts of aligned reads as values.
    """
    # Initialize a dictionary to hold the count of reads per reference
    reference_counts = {}
    notseen = set()
    unaligned = 0
    aligned_reads = 0
    total_reads = 0
    amount_pre_read = dict()
    # Open the BAM file for reading
    with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
        for read in bam_file.fetch(until_eof=True):
            total_reads += 1
            if not read.is_unmapped:  # Check if the read is aligned
                reference_name = bam_file.get_reference_name(read.reference_id)
                aligned_reads +=1
                # get read id
                readid = read.query_name
                if amount_pre_read.get(readid):
                    amount_pre_read[readid] += 1
                else:
                    amount_pre_read[readid] = 1
                if reference_name in matchdct:
                    reference_name = matchdct[reference_name]
                else:
                    if reference_name not in notseen:
                        notseen.add(reference_name)
                if reference_name in reference_counts:
                    reference_counts[reference_name] += 1
                else:
                    reference_counts[reference_name] = 1
            else:
                unaligned += 1
    bam_file.close()
    for ref in notseen:
        print(f"Reference {ref} not found in match file")
    # iterate through amount_pre_read and print out the ones that are greater than 1
    summation = 0
    for readid, amount in amount_pre_read.items():
        if amount > 1:
            # print(f"Read {readid} has {amount} alignments")
            summation += amount - 1
    print(f"Total reads with multiple alignments: {summation}")
    print(f"\nUnaligned reads: {unaligned}\nAligned Reads: {aligned_reads}\nTotal Reads: {total_reads}\nPercent Reads Aligned: {100*(1-(unaligned/total_reads))}\n\n")
    # make a new dict which is name then count, perentage across total reads and percentage across aligned reads
    newdct = dict()
    for ref, count in reference_counts.items():
        newdct[ref] = {
            'count': count,
            'percent_of_total': round(100*(count/total_reads), 2), # make it 2 decimal places
            'percent_of_aligned': round(100*(count/aligned_reads),2)
        }
    return newdct

def main():
    args = parse_args()
    inputfile = args.input
    pathogenfile = args.pathogens
    output = args.output
    matcher = args.match
    matchdct = dict()
    if args.match:
        # open the match file and import the match file

        with open (matcher, 'r') as f:
            accindex = args.accessioncol
            nameindex = args.namecol
            for line in f:
                splitline = line.split('\t')
                if len(splitline) > 0:
                    accession = splitline[accindex]
                else:
                    accession = None
                if len(splitline) > 1:
                    name = splitline[nameindex]
                else:
                    name = None
                matchdct[accession] = name

        f.close()
    pathogens = import_pathogens(pathogenfile)
    # Next go through the BAM file (inputfile) and see what pathogens match to the reference, use biopython
    # to do this
    reference_hits = count_reference_hits(inputfile, matchdct)

    write_to_tsv(
        reference_hits=reference_hits,
        pathogens=pathogens,
        output_file_path=output,
        sample_name=args.samplename,
        sample_type = args.sampletype
    )


def write_to_tsv(reference_hits, pathogens, output_file_path, sample_name="No_Name", sample_type="Unknown"):
    """
    Write reference hits and pathogen information to a TSV file.

    Args:
    reference_hits (dict): Dictionary with reference names as keys and counts as values.
    pathogens (dict): Dictionary with reference names as keys and dictionaries of additional attributes as values.
    output_file_path (str): Path to the output TSV file.
    """
    with open(output_file_path, 'w') as file:
        # Write the header row

        header =  "Name\tSample\tSample Type\t% Aligned\t% Total Reads\t# Aligned\tIsAnnotated\tSites\tType\tTaxid\tStatus"
        file.write(f"{header}\n")
        for ref, count in reference_hits.items():
            if ref in pathogens:
                if pathogens[ref]['callclass'] == "commensal":
                    is_pathogen = "Commensal"
                else :
                    is_pathogen = "Pathogen"
                taxid = pathogens[ref]['taxid']
                is_annotated = "Yes"
                callclass = pathogens[ref]['callclass']
                sites = pathogens[ref]['sites']
                status = pathogens[ref]['status']

            else:
                is_pathogen = "N/A"
                is_annotated = "No"
                taxid = ""
                callclass = ""
                sites = ""
                status = ""

            countreads = count['count']
            percent_aligned = count['percent_of_aligned']
            percent_total = count['percent_of_total']
            # Assuming 'count' is a simple value; if it's a dictionary or complex structure, adjust accordingly.
            file.write(f"{ref}\t{sample_name}\t{sample_type}\t{percent_aligned}\t{percent_total}\t{countreads}\t{is_annotated}\t{sites}\t{is_pathogen}\t{taxid}\t{status}\n")

if __name__ == "__main__":
    sys.exit(main())


