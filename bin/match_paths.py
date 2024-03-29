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
from collections import defaultdict
import sys
import os
import gzip
import matplotlib.pyplot as plt
import argparse
import numpy as np
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
        required=True,
        metavar="INPUT",
        help="BAM File to process",
    )
    parser.add_argument(
        "-d",
        "--depth",
        metavar="DEPTH",

        help="Depth File (from samtools) corresponding to your sample",
    )
    parser.add_argument(
        "-m",
        "--match",
        metavar="MATCH",
        help="Reference accessions match, 2 cols min: accession (NZ/NC), name",
    )
    parser.add_argument(
        "-r",
        "--min_reads_align",
        metavar="MINREADSALIGN",
        default=3,
        help="Filter for minimum reads aligned to reference per organism. Default is 1",
    )
    parser.add_argument(
        "-v",
        "--mincoverage",
        metavar="Minimum coverage value to consider acceptable cutoff for confidence. Anything above == confidence",
        default=1,
        type=int,
    )
    parser.add_argument(
        "-c",
        "--capval",
        metavar="Cap val for depth for entropy calculation. Default: 15",
        default=15,
        type=int,
        help="At what threshold to cutoff for determining shannon entropy stats for organism depth of cov. Default is 15, -1 for no cap value",
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
        "-x",
        "--coverage",
        metavar="COVERAGEFILE",
        default=None,
        help="Samtools coverage file",
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

def calculate_entropy(values):
    """Calculate the Shannon entropy of the given values."""
    probabilities = np.array(values) / sum(values)
    return -sum(p * np.log2(p) for p in probabilities if p > 0)

def calculate_gini(array):
    """Calculate the Gini coefficient of a numpy array."""
    # Check for zero-size array to avoid ValueError
    if array.size == 0:
        return 0  # Define behavior for empty arrays, perhaps Gini = 0

    # Continue with Gini calculation
    array = array.flatten()  # All values must be treated equally, arrays must be 1D.
    if np.amin(array) < 0:
        raise ValueError("Array cannot contain negative values.")
    # Ensure the array does not sum to zero
    if np.sum(array) == 0:
        return 0  # No inequality if there is no coverage

    # Sort the array
    array_sorted = np.sort(array)
    n = array.shape[0]
    index = np.arange(1, n+1)
    # Gini coefficient calculation
    return (np.sum((2 * index - n - 1) * array_sorted)) / (n * np.sum(array_sorted))

def make_plot(reference_coverage, plotname="test.png"):
    # Sample data: Replace these lists with your actual data
    mean_coverages = []
    gini_coefficients = []
    for ref, data in reference_coverage.items():
        mean_coverages.append(data['depth_of_coverage'])
        gini_coefficients.append(data['gini_coefficient'])
    # convert mean coverage to np list


    adjusted_mean_coverages = np.array(mean_coverages)
    # adjusted_mean_coverages = np.log([x + 0.01 for x in mean_coverages])
    # adjusted_mean_coverages = np.sqrt([x for x in mean_coverages])
    # adjusted_mean_coverages = np.log([x + 0.01 for x in mean_coverages])


    # Fit a trendline in the log-transformed space
    z = np.polyfit(adjusted_mean_coverages, gini_coefficients, 1)
    p = np.poly1d(z)
    x_for_plot = np.linspace(adjusted_mean_coverages.min(), adjusted_mean_coverages.max(), 100)
    y_for_plot = p(x_for_plot)

    # Generate x values from the minimum to the maximum log-transformed coverage for plotting the trendline
   # Create a scatter plot in the log-transformed space
    plt.scatter(adjusted_mean_coverages, gini_coefficients, color='blue', label='Data Points')
    plt.plot(x_for_plot, y_for_plot, "r--", label='Trendline in Space')

    # Label the axes and provide a title
    plt.xlabel('Log of Mean Coverage (adjusted)')
    plt.ylabel('Gini Coefficient')
    plt.title('Coverage Uniformity Analysis')
    plt.legend()

    # Display the plot
    plt.show()

# Function to calculate weighted mean
def calculate_weighted_mean(data):
    total_weight = sum(pair[0] for pair in data)  # Sum of all weights (numreads)
    weighted_sum = sum(pair[0] * pair[1] for pair in data)  # Sum of weight*value
    weighted_mean = weighted_sum / total_weight
    return weighted_mean


def count_reference_hits(bam_file_path, depthfile, covfile, matchdct):
    """
    Count the number of reads aligned to each reference in a BAM file.

    Args:
    bam_file_path (str): Path to the BAM file.

    Returns:
    dict: A dictionary with reference names as keys and counts of aligned reads as values.
    """
    # Initialize a dictionary to hold the count of reads per reference
    reference_counts = {}
    reference_coverage = defaultdict(lambda: defaultdict(dict))
    reference_lengths = {}
    notseen = set()
    unaligned = 0
    aligned_reads = 0
    total_reads = 0
    amount_pre_read = dict()
    # Open the BAM file for reading


    with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
        for ref in bam_file.header.references:
            reference_lengths[ref] = bam_file.get_reference_length(ref)
            reference_coverage[ref] = dict(
                length = reference_lengths[ref],
                depths =  defaultdict(int),
                gini_coefficient = 0,
                breadth_of_coverage = 0,
                depth_of_coverage = 0,
                mapqs = [],
                baseqs = [],
                meanmapq = 0,
                numreads = 0,
                meanbaseq = 0,
                coverage = 0,
                meandepth = 0
            )
        if covfile:
            print("Reading coverage information from coverage file")
            with open(covfile, 'r') as f:
                for line in f:
                    line = line.strip()
                    splitline = line.split('\t')
                    # if not header line
                    if "#rname"  in splitline:
                        continue
                    reference_name = splitline[0]
                    numreads = int(splitline[3])
                    covbases = float(splitline[4])
                    coverage = float(splitline[5])
                    meandepth = float(splitline[6])
                    meanbaseq = float(splitline[7])
                    meanmapq = float(splitline[8])
                    if reference_name in reference_coverage:
                        reference_coverage[reference_name]['numreads'] = numreads
                        reference_coverage[reference_name]['coverage'] = coverage
                        reference_coverage[reference_name]['meandepth'] = meandepth
                        reference_coverage[reference_name]['meanbaseq'] = meanbaseq
                        reference_coverage[reference_name]['meanmapq'] = meanmapq
            f.close()
        if not depthfile or not covfile:
            print("No depthfile or covfile supplied, reading input from bam file")
            for read in bam_file.fetch(until_eof=True):
                total_reads += 1
                if not read.is_unmapped:  # Check if the read is aligned
                    reference_name = bam_file.get_reference_name(read.reference_id)
                    aligned_reads +=1
                    # Accumulate coverage data
                    if not covfile:
                        # get the baseq of the read
                        # Get the base qualities
                        base_qualities = read.query_qualities
                        # Convert base qualities to ASCII
                        base_qualities_ascii = [chr(qual + 33) for qual in base_qualities]
                        # Convert ASCII to integers
                        base_qualities_int = [ord(qual) - 33 for qual in base_qualities_ascii]
                        baseq = np.mean(base_qualities)
                        # get the mapq of the read
                        mapq = float(read.mapping_quality)
                        reference_coverage[reference_name]['mapqs'].append(mapq)
                        reference_coverage[reference_name]['baseqs'].append(baseq)
                        reference_coverage[reference_name]['numreads']+=1
                    if not depthfile:
                        for pos in range(read.reference_start, read.reference_end):
                            reference_coverage[reference_name]['depths'][pos] += 1
                else:
                    unaligned += 1
        if not covfile:
            # make meanbaseq and meanmapq
            for ref, data in reference_coverage.items():
                if data['numreads'] > 0:
                    data['meanbaseq'] = np.mean(data['baseqs'])
                    data['meanmapq'] = np.mean(data['mapqs'])

    bam_file.close()
    if depthfile:
        print("Reading depth information from depth file")
        # read in depthfile and reference_coverage[reference_name][pos]  as value
        with open(depthfile, 'r') as f:
            for line in f:
                splitline = line.split('\t')
                reference_name = splitline[0]
                pos = int(splitline[1])
                depth = int(splitline[2])
                reference_coverage[reference_name]['depths'][pos] = depth
        f.close()
     # make a new dict which is name then count, percentage across total reads and percentage across aligned reads
    i=0
    # Aggregate coverage by organism
    organism_coverage = {}
    for ref, data in reference_coverage.items():
        organism = matchdct.get(ref)
        if organism:
            if organism not in organism_coverage:
                organism_coverage[organism] = {
                    'total_length': 0,
                    'depths': defaultdict(int),
                    "breadth_of_coverage": 0,
                    "mean_coverage": 0,
                    "gini_coefficient": 0,
                    "depth_of_coverage": 0,
                    "numreads": 0,
                    "accessions": [],
                    "mapqs": [],
                    "baseqs": [],
                    "meanmapq": 0,
                    "meanbaseq": 0

                }
            organism_info = organism_coverage[organism]
            organism_info['numreads'] += data['numreads']
            organism_info['total_length'] += data['length']
            organism_info['mapqs'].append([data['length'], data['meanmapq']])
            organism_info['baseqs'].append([data['length'], data['meanbaseq']])
            organism_info['accessions'].append(ref)
            for pos, depth in data['depths'].items():
                # Offset positions for each accession to ensure uniqueness
                adjusted_pos = pos + organism_info['total_length'] - data['length']
                organism_info['depths'][adjusted_pos] += depth
    # Now, calculate metrics for each organism based on aggregated coverage

    for organism, data in organism_coverage.items():
        # calulcate mean baseq and meanmapq for 2 index list, relative to first index length with baseq
        # and mapq as second index
        # calculate mean baseq and meanmapq
        if data['numreads'] > 0:
            weighted_baseqs_mean = calculate_weighted_mean(data['baseqs'])
            weighted_mapqs_mean = calculate_weighted_mean(data['mapqs'])
            data['meanbaseq'] = weighted_baseqs_mean
            data['meanmapq'] = weighted_mapqs_mean

        # Convert depths to a numpy array or ensure it has content before operations
        depths_array = np.array(list(data['depths'].values())) if data['depths'].values() else np.array([0])
        # Now, since depths_array is guaranteed to be non-empty (at least containing [0]),
        # you can safely perform numpy operations without encountering the zero-size array error.
        breadth_of_coverage = np.count_nonzero(depths_array) / data['total_length'] * 100 if data['total_length'] > 0 else 0
        mean_coverage = np.mean(depths_array)  # Safe due to the default value
        gini_coefficient = calculate_gini(depths_array)  # calculate_gini handles empty arrays as shown above
        # Calculate the mean depth of coverage
        mean_depth = np.sum(depths_array) / data['total_length'] if data['total_length'] > 0 else 0
        # Update the organism info with the calculated metrics
        organism_coverage[organism]['breadth_of_coverage'] = breadth_of_coverage
        organism_coverage[organism]['mean_coverage'] = mean_coverage
        organism_coverage[organism]['gini_coefficient'] = gini_coefficient
        organism_coverage[organism]['depth_of_coverage'] = mean_depth

    # make_plot(organism_coverage, "test.png")
    return organism_coverage, total_reads

def main():
    args = parse_args()
    inputfile = args.input
    pathogenfile = args.pathogens
    covfile = args.coverage
    output = args.output
    matcher = args.match
    matchdct = dict()
    organisms = set()
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
    capval  = args.capval
    mincov = args.mincoverage
    reference_hits, total_reads = count_reference_hits(
        inputfile,
        args.depth,
        covfile,
        matchdct
    )

    if args.min_reads_align:
        # filter the reference_hits based on the minimum number of reads aligned
        print(f"Filtering for minimum reads aligned: {args.min_reads_align}")
        reference_hits = {k: v for k, v in reference_hits.items() if v['numreads'] >= int(args.min_reads_align)}
    for ref, data in reference_hits.items():
        print(f"Reference: {ref}")
        print(f"\tNumber of Reads: {data['numreads']}")
        print(f"\tMean Coverage: {data['mean_coverage']}")
        print(f"\tGini Coefficient: {data['gini_coefficient']}")
        print(f"\tBreadth of Coverage: {data['breadth_of_coverage']}")
        print(f"\tDepth of Coverage: {data['depth_of_coverage']}")
        print(f"\tMean BaseQ: {data['meanbaseq']}")
        print(f"\tMean MapQ: {data['meanmapq']}")
        print()
    write_to_tsv(
        reference_hits=reference_hits,
        pathogens=pathogens,
        output_file_path=output,
        sample_name=args.samplename,
        sample_type = args.sampletype,
        total_reads = total_reads
    )
def format_non_zero_decimals(number):
    # Convert the number to a string
    num_str = str(number)
    if '.' not in num_str:
        # If there's no decimal point, return the number as is
        return number
    else:
        # Split into integer and decimal parts
        integer_part, decimal_part = num_str.split('.')
        non_zero_decimals = ''.join([d for d in decimal_part if d != '0'][:2])
        # Count leading zeros in the decimal part
        leading_zeros = len(decimal_part) - len(decimal_part.lstrip('0'))
        # Construct the new number with two significant decimal digits
        formatted_number = f"{integer_part}.{('0' * leading_zeros) + non_zero_decimals}"
        # Convert back to float, then to string to remove trailing zeros
        return str(float(formatted_number))



def write_to_tsv(reference_hits, pathogens, output_file_path, sample_name="No_Name", sample_type="Unknown", total_reads=0):
    """
    Write reference hits and pathogen information to a TSV file.

    Args:
    reference_hits (dict): Dictionary with reference names as keys and counts as values.
    pathogens (dict): Dictionary with reference names as keys and dictionaries of additional attributes as values.
    output_file_path (str): Path to the output TSV file.
    """
    with open(output_file_path, 'w') as file:
        # Write the header row
        total_reads_aligned = 0
        for ref, data in reference_hits.items():
            total_reads_aligned += data['numreads']

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

            countreads = count['numreads']
            # if total_reads_aligned is 0 then set percent_aligned to 0
            if total_reads_aligned == 0:
                percent_aligned = 0
            else:
                percent_aligned = format_non_zero_decimals(countreads / total_reads_aligned)
            if total_reads == 0:
                percent_total = 0
            else:
                percent_total = format_non_zero_decimals(countreads / total_reads)
            # Assuming 'count' is a simple value; if it's a dictionary or complex structure, adjust accordingly.
            file.write(f"{ref}\t{sample_name}\t{sample_type}\t{percent_aligned}\t{percent_total}\t{countreads}\t{is_annotated}\t{sites}\t{is_pathogen}\t{taxid}\t{status}\n")

if __name__ == "__main__":
    sys.exit(main())


