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
import statistics
import math as Math
import argparse
import re
import csv
import math
import pysam
from math import log2



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
        type=int,
        help="Filter for minimum reads aligned to reference per organism. Default is 1",
    )
    parser.add_argument(
        "--min_cds_found",
        metavar="MINCDSFOUND",
        default=3,
        type=int,
        help="Filter for minimum coding regions found per organism post denovo assembly. Default is 3",
    )
    parser.add_argument(
        "--k2",
        metavar="K2FILE",
        default=None,
        help="Provide kraken2 output report to match taxid to reference",
    )
    parser.add_argument(
        "-v",
        "--mincoverage",
        metavar="Minimum coverage value to consider acceptable cutoff for confidence. Anything above == confidence",
        default=1,
        type=int,
    )
    parser.add_argument(
        "--diamond",
        metavar="Output of diamond txt file (outfmt:6)",
        default=None,
        type=str,
        help="DIAMOND BLASTX Output file",
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
        "--parent_k2_match",
        metavar="PARENTK2MATCH",
        default="G",
        type=str,
        choices=["G", "F", "O", "C", "P", "K", "D", "R"],
        help="Which parent to match for k2 hierarchy of taxonomy for disparity across genus, family, etc calls",
    )
    parser.add_argument(
        "-a",
        "--accessioncol",
        metavar="ACCCOL",
        default=0,
        help="Index of the column in mapfile (if specified) to match to the reference accession. 0 index start",
    )
    parser.add_argument(
        "-q",
        "--taxcol",
        metavar="TAXCOL",
        default=4, required =False,
        help="Index of the column in mapfile (if specified) to add taxid col",
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
        "-j",
        "--assembly",
        metavar="ASSEMBLY",
        default=None, required=False,
        help="Assembly refseq file",
    )
    parser.add_argument(
        "-k",  "--compress_species", default=False,  help="Compress species to species level",  action='store_true'
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

def lorenz_curve(depths):
    """Compute the Lorenz curve for a list of depths."""
    sorted_depths = sorted(depths)
    cumulative_depths = [0]
    total = sum(sorted_depths)
    for depth in sorted_depths:
        cumulative_depths.append(cumulative_depths[-1] + depth)
    lorenz_curve = [x / total if x >0 else 0 for x in cumulative_depths]
    return lorenz_curve


def gini_coefficient(depths):
    """Calculate the Gini coefficient for a list of depths."""
    lorenz = lorenz_curve(depths)
    n = len(lorenz)
    area_under_lorenz = sum((lorenz[i] + lorenz[i + 1]) / 2 for i in range(n - 1)) / (n - 1)
    gini = 1 - 2 * area_under_lorenz
    return gini

def breadth_of_coverage(depths, genome_length):
    """Calculate the breadth of coverage for a list of depths."""
    non_zero_positions = len([depth for depth in depths if depth > 0])
    return non_zero_positions / genome_length if genome_length > 0 else 0

def adjusted_fair_distribution_score(depths, genome_length):
    """Calculate the adjusted score for fair distribution of depths."""
    gini = gini_coefficient(depths)
    breadth = breadth_of_coverage(depths, genome_length)

    total_depth = sum(depths)
    if total_depth == 0:
        return 0  # If there's no coverage, we consider it the lowest score

    # Compute penalty based on the distribution of depths
    max_depth = max(depths)
    avg_depth = sum(depths) / len(depths) if len(depths) > 0 else 0

    # Define a penalty factor that reduces the impact of high depths
    penalty_factor = 0.5  # You can adjust this factor based on desired sensitivity

    if avg_depth > 0:
        penalty = penalty_factor * (max_depth / avg_depth - 1) / (max_depth / avg_depth + 1)
    else:
        penalty = penalty_factor * max_depth

    # Adjust the score: 1 - Gini coefficient, penalized by the penalty factor
    gini_score = (1 - gini) * (1 - penalty)

    # Ensure both scores are between 0 and 1
    gini_score = max(0, min(1, gini_score))
    breadth = max(0, min(1, breadth))
    # log transform breadth to 0 and 1, more weight closer to 1
    if breadth > 0:
        breadth = log2(breadth + 1) / log2(2)
    # Combine the Gini score and breadth of coverage with equal weight
    final_score = 0.1 * (1-gini_score) + 0.9 * breadth
    return final_score


def get_fair_distribution_score(data):
    # Assuming 'data' is a dictionary with 'depths' as a list and 'total_length' as an int
    depths = data.get('depths', [])
    genome_length = data.get('total_length', len(depths))  # Use provided total length or length of depths

    # Calculate the adjusted fair distribution score
    fair_distribution_score = adjusted_fair_distribution_score(depths, genome_length)

    return fair_distribution_score


def import_pathogens(pathogens_file):
    """Import the pathogens from the input CSV file, correctly handling commas in quoted fields."""
    pathogens_dict = {}
    # Open the file using the `with` statement
    with open(pathogens_file, 'r', newline='', encoding='ISO-8859-1') as file:
        # Create a CSV reader object that handles commas inside quotes automatically
        reader = csv.reader(file, delimiter=',', quotechar='"')

        # Iterate over each row in the CSV file
        for row in reader:
            # Assign each part of the row to variables if available
            pathogen_name = row[0] if len(row) > 0 else None
            taxid = row[1] if len(row) > 1 else None
            call_class = row[2] if len(row) > 2 else None

            pathogenic_sites = row[3] if len(row) > 3 else None
            commensal_sites = row[4] if len(row) > 4 else None
            status = row[5] if len(row) > 5 else None
            pathology = row[6] if len(row) > 6 else None

            # Store the data in the dictionary, keyed by pathogen name
            pathogens_dict[pathogen_name] = {
                'taxid': taxid,
                'callclass': call_class,
                'pathogenic_sites': pathogenic_sites,
                'name': pathogen_name,
                'commensal_sites': commensal_sites,
                'status': status,
                'pathology': pathology
            }
            pathogens_dict[taxid] = pathogens_dict[pathogen_name]


    # No need to explicitly close the file, `with` statement handles it.
    return pathogens_dict

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
    """Calculate the Shannon entropy of the given values without NumPy."""
    total = sum(values)
    probabilities = [value / total for value in values]
    return -sum(p * log2(p) for p in probabilities if p > 0)

def calculate_mean(diamond_list, key):
    """
    Calculate the weighted mean of a key, weighted by 'cds' values.

    Parameters:
    diamond_list (list of dict): List of dictionaries containing the data.
    key (str): The key for which the weighted mean is to be calculated.

    Returns:
    float: The weighted mean of the values.
    """
    total_weighted_value = 0
    total_cds = 0

    # Loop through each item in the list and calculate the weighted value
    for item in diamond_list:
        value = float(item[key])  # Convert the key value to a float
        cds = int(item['cds'])    # Get the cds value (as the weight)

        # Add to the weighted sum
        total_weighted_value += value * cds
        total_cds += cds

    # Return the weighted mean
    return total_weighted_value / total_cds if total_cds != 0 else 0
# Function to calculate disparity based on the sum of numreads and k2_reads
def calculate_disparity_siblings(numreads, k2_reads):
    sumreads = sum(numreads)

    if sumreads + k2_reads > 0:
        # Calculate the harmonic mean
        harmonic_mean = 2 * (sumreads * k2_reads) / (sumreads + k2_reads)

        # Calculate the arithmetic mean for normalization
        arithmetic_mean = (sumreads + k2_reads) / 2

        # Normalize the harmonic mean by dividing by the arithmetic mean
        normalized_harmonic_mean = harmonic_mean / arithmetic_mean
    else:
        normalized_harmonic_mean = 0  # No disparity if both reads are 0

    return normalized_harmonic_mean




# Function to calculate weighted mean
def calculate_weighted_mean(data, numreads):
    if len(numreads) > 0:
        total_weight = sum(numreads)  # Sum of all weights (numreads)
        weighted_sum = []
        for x in range(0, len(numreads)):
            if x >= len(data):
                break
            if numreads[x] > 0:
                weighted_sum.append(data[x] * numreads[x]) # Sum of weight*value
            else:
                weighted_sum.append(0)
        weighted_sum = sum(weighted_sum)
        if total_weight > 0:
            weighted_mean = weighted_sum / total_weight
        else:
            return 0
    else:
        weighted_mean = 0
    return weighted_mean

def import_k2_file(filename):
    """Import the Kraken2 output file with parent-child relationships based on name indentation"""

    # Open the Kraken2 file
    tsv_file = open(filename, newline='')
    read_tsv = csv.reader(tsv_file, delimiter="\t")

    # Regex to capture the indentation in the name
    k2_regex = re.compile(r"^(\s+)(.+)")

    # Initialize variables
    taxids = dict()
    depth = dict()  # To store the depth (level of indentation) of each taxid
    lastparents = dict()  # To track the last parent at each depth level

    # Define header columns
    header = ['abundance', 'clade_fragments_covered', 'number_fragments_assigned', 'rank', 'taxid', 'name', 'parents']

    # Iterate through each row in the Kraken2 file
    for row in read_tsv:
        entry = dict()
        for x in range(0, len(header)):
            if (header[x] != 'name' and header[x] != 'rank' and header[x] != 'taxid' and x < len(row)):
                entry[header[x]] = float(row[x])
            elif x < len(row):
                if header[x] == 'taxid':
                    entry[header[x]] = str(row[x])
                else:
                    entry[header[x]] = row[x]
            else:
                entry[header[x]] = ""

        # Use regex to adjust entry['name'] and determine depth (indentation level)
        match = k2_regex.match(row[5])  # row[5] corresponds to the 'name' field
        if match:
            indent = len(match.group(1)) // 2  # Determine depth based on the number of spaces (2 spaces per level)
            entry['name'] = match.group(2)  # Remove indentation from the name
        else:
            indent = 0  # If no match, it is at the top level

        # Set depth for the current entry
        depth[entry['taxid']] = indent
        entry['depth'] = indent

        # Determine parents for the current entry (2D list of [taxid, rank])
        parents = []
        for i in range(indent - 1, 0, -1):
            parents.append([lastparents[i], taxids.get(lastparents[i], {}).get('rank', '')])

        # Store the parent taxid at the current depth
        lastparents[indent] = entry['taxid']

        # Add parents to the entry
        entry['parents'] = parents

        # Add taxid information to the taxids dictionary
        taxids[str(entry['taxid'])] = dict(
            name=entry['name'],
            abundance=entry['abundance'],
            taxid=str(entry['taxid']),
            number_fragments_assigned=entry['number_fragments_assigned'],
            rank=entry['rank'],
            clades_covered=entry['clade_fragments_covered'],
            depth=entry['depth'],
            parents=entry['parents']
        )

    # Return the mapping with parent-child relationships
    return taxids



def detect_regions_from_depth(reference_coverage, depthfile, avg_read_length):
    """
    Detect regions based on gaps in depth and the average read length.

    Args:
    - reference_coverage (dict): Dictionary holding reference coverage information.
    - depthfile (str): Path to the depth file.
    - avg_read_length (int): The average read length for aligning reads.

    Returns:
    - reference_coverage (dict): Updated reference coverage with detected regions.
    """

    # Read the depth information from the depth file
    print("Reading depth information from depth file")
    with open(depthfile, 'r') as f:
        for line in f:
            splitline = line.split('\t')
            reference_name = splitline[0]
            pos = int(splitline[1])
            depth = int(splitline[2])
            reference_coverage[reference_name]['depths'][pos] = depth

    # Process each reference to detect regions
    for ref, data in reference_coverage.items():
        depths = data['depths']
        regions = []
        current_region_start = None
        previous_depth = None

        # Traverse through all positions
        for pos in range(1, data['length'] + 1):
            depth = depths.get(pos, 0)  # Default to 0 if no depth recorded

            # Check for the start of a new region
            if current_region_start is None and depth >= 1:
                current_region_start = pos

            # Check for gaps or the end of a region
            if depth == 0 and current_region_start is not None:
                regions.append((current_region_start, pos - 1))
                current_region_start = None

            # Update the previous depth
            previous_depth = depth

        # Handle the last region (if there was no gap at the end)
        if current_region_start is not None:
            regions.append((current_region_start, data['length']))

        # Adjust regions based on the average read length
        adjusted_regions = []
        for region_start, region_end in regions:
            region_length = region_end - region_start + 1
            if region_length > avg_read_length:
                # Split region into multiple parts if it spans more than the average read length
                num_subregions = max(1, region_length // avg_read_length)
                subregion_size = region_length // num_subregions
                for i in range(num_subregions):
                    subregion_start = region_start + i * subregion_size
                    subregion_end = min(region_start + (i + 1) * subregion_size - 1, region_end)
                    adjusted_regions.append((subregion_start, subregion_end))
            else:
                # Keep the region as is
                adjusted_regions.append((region_start, region_end))

        # Update the reference coverage with the detected regions
        reference_coverage[ref]['covered_regions'] = adjusted_regions

    return reference_coverage



def count_reference_hits(bam_file_path, depthfile, covfile, matchdct, min_reads_align):
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
        # get total reads


        for ref in bam_file.header.references:
            # get average read length
            reference_lengths[ref] = bam_file.get_reference_length(ref)
            reference_coverage[ref] = dict(
                length = reference_lengths[ref],
                depths =  defaultdict(int),
                gini_coefficient = 0,
                mapqs = [],
                baseqs = [],
                meanmapq = 0,
                numreads = 0,
                meanbaseq = 0,
                coverage = 0,
                meandepth = 0,
                name = ref,
                accession = ref,
                isSpecies = False,
                covered_regions = [],  # Store regions covered by reads
            )
        # Open BAM file
        unique_read_ids = set()  # To track unique read IDs
        total_length = 0
        total_reads = 0
        readlengths = []

        for read in bam_file.fetch():
            read_id = read.query_name

            # Check if the read is already processed (to handle multiple alignments)
            if read_id not in unique_read_ids:
                unique_read_ids.add(read_id)
                total_reads += 1
                readlengths.append(read.query_length)
                total_length += read.query_length  # Add the length of the read
            if not read.is_unmapped:  # Check if the read is aligned
                aligned_reads += 1
        # Calculate average read length
        average_read_length = Math.ceil(total_length / total_reads if total_reads > 0 else 0)

        print(f"Total unique reads: {total_reads}")
        print(f"Average read length: {average_read_length}")

        if depthfile:
            # Detect regions based on the depth file and average read length
            reference_coverage = detect_regions_from_depth(reference_coverage, depthfile, average_read_length)

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
            for read in bam_file.fetch():
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
                        baseq = sum(base_qualities) / len(base_qualities)
                        # get the mapq of the read
                        mapq = float(read.mapping_quality)
                        reference_coverage[reference_name]['mapqs'].append(mapq)
                        reference_coverage[reference_name]['baseqs'].append(baseq)
                        reference_coverage[reference_name]['numreads']+=1
                    ### This needs updating.
                    if not depthfile:
                        for pos in range(read.reference_start, read.reference_end):
                            reference_coverage[reference_name]['depths'][pos] += 1
                    reference_coverage[reference_name]['covered_regions'].append((read.reference_start, read.reference_end))
                else:
                    unaligned += 1

        if not covfile:
            # make meanbaseq and meanmapq
            for ref, data in reference_coverage.items():
                if data['numreads'] > 0:
                    data['meanbaseq'] = sum(data['baseqs']) / len(data['baseqs'])
                    data['meanmapq'] = sum(data['mapqs']) / len(data['mapqs'])

    bam_file.close()
    i=0
    return reference_coverage, aligned_reads, total_reads
def main():
    args = parse_args()
    inputfile = args.input
    pathogenfile = args.pathogens
    covfile = args.coverage
    """
    # Final Score Calculation

    The final score is calculated as the average of four tests, each contributing 25% (or 0.25 proportion) to the overall score. The four tests are based on the following criteria:

    1. **Gini Coefficient (25% of final score)**:
        - The Gini coefficient is a measure of inequality and is expected to be between 0 and 1. It directly contributes to 25% of the final score.
        - Expected Range: [0, 1]

    2. **MAPQ Score (25% of final score)**:
        - The MAPQ (Mapping Quality) score is derived from the Phred score using the formula:
            \[
            \text{MAPQ} = 10 \times \log_{10}(\text{Phred score})
            \]
        - This score represents the average quality of the alignment and contributes 25% to the final score.
        - Handle cases where the Phred score is 0 to avoid logarithmic errors.
        - Expected Range: [0, ~60] depending on the Phred score.

    3. **Disparity between k2 reads (25% of final score)**:
        - The disparity represents the difference between two sets of reads. This value is provided or calculated externally.
        - Expected Range: [0, 1] (assuming the disparity is normalized for the calculation).

    4. **CD Limit (25% of final score)**:
        - The CD (Cutoff Determination) limit is a binary value (0 or 1) indicating whether a certain threshold or limit has been reached.
        - It directly contributes to 25% of the final score.
        - Expected Range: {0, 1} (binary)

    # Formula to Calculate Final Score:

        final_score = gini_coeff * gini_coeff_weight \
            + mapq_score * mapq_score_weight \
            + disparity * disparity_weight \
            + CD_limit * CD_limit_weight

    Each component contributes equally (25%) to the final score.
    """
    accumulated_scores = defaultdict(dict)
    output = args.output
    matcher = args.match
    k2_mapping = dict()
    matchdct = dict()
    header = True
    i =0
    capval  = args.capval
    mincov = args.mincoverage
    if args.k2:
        k2_mapping = import_k2_file(args.k2)
    dmnd = defaultdict()
    if args.diamond:
        # read in diamond file as dataframe
        with open(args.diamond, 'r') as f:
            # use csv reader
            reader = csv.reader(f, delimiter="\t")
            next(reader)
            for line in reader:
                if len(line) > 0:
                    try:
                        dmnd[line[0]] = dict(
                            contigs= int(line[1]),
                            cds = int(line[2]),
                            identity= float(line[3])/100,
                            lengthmedian = line[4],
                            mismatchedmedian= line[5],
                            medianevalue= line[6],
                        )
                    except Exception as e:
                        print(f"Error: {e}")
                i+=1
    reference_hits, aligned_total, total_reads = count_reference_hits(
        inputfile,
        args.depth,
        covfile,
        matchdct,
        args.min_reads_align
    )

    assembly_to_accession = defaultdict(set)

    if args.match:
        # open the match file and import the match file
        header = True
        i=0
        with open (matcher, 'r') as f:
            accindex = args.accessioncol
            nameindex = args.namecol
            taxcol = args.taxcol
            for line in f:
                line = line.strip()
                if header and i==0:
                    header = False
                    continue
                taxid=None
                assembly = None
                splitline = line.split('\t')
                if len(splitline) > 0:
                    accession = splitline[accindex]
                else:
                    accession = None
                if len(splitline) > 0:
                    assembly = splitline[1]
                if len(splitline) > taxcol:
                    taxid = splitline[taxcol]
                if len(splitline) > nameindex:
                    name = splitline[nameindex]
                else:
                    name = None
                if accession in reference_hits:
                    reference_hits[accession]['taxid'] = taxid
                    reference_hits[accession]['name'] = name
                    reference_hits[accession]['assembly'] = assembly
                if accession not in assembly_to_accession[assembly]:
                    assembly_to_accession[assembly].add(accession)
                i+=1

        f.close()

    if args.assembly:
        i=0
        with open(args.assembly, 'r') as f:
            for line in f:
                line = line.strip()
                i+=1
                if i<2:
                    continue
                splitline = line.split('\t')
                if len(splitline) > 0:
                    accession = splitline[0]
                    taxid = splitline[5]
                    species_taxid = splitline[6]
                    name = splitline[7]
                    strain = splitline[8].replace("strain=", "")
                    if strain == "na":
                        strain = f"≡"

                    isSpecies = False if species_taxid != taxid else True
                    # fine value where assembly == accession from reference_hits
                    if accession in assembly_to_accession:
                        for acc in assembly_to_accession[accession]:
                            if acc in reference_hits:
                                reference_hits[acc]['isSpecies'] = isSpecies
                                reference_hits[acc]['toplevelkey'] = species_taxid
                                if strain == "na":
                                    strain = "≡"
                                reference_hits[acc]['strain'] = strain
                                reference_hits[acc]['assemblyname'] = name
                                reference_hits[acc]['name'] = name

        f.close()
    final_format = defaultdict(dict)


    if args.compress_species:
        # Need to convert reference_hits to only species species level in a new dictionary
        for key, value in reference_hits.items():
            value['accession'] = key
            if 'toplevelkey' in value and value['toplevelkey']:
                valtoplevel = value['toplevelkey']
            else:
                valtoplevel = key
            valkey = value['accession']
            if value['numreads'] > 0 and value['meandepth'] > 0:
                final_format[valtoplevel][valkey]= value
    else:
        # We don't aggregate, so do final format on the organism name only
        for key, value in reference_hits.items():
            if 'toplevelkey' in value and value['toplevelkey']:
                valtoplevel = value['toplevelkey']
            else:
                valtoplevel = key
            valkey = key
            if value['numreads'] > 0 and value['meandepth'] > 0:
                final_format[valtoplevel][valkey] = value
    # Dictionary to store aggregated species-level data
    species_aggregated = {}
    def getGiniCoeff(data, acc_length=0):
        # Assuming 'data' is a dictionary with 'depths' as a list and 'total_length' as an int.
        depths = list(data.values()) if data.values() else [0]
        # Get get proportion of depths in increments 10 or higher like 0-10 x 11-20 x etc


        # Gini Coefficient
        gini_coefficient = get_fair_distribution_score({"depths": depths, 'total_length': acc_length})
        return gini_coefficient

    # Step 2: Define a function to calculate disparity for each organism
    # Define a function to calculate disparity with softer variance influence
    # Step 2: Define a function to dynamically dampen variance based on the proportion of reads
    def calculate_disparity(numreads, total_reads, variance_reads, k=10):
        """
        Dynamically dampens the variance effect based on the proportion of reads.
        numreads: Total number of reads aligned to the organism (sum of reads)
        total_reads: Total number of reads aligned in the sample
        variance_reads: Variance of the aligned reads across all organisms
        k: Damping factor to control the influence of the proportion on the penalty
        """
        if total_reads == 0:
            return 0  # Avoid division by zero

        # Calculate the proportion of aligned reads
        proportion = numreads / total_reads

        # Dynamically adjust the variance penalty based on the proportion of reads
        dampened_variance = variance_reads / (1 + k * proportion)

        # Calculate disparity based on the proportion and the dynamically dampened variance
        disparity = proportion * (1 + dampened_variance)

        return disparity

    # Aggregate data at the species level
    for top_level_key, entries in final_format.items():

        # all_assemblies = [[x['assemblyname'], x['meanmapq'], x['numreads'], x['taxid']] for x in entries.values()]
        for val_key, data in entries.items():
            # get all organisms that are NOT the same assemblyname from all_assemblies
            # filtered_assemblies = [x for x in all_assemblies if x[0] != data['assemblyname']]
            # if len(filtered_assemblies) > 0:
            #     print(data['assemblyname'], data['meanmapq'], data['taxid'], data['numreads'], filtered_assemblies)
            if top_level_key not in species_aggregated:
                species_aggregated[top_level_key] = {
                    'key': top_level_key,
                    'numreads': [],
                    'mapqs': [],
                    'depths': [],
                    "taxids": [],
                    "accs": [],
                    'assemblies': [],
                    "coverages": [],
                    "prevalence_disparity": 0,
                    "coeffs": [],
                    'baseqs': [],
                    "isSpecies": True if  args.compress_species  else data['isSpecies'],
                    'strainslist': [],
                    'covered_regions': 0,
                    'name': data['name'],  # Assuming the species name is the same for all strains
            }

            try:
                gini_strain = getGiniCoeff(data['depths'], data['length'])
                species_aggregated[top_level_key]['coeffs'].append(gini_strain)
                species_aggregated[top_level_key]['taxids'].append(data['taxid'])
                species_aggregated[top_level_key]['numreads'].append(data['numreads'])
                species_aggregated[top_level_key]['covered_regions'] += len(data['covered_regions'])
                species_aggregated[top_level_key]['coverages'].append(data['coverage'])
                species_aggregated[top_level_key]['baseqs'].append(data['meanbaseq'])
                species_aggregated[top_level_key]['accs'].append(data['accession'])
                species_aggregated[top_level_key]['mapqs'].append(data['meanmapq'])
                species_aggregated[top_level_key]['depths'].append(data['meandepth'])
                name = data['name']
                if 'strain' in data:
                    species_aggregated[top_level_key]['strainslist'].append({
                        "strainname": data['strain'],
                        "fullname":data['name'],
                        "subkey": val_key,
                        "numreads": data['numreads'],
                        "taxid": data['taxid'] if "taxid" in data else None,
                    })
                else:
                    species_aggregated[top_level_key]['strainslist'].append({
                        "strainname":val_key,
                        "fullname":val_key,
                        "subkey": val_key,
                        "numreads": data['numreads'],
                        "taxid": data['taxid'] if "taxid" in data else None,
                    })
            except Exception as e:
                print(top_level_key, val_key,"___")
                print(f"Error: {e}")
    all_readscounts = [sum(x['numreads']) for x in species_aggregated.values()]
    def calculate_var(read_counts):
        """
        Manually calculate variance for a list of read counts.
        read_counts: A list of aligned reads for each organism
        """
        n = len(read_counts)
        if n == 0:
            return 0  # Avoid division by zero if no organisms

        # Step 1: Calculate the mean of the reads
        mean_reads = sum(read_counts) / n

        # Step 2: Calculate the squared differences
        squared_diffs = [(x - mean_reads) ** 2 for x in read_counts]

        # Step 3: Calculate variance
        variance = sum(squared_diffs) / n
        return variance

    variance_reads = calculate_var(all_readscounts)
    print(f"Variance of reads: {variance_reads}")
    for top_level_key, aggregated_data in species_aggregated.items():
        numreads = aggregated_data['numreads']
        aggregated_data['meanmapq'] = calculate_weighted_mean(aggregated_data['mapqs'], numreads)
        aggregated_data['meanbaseq'] = calculate_weighted_mean(aggregated_data['baseqs'],numreads)
        aggregated_data['meandepth'] = calculate_weighted_mean(aggregated_data['depths'],numreads)
        aggregated_data['meancoverage'] = calculate_weighted_mean(aggregated_data['coverages'],numreads)
        aggregated_data['meangini'] = calculate_weighted_mean(aggregated_data['coeffs'],numreads)
        # Step 4: Calculate the disparity for this organism
        aggregated_data['disparity'] = calculate_disparity(sum(numreads), total_reads, variance_reads)

        # if taxid is in taxid_mapping, use that, otherwise use the strainname
        k2_reads = 0
        if args.k2:
            taxid = aggregated_data.get('key', None)
            name = aggregated_data.get('name', None)
            if taxid and taxid in k2_mapping:
                k2_reads = k2_mapping[taxid].get('clades_covered', 0)
            elif name and name in k2_mapping:
                k2_reads = k2_mapping[name].get('clades_covered', 0)
            else:
                # get all taxids from strainslist
                taxids = [x.get('taxid', None) for x in aggregated_data['strainslist']]
                # remove None from taxids
                taxids = [x for x in taxids if x]
                k2_reads = 0
                for taxid in taxids:
                    if taxid in k2_mapping:
                        k2_reads += k2_mapping[taxid].get('clades_covered', 0)
        aggregated_data['k2_numreads'] = k2_reads

    # Step 4: Find the min and max disparity values
    disparities = [aggregated_data['disparity'] for aggregated_data in species_aggregated.values()]
    min_disparity = min(disparities)
    max_disparity = max(disparities)

    for top_level_key, aggregated_data in species_aggregated.items():
        disparity = aggregated_data['disparity']
        if max_disparity > min_disparity:
            normalized_disparity = (disparity - min_disparity) / (max_disparity - min_disparity)
        else:
            normalized_disparity = 0  # If all disparities are the same, set to 0

        # Store the normalized disparity
        aggregated_data['normalized_disparity'] = normalized_disparity
        print(f"Entry Top Key: {top_level_key}")
        print(f"\tName: {aggregated_data['name']}")
        print(f"\tNum Reads: {aggregated_data['numreads']}")
        print(f"\tK2 Reads: {aggregated_data['k2_numreads']}")
        print(f"\tPrev. Disparity: {aggregated_data['disparity']}")
        print(f"\t^Norm. Disparity: {aggregated_data['normalized_disparity']}")
    # exit()

    # Function to normalize the MAPQ score to 0-1 based on a maximum MAPQ value
    def normalize_mapq(mapq_score, max_mapq=60):
        # Normalize the MAPQ score to be between 0 and 1
        probability = 10 ** (-mapq_score / 10)
        return 1-probability  # Ensure values stay between 0 and 1

    def calculate_disparity_cv(numreads):
        if len(numreads) > 1:
            # Ensure the values are positive
            numreads = [max(0, x) for x in numreads]

            mean_reads = statistics.mean(numreads)
            stddev_reads = statistics.stdev(numreads)

            # Calculate CV
            if mean_reads != 0:
                cv = stddev_reads / mean_reads
            else:
                cv = 0

            # Clamp the CV value to be between 0 and 1
            return max(0, min(cv, 1))

        return 0
    def get_species_by_parent_rank(species_aggregated, parent_taxid, my_rank, rank_mapping_match):
        """
        Get all species that have the specified rank (e.g., 'S' for species) and belong to a parent with the
        specified taxonomic rank (e.g., 'F' for family).

        Parameters:
            species_aggregated (dict): Dictionary of species.
            parent_taxid (str): The parent taxid to match (e.g., for family 'F').
            my_rank (str): The rank to filter species by (e.g., 'S').
            rank_mapping_match (str): The taxonomic rank to check for the parent (e.g., 'F').

        Returns:
            List of species that match the specified rank and belong to the parent with the given taxonomic rank.
        """
        sibling_species = []

        # Loop through all species in species_aggregated
        for species in species_aggregated.values():
            parents = species.get('parents', [])
            rank = species.get('rank', '')

            # Check if the species rank matches `my_rank` (e.g., 'S' for species)
            if rank == my_rank:
                # Look through the parents to find the one that matches the `rank_mapping_match` (e.g., 'F')
                for parent in parents:
                    if parent[1] == rank_mapping_match and parent[0] == parent_taxid:
                        sibling_species.append(species)
                        break  # Exit loop after finding the match

        return sibling_species

    for key, value in species_aggregated.items():
        ## Test: Disaprity across genus and its species

        if args.k2 and key in k2_mapping:
            # Define the rank code you're looking for, for example "G" for Genus
            rank_mapping_match = None
            if args.parent_k2_match:
                rank_mapping_match = args.parent_k2_match

            # Step 1: Get the 'parents' list for the current key from k2_mapping
            parents = k2_mapping.get(key, {}).get('parents', [])
            my_rank = k2_mapping.get(key, {}).get('rank', 'U')

            # Step 2: Find the parent taxid where the rank matches rank_mapping_match (e.g., "F" for Family)
            parent_taxid = None
            if rank_mapping_match:
                for parent in parents:
                    if parent[1] == rank_mapping_match:
                        parent_taxid = parent[0]
                        break
            else:
                parent_taxid = parents[0][0]  # Use the first parent if no rank_mapping_match is provided

            # If no parent_taxid was found, skip the rest of the steps for this key
            if not parent_taxid and parent_taxid != 0:
                continue

            # Step 3: Get all species ('S') that belong to the same family ('F') using the found parent taxid
            related_species = get_species_by_parent_rank(k2_mapping, parent_taxid, my_rank, rank_mapping_match)

            # Step 4: Aggregate k2_numreads for all species in the same genus
            sibling_k2_reads = [[species.get('clades_covered', 0), species.get('taxid')] for species in related_species]
            # Save the aggregated k2 reads as 'parent_k2_reads'
            value['parent_k2_reads'] = sum([species.get('clades_covered', 0) for species in related_species])
            value['sibling_k2_reads'] = sibling_k2_reads

            # Step 5: Run disparity calculations between the species-specific numreads and genus_k2_reads
            # Get the numreads from k2_mapping for this key
            numreads = value.get('k2_numreads', 0)

            # Disparity calculations
            value['disparity_cv'] = 1-calculate_disparity_cv([x[0] for x in sibling_k2_reads])
            # Optionally calculate disparity_ratio or disparity_harmonic
            # k2_mapping[key]['disparity_ratio'] = calculate_disparity_ratio(numreads, genus_k2_reads)
            # k2_mapping[key]['disparity_harmonic'] = calculate_normalized_harmonic_mean(numreads, genus_k2_reads)

        ## Test: Disparity of k2 or other classifier to alignment stats of NTs
        if value.get('k2_numreads', None):
            value['raw_disparity'] = abs(calculate_disparity_siblings(value['numreads'], value['k2_numreads']))
        else:
            value['raw_disparity'] = 0
        ## Test: Mapq of NT,
        normalized_mapq = normalize_mapq(value.get('meanmapq', 0))
        value['alignment_score'] = normalized_mapq


        ## Test: Min threshold of CDs found  - are there coding regions at the min amount?
        species_aggregated[key]['diamond'] = {'identity': 0, 'maxvalereached': False}
        if args.diamond:

            if key in dmnd:
                dmnd[key]['maxvalereached'] = dmnd[key].get('cds', 0) > args.min_cds_found
                species_aggregated[key]['diamond'] = dmnd[key]
            elif species_aggregated[key].get('strainslist', None):
                # get all taxids in strainslist
                taxids = set([x.get('taxid', None) for x in value['strainslist']])
                # check if any of the taxids are in diamond
                allstrains = []
                for taxid in taxids:
                    if taxid in dmnd:
                        allstrains.append(dmnd[taxid])
                if len(allstrains) > 0:
                    # aggregate all as medians
                    # Calculate medians for 'identity', 'lengthmedian', 'mismatchedmedian', and 'medianevalue'
                    mean_identity = calculate_mean(allstrains, 'identity')
                    mean_length = calculate_mean(allstrains, 'lengthmedian')
                    mean_mismatched = calculate_mean(allstrains, 'mismatchedmedian')
                    mean_evalue = calculate_mean(allstrains, 'medianevalue')
                    mean_contigs = calculate_mean(allstrains, 'contigs')
                    mean_cds = calculate_mean(allstrains, 'cds')
                    max_val_reached =  mean_cds > args.min_cds_found
                    dmnd_results = {
                        'identity': mean_identity,
                        'lengthmean': mean_length,
                        'mismatchedmean': mean_mismatched,
                        'meanevalue': mean_evalue,
                        'contigs': mean_contigs,
                        'cds': mean_cds,
                        'maxvalereached': max_val_reached
                    }

                    species_aggregated[key]['diamond'] = dmnd_results
            else:
                print(f"Key {key} not found in diamond file")


    #  # Print the final aggregated data
    # for top_level_key, aggregated_data in species_aggregated.items():
    #     print(f"Entry Top Key: {top_level_key}")
    #     print(f"\tNum Reads: {sum(aggregated_data['numreads'])}, {aggregated_data['numreads']}")
    #     print(f"\tMean MapQ: {aggregated_data.get('meanmapq', 0)}")
    #     print(f"\tMean BaseQ: {aggregated_data.get('meanbaseq', 0)}")
    #     print(f"\tStrains List: {[x.get('fullname', '') for x in aggregated_data['strainslist']]}")
    #     print(f"\tK2 Reads: {aggregated_data.get('k2_numreads', 0)}")
    #     print(f"\tParent K2 Reads: {aggregated_data.get('parent_k2_reads', 0)}")
    #     print(f"\tK2 Disparity: {aggregated_data.get('raw_disparity',0)}")
    #     print(f"\tK2 Parent-Level Disparity: {aggregated_data.get('disparity_cv', 0)}")
    #     print(f"\tPrev. Disparity: {aggregated_data.get('disparity', 0)}")
    #     print(f"\tNorma. Prev. Disparity: {aggregated_data.get('normalized_disparity', 0)}")
    #     print(f"\tDiamond Identity: {aggregated_data.get('diamond', {}).get('identity',0)}")
    #     print(f"\tScore MapQ: {aggregated_data.get('alignment_score', 0)}")
    #     print(f"\tGini Score: {aggregated_data.get('meangini', 0)}")
    #     print()

    pathogens = import_pathogens(pathogenfile)
    # Next go through the BAM file (inputfile) and see what pathogens match to the reference, use biopython
    # to do this

    if args.min_reads_align:
        # filter the reference_hits based on the minimum number of reads aligned
        print(f"Filtering for minimum reads aligned: {args.min_reads_align}")
        species_aggregated = {k: v for k, v in species_aggregated.items() if ((v['covered_regions'])) >= int(args.min_reads_align)}
        species_aggregated = {k: v for k, v in species_aggregated.items() if sum(v['numreads']) >= int(args.min_reads_align)}
    # Create a new dictionary to store aggregated species-level data
    for top_level_key, data in species_aggregated.items():
        print(f"Reference: {data['name']} (Taxid: {data['key']})")
        print(f"\tNumber of Reads: {data['numreads']}")
        print(f"\tMean Coverage: {data['meancoverage']}")
        print(f"\tAlignment Conf: {data['meangini']}")
        print(f"\tDepth of Coverage: {data['meandepth']}")
        print(f"\tMean BaseQ: {data['meanbaseq']}")
        print(f"\tMean MapQ: {data['meanmapq']}")
        print(f"\tisSpecies: {data['isSpecies']}")
        print(f"\tPathogenic Strains: {len(data['strainslist'])}")
        print(f"\tRegions: {data['covered_regions']}")
        print()
    write_to_tsv(
        aggregated_stats=species_aggregated,
        pathogens=pathogens,
        output_file_path=output,
        sample_name=args.samplename,
        sample_type = args.sampletype,
        total_reads = total_reads,
        aligned_total = aligned_total
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



def write_to_tsv(aggregated_stats, pathogens, output_file_path, sample_name="No_Name", sample_type="Unknown", total_reads=0, aligned_total = 0):
    """
    Write reference hits and pathogen information to a TSV file.

    Args:
    reference_hits (dict): Dictionary with reference names as keys and counts as values.
    pathogens (dict): Dictionary with reference names as keys and dictionaries of additional attributes as values.
    output_file_path (str): Path to the output TSV file.
    """
    with open(output_file_path, 'w') as file:
        # Write the header row

        header = "Detected Organism\tSpecimen ID\tSample Type\t% Reads\t% Aligned Reads\t# Reads Aligned\tIsAnnotated\tPathogenic Sites\tMicrobial Category\tTaxonomic ID #\tStatus\tGini Coefficient\tMean BaseQ\tMean MapQ\tMean Coverage\tMean Depth\tAnnClass\tisSpecies\tPathogenic Subsp/Strains\tK2 Reads\tParent K2 Reads\tMapQ Score\tDisparity Score\tProtein Identity Score\tSiblings score\tTASS Score\n"
        file.write(f"{header}")
        print("________________________________________")

        for ref, count in aggregated_stats.items():
            strainlist = count.get('strainslist', [])
            isSpecies = count.get('isSpecies', False)
            is_pathogen = "Unknown"
            callfamclass = ""
            derived_pathogen = False
            isPath = False
            annClass = "None"
            pathogenic_sites = []
            refpath = pathogens.get(ref)
            pathogenic_sites = refpath.get('pathogenic_sites', []) if refpath else []
            def pathogen_label(ref):
                is_pathogen = "Unknown"
                isPathi = False
                callclass = ref.get('callclass', "N/A")
                pathogenic_sites = ref.get('pathogenic_sites', [])
                commensal_sites = ref.get('commensal_sites', [])
                if sample_type in pathogenic_sites:
                    if callclass != "commensal":
                        is_pathogen = callclass.capitalize()
                        isPathi = True
                    else:
                        is_pathogen = "Potential"
                elif sample_type in commensal_sites:
                    is_pathogen = "Commensal"
                elif callclass and callclass != "":
                    is_pathogen = callclass.capitalize() if callclass else "Unknown"
                    isPathi = True
                return is_pathogen, isPathi

            formatname = count.get('name', "N/A")

            if refpath:
                is_pathogen, isPathi = pathogen_label(refpath)
                isPath = isPathi
                if isPathi:
                    annClass = "Direct"
                taxid = refpath.get(ref, count.get(ref, ""))
                is_annotated = "Yes"
                commsites = refpath.get('commensal_sites', [])
                status = refpath.get('status', "N/A")
                formatname = refpath.get('name', formatname)
                if is_pathogen == "Commensal":
                    callfamclass = "Commensal Listing"
            else:
                is_annotated = "No"
                taxid = count.get(ref, "")
                status = ""

            listpathogensstrains = []
            fullstrains = []

            if strainlist:
                pathogenic_reads = 0
                merged_strains = defaultdict(dict)

                for x in strainlist:
                    keyx = x.get('strainname', x.get('taxid', ""))
                    # strainanme = x.get('strainname', None)

                    # if formatname != strainanme and strainanme:
                    #     formatname = formatname.replace(strainanme, "")
                    # elif not strainanme or not formatname:
                    #     formatname = x.get('fullname', "Unnamed Organism")
                    if keyx in merged_strains:
                        merged_strains[keyx]['numreads'] += x.get('numreads', 0)
                        merged_strains[keyx]['subkeys'].append(x.get('subkey', ""))
                    else:
                        merged_strains[keyx] = x
                        merged_strains[keyx]['numreads'] = x.get('numreads', 0)
                        merged_strains[keyx]['subkeys'] = [x.get('subkey', "")]

                for xref, x in merged_strains.items():
                    pathstrain = None
                    if x.get('taxid'):
                        fullstrains.append(f"{x.get('strainname', 'N/A')} ({x.get('taxid', '')}: {x.get('numreads', 0)} reads)")
                    else:
                        fullstrains.append(f"{x.get('strainname', 'N/A')} ({x.get('numreads', 0)} reads)")

                    if x.get('taxid') in pathogens:
                        pathstrain = pathogens.get(x.get('taxid'))
                    elif x.get('fullname') in pathogens:
                        pathstrain = pathogens.get(x.get('fullname'))

                    if pathstrain:
                        taxx = x.get('taxid', "")
                        if sample_type in pathstrain.get('pathogenic_sites', []) or pathstrain.get('callclass') != "commensal":
                            pathogenic_reads += x.get('numreads', 0)
                            percentreads = f"{x.get('numreads', 0)*100/aligned_total:.1f}" if aligned_total > 0 and x.get('numreads', 0) > 0 else "0"
                            listpathogensstrains.append(f"{x.get('strainname', 'N/A')} ({percentreads}%)")

                if callfamclass == "" or len(listpathogensstrains) > 0:
                    callfamclass = f"{', '.join(listpathogensstrains)}" if listpathogensstrains else ""
                if (is_pathogen == "N/A" or is_pathogen == "Unknown" or is_pathogen == "Commensal" or is_pathogen=="Potential") and listpathogensstrains:
                    is_pathogen = "Primary"
                    annClass = "Derived"

            meanbaseq = format_non_zero_decimals(count.get('meanbaseq', 0))
            gini_coefficient = format_non_zero_decimals(count.get('meangini', 0))
            meanmapq = format_non_zero_decimals(count.get('meanmapq', 0))
            meancoverage = format_non_zero_decimals(count.get('meancoverage', 0))
            meandepth = format_non_zero_decimals(count.get('meandepth', 0))
            k2_reads = format_non_zero_decimals(count.get('k2_numreads', 0))
            k2_parent_reads = format_non_zero_decimals(count.get('parent_k2_reads', 0))
            k2_disparity_score = format_non_zero_decimals(count.get('disparity_cv', 0))
            disparity_score = format_non_zero_decimals(count.get('normalized_disparity', 0))
            mapq_score = format_non_zero_decimals(count.get('alignment_score', 0))
            diamond_identity = format_non_zero_decimals(count.get('diamond', {}).get('identity', 0))
            siblings_score = format_non_zero_decimals(count.get('raw_disparity', 0))



            # if total_reads_aligned is 0 then set percent_aligned to 0
            countreads = sum(count['numreads'])
            if aligned_total == 0:
                percent_aligned = 0
            else:
                percent_aligned = format_non_zero_decimals(100*countreads / aligned_total)
            if total_reads == 0:
                percent_total = 0
            else:
                percent_total = format_non_zero_decimals(100*countreads / total_reads)
            if len(pathogenic_sites) == 0:
                pathogenic_sites = ""
            # Write to file in a newline format
            weights = {
                'mapq_score': 0.05,
                'diamond_identity': 0.2,
                'disparity_score': 0.5,
                'gini_coefficient': 0.2,
                "k2_disparity": 0.05,
                'siblings_score': 0
            }
            # if sum of vals in weights isnt 1 then normalize to 1
            total_weight = sum(weights.values())
            if total_weight != 1:
                for key in weights:
                    weights[key] = weights[key] / total_weight

            # Function to apply weights and format the result (assuming `format_non_zero_decimals` is defined)
            def apply_weight(value, weight):
                try:
                    # Convert the value to float before multiplying by the weight
                    value = float(value)
                    return value * weight
                except (TypeError, ValueError):
                    # If value cannot be converted to a float, return 0 or handle as needed
                    return 0


            # Apply weights to the relevant scores
            tass_score = format_non_zero_decimals(sum([
                apply_weight(count.get('normalized_disparity', 0), weights.get('disparity_score',0)),
                apply_weight(count.get('alignment_score', 0), weights.get('mapq_score',0)),
                apply_weight(count.get('meangini', 0), weights.get('gini_coefficient',0)),
                apply_weight(count.get('diamond', {}).get('identity', 0), weights.get('diamond_identity', 0)),
                apply_weight(count.get('k2_disparity', 0), weights.get('k2_disparity',0))
            ]))
            if is_pathogen == "Primary" or is_pathogen=="Potential":
                print(f"Reference: {ref} - {formatname}")
                print(f"\tIsPathogen: {is_pathogen}")
                print(f"\tCallClass: {callfamclass}")
                print(f"\tPathogenic Strains: {listpathogensstrains}")
                percentreads = f"{100*pathogenic_reads/aligned_total:.1f}" if aligned_total > 0 and pathogenic_reads>0 else 0
                print(f"\tPathogenic SubStrain Reads: {pathogenic_reads} - {percentreads}%")
                print(f"\tAligned Strains: {fullstrains}")
                print(f"\tTotal reads: {sum(count['numreads'])}")
                print(f"\tGini Conf: {count.get('meangini', 0)}")
                print(f"\tAlignment Score: {count.get('alignment_score', 0)}")
                print(f"\tK2 Disparity Score: {count.get('disparity_cv', 0)}")
                print(f"\tDisparity Score: {count.get('normalized_disparity', 0)}")
                print(f"\tDiamond Identity: {count.get('diamond', {}).get('identity', 0)}")
                print(f"\tFinal Score: {tass_score}")
                print(f"\tCovered Regions: {count.get('covered_regions', 0)}")
                print(f"\tK2 Reads: {count['k2_numreads']}")
                print()

            file.write(
                f"{formatname}\t{sample_name}\t{sample_type}\t{percent_total}\t{percent_aligned}\t{countreads}\t"
                f"{is_annotated}\t{pathogenic_sites}\t{is_pathogen}\t{ref}\t{status}\t{gini_coefficient}\t"
                f"{meanbaseq}\t{meanmapq}\t{meancoverage}\t{meandepth}\t{annClass}\t{isSpecies}\t{callfamclass}\t"
                f"{k2_reads}\t{k2_parent_reads}\t{mapq_score}\t{disparity_score}\t{diamond_identity}\t"
                f"{siblings_score}\t{tass_score}\n"
            )
if __name__ == "__main__":
    sys.exit(main())
