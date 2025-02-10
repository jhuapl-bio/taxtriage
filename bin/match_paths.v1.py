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
import random

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
        "-b",
        "--bedgraph",
        metavar="BEDGRAPH",
        default=None,
        help="Depth BEDGRAPH rather than 1 based per positions",
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
        default=2,
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
        type=int,
        help="Index of the column in mapfile (if specified) to match to the reference accession. 0 index start",
    )
    parser.add_argument(
        "-q",
        "--taxcol",
        metavar="TAXCOL",
        type=int,
        default=4, required =False,
        help="Index of the column in mapfile (if specified) to add taxid col",
    )
    parser.add_argument(
        "-n",
        "--namecol",
        default=2,
        type=int,
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
        "--ignore_missing_inputs", default=False,  help="If K2 or Dimaond output is not provided, dont reduce confidence",  action='store_true'
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
    parser.add_argument(
        "--disparity_score_weight",
        metavar="DISPARITYSCOREWEIGHT",
        type=float,
        default=0.05,
        help="value of weight for disparity reads vs other organisms in final TASS Score",
    )
    parser.add_argument(
        "--k2_disparity_score_weight",
        metavar="DISPARITYSCOREWEIGHT",
        type=float,
        default=0.05,
        help="value of weight for disparity of k2 and alignment in final TASS Score",
    )
    parser.add_argument(
        "--diamond_identity_weight",
        metavar="DISPARITYSCOREWEIGHT",
        type=float,
        default=0.05,
        help="value of weight for disparity of diamond_identity in final TASS Score",
    )
    parser.add_argument(
        "--mapq_weight",
        metavar="MAPQWEIGHT",
        type=float,
        default=0.05,
        help="value of weight for disparity ofmapq in final TASS Score",
    )
    parser.add_argument(
        "--gt",
        metavar="GTFILE",
        type=str,
        default=None,
        help="You want to provide a gt file to see how well you fare with the detections",
    )
    parser.add_argument(
        "--optimize",
        default=False,
        help="Optimize the gt file to get the best possible score",  action='store_true'
    )
    parser.add_argument("--max_iterations",
        type=int,
        default=20,
        help="Maximum number of random tweaks to explore if you provide a gt file and want to optimize"
    )
    parser.add_argument(
        '--breadth_weight',
        metavar="BREADTHSCORE",
        type=float,
        default=0.15,
        help="value of weight for breadth of coverage in final TASS Score",
    )
    parser.add_argument(
        "--gini_weight",
        metavar="GINIWEIGHT",
        type=float,
        default=0.85,
        help="value of weight for gini coefficient in final TASS Score",
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



def transform_func(d):
        return math.log2(1 + d)



def read_ground_truth(gt_file):
    """
    Reads a ground-truth file of the form:
        <accession> <coverage>
    Returns a dict: { accession: float(coverage), ... }
    """
    coverage_dict = {}
    with open(gt_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            acc = parts[0]
            cov = float(parts[1])
            coverage_dict[acc] = cov
    return coverage_dict

def build_transformed_coverage_hist(regions, genome_length):
    """
    Build a histogram of 'transformed' coverage for the entire genome.
    transform_func(depth) should return the transformed coverage value.
    """
    from collections import defaultdict
    coverage_hist = defaultdict(int)
    total_covered_bases = 0
    for (start, end, depth) in regions:
        length = end - start + 1
        transformed_depth = transform_func(depth)
        coverage_hist[transformed_depth] += length
        total_covered_bases += length

    # Add zeros (transformed) for uncovered portion
    uncovered = genome_length - total_covered_bases
    if uncovered > 0:
        coverage_hist[transform_func(0)] += uncovered

    return coverage_hist


def build_coverage_hist(regions, genome_length):
    """
    Build a histogram of coverage depth for the entire genome.
    coverage_hist[d] = number of bases with coverage d
    """
    coverage_hist = defaultdict(int)
    total_covered_bases = 0

    for (start, end, depth) in regions:
        length = end - start + 1   # how many bases in this interval
        coverage_hist[depth] += length
        total_covered_bases += length

    # All remaining bases have coverage = 0
    uncovered_bases = genome_length - total_covered_bases
    if uncovered_bases > 0:
        coverage_hist[0] += uncovered_bases

    return coverage_hist

def gini_coefficient_from_hist(coverage_hist):
    """
    Calculate the Gini coefficient from a histogram mapping coverage->number_of_bases.
    """
    if not coverage_hist:
        return 0.0  # no data, gini is 0 by convention

    # Sort coverage values in ascending order
    coverage_values = sorted(coverage_hist.keys())

    # N = total number of bases
    N = sum(coverage_hist[c] for c in coverage_values)
    # If somehow N is 0, return 0 to avoid division by zero
    if N == 0:
        return 0.0

    # total coverage across all bases
    total_coverage = sum(c * coverage_hist[c] for c in coverage_values)
    if total_coverage == 0:
        # All depths are zero => distribution is uniform(=0) => Gini = 0
        return 0.0

    # Build the Lorenz points (x_i, y_i)
    # x_i = cumulative fraction of bases
    # y_i = cumulative fraction of coverage
    lorenz_points = []
    pop_cum = 0
    coverage_cum = 0

    # Start from (0, 0) on the Lorenz curve
    lorenz_points.append((0.0, 0.0))

    for c in coverage_values:
        freq = coverage_hist[c]
        pop_cum_prev = pop_cum
        coverage_cum_prev = coverage_cum

        pop_cum += freq
        coverage_cum += c * freq

        x1 = pop_cum_prev / N
        y1 = coverage_cum_prev / total_coverage
        x2 = pop_cum / N
        y2 = coverage_cum / total_coverage

        # Append the new point (x2, y2)
        lorenz_points.append((x2, y2))

    # Now, approximate the area under the Lorenz curve via trapezoids
    area_under_lorenz = 0.0
    for i in range(len(lorenz_points) - 1):
        x1, y1 = lorenz_points[i]
        x2, y2 = lorenz_points[i + 1]
        base = x2 - x1           # horizontal distance
        avg_height = (y1 + y2) / 2
        area_under_lorenz += base * avg_height

    # Finally, Gini = 1 - 2 * area under the Lorenz curve
    gini = 1 - 2 * area_under_lorenz
    return gini

def getGiniCoeff(regions, genome_length, alpha=1.2):
    """Calculate the adjusted score for fair distribution of depths."""
    # 1) Build a histogram of coverage
    # start = time.time()
    coverage_hist = build_coverage_hist(regions, genome_length)
    coverage_hist_transformed = build_transformed_coverage_hist(regions, genome_length)
    # 2) Compute the Gini coefficient from that histogram
    # gini_hist = gini_coefficient_from_hist(coverage_hist)
    gini_hist_transformed = gini_coefficient_from_hist(coverage_hist_transformed)
    # print("\tTime to build gini", time.time() - start)
    gini = gini_hist_transformed
    # depths = []
    # deths = defaultdict(int)
    # start = time.time()
    # for region in regions:
    #     for i in range(region[0], region[1] + 1):
    #         depths.append(region[2])
    #         deths[i] = region[2]
    # gini = gini_coefficient(depths)
    # print("Time to build nonzero gini", time.time() - start)
    # depths_wzero = []
    # start = time.time()
    # for ix in range(0, (genome_length)):
    #     if ix in deths:
    #         depths_wzero.append(deths[ix])
    #     else:
    #         depths_wzero.append(0)
    # gini_with_zeros = gini_coefficient(depths_wzero)
    # print("Time to build zero gini", time.time() - start)
    # print(f"Histogram gini: {gini_hist:.4f}, regular gini: {gini:.4f}, with 0s: {gini_with_zeros:.4f}")

    #
    # Option A: log-transform Gini directly, so that 0 => 1, 1 => 0:
    #
    gini_log = 0.0
    if gini <= 1.0:
        # gini_log = math.log2(2 - gini)  # range ~ [0..1]
        gini_log = alpha* math.sqrt(1 - gini)
        # gini_log = 1-gini
        # clamp
        gini_log = max(0.0, min(1.0, gini_log))
        # multiply by alpha
    #
    # Option B: also log-transform the breadth
    # breadth = fraction in [0..1]
    #

    #
    # Now combine them.  Suppose we do an equal-weighted average:
    #
    return gini_log

def getBreadthOfCoverage(breadth, genome_length):
    breadth = max(0.0, min(1.0, breadth))
    if breadth > 0:
        # log2(1 + breadth) maps [0..1] -> [0..1] but “compresses” near 1
        breadth_log = math.log2(1 + breadth)
        # For breadth=1 => log2(2) => 1.0
        # For breadth=0 => log2(1) => 0.0
    else:
        breadth_log = 0.0

    # clamp to [0..1]
    breadth_log = max(0.0, min(1.0, breadth_log))
    return breadth_log

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


def detect_regions_from_runlength(ref_runs, avg_read_length):
    """
    ref_runs: dict { refName: [(start, end, cov), ...] }
    Returns a dict { refName: [(region_start, region_end), ...] }
    where each region is coverage >= 1
    """

    coverage_regions = {}

    for ref, runs in ref_runs.items():
        # We'll find contiguous runs where coverage>=1
        covered_regions = []
        current_start = None

        for (start, end, cov) in runs:
            if cov >= 1:
                # This run is covered
                if current_start is None:
                    # start new region
                    current_start = start
                # else keep extending
            else:
                # coverage=0 => close any open region
                if current_start is not None:
                    covered_regions.append((current_start, start - 1))
                    current_start = None

        # if we ended with a region open
        if current_start is not None:
            # the last run had cov≥1, so close at that run's end
            # but we need the 'end' from the run if we had multiple consecutive ≥1 runs
            # so we can track that in a separate variable or just do:
            # Actually, let's handle it properly:
            last_run = runs[-1]
            covered_regions.append((current_start, last_run[1]))
            current_start = None

        # Now optionally subdivide big regions by avg_read_length
        adjusted_regions = []
        for (rstart, rend) in covered_regions:
            region_len = rend - rstart + 1
            if region_len > avg_read_length:
                # subdivide
                n = max(1, region_len // avg_read_length)
                chunk_size = region_len // n
                for i in range(n):
                    sub_start = rstart + i*chunk_size
                    # be careful with the last chunk
                    sub_end = rstart + (i+1)*chunk_size - 1
                    if i == n-1:
                        # last chunk goes to rend
                        sub_end = rend
                    adjusted_regions.append((sub_start, sub_end))
            else:
                adjusted_regions.append((rstart, rend))

        coverage_regions[ref] = adjusted_regions

    return coverage_regions




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
    print("Reading depth information from depth file, this can take quite some time....")
    regions_by_chr = defaultdict(list)

    last_pos = None
    current_chrom = None
    last_depth = None
    current_region = []

    def add_region(positions, depth):
        """
        Convert a list of positions and their shared depth into a (start, end, depth) tuple.
        Skips regions where depth is 0.
        """
        if not positions or depth == 0:
            return None
        start = min(positions)
        end = max(positions)
        return (start, end, depth)

    with open(depthfile, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue  # skip empty lines if any

            splitline = line.split('\t')
            if len(splitline) < 3:
                continue
            reference_name = splitline[0]
            pos = int(splitline[1])
            depth = int(splitline[2])

            # If this is the first line or the chromosome changed:
            if current_chrom != reference_name:
                # If there was a region in progress for the previous chromosome, close it out.
                if current_chrom and current_region:
                    region_tuple = add_region(current_region, last_depth)
                    if region_tuple:
                        regions_by_chr[current_chrom].append(region_tuple)

                # Reset for this new chromosome
                current_chrom = reference_name
                last_pos = pos
                last_depth = depth
                current_region = [pos]
                continue

            # If on the same chromosome, check continuity and depth
            if (pos == last_pos + 1) and (depth == last_depth):
                # Same depth and consecutive position -> extend current region
                current_region.append(pos)
            else:
                # Finalize the old region first
                region_tuple = add_region(current_region, last_depth)
                if region_tuple:
                    regions_by_chr[current_chrom].append(region_tuple)

                # Start a new region
                current_region = [pos]

            # Update trackers
            last_pos = pos
            last_depth = depth

        # After the loop, if there's a region that hasn't been saved yet, finalize it:
        if current_region:
            region_tuple = add_region(current_region, last_depth)
            if region_tuple:
                regions_by_chr[current_chrom].append(region_tuple)
    return regions_by_chr


def read_depth_as_runs(depthfile, reference_lengths):
    """
    Reads a depth file that only has lines for depth >= 1 and
    constructs a run-length representation for each reference.

    :param depthfile: Path to the depth file.
    :param reference_lengths: dict { refname: length_of_ref, ... }
    :return: dict of { refname: [(start, end, depth), ...], ... }
    """

    # We'll store lines keyed by reference, each a list of (position, depth),
    # so we can sort them and convert to runs.
    coverage_positions = defaultdict(list)

    # 1) Read file lines
    print("Reading sparse depth file...")
    with open(depthfile, 'r') as fh:
        for line in fh:
            # typical format:  refName pos depth
            ref, p, d = line.strip().split('\t')
            p = int(p)
            d = int(d)
            coverage_positions[ref].append((p, d))

    # 2) Build run-length coverage for each ref
    ref_runs = {}
    for ref, pos_depth_list in coverage_positions.items():
        # Sort by position
        pos_depth_list.sort(key=lambda x: x[0])

        # total length of the reference from your dictionary
        ref_len = reference_lengths.get(ref, 0)
        if ref_len == 0:
            print(f"Warning: reference {ref} not in reference_lengths, skipping.")
            continue

        runs = []
        last_pos = 0

        # Go line-by-line in ascending position order
        for (pos, depth) in pos_depth_list:
            # If there's a gap from last_pos+1 to pos-1 => coverage=0
            if pos > last_pos + 1:
                runs.append((last_pos + 1, pos - 1, 0))

            # This position has coverage≥1
            runs.append((pos, pos, depth))
            last_pos = pos

        # If we still haven't covered the end of the reference, fill with 0
        if last_pos < ref_len:
            runs.append((last_pos + 1, ref_len, 0))

        # 3) Merge consecutive runs with the same coverage
        merged = []
        for run in runs:
            if not merged:
                merged.append(run)
            else:
                # compare to the last run
                last_start, last_end, last_cov = merged[-1]
                cur_start,  cur_end,  cur_cov = run

                # If coverage is the same and they are contiguous
                if last_cov == cur_cov and last_end + 1 == cur_start:
                    # merge
                    merged[-1] = (last_start, cur_end, last_cov)
                else:
                    merged.append(run)

        ref_runs[ref] = merged

    # Also handle references that had no positions with coverage≥1
    # (they won't appear in coverage_positions at all)
    for ref, length in reference_lengths.items():
        if ref not in ref_runs:
            # entire ref is coverage=0
            ref_runs[ref] = [(1, length, 0)]

    return ref_runs
# Function to apply weights and format the result (assuming `format_non_zero_decimals` is defined)
def apply_weight(value, weight):
    try:
        # Convert the value to float before multiplying by the weight
        value = float(value)
        return value * weight
    except (TypeError, ValueError):
        # If value cannot be converted to a float, return 0 or handle as needed
        return 0

def compute_tass_score(count, weights):
    """
    count is a dictionary that might look like:
    {
      'normalized_disparity': <some_value>,
      'alignment_score': <some_value>,
      'meangini': <some_value>,
      'diamond': {'identity': <some_value>},
      'k2_disparity_score': <some_value>,
      ...
    }
    We apply the known formula for TASS Score using the provided weights.
    """
    # Summation of each sub-score * weight
    tass_score = sum([
        apply_weight(count.get('normalized_disparity', 0), weights.get('disparity_score', 0)),
        apply_weight(count.get('alignment_score', 0),       weights.get('mapq_score', 0)),
        apply_weight(count.get('meangini', 0),             weights.get('gini_coefficient', 0)),
        apply_weight(count.get('breadth_total', 0),             weights.get('breadth_total', 0)),
        apply_weight(count.get('diamond', {}).get('identity', 0),
                     weights.get('diamond_identity', 0)),
        apply_weight(count.get('k2_disparity_score', 0),         weights.get('k2_disparity_score', 0))
    ])
    return tass_score

def count_reference_hits(bam_file_path, depthfile, covfile, bedgraph = None):
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
                length = reference_lengths.get(ref, 0),
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
        if bedgraph:
            # read in the bedgraph
            with open(bedgraph, 'r') as bed_file:
                reader = csv.reader(bed_file, delimiter='\t')

                for row in reader:
                    chrom, start, end, depth = row[0], int(row[1]), int(row[2]), int(row[3])
                    # add the regions to the reference_coverage at the chrom

                    if "covered_regions" not in reference_coverage[chrom]:
                        reference_coverage[chrom]['covered_regions'] = []
                    if chrom in reference_coverage:
                        #add one for start since bedgraph is 0 based
                        reference_coverage[chrom]['covered_regions'].append((start+1, end, depth))
        elif depthfile:
            # Detect regions based on the depth file and average read length
            print("Determining regions from a depthfile")
            ref_regions = detect_regions_from_depth(
                reference_coverage,
                depthfile,
                average_read_length
            )
            for k, v in ref_regions.items():
                if k in reference_coverage:
                    reference_coverage[k]['covered_regions'] = v
        else:
            print(f"Please provide a bedgraph from bedtools genomecov -ibam OR a depth file from samtools")

        for reference_name, refd in reference_coverage.items():
            # get the % of positions > 0 vs. the length
            total_positions = refd['length']
            sum_of_positions_covered = 0
            for region in refd['covered_regions']:
                start, end, depth = region
                sum_of_positions_covered += end - start
            coverage = (100*sum_of_positions_covered) / total_positions if total_positions > 0 else 0
            refd['coverage'] = f"{coverage:.2f}"

        if not covfile:
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
                else:
                    unaligned += 1
        else:
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
    # Write to file in a newline format
    weights = {
        'mapq_score': args.mapq_weight,
        'diamond_identity': args.diamond_identity_weight,
        'disparity_score': args.disparity_score_weight,
        'gini_coefficient': args.gini_weight,
        "breadth_weight": args.breadth_weight,
        "k2_disparity_score": args.k2_disparity_score_weight,
        'siblings_score': 0
    }

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
            + breadth * breadth_weight \
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
                        print(f"Error in: {e}")
                i+=1
    reference_hits, aligned_total, total_reads = count_reference_hits(
        inputfile,
        args.depth,
        covfile,
        args.bedgraph
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
                    isolate = splitline[9]

                    isSpecies = False if species_taxid != taxid else True
                    # fine value where assembly == accession from reference_hits
                    if accession in assembly_to_accession:
                        for acc in assembly_to_accession[accession]:
                            if acc in reference_hits:
                                reference_hits[acc]['isSpecies'] = isSpecies
                                if args.compress_species:
                                    reference_hits[acc]['toplevelkey'] = species_taxid
                                else:
                                    reference_hits[acc]['toplevelkey'] = taxid
                                current_strain = strain  # Start with the original strain value
                                if current_strain == "na":
                                    current_strain = "≡"
                                if isolate and isolate != "" and isolate != "na":
                                    current_strain = f"{current_strain} {isolate}"
                                reference_hits[acc]['strain'] = current_strain
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

    # Step 2: Define a function to calculate disparity for each organism
    # Define a function to calculate disparity with softer variance influence
    # Step 2: Define a function to dynamically dampen variance based on the proportion of reads
    def calculate_disparity(numreads, total_reads, variance_reads, k=1000):
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
                    'lengths': [],
                    'depths': [],
                    "taxids": [],
                    "accs": [],
                    'assemblies': [],
                    "coverages": [],
                    "breadth_total": 0,
                    "prevalence_disparity": 0,
                    "coeffs": [],
                    'baseqs': [],
                    "isSpecies": True if  args.compress_species  else data['isSpecies'],
                    'strainslist': [],
                    'covered_regions': 0,
                    'name': data['name'],  # Assuming the species name is the same for all strains
            }

            try:
                # if accession isn't NC_042114.1 skip
                if len(data.get('covered_regions', [])) > 0 :
                    gini_strain = getGiniCoeff(data['covered_regions'], data['length'])
                else:
                    gini_strain = 0
                species_aggregated[top_level_key]['coeffs'].append(gini_strain)
                species_aggregated[top_level_key]['taxids'].append(data.get('taxid', "None"))
                species_aggregated[top_level_key]['numreads'].append(data['numreads'])
                species_aggregated[top_level_key]['covered_regions'] += len(data['covered_regions'])
                species_aggregated[top_level_key]['coverages'].append(data['coverage'])
                species_aggregated[top_level_key]['baseqs'].append(data['meanbaseq'])
                species_aggregated[top_level_key]['lengths'].append(data['length'])
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
                print(f"Error in top level: {e}")
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
        # calcualte the total coverage by summing all lengths and getting covered bases
        total_length = sum(aggregated_data['lengths'])
        covered_bases = sum([x * y for x, y in zip(aggregated_data['lengths'], aggregated_data['coverages'])])
        aggregated_data['breadth_total'] = covered_bases / total_length if total_length > 0 else 0
        aggregated_data['total_length'] = total_length
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
        else:
            if args.ignore_missing_inputs:
                weights['k2_disparity_score'] = 0
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
            if len(disparities) > 1:
                normalized_disparity = 0  # If all disparities are the same, set to 0
            else:
                normalized_disparity = 1
        # Store the normalized disparity
        aggregated_data['normalized_disparity'] = normalized_disparity
        # print(f"Entry Top Key: {top_level_key}")
        # print(f"\tName: {aggregated_data['name']}")
        # print(f"\tNum Reads: {aggregated_data['numreads']}")
        # print(f"\tK2 Reads: {aggregated_data['k2_numreads']}")
        # print(f"\tPrev. Disparity: {aggregated_data['disparity']}")
        # print(f"\t^Norm. Disparity: {aggregated_data['normalized_disparity']}")

    # Function to normalize the MAPQ score to 0-1 based on a maximum MAPQ value
    def normalize_mapq(mapq_score, max_mapq=60, min_mapq=0):
        # Normalize the MAPQ score to be between 0 and 1
        probability = 1-(10 ** (-mapq_score / 10))
        # convert the min and max mapq
        # Normalize the MAPQ score to be between 0 and 1 based on the min and max
        return probability  # Ensure values stay between 0 and 1

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
    all_mapqs = [x['meanmapq'] for x in species_aggregated.values()]
    # get min, max of all mapqs
    min_mapq = min(all_mapqs)
    max_mapq = max(all_mapqs)

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
        # get all mapq scores
        normalized_mapq = normalize_mapq(value.get('meanmapq', 0), max_mapq, min_mapq)
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
        else:
            if args.ignore_missing_inputs:
                weights['diamond_identity'] = 0


    pathogens = import_pathogens(pathogenfile)
    # Next go through the BAM file (inputfile) and see what pathogens match to the reference, use biopython
    # to do this

    if args.min_reads_align:
        # filter the reference_hits based on the minimum number of reads aligned
        print(f"Filtering for minimum reads aligned: {args.min_reads_align}")
        species_aggregated = {k: v for k, v in species_aggregated.items() if ((v['covered_regions'])) >= int(args.min_reads_align)}
        species_aggregated = {k: v for k, v in species_aggregated.items() if sum(v['numreads']) >= int(args.min_reads_align)}
    # if sum of vals in weights isnt 1 then normalize to 1
    total_weight = sum(weights.values())
    if total_weight != 1:
        for key in weights:
            if total_weight != 0:
                weights[key] = weights[key] / total_weight
            else:
                weights[key] = 0
    if args.gt:
        # import the gtdata and apply to a dict
        with open(args.gt, 'r') as f:
            gtdatareader = csv.reader(f, delimiter="\t")
            gtdata = {}
            for line in gtdatareader:
                if len(line) > 0:
                    gtdata[line[0]] = float(line[1])
            f.close()
        # iterate through the aggregated_stats and add missing accessions
        for key, value in species_aggregated.items():
            gtcov = gtdata.get(key, 0)
            species_aggregated[key]["gtcov"]=gtcov
    final_scores = calculate_scores(
        aggregated_stats=species_aggregated,
        pathogens=pathogens,
        sample_name=args.samplename,
        sample_type = args.sampletype,
        total_reads = total_reads,
        aligned_total = aligned_total,
        weights = weights
    )
    # print("Final Scores:")
    # for entry in final_scores:
    #     if entry['tass_score'] < 0.4:
    #         continue
    #     print("\t",entry['formatname'])
    #     print("\t\tTASS Score:",entry['tass_score'],
    #           "\n\t\tGini:", entry['gini_coefficient'],
    #           "\n\t\tAlignment:", entry['alignment_score'],
    #           "\n\t\tSiblings:", entry['siblings_score'],
    #           "\n\t\tDMND:", entry['diamond_identity'],
    #           "\n\t\tDisparity:", entry['disparity_score'],
    #           "\n\t\tK2 Disparity:", entry['k2_disparity_score'],
    #           "\n\t\tBreadth:", entry['breadth_total'],
    #           "\n\t\tBreadth2:", ((entry['breadth_total'])**2),
    #           "\n\t\tBreadthLog:", math.log2(2-entry['breadth_total']),
    #           "\n\t\tBreadthsqrt:", math.sqrt(1-entry['breadth_total'])
    #     )
    write_to_tsv(output, final_scores)
    if (args.gt and args.optimize):
        from scipy.optimize import minimize
        # Usage
        best_weights = optimize_weights(
            aggregated_stats=species_aggregated,
            pathogens=pathogens,
            sample_name="Test_Sample",
            sample_type="Unknown",
            total_reads=total_reads,
            aligned_total=aligned_total,
            initial_weights=[1, 1, 1, 1, 1, 1, 0],
            max_iterations=args.max_iterations
        )
        print("Optimized Weights:")
        for key, value in best_weights.items():
            print(f"\t{key}: {value:.3f}")
        # convert all np floats





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


def compute_cost(
    aggregated_stats,
    pathogens,
    sample_name,
    sample_type,
    total_reads,
    aligned_total,
    weights
):
    """
    Runs calculate_scores(...), then sums up penalties based on meancoverage & tass_score.
    """
    total_weight = sum(weights.values())
    if total_weight != 1:
        for key in weights:
            if total_weight != 0:
                weights[key] = weights[key] / total_weight
            else:
                weights[key] = 0
    # First, get final_scores from your existing pipeline:

    final_scores = calculate_scores(
        aggregated_stats=aggregated_stats,
        pathogens=pathogens,
        sample_name=sample_name,
        sample_type=sample_type,
        total_reads=total_reads,
        aligned_total=aligned_total,
        weights=weights,
    )

    cost = 0.0
    # print("__________")
    perrefcost = {}
    for rec in final_scores:
        ref = rec.get('ref')
        tass_score = rec['tass_score']
        actual_coverage = rec.get('meancoverage', 0.0)

        # The "ground-truth" coverage for this ref (already placed in aggregated_stats)
        # e.g. aggregated_stats[ref]['gtcov'] may be 0 or >0
        coverage = aggregated_stats.get(ref, {}).get('gtcov', 0.0)
        # CASE 1: coverage == 0 => we want TASS ~ 0,
        #         also penalize if actual_coverage is large (contradiction).
        if coverage <= 0:
            # The bigger the TASS (wrongly claiming positivity),
            # or the bigger the actual coverage (contradiction),
            # the larger this penalty.
            cost += (1.0 + actual_coverage) * (tass_score ** 2)
            perrefcost[ref] = (1.0 + actual_coverage) * (tass_score ** 2)
        else:
            # CASE 2: coverage > 0 => we want TASS near 1 for big coverage.
            #  (A) penalize if TASS is low => cost ~ coverage * (1 - TASS)^2
            #  (B) penalize if TASS is high but actual coverage is too small
            #      => e.g. coverage_short = (coverage - actual_coverage)
            #             cost2 = (tass_score^2) * coverage_short
            #      meaning if ground-truth coverage is large but the pipeline
            #      found less, it’s contradictory to have TASS=high.
            cost1 = coverage * ((1.0 - tass_score) ** 2)

            coverage_short = max(0.0, coverage - actual_coverage)
            cost2 = (tass_score ** 2) * coverage_short

            cost += (cost1 + cost2)
            perrefcost[ref] =  (cost1 + cost2)
        # print("\t",ref, actual_coverage, coverage, tass_score, rec['gini_coefficient'])
    return cost, perrefcost

def random_tweak(weights, scale=0.2):
    """
    Returns a new dict, each weight randomly offset by up to ±scale.
    Clamps at 0 (no negative).
    """
    new_w = {}
    for k, v in weights.items():
        delta = random.uniform(-scale, scale)
        candidate = v + delta
        if candidate < 0 :
            candidate = 0.0

        new_w[k] = candidate
    return new_w

def optimize_weights(
    aggregated_stats,
    pathogens,
    sample_name="No_Name",
    sample_type="Unknown",
    aligned_total=0,
    total_reads=0,
    initial_weights=None,
    max_iterations=1000
):
    if initial_weights is None:
        # Initialize with equal weights
        initial_weights = [1, 1, 1, 1, 1, 1, 0]  # Using dtype=object for mixed types


    # Normalize initial weights
    total_sum = 0
    for weight in initial_weights:
        if isinstance(weight, (list)):
            total_sum += sum(weight)
        else:
            total_sum += weight
    normalized_weights = []
    for weight in initial_weights:
        if isinstance(weight, list):
            # Normalize each element in the list
            normalized_sub_weights = [w / total_sum for w in weight]
            normalized_weights.append(normalized_sub_weights)
        else:
            # Normalize the single weight
            normalized_weight = weight / total_sum
            normalized_weights.append(normalized_weight)
    # Define bounds for each weight
    bounds = [(0, 1) for _ in initial_weights]

    # Define the constraint that weights sum to 1
    constraints = {'type': 'eq', 'fun': lambda w: sum(w) - 1}

    # Objective function to minimize
    def objective(w):
        weight_dict = {
            'mapq_score': w[0],
            'diamond_identity': w[1],
            'disparity_score': w[2],
            'gini_coefficient': w[3],
            'breadth_score': w[4],
            'k2_disparity_score': w[5],
            'siblings_score': w[6],
        }
        cost, perrefcost = compute_cost(
            aggregated_stats,
            pathogens,
            sample_name,
            sample_type,
            total_reads,
            aligned_total,
            weight_dict
        )
        return -cost  # Negate if you want to maximize the cost

    # Perform optimization
    result = minimize(
        objective,
        initial_weights,
        method='SLSQP',
        bounds=bounds,
        constraints=constraints,
        options={'maxiter': max_iterations, 'disp': True}
    )

    if result.success:
        optimized_weights = {
            'mapq_score': result.x[0],
            'diamond_identity': result.x[1],
            'disparity_score': result.x[2],
            'gini_coefficient': result.x[3],
            'breadth_score': result.x[4],
            'k2_disparity_score': result.x[5],
            'siblings_score': result.x[6],
        }
        # Extract the optimized weights
        optimized_weights = {
            'mapq_score': result.x[0],
            'diamond_identity': result.x[1],
            'disparity_score': result.x[2],
            'gini_coefficient': result.x[3],
            'breadth_score': result.x[4],
            'k2_disparity_score': result.x[5],
            'siblings_score': result.x[6],
        }

        # Compute final cost and per-reference scores with optimized weights
        final_cost, final_per_reference_scores = compute_cost(
            aggregated_stats,
            pathogens,
            sample_name,
            sample_type,
            total_reads,
            aligned_total,
            optimized_weights
        )

        for k, v in final_per_reference_scores.items():
            print(f"\t{k}: {v:.4f}")
        return optimized_weights
    else:
        raise RuntimeError(f"Optimization failed: {result.message}")


def calculate_scores(
        aggregated_stats,
        pathogens,
        sample_name="No_Name",
        sample_type="Unknown",
        total_reads=0,
        aligned_total = 0,
        weights={}
    ):
    """
    Write reference hits and pathogen information to a TSV file.

    Args:
    reference_hits (dict): Dictionary with reference names as keys and counts as values.
    pathogens (dict): Dictionary with reference names as keys and dictionaries of additional attributes as values.
    output_file_path (str): Path to the output TSV file.
    """
    final_scores = []


    # Write the header row
    sumfgini = 0
    totalcounts = 0

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
        breadth_total = count.get('breadth_total', 0) / 100 if count.get('breadth_total', 0) else 0
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

        # Apply weights to the relevant scores
        tass_score = compute_tass_score(
            count,
            weights,
        )
        final_scores.append(
            dict(
                tass_score=tass_score,
                formatname=formatname,
                sample_name=sample_name,
                sample_type=sample_type,
                fullstrains=fullstrains,
                listpathogensstrains=listpathogensstrains,
                percent_total=percent_total,
                percent_aligned=percent_aligned,
                callfamclass=callfamclass,
                total_reads=sum(count['numreads']),
                reads_aligned=countreads,
                meancoverage=count.get('meancoverage', 0),
                breadth_total = breadth_total,
                total_length = count.get('total_length', 0),
                is_annotated=is_annotated,
                pathogenic_sites=pathogenic_sites,
                is_pathogen=is_pathogen,
                ref=ref,
                status=status,
                pathogenic_reads = pathogenic_reads,
                gini_coefficient=count.get('meangini', 0),
                meanbaseq=count.get('meanbaseq', 0),
                meanmapq=count.get('meanmapq', 0),
                alignment_score=count.get('alignment_score', 0),
                k2_disparity_score=count.get('disparity_cv', 0),
                disparity_score=count.get('normalized_disparity', 0),
                diamond_identity=count.get('diamond', {}).get('identity', 0),
                covered_regions=count.get('covered_regions', 0),
                siblings_score=count.get('raw_disparity', 0),
                k2_reads=count.get('k2_numreads', 0),
                gtcov=count.get('gtcov', 0),


            )
        )

    return final_scores

def write_to_tsv(output_path, final_scores):
    total=0
    fulltotal=0
    with open (output_path, 'w') as file:
        header = "Detected Organism\tSpecimen ID\tSample Type\t% Reads\t# Reads Aligned\t% Aligned Reads\tCoverage\tIsAnnotated\tPathogenic Sites\tMicrobial Category\tTaxonomic ID #\tStatus\tGini Coefficient\tMean BaseQ\tMean MapQ\tMean Coverage\tMean Depth\tAnnClass\tisSpecies\tPathogenic Subsp/Strains\tK2 Reads\tParent K2 Reads\tMapQ Score\tDisparity Score\tProtein Identity Score\tSiblings score\tTASS Score\n"
        file.write(f"{header}")
        for entry in final_scores:
            is_pathogen = entry.get('is_pathogen', "Unknown")
            formatname = entry.get('formatname', "N/A")
            sample_name = entry.get('sample_name', "N/A")
            ref = entry.get('ref', "N/A")
            percent_aligned = entry.get('percent_aligned', 0)
            is_pathogen = entry.get('is_pathogen', "Unknown")
            status = entry.get('status', "N/A")
            is_annotated= entry.get('is_annotated', "N/A")
            sample_type = entry.get('sample_type', "N/A")
            callfamclass = entry.get('callfamclass', "N/A")
            gini_coefficient = entry.get('gini_coefficient', 0)
            meanbaseq = entry.get('meanbaseq', 0)
            meanmapq = entry.get('meanmapq', 0)
            meancoverage = entry.get('meancoverage', 0)
            meandepth = entry.get('meandepth', 0)
            annClass = entry.get('annClass', "N/A")
            isSpecies = entry.get('isSpecies', False)
            k2_reads = entry.get('k2_reads', 0)
            k2_parent_reads = entry.get('k2_parent_reads', 0)
            mapq_score = entry.get('mapq_score', 0)
            disparity_score = entry.get('disparity_score', 0)
            diamond_identity = entry.get('diamond_identity', 0)
            siblings_score = entry.get('siblings_score', 0)
            tass_score = entry.get('tass_score', 0)
            pathogenic_sites = entry.get('pathogenic_sites', "")
            listpathogensstrains = entry.get('listpathogensstrains', [])
            countreads = entry.get('reads_aligned', 0)
            breadth_of_coverage = entry.get('breadth_total', 0)
            aligned_total = entry.get('total_reads', 0)
            pathogenic_reads = entry.get('pathogenic_reads', 0)
            percent_total = entry.get('percent_total', 0)
            if  (is_pathogen == "Primary" or is_pathogen=="Potential") and tass_score >= 0.0  :
                print(f"Reference: {ref} - {formatname}")
                print(f"\tIsPathogen: {is_pathogen}")
                print(f"\tPathogenic Strains: {listpathogensstrains}")
                percentreads = f"{100*pathogenic_reads/aligned_total:.1f}" if aligned_total > 0 and pathogenic_reads>0 else 0
                print(f"\tPathogenic SubStrain Reads: {pathogenic_reads} - {percentreads}%")
                print(f"\tAligned Strains:")
                for f in entry.get('fullstrains', []):
                    print(f"\t\t{f}")
                print(f"\tTotal reads: {entry.get('total_reads', 0)}")
                print(f"\tGini Conf: {entry.get('gini_coefficient', 0):.4f}")
                print(f"\tAlignment Score: {entry.get('alignment_score', 0):.2f}")
                print(f"\tK2 Disparity Score: {entry.get('k2_disparity_score', 0):.2f}")
                print(f"\t# Reads Aligned: {entry.get('reads_aligned', 0)}")
                print(f"\tDisparity Score: {entry.get('disparity_score', 0):.2f}")
                print(f"\tDiamond Identity: {entry.get('diamond_identity', 0):.2f}")
                print(f"\tCoverage (Mean%): {entry.get('meancoverage', 0):.2f}")
                print(f"\tBreadth: {entry.get('breadth_total', 0):.2f}")
                print(f"\tDepth (Mean): {entry.get('meandepth', 0):.2f}")
                print(f"\tGT Coverage: {entry.get('gtcov', 0):.2f}")
                print(f"\tFinal Score: {entry.get('tass_score', 0):.2f}")
                print(f"\tCovered Regions: {entry.get('covered_regions', 0)}")
                print(f"\tK2 Reads: {entry.get('k2_reads', 0)}")
                print()
                total+=1
            file.write(
                f"{formatname}\t{sample_name}\t{sample_type}\t{percent_total}\t{countreads}\t{percent_aligned}\t{breadth_of_coverage:.2f}\t"
                f"{is_annotated}\t{entry.get('pathogenic_sites')}\t{is_pathogen}\t{ref}\t{status}\t{gini_coefficient:.2f}\t"
                f"{meanbaseq:.2f}\t{meanmapq:.2f}\t{meancoverage:.2f}\t{meandepth:.2f}\t{annClass}\t{isSpecies}\t{callfamclass}\t"
                f"{k2_reads}\t{k2_parent_reads}\t{mapq_score:.2f}\t{disparity_score:.2f}\t{diamond_identity:.2f}\t"
                f"{siblings_score:.2f}\t{tass_score:.2f}\n"
            )
            fulltotal+=1
    print(f"Total pathogenic orgs: {total}, Total entire: {fulltotal}")
if __name__ == "__main__":
    sys.exit(main())
