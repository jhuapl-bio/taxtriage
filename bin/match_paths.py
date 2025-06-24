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
from scipy.stats import norm
import time
from intervaltree import Interval, IntervalTree
import statistics
import math as Math
from distributions import import_distributions, body_site_map
import argparse
import re
import csv
import math
import os
from conflict_regions import determine_conflicts
import pysam
from scipy.optimize import minimize

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
        "-f",
        "--fasta",
        metavar="FASTAFILE",
        nargs='+',
        default=[],
        help="The fasta file(s) used for alignment originally. Optional but used for the signature creation step if you dont want to make signatures of all reads. Dramatically improves memory and speed",
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
        "--sensitive", default=False,  help="Use sensitive mode to detect greater array of variants",  action='store_true'
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
    parser.add_argument("--output_dir", required=False, help="Output directory for results.")
    parser.add_argument(
        "--disparity_score_weight",
        metavar="DISPARITYSCOREWEIGHT",
        type=float,
        default=0.05,
        help="value of weight for disparity reads vs other organisms in final TASS Score",
    )
    parser.add_argument(
        "--alpha",
        metavar="MAPQWEIGHT",
        type=float,
        default=1.2,
        help="alpha value for lorenz curve with gini calculation",
    )
    parser.add_argument(
        "-X",
        '--cpu_count',
        metavar="CPUCOUNT",
        type=int,
        default=None,
        help="Overwrite the number of CPUs to use.",
    )
    parser.add_argument(
        "--hmp",
        metavar="HMP",
        type=str,
        default=None,
        help="HMP distribution of abundances of organisms on a per-sample basis. Optional but useful for some sample types that are highly diverse in content",
    )
    parser.add_argument(
        "--readcount",
        default=None,
        metavar="READCOUNT",
        type=float,
        help="What is the full read count, overwrites the total count from the BAM alignment.",
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
        default=0.1,
        help="value of weight for breadth of coverage in final TASS Score",
    )
    parser.add_argument(
        "--minhash_weight",
        metavar="MINHASHSCORE",
        type=float,
        default=0.5,
        help="value of weight for minhash signature reduction in final TASS Score",
    )
    parser.add_argument(
        "--k2_disparity_weight",
        metavar="DISPARITYSCOREWEIGHT",
        type=float,
        default=0.0,
        help="value of weight for disparity of k2 and alignment in final TASS Score",
    )
    parser.add_argument(
        "--diamond_identity_weight",
        metavar="DISPARITYSCOREWEIGHT",
        type=float,
        default=0.0,
        help="value of weight for disparity of diamond_identity in final TASS Score",
    )
    parser.add_argument(
        "--hmp_weight",
        metavar="HMPWEIGHT",
        type=float,
        default=0.0,
        help="value of weight for hmp abundance in final TASS Score",
    )

    parser.add_argument(
        "--gini_weight",
        metavar="GINIWEIGHT",
        type=float,
        default=0.50,
        help="value of weight for gini coefficient in final TASS Score",
    )
    parser.add_argument(
        "--reference_signatures",
        help="Path to the CSV file generated by `sourmash compare --csv ...`",
        required=False
    )
    parser.add_argument(
        "--dispersion_factor",
        metavar="DISP_FACTOR",
        type=float,
        default=0.50,
        help="The disperion factor boost for determining the positve nature of spread out regions of coverage, used during the gini coeff. calculations. ",
    )
    parser.add_argument(
        "--reward_factor",
        metavar="REWARD_FACTOR",
        type=float,
        default=2,
        help="Factor the reward used for the overall confidence in gini. Boosts it, essentially",
    )
    parser.add_argument(
        "--min_similarity_comparable",
        type=float,
        default=0.0,
        help="Minimum ANI threshold to consider references comparable (default=0.7)"
    )
    parser.add_argument("--only_filter", required=False, action='store_true', help="Stop after creating a filtered bamfile")
    parser.add_argument("--kmer_size", type=int, default=51, help="k-mer size for MinHash.")
    parser.add_argument("--matrix", required=False, help="A Matrix file for ANI in long format from fastANI")
    parser.add_argument("--comparisons", required=False, help="Skip comparison metrics if present, can be either csv, tsv, or xlsx")
    parser.add_argument("--failed_reads", required=False, help="Load a 2 col tsv of reference   read_id that is to be the passed reads. Remove all others and update bedgraph and cov file(s)")
    parser.add_argument("--scaled", type=int, default=2000, help="scaled factor for MinHash.")
    parser.add_argument("--coverage_length", type=int, default=500, help="Length of each coverage chunk.")
    parser.add_argument("--coverage_spacing", type=int, default=4, help="Allowed gap of zero coverage in a chunk.")
    parser.add_argument("--min_threshold", type=float, default=0.2, help="Min Jaccard similarity to report.")
    parser.add_argument("--abu_file", required=False, help="Path to ground truth coverage file (abu.txt).")
    parser.add_argument("--use_variance", required=False, action='store_true', help="Use variance instead of Gini index.")
    parser.add_argument("--apply_ground_truth", required=False, action='store_true', help="Overwrites abu file for F1 score metrics. Parses the read name for the actual reference.")
    parser.add_argument("--sigfile", required=False, type=str, help="Skip signatures comparison if provided ad load it instead")
    parser.add_argument("--config", required=False, type=str, help="Configuration file for generating a minimal output file for LIMS integration or other data import system. ")
    parser.add_argument("--fast", required=False, action='store_true', help="FAST Mode enabled. Uses Sourmash's SBT bloom factory for querying similarity of jaccard scores per signature per region. This is much faster than the original method but requires a pre-built SBT file which takes time and can lead to false positive region matches.")
    parser.add_argument('--gap_allowance', type=float, default=0.1, help="Gap allowance for determining merging of regions")
    parser.add_argument('--jump_threshold', type=float, default=None, help="Gap allowance for determining merging of regions")
    parser.add_argument(
        "--filtered_bam", default=False,  help="Create a filtered bam file of a certain name post sourmash sigfile matching..", type=str
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
        return math.log10(1 + d)



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
        length = end - start
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
        length = end - start    # how many bases in this interval
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

def get_dynamic_reward_factor(genome_length, baseline=5e4, max_length=1e7, max_reward=2):
    """
    Computes a dynamic reward factor:
      - Returns 1 when genome_length is at or below baseline.
      - Returns max_reward when genome_length is at or above max_length.
      - Otherwise, linearly interpolates between 1 and max_reward.
    """
    if genome_length <= baseline:
        return 1.0
    elif genome_length >= max_length:
        return max_reward
    else:
        # Linear interpolation between 1 and max_reward
        return 1.0 + (max_reward - 1.0) * ((genome_length - baseline) / (max_length - baseline))


def getGiniCoeff(regions, genome_length, alpha=1.8, baseline=5e5, max_length=1e9, reward_factor=2, beta=0.5):
    """
    Calculate an adjusted 'Gini-based' score for the fair distribution of coverage,
    and then penalize (or boost) it according to the disparity in positions of the regions.

    Parameters:
      - regions: list of tuples (start, end, depth)
      - genome_length: total genome length
      - alpha: parameter for transforming the raw Gini
      - baseline, max_length, reward_factor: parameters for length-based scaling
      - beta: weight for the positional dispersion factor

    The final score is a product of the (transformed) Gini measure, a scaling factor
    based on genome length, and a term (1 + beta * dispersion) where dispersion is higher
    when the regions are more spread out.
    """
    # 1) Build Histograms
    coverage_hist = build_coverage_hist(regions, genome_length)
    coverage_hist_transformed = build_transformed_coverage_hist(regions, genome_length)

    # 2) Compute raw Gini from the transformed histogram
    gini = gini_coefficient_from_hist(coverage_hist_transformed)

    # 3) Transform the raw Gini (ensuring the result is in [0,1])
    if 0.0 <= gini <= 1.0:
        gini_log = alpha * math.sqrt(1 - gini)
        gini_log = max(0.0, min(1.0, gini_log))
    else:
        gini_log = 0.0

    # 4) Compute length-based scaling (using a log scale)
    gl_capped = min(genome_length, max_length)
    ratio = gl_capped / baseline
    ratio = max(ratio, 1.0)
    scaling_factor = 1.0 + reward_factor * math.log10(ratio)

    # 5) Compute the positional dispersion factor
    dispersion = position_dispersion_factor(regions, genome_length)

    # 6) Combine the measures.
    # The idea is to boost the score if the covered regions are spread out.
    final_score = gini_log * scaling_factor * (1 + beta * dispersion)
    final_score = min(1.0, final_score)
    return final_score

def position_dispersion_factor(regions, genome_length):
    """
    Compute a dispersion factor based on the positions of the regions.
    We use the midpoints of each region and compute their variance.
    For a uniformly distributed set of midpoints on [0, genome_length],
    the maximum variance is (genome_length^2)/12.
    We then take the square root of the normalized variance to obtain a
    factor between 0 and 1.
    """
    if not regions:
        return 0.0

    # Compute midpoints for each region
    midpoints = [(start + end) / 2.0 for (start, end, depth) in regions]
    mean_mid = sum(midpoints) / len(midpoints)
    variance = sum((m - mean_mid)**2 for m in midpoints) / len(midpoints)

    # Maximum variance for a uniform distribution in [0, genome_length]
    max_variance = (genome_length**2) / 12.0
    normalized_variance = variance / max_variance  # in [0,1]

    # Taking square root to keep the metric in a similar scale as a coefficient of variation
    dispersion = math.sqrt(normalized_variance)
    return dispersion

def gini_coefficient(values):
    """
    Compute the Gini coefficient using a sorted approach.
    """
    n = len(values)
    if n == 0:
        return 0.0

    # Sort the values in ascending order
    sorted_values = sorted(values)
    total = sum(sorted_values)

    if total == 0:
        return 0.0

    # Compute the weighted sum: sum(i * x_i) for i=1..n (using 1-indexing)
    weighted_sum = sum((i + 1) * x for i, x in enumerate(sorted_values))

    # Apply the formula:
    gini = (2 * weighted_sum) / (n * total) - (n + 1) / n
    return gini


def gini_coverage_spread(regions, genome_length):
    """
    Compute a spatial Gini coefficient that reflects how spread out the
    covered regions are across the genome.

    The idea:
      1. Assume each region in 'regions' is given as a tuple (start, end, depth).
         (For now, we ignore 'depth' and assume any region has a binary coverage = 1.)
      2. Calculate the gaps (uncovered segments) across the genome:
           - From position 0 to the start of the first region.
           - Between the end of one region and the start of the next.
           - From the end of the last region to the genome end.
      3. Compute the Gini coefficient on these gap lengths.

    A low Gini value means the gaps are fairly equal in length (i.e. the covered regions
    are clumped together), while a high Gini value means there is great inequality among gap
    sizes (i.e. the covered regions are very spread out).
    """
    # Ensure the regions are sorted by start coordinate
    # regions = sorted(regions, key=lambda r: r[0])

    gaps = []

    if regions:
        # Gap from genome start to first region start
        first_gap = regions[0][0]  # since genome starts at 0
        gaps.append(first_gap)

        # Gaps between successive regions
        for i in range(len(regions) - 1):
            current_end = regions[i][1]
            next_start = regions[i+1][0]
            gap = next_start - current_end
            gaps.append(gap)

        # Gap from the end of the last region to the end of the genome
        last_gap = genome_length - regions[-1][1]
        gaps.append(last_gap)
    else:
        # No covered regions: one gap equals the entire genome
        gaps.append(genome_length)
    return gini_coefficient(gaps)



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
def import_pathogens(pathogens_file):
    """Import the pathogens from the input CSV file, correctly handling commas in quoted fields."""
    pathogens_dict = {}
    # Open the file using the `with` statement
    with open(pathogens_file, 'r', newline='', encoding='utf-8') as file:
        # Create a CSV reader object that handles commas inside quotes automatically
        reader = csv.reader(file, delimiter=',', quotechar='"')

        # Iterate over each row in the CSV file
        for row in reader:
            # Assign each part of the row to variables if available
            pathogen_name = row[0] if len(row) > 0 else None
            taxid = row[1] if len(row) > 1 else None
            call_class = row[2] if len(row) > 2 else None

            pathogenic_sites = row[4] if len(row) > 4 else None
            commensal_sites = row[5] if len(row) > 5 else None
            # split and strip commensal sites spaces
            commensal_sites = [x.strip() for x in commensal_sites.split(',')] if commensal_sites else []
            pathogenic_sites = [x.strip() for x in pathogenic_sites.split(',')] if pathogenic_sites else []
            # body_site map to the commensal sites and pathogenic_stites
            commensal_sites = [body_site_map(x.lower()) for x in commensal_sites]
            pathogenic_sites = [body_site_map(x.lower()) for x in pathogenic_sites ]
            status = row[6] if len(row) > 6 else None
            pathology = row[8] if len(row) > 8 else None

            high_cons = row[7] if len(row) > 7 else False
            # denote that if FALSE set to None, else set to True
            if high_cons and high_cons != "'":
                if high_cons.lower() == "false":
                    high_cons = False
                else:
                    high_cons = True
            # Store the data in the dictionary, keyed by pathogen name
            pathogens_dict[pathogen_name] = {
                'taxid': taxid,
                'callclass': call_class,
                'pathogenic_sites': pathogenic_sites,
                'name': pathogen_name,
                'commensal_sites': commensal_sites,
                'status': status,
                'pathology': pathology,
                'high_cons': high_cons
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
    return sum(data) / len(data) if sum(numreads) > 0 else 0
    # if len(numreads) > 0:
    #     total_weight = sum(numreads)  # Sum of all weights (numreads)
    #     weighted_sum = []
    #     for x in range(0, len(numreads)):
    #         if x >= len(data):
    #             break
    #         if numreads[x] > 0:
    #             weighted_sum.append(data[x] * numreads[x]) # Sum of weight*value
    #         else:
    #             weighted_sum.append(0)
    #     weighted_sum = sum(weighted_sum)
    #     if total_weight > 0:
    #         weighted_mean = weighted_sum / total_weight
    #     else:
    #         return 0
    # else:
    #     weighted_mean = 0
    # return weighted_mean

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
            region_len = rend - rstart
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
def logarithmic_weight(breadth, min_breadth=1e-3):
    """
    Returns a weight between 0 and 1 for a given breadth value.
    - When breadth == min_breadth, the weight is 1.
    - When breadth == 1, the weight is 0.
    - Values in between are mapped logarithmically.

    Parameters:
    breadth (float): The raw breadth value (expected to be between 0 and 1).
    min_breadth (float): The minimum expected breadth (should be > 0 to avoid log(0)).
    """
    # Clamp breadth to the [min_breadth, 1] range
    breadth = max(min_breadth, min(breadth, 1.0))

    # Normalize the log value
    normalized = (math.log(breadth) - math.log(min_breadth)) / (0 - math.log(min_breadth))
    weight = 1 - normalized
    return weight
def compute_tass_score(count, weights):
    """
    count is a dictionary that might look like:
    {
      'normalized_disparity': <some_value>,
      'alignment_score': <some_value>,
      'meangini': <some_value>,
      ...
    }
    We apply the known formula for TASS Score using the provided weights.
    """
    # normalize score of count to 0-1, min max range is 3, if larger than 3 set to 3 first
    # convert z score to percentile

    # Summation of each sub-score * weight



    tass_score = sum([
        apply_weight(count.get('normalized_disparity', 0), weights.get('disparity_score', 0)),
        apply_weight(count.get('alignment_score', 0),       weights.get('mapq_score', 0)),
        apply_weight(count.get('meanminhash_reduction', 0),       weights.get('minhash_weight', 0)),
        apply_weight(count.get('meangini', 0),             weights.get('gini_coefficient', 0)),
        apply_weight(count.get('hmp_percentile', 0),             weights.get('hmp_weight', 0)),
        apply_weight(count.get('log_weight_breadth', 0),             weights.get('breadth_weight', 0)),
        apply_weight(count.get('k2_disparity', 0),         weights.get('k2_disparity_weight', 0)),
        apply_weight(count.get('diamond', {}).get('identity', 0),
                     weights.get('diamond_identity', 0))  ])

    return tass_score


def adjust_bedgraph_coverage(bedgraph_path, excluded_interval_trees, output_bedgraph=None):
    """
    Adjust bedgraph coverage by removing 1 depth value from excluded regions.

    Parameters:
    - bedgraph_path (str): Path to the input bedgraph file.
    - excluded_interval_trees (dict): Dictionary of IntervalTrees per reference.
    - output_bedgraph (str, optional): Path to save the adjusted bedgraph. If None, data is not saved.

    Returns:
    - list of tuples: Adjusted bedgraph entries as (chrom, start, end, adjusted_depth).
    """
    adjusted_bedgraph = []

    with open(bedgraph_path, 'r') as bed_file:
        reader = csv.reader(bed_file, delimiter='\t')
        for row in reader:
            if len(row) < 4:
                continue  # Skip malformed lines
            chrom, start, end, depth = row[0], int(row[1]), int(row[2]), int(row[3])
            # Check if chrom has excluded regions
            if chrom in excluded_interval_trees:
                tree = excluded_interval_trees[chrom]
                overlaps = sorted(tree.overlap(start, end))

                if overlaps:
                    # Initialize current position
                    current_start = start
                    for interval in sorted(overlaps, key=lambda x: x.begin):
                        exclude_start, exclude_end = interval.begin, interval.end

                        # No overlap if exclude_end <= current_start or exclude_start >= end
                        if exclude_end <= current_start or exclude_start >= end:
                            continue

                        # Determine overlapping region
                        overlap_start = max(current_start, exclude_start)
                        overlap_end = min(end, exclude_end)

                        # Add non-overlapping region before the excluded region
                        if current_start < overlap_start:
                            adjusted_bedgraph.append((chrom, current_start, overlap_start, depth))

                        # Add overlapping region with depth decremented by 1
                        adjusted_depth = max(depth - 1, 0)  # Ensure depth doesn't go below 0
                        adjusted_bedgraph.append((chrom, overlap_start, overlap_end, adjusted_depth))

                        # Update current_start
                        current_start = overlap_end

                    # Add any remaining non-overlapping region after the last excluded region
                    if current_start < end:
                        adjusted_bedgraph.append((chrom, current_start, end, depth))
                else:
                    # No overlaps; retain the original bedgraph entry
                    adjusted_bedgraph.append((chrom, start, end, depth))
            else:
                # No excluded regions for this chrom; retain the original bedgraph entry
                adjusted_bedgraph.append((chrom, start, end, depth))

    # Optionally, write the adjusted bedgraph to a new file
    if output_bedgraph:
        with open(output_bedgraph, 'w', newline='') as outfile:
            writer = csv.writer(outfile, delimiter='\t')
            for entry in adjusted_bedgraph:
                writer.writerow(entry)

    return adjusted_bedgraph

def count_reference_hits(bam_file_path,alignments_to_remove=None):
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

        reference_stats = defaultdict(dict)
        for ref in bam_file.header.references:
            # get average read length
            reference_lengths[ref] = bam_file.get_reference_length(ref)
            reference_stats[ref] = dict(
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
                sum_mapq = 0,
                sum_baseq = 0,
                count_baseq = 0,
                count_mapq = 0,
                total_reads = 0,
                read_positions = [],
                total_length = 0,
                unique_read_ids = set(),
            accession = ref,
                isSpecies = False,
                covered_regions = [],  # Store regions covered by reads
            )
        # Open BAM file
        unique_read_ids = set()  # To track unique read IDs
        total_length = 0
        total_reads = 0
        readlengths = []
        removed_reads = defaultdict(dict )
        read_stats = defaultdict(list)
        excluded_regions = defaultdict(list)
        # check if the alignment was paired end or single end
        start_time = time.time()
        for read in bam_file.fetch():
            # Only process paired reads
            # if not read.is_paired:
            #     continue
            # For paired-end reads, only count the first mate to avoid double-counting

            # if not read.is_read1 and read.is_paired:
            #     continue

            ref = read.reference_name
            total_reads += 1

            # Skip reads that are in alignments_to_remove if provided
            if alignments_to_remove and read.query_name in alignments_to_remove and ref in alignments_to_remove[read.query_name]:
                continue

            # Create a unique key for the read based on its query name and strand
            read_id_key = f"{read.query_name}:{read.is_reverse}"

            # Handle unique reads
            if read_id_key not in reference_stats[ref]["unique_read_ids"]:
                reference_stats[ref]["unique_read_ids"].add(read_id_key)
                reference_stats[ref]["total_reads"] += 1
                reference_stats[ref]["total_length"] += read.query_length

            # Process only mapped reads
            if not read.is_unmapped:
                aligned_reads += 1
                # Accumulate base quality scores
                if read.query_qualities:
                    reference_stats[ref]["sum_baseq"] += sum(read.query_qualities)
                    reference_stats[ref]["count_baseq"] += len(read.query_qualities)

                # Accumulate mapping quality scores
                reference_stats[ref]["sum_mapq"] += read.mapping_quality
                reference_stats[ref]["count_mapq"] += 1

                # Record read positions for coverage calculation
                reference_stats[ref]["read_positions"].append((read.reference_start, read.reference_end))
        # Calculate statistics for each reference
        for ref, length in reference_lengths.items():
            stats = reference_stats.get(ref, None)
            if stats is None:
                reference_stats[ref]["avg_read_length"] = 0
                reference_stats[ref]["meanbaseq"] =0
                reference_stats[ref]["length"] = reference_lengths.get(ref, 0)
                reference_stats[ref]["meanmapq"] = 0
                reference_stats[ref]['meandepth'] = 0
                reference_stats[ref]["coverage"] = 0
                reference_stats[ref]["covered_regions"] = []
                reference_stats[ref]['numreads'] = 0
                reference_stats[ref]['accession'] = ref
            else:

                # Calculate average read length
                avg_read_length = math.ceil(stats["total_length"] / stats["total_reads"]) if stats["total_reads"] > 0 else 0
                # Calculate average base quality
                avg_baseq = stats["sum_baseq"] / stats["count_baseq"] if stats["count_baseq"] > 0 else 0
                # Calculate average mapping quality
                avg_mapq = stats.get("sum_mapq", 0) / stats.get("count_mapq", 0) if stats.get("count_mapq", 0) > 0 else 0



                # Compute coverage regions
                events = []
                for start, end in stats.get("read_positions", 0):
                    events.append( (start, 1) )   # Read start
                    events.append( (end, -1) )    # Read end

                # Sort events by position
                events.sort()
                coverage_regions = []
                current_depth = 0
                last_pos = 0

                for pos, delta in events:
                    if pos > last_pos and current_depth > 0:
                        coverage_regions.append( (last_pos, pos, current_depth) )
                    current_depth += delta
                    last_pos = pos

                # Calculate average depth
                sum_depth_bases = sum(depth * (end - start) for start, end, depth in coverage_regions)
                avg_depth = sum_depth_bases / reference_lengths.get(ref,0) if reference_lengths.get(ref,0)   > 0 else 0

                # Calculate coverage breadth
                covered_length = sum(end - start for start, end, depth in coverage_regions)
                coverage_breadth = covered_length / reference_lengths.get(ref, 0) if reference_lengths.get(ref, 0) > 0 else 0

                # Assign computed statistics
                reference_stats[ref]["avg_read_length"] = avg_read_length
                reference_stats[ref]["meanbaseq"] = avg_baseq
                reference_stats[ref]["meanmapq"] = avg_mapq
                reference_stats[ref]['meandepth'] = avg_depth
                reference_stats[ref]["coverage"] = coverage_breadth
                reference_stats[ref]["covered_regions"] = coverage_regions
                reference_stats[ref]['numreads'] = stats["total_reads"]
                reference_stats[ref]['accession'] = ref

                # Clean up intermediate fields
                del reference_stats[ref]["sum_baseq"]
                del reference_stats[ref]["count_baseq"]
                del reference_stats[ref]["sum_mapq"]
                del reference_stats[ref]["count_mapq"]
                del reference_stats[ref]["read_positions"]
                del reference_stats[ref]["unique_read_ids"]
                del reference_stats[ref]["total_reads"]
                del reference_stats[ref]["total_length"]

    print(f"Processed {total_reads} reads from {len(reference_lengths)} references in {time.time()-start_time} seconds.")
    # for ref, stats in reference_stats.items():
    #     print(ref)
    #     print("\t", stats.get('meanmapq'))
    #     print("\t", stats.get('meanbaseq'))
    #     print("\t", stats.get('meandepth'))
    #     print("\t", stats.get('coverage'))
    #     print("\t", stats.get('numreads'))
    bam_file.close()
    i=0
    return reference_stats, aligned_reads, total_reads
def main():
    args = parse_args()
    inputfile = args.input
    pathogenfile = args.pathogens
    cfig = None

    # Write to file in a newline format
    weights = {
        'mapq_score': args.mapq_weight,
        'disparity_score': args.disparity_score_weight,
        'hmp_weight': args.hmp_weight,
        'gini_coefficient': args.gini_weight,
        "breadth_weight": args.breadth_weight,
        "minhash_weight": args.minhash_weight,
        'siblings_score': 0,
        'diamond_identity': args.diamond_identity_weight,
        "k2_disparity_weight": args.k2_disparity_weight,
    }

    """
    # Final Score Calculation

    The final score is calculated as the average of four tests, each contributing 25% (or 0.25 proportion) to the overall score. The four tests are based on the following criteria:

    1. **Gini Coefficient (25% of final score)**:
        - The Gini coefficient is a measure of inequality and is expected to be between 0 and 1. It directly contributes to 25% of the final score.
        - Expected Range: [0, 1]

    2. **MAPQ Score (25% of final score)**:
        - The MAPQ (Mapping Quality) score is derived from the Phred score using the formula:
            [
            text{MAPQ} = 10 times log_{10}(text{Phred score})
            ]
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
        print("Importing k2 file")
        k2_mapping = import_k2_file(args.k2)
    import pandas as pd
    comparison_df =  pd.DataFrame()
    alignments_to_remove = defaultdict(set)
    if args.minhash_weight > 0:



        if args.comparisons and os.path.exists(args.comparisons):

            # if ends with csv
            if args.comparisons.endswith('.csv'):
                comparison_df = pd.read_csv(args.comparisons, sep=',')
            elif args.comparisons.endswith('.tsv'):
                comparison_df = pd.read_csv(args.comparisons, sep='\t')
            else:
                comparison_df = pd.read_excel(args.comparisons)
        else:
            print("Generating conflict regions info")
            if not args.output_dir:
                args.output_dir = os.path.dirname(args.output)
            alignments_to_remove, comparison_df = determine_conflicts(
                output_dir = args.output_dir,
                input_bam = args.input,
                matrix = args.matrix,
                min_threshold = args.min_threshold,
                abu_file = args.abu_file,
                fasta_files = args.fasta,
                min_similarity_comparable = args.min_similarity_comparable,
                use_variance = args.use_variance,
                apply_ground_truth = args.apply_ground_truth,
                sigfile = args.sigfile,
                bedfile = args.bedgraph,
                reference_signatures = args.reference_signatures,
                scaled = args.scaled,
                kmer_size = args.kmer_size,
                FAST_MODE=args.fast,
                filtered_bam_create=args.filtered_bam,
                sensitive=args.sensitive,
                cpu_count=args.cpu_count,
                jump_threshold = args.jump_threshold,
                gap_allowance=args.gap_allowance
            )
        if args.failed_reads:
            alignments_to_remove = defaultdict(set)
            with open(args.failed_reads, 'r') as f:
                for line in f:
                    lines = line.strip().split("\t")
                    ref, read_id = lines
                    alignments_to_remove[ref].add((read_id))
                f.close()
    if args.only_filter:
        print(f"Only filtering, exiting...")
        return

    if args.config:
        import yaml
        with open(args.config, "r") as file:
            cfig = yaml.safe_load(file)
        file.close()

    if comparison_df is not None and not comparison_df.empty:
        # set index to Reference column
        comparison_df.set_index('Reference', inplace=True)
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
        alignments_to_remove=alignments_to_remove,
    )


    assembly_to_accession = defaultdict(set)
    taxid_to_accession = defaultdict(int)
    if args.match and os.path.exists(matcher):
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
                if accession in reference_hits and not reference_hits[accession].get('taxid', None):
                    reference_hits[accession]['taxid'] = taxid
                reference_hits[accession]['name'] = name
                reference_hits[accession]['assembly'] = assembly
                taxid_to_accession[accession] = taxid
                if accession not in assembly_to_accession[assembly]:
                    assembly_to_accession[assembly].add(accession)
                i+=1

        f.close()
    taxid_to_parent = defaultdict(int)
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
                    taxid_to_parent[taxid] = species_taxid
                    # fine value where assembly == accession from reference_hits
                    if accession in assembly_to_accession:
                        for acc in assembly_to_accession[accession]:
                            if acc in reference_hits:

                                reference_hits[acc]['isSpecies'] = isSpecies
                                if args.compress_species:
                                    reference_hits[acc]['toplevelkey'] = species_taxid
                                else:
                                    reference_hits[acc]['toplevelkey'] = taxid
                                reference_hits[acc]['species_taxid'] = species_taxid
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
    seen = dict()
    for k, v in reference_hits.items():
        v['species_taxid'] = taxid_to_parent.get(v.get('taxid'), v.get('taxid'))
        if not v.get('toplevelkey'):
            # check if taxid is present, and if so then set it to that, and if that isnt then set it to the accession
            if v.get('taxid', None):
                v['toplevelkey'] = v['taxid']
            else:
                v['toplevelkey'] = k

    if args.compress_species:
        # Need to convert reference_hits to only species species level in a new dictionary
        for key, value in reference_hits.items():
            key = value.get('accession')
            valtoplevel = value.get('toplevelkey', key)
            valkey = key
            if value.get('numreads', 0)> 0 and  value.get('meandepth', 0) > 0:
                final_format[valtoplevel][valkey]= value
    else:
        # We don't aggregate, so do final format on the organism name only
        for key, value in reference_hits.items():
            valtoplevel = value.get('toplevelkey', key)
            valkey = key
            # if value['numreads'] > 0 and value['meandepth'] > 0:
            final_format[valtoplevel][valkey] = value
    # Dictionary to store aggregated species-level data
    species_aggregated = {}


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
    # Step 2: Define a function to calculate disparity for each organism
    # Define a function to calculate disparity with softer variance influence
    # Step 2: Define a function to dynamically dampen variance based on the proportion of reads
    def z(numreads, total_reads, variance_reads, k=1000):
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
    i=0
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
                    'species_taxid': data.get('species_taxid', None),
                    'mapqs': [],
                    'lengths': [],
                    'depths': [],
                    "taxids": [],
                    'minhash_reductions': [],
                    "accs": [],
                    'assemblies': [],
                    "coverages": [],
                    "breadth_total": 0,
                    "prevalence_disparity": 0,
                    "coeffs": [],
                    'baseqs': [],
                    "isSpecies": True if  args.compress_species  else data.get('isSpecies', False),
                    'strainslist': [],
                    'covered_regions': 0,
                    'name': data.get('name', ""),  # Assuming the species name is the same for all strains
            }



            try:

                # if accession isn't NC_042114.1 skip
                if len(data.get('covered_regions', [])) > 0 :
                    # gini_strain2 = gini_coverage_spread(data.get('covered_regions', []),  data['length'])
                    gini_strain = getGiniCoeff(data['covered_regions'],
                                               data['length'], alpha=args.alpha,
                                               reward_factor=args.reward_factor,
                                               beta=args.dispersion_factor
                                    )
                    # if "Vaccinia" in data['name'] or "Monkey" in data['name']:
                    #     print("\t",data.get('name'), gini_strain, gini_strain2)

                else:
                    gini_strain = 0
                # get Δ All% column value if it exists, else set to 0
                # col_stat = 'Sum Rs %'
                col_stat2 = 'Δ All%'
                col_stat = 'Δ^-1 Breadth'
                if not comparison_df.empty:
                    # Ensure 'accession' is a string and matches the index type
                    accession = str(data['accession']).strip()
                    if accession in comparison_df.index:
                        c1 = float(comparison_df.loc[accession, col_stat])
                        c2 = 1+(float(comparison_df.loc[accession, col_stat2]) / 100)
                        comparison_value = min(
                            1,
                            (c1+c2) / 2
                        )
                        # weight it so that any value less than 0.9 is even lower by getting the log value
                        # if comparison_value < 0.75:
                        #     # Using an exponent of 3.3 will reduce 0.64 to roughly 0.23.
                        #     comparison_value = comparison_value ** 3.3

                        data['comparison'] =  ( comparison_value  )
                    else:
                        data['comparison'] = 1
                species_aggregated[top_level_key]['minhash_reductions'].append(data.get('comparison', 1))
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
        aggregated_data['meanminhash_reduction'] = calculate_weighted_mean(aggregated_data['minhash_reductions'],numreads)
        # Step 4: Calculate the disparity for this organism
        aggregated_data['disparity'] = calculate_disparity(sum(numreads), total_reads, variance_reads)
        # calcualte the total coverage by summing all lengths and getting covered bases
        total_length = sum(aggregated_data['lengths'])
        covered_bases = sum([float(x) * float(y) for x, y in zip(aggregated_data['lengths'], aggregated_data['coverages'])])
        aggregated_data['breadth_total'] = covered_bases / total_length if total_length > 0 else 0
        aggregated_data['total_length'] = total_length
        k2_reads = 0
        if args.k2:
            taxid = aggregated_data.get('taxid', None)
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
        taxids = value.get('taxids', None)
        if args.k2:
            # Define the rank code you're looking for, for example "G" for Genus
            rank_mapping_match = None
            if args.parent_k2_match:
                rank_mapping_match = args.parent_k2_match
            parent_k2_reads_total  = 0
            sibling_k2_reads_total = []
            k2_numreads_total  = 0
            for taxid in taxids:
                # Step 1: Get the 'parents' list for the current key from k2_mapping
                parents = k2_mapping.get(taxid, {}).get('parents', [])
                my_rank = k2_mapping.get(taxid, {}).get('rank', 'U')
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
                parent_k2_reads_total += sum([species.get('clades_covered', 0) for species in related_species])
                sibling_k2_reads_total.extend(sibling_k2_reads)

                # Step 5: Run disparity calculations between the species-specific numreads and genus_k2_reads
                # Get the numreads from k2_mapping for this key
                k2_numreads_total += value.get('k2_numreads', 0)
            # get index where taxid is in k2_mapping
            # sort sibling_k2_reads_total by the taxid
            # Ensure the list is sorted based on the second element of each tuple
            sibling_k2_reads_total = sorted(sibling_k2_reads_total, key=lambda x: x[0], reverse=True)
            # Get the index of the value
            try:
                idx = [x[1] for x in sibling_k2_reads_total].index(value['taxids'][0])
            except (ValueError, IndexError):
                # If value['taxids'][0] is not found, handle it appropriately
                idx = -1

            # Ensure the list has more than 1 element to avoid division by zero
            if len(sibling_k2_reads_total) > 1 and idx >= 0:
                value['k2_disparity'] = 1 - (idx / (len(sibling_k2_reads_total) - 1))
            else:
                # Default to 0 if the list is empty, has one element, or idx is invalid
                value['k2_disparity'] = 0


        ## Test: Mapq of NT,
        # get all mapq scores
        normalized_mapq = normalize_mapq(value.get('meanmapq', 0), max_mapq, min_mapq)
        value['alignment_score'] = normalized_mapq


        ## Test: Min threshold of CDs found  - are there coding regions at the min amount?

    pathogens = import_pathogens(pathogenfile)

    # for values of pathogens, klust the ones with high_cons != ''
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
    if args.readcount:
        total_reads = float(args.readcount)

    print(f"Total Read Count in Entire Sample pre-filter: {total_reads}")
    if args.sampletype:
        sampletype = body_site_map(args.sampletype.lower())
    else:
        sampletype = "Unknown"
    if args.hmp:
        if sampletype == "Unknown":
            body_sites = []
        else:
            body_sites = [sampletype]
        dists, site_counts = import_distributions(
            args.hmp,
            "tax_id",
            body_sites
        )
        # for each entry in species_aggregated, get the taxid
        for key, value in species_aggregated.items():
            abus = []
            for body_site in body_sites:
                taxids = value.get('taxids', [])
                if args.compress_species:
                    taxids = [key]
                for taxid in taxids:
                    if not taxid:
                        continue
                    try:
                        # get the key where it is the (taxid, body_site)
                        if (int(taxid), body_site) in dists:
                            abus.append(
                                dict(
                                    norm_abundance = dists[(int(taxid), body_site)].get('norm_abundance', 0),
                                    std = dists[(int(taxid), body_site)].get('std', 0),
                                    mean = dists[(int(taxid), body_site)].get('mean', 0),
                                )
                            )
                    except Exception as e:
                        print(f"Error in taxid lookup for hmp: {e}")
            percent_total_reads_observed = 100*(sum(value.get('numreads', [])) / total_reads) if total_reads > 0 else 0
            sum_abus_expected = sum([x.get('mean',0)/100 for x in abus])
            sum_norm_abu = sum([x.get('norm_abundance',0) for x in abus])
            percent_total_reads_expected =  sum_abus_expected * total_reads
            stdsum = sum([x.get('std',0) for x in abus])
            zscore = ((percent_total_reads_observed -sum_norm_abu)/ stdsum)  if stdsum > 0 else 3
            # set zscore max 3 if greater
            # if zscore > 3:
            #     zscore = 3.1
            value['zscore'] = zscore
            percentile = norm.cdf(zscore)

            value['hmp_percentile'] = percentile

    else:
        for key, value in species_aggregated.items():
            value['zscore'] = 3
            value['hmp_percentile'] = 100
    final_scores = calculate_scores(
        aggregated_stats=species_aggregated,
        pathogens=pathogens,
        sample_name=args.samplename,
        sample_type = sampletype,
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
    #           "\n\t\tDisparity:", entry['disparity_score'],
    #           "\n\t\tBreadth:", entry['breadth_total'],
    #           "\n\t\tBreadth2:", ((entry['breadth_total'])**2),
    #           "\n\t\tBreadthLog:", math.log2(2-entry['breadth_total']),
    #           "\n\t\tBreadthsqrt:", math.sqrt(1-entry['breadth_total'])
    #     )
    header = [
        "Detected Organism",
        "Specimen ID",
        "Sample Type",
        "% Reads",
        "# Reads Aligned",
        "% Aligned Reads",
        "Coverage",
        'HHS Percentile',
        "IsAnnotated",
        "AnnClass",
        "Microbial Category",
        'High Consequence',
        "Taxonomic ID #",
        "Status",
        "Gini Coefficient",
        "Mean BaseQ",
        "Mean MapQ",
        "Mean Coverage",
        "Mean Depth",
        "isSpecies",
        "Pathogenic Subsp/Strains",
        "K2 Reads",
        "Parent K2 Reads",
        "MapQ Score",
        "Disparity Score",
        "Minhash Score",
        "Diamond Identity",
        "K2 Disparity Score",
        "Siblings score",
        "Breadth Weight Score",
        "TASS Score"
    ]
    write_to_tsv(output, final_scores, "\t".join(header))
    if cfig:
        idxes_header = []
        # retrive each of the attributes from the items in the imported dict. Assign to a new list with header
        for item in cfig['items']:
            header = item.get('label')
            key = item.get('key')
            fr = item.get('from')
            if fr == "report": # pull it from final_scores variable
                idx_header = headers.index(header)
                idxes_header.append(idx_header)



    if (args.gt and args.optimize):

        # Usage
        best_weights = optimize_weights(
            aggregated_stats=species_aggregated,
            pathogens=pathogens,
            sample_name="Test_Sample",
            sample_type="Unknown",
            total_reads=total_reads,
            aligned_total=aligned_total,
            initial_weights=[1, 1, 1, 1, 1, 1, 1, 0],
            max_iterations=args.max_iterations
        )
        print("Optimized Weights:")
        for key, value in best_weights.items():
            print(f"\t{key}: {value:.3f}")
        # convert all np floats





def format_non_zero_decimals(number):
    # Convert number from the scientific notation to the tenth digit
    number = "{:.10f}".format(number).rstrip('0').rstrip('.')
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
        initial_weights =[1, 1, 1, 1, 1, 1, 0],



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
            'disparity_score': w[1],
            'gini_coefficient': w[2],
            'breadth_score': w[3],
            'minhash_score': w[4],
            'siblings_score': w[5],
            'k2_disparity_weight': w[6],
            'diamond_identity': w[7],
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
            'disparity_score': result.x[1],
            'gini_coefficient': result.x[2],
            'breadth_score': result.x[3],
            'minhash_score': result.x[4],
            'siblings_score': result.x[5],
            'k2_disparity_weight': result.x[6],
            'diamond_identity': result.x[7],
        }
        # Extract the optimized weights
        optimized_weights = {
            'mapq_score': result.x[0],
            'disparity_score': result.x[1],
            'gini_coefficient': result.x[2],
            'breadth_score': result.x[3],
            'minhash_score': result.x[4],
            'siblings_score': result.x[5],
            'k2_disparity_weight': result.x[6],
            'diamond_identity': result.x[7],
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

        # for k, v in final_per_reference_scores.items():
        #     print(f"\t{k}: {v:.4f}")
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
    for ref, count in aggregated_stats.items():
        strainlist = count.get('strainslist', [])
        is_pathogen = "Unknown"
        callfamclass = ""
        annClass = "None"
        refpath = pathogens.get(ref)
        pathogenic_sites = refpath.get('pathogenic_sites', []) if refpath else []
        # check if the sample type is in the pathogenic sites
        direct_match = False
        high_cons = False
        def pathogen_label(ref):
            is_pathogen = "Unknown"
            isPathi = False
            direct_match = False
            callclass = ref.get('callclass', "N/A")
            high_cons = ref.get('high_cons', False)
            # pathogenic_sites = ref.get('pathogenic_sites', [])
            commensal_sites = ref.get('commensal_sites', [])
            if sample_type in pathogenic_sites:
                if callclass != "commensal":
                    is_pathogen = callclass.capitalize()
                    isPathi = True
                else:
                    is_pathogen = "Potential"
                direct_match = True
            elif sample_type in commensal_sites:
                is_pathogen = "Commensal"
                direct_match = True
            elif callclass and callclass != "":
                is_pathogen = callclass.capitalize() if callclass else "Unknown"
                isPathi = True
            return is_pathogen, isPathi, direct_match, high_cons

        formatname = count.get('name', "N/A")

        if refpath:
            is_pathogen, isPathi, direct_match, high_cons = pathogen_label(refpath)
            is_annotated = "Yes"
            status = refpath.get('status', "N/A")
            formatname = refpath.get('name', formatname)
            if is_pathogen == "Commensal":
                callfamclass = "Commensal Listing"
        else:
            is_annotated = "No"
            taxid = count.get(ref, "")
            status = ""
        if direct_match:
            annClass = "Direct"
        else:
            annClass = "Derived"
            # take the species_taxid and see if it is a pathogen
            if ref != count.get('species_taxid'):
                ref_spec = pathogens.get(count.get('species_taxid'), None)
                if ref_spec:
                    is_pathogen_spec, _,_, high_cons_spec = pathogen_label(ref_spec)
                    is_pathogen = is_pathogen_spec
                    high_cons = high_cons_spec
        listpathogensstrains = []
        fullstrains = []

        callclasses = set()
        if strainlist:
            pathogenic_reads = 0
            merged_strains = defaultdict(dict)

            for x in strainlist:
                strainname = x.get('strainname', None)
                taxid = x.get('taxid', None)
                if taxid:
                    keyx = taxid
                elif not taxid and strainname:
                    keyx = strainname
                else:
                    keyx = None
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
                    if pathstrain.get('callclass') not in ["commensal", "Unknown", 'unknown', '', None]:
                        callclasses.add(pathstrain.get('callclass').capitalize())
                    # annClassN = pathstrain.get('callclass', "Unknown")
                    # if the pathstrain is high consequence set high_cons to True
                    # callclasses.add(annClassN.capitalize())
                    if pathstrain.get('high_cons', False):
                        high_cons = True
                    if sample_type in pathstrain.get('pathogenic_sites', []) or sample_type == pathstrain.get('general_classification', ''):
                        pathogenic_reads += x.get('numreads', 0)
                        percentreads = f"{x.get('numreads', 0)*100/aligned_total:.1f}" if aligned_total > 0 and x.get('numreads', 0) > 0 else "0"
                        listpathogensstrains.append(f"{x.get('strainname', 'N/A')} ({percentreads}%)")

            if callfamclass == "" or len(listpathogensstrains) > 0:
                callfamclass = f"{', '.join(listpathogensstrains)}" if listpathogensstrains else ""
        if len(callclasses) > 0:
            # if Primary set is_pathogen to primary, if opposite set to opportunistic if potential set to potential
            if "Primary" in callclasses:
                is_pathogen = "Primary"
            elif "Opportunistic" in callclasses:
                is_pathogen = "Opportunistic"
            elif "Potential" in callclasses:
                is_pathogen = "Potential"
        breadth_total = count.get('breadth_total', 0)
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
        # Example usage:
        breadth = count.get('breadth_total')  # your raw value between 0 and 1
        log_weight_breadth = 1-logarithmic_weight(breadth)
        # get log weight of the breadth_total
        count['log_weight_breadth'] = log_weight_breadth
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
                log_breadth_weight = count.get('log_weight_breadth', 0),
                total_reads=sum(count['numreads']),
                reads_aligned=countreads,
                hmp_percentile = count.get('hmp_percentile', 0),
                zscore= count.get('zscore', 0),
                meancoverage=count.get('meancoverage', 0),
                breadth_total = breadth_total,
                total_length = count.get('total_length', 0),
                is_annotated=is_annotated,
                is_pathogen=is_pathogen,
                ref=ref,
                status=status,
                annClass=annClass,
                high_cons = high_cons,
                pathogenic_reads = pathogenic_reads,
                gini_coefficient=count.get('meangini', 0),
                meanbaseq=count.get('meanbaseq', 0),
                meanmapq=count.get('meanmapq', 0),
                alignment_score=count.get('alignment_score', 0),
                disparity_score=count.get('normalized_disparity', 0),
                covered_regions=count.get('covered_regions', 0),
                siblings_score=count.get('raw_disparity', 0),
                k2_reads=count.get('k2_numreads', 0),
                k2_disparity=count.get('k2_disparity', 0),
                gtcov=count.get('gtcov', 0),
                minhash_score=count.get('meanminhash_reduction', 0),

            )
        )
    return final_scores

def write_to_tsv(output_path, final_scores, header):
    total=0
    fulltotal=0
    with open (output_path, 'w') as file:
        file.write(f"{header}\n")
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
            diamond_identity = entry.get('diamond_identity', 0)
            k2_reads = entry.get('k2_reads', 0)
            k2_parent_reads = entry.get('k2_parent_reads', 0)
            mapq_score = entry.get('alignment_score', 0)
            disparity_score = entry.get('disparity_score', 0)
            siblings_score = entry.get('siblings_score', 0)
            tass_score = entry.get('tass_score', 0)
            minhash_score = entry.get('minhash_score', 0)
            listpathogensstrains = entry.get('listpathogensstrains', [])
            countreads = entry.get('reads_aligned', 0)
            breadth_of_coverage = entry.get('breadth_total', 0)
            aligned_total = entry.get('total_reads', 0)
            pathogenic_reads = entry.get('pathogenic_reads', 0)
            percent_total = entry.get('percent_total', 0)
            k2_disparity_score = entry.get('k2_disparity', 0)
            hmp_percentile = entry.get('hmp_percentile', 0)
            log_breadth_weight = entry.get('log_breadth_weight', 0)
            high_conse = entry.get('high_cons', False)
            # if plasmid uper or lower case doesnt matter matches then skip
            if "plasmid" in formatname.lower():
                continue
            if  (is_pathogen == "Primary" or is_pathogen=="Potential" or is_pathogen=="Opportunistic") and tass_score >= 0.990  :
                print(f"Reference: {ref} - {formatname}, High Cons?: {high_conse}")
                print(f"\tIsPathogen: {is_pathogen}")
                print(f"\tPathogenic Strains: {listpathogensstrains}")
                percentreads = f"{100*pathogenic_reads/aligned_total:.1f}" if aligned_total > 0 and pathogenic_reads>0 else 0
                print(f"\tPathogenic SubStrain Reads: {pathogenic_reads} - {percentreads}%")
                print(f"\tAligned Strains:")
                for f in entry.get('fullstrains', []):
                    print(f"\t\t{f}")
                print(f"\tTotal reads: {entry.get('total_reads', 0)}")
                print(f"\tGini Conf: {entry.get('gini_coefficient', 0):.4f}")
                print(f"\tDiamond Identity: {entry.get('diamond_identity', 0):.2f}")
                print(f"\tAlignment Score: {entry.get('alignment_score', 0):.2f}")
                print(f"\t# Reads Aligned: {entry.get('reads_aligned', 0)}")
                print(f"\tDisparity Score: {entry.get('disparity_score', 0):.2f}")
                print(f"\tK2 Disparity Score: {entry.get('k2_disparity', 0):.2f}")
                print(f"\tCoverage (Mean%): {entry.get('meancoverage', 0):.2f}")
                print(f"\tBreadth: {entry.get('breadth_total', 0):.2f}")
                print(f"\tMinhash score: {entry.get('minhash_score', 0):.2f}")
                print(f"\tDepth (Mean): {entry.get('meandepth', 0):.2f}")
                print(f"\tGT Coverage: {entry.get('gtcov', 0):.2f}")
                print(f"\tFinal Score: {entry.get('tass_score', 0):.2f}")
                print(f"\tCovered Regions: {entry.get('covered_regions', 0)}")
                print(f"\tK2 Reads: {entry.get('k2_reads', 0)}")
                print(f"\tHMP ZScore: {entry.get('zscore', 0)}")
                print(f"\tHMP Percentile: {entry.get('hmp_percentile', 0)}")
                print(f"\tLogWeightBreath: {entry.get('log_breadth_weight', 0)}")
                print()
                total+=1
        # header = "Detected Organism\tSpecimen ID\tSample Type\t% Reads\t# Reads Aligned\t% Aligned Reads\tCoverage\tIsAnnotated\tPathogenic Sites\tMicrobial Category\tTaxonomic ID #\tStatus\tGini Coefficient\tMean BaseQ\tMean MapQ\tMean Coverage\tMean Depth\tAnnClass\tisSpecies\tPathogenic Subsp/Strains\tK2 Reads\tParent K2 Reads\tMapQ Score\tDisparity Score\tMinhash Score\tDiamond Identity\tK2 Disparity Score\tSiblings score\tTASS Score\n"
            file.write(
                f"{formatname}\t{sample_name}\t{sample_type}\t{percent_total}\t{countreads}\t{percent_aligned}\t{breadth_of_coverage:.2f}\t{hmp_percentile:.2f}\t"
                f"{is_annotated}\t{annClass}\t{is_pathogen}\t{high_conse}\t{ref}\t{status}\t{gini_coefficient:.2f}\t"
                f"{meanbaseq:.2f}\t{meanmapq:.2f}\t{meancoverage:.2f}\t{meandepth:.2f}\t{isSpecies}\t{callfamclass}\t"
                f"{k2_reads}\t{k2_parent_reads}\t{mapq_score:.2f}\t{disparity_score:.2f}\t{minhash_score:.2f}\t"
                f"{diamond_identity:.2f}\t{k2_disparity_score:.2f}\t{siblings_score:.2f}\t{log_breadth_weight}\t{tass_score:.2f}\n"
            )
            fulltotal+=1
    print(f"Total pathogenic orgs: {total}, Total entire: {fulltotal}")
if __name__ == "__main__":
    sys.exit(main())
