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
import pandas as pd

import time
from distributions import import_distributions, body_site_map
import argparse
import re
import csv
import math
import os
from conflict_regions import determine_conflicts, generate_ani_matrix
import pysam
import random
from ground_truth import build_ground_truth_metrics_df, optimize_weights
from optimize_weights import annotate_aggregate_dict, compute_scores_per, calculate_aggregate_scores, calculate_classes, calculate_normalized_groups, compute_tass_score, pathogen_label, normalize_category
from map_taxid import load_taxdump, load_names
from utils import taxid_to_rank, calculate_var

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
        "--orgcol",
        metavar="ORGCOL",
        type=int,
        default=3, required =False,
        help="Index of the organism name column, used as a fallback to match organism name if taxid is missing for merged accessions & taxid no present",
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
        "--taxdump",
        metavar="TAXDUMP",
        default=None, required=False,
        help="Taxdump directory containing nodes.dmp and names.dmp files",
    )
    parser.add_argument(
        "--compare_references", default=False,  help="Compress species to species level",  action='store_true'
    )
    parser.add_argument(
        "-k",  "--rank",  help="Compress species to species level",  default=None, type=str
    )
    parser.add_argument(
        "--sensitive", default=False,  help="Use sensitive mode to detect greater array of variants",  action='store_true'
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
        default=0.0,
        help="value of weight for disparity reads vs other organisms in final TASS Score",
    )
    parser.add_argument(
        "--alpha",
        metavar="MAPQWEIGHT",
        type=float,
        default=1,
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
        type=bool,
        default=None,
        help="You want to provide a gt file to see how well you fare with the detections",
    )
    parser.add_argument("--optimize", action="store_true",
                        help="Optimize breadth/minhash/disparity/gini weights using BAM-derived TP/FP labels.")
    parser.add_argument("--optimize_maxiter", type=int, default=80,
                        help="Max iterations for the global stage (differential evolution).")
    parser.add_argument("--optimize_local_maxiter", type=int, default=200,
                        help="Max iterations for the local refinement stage (SLSQP).")
    parser.add_argument("--optimize_seed", type=int, default=1,
                        help="RNG seed for optimizer reproducibility.")
    parser.add_argument("--optimize_pos_weight", type=float, default=0.5,
                        help="Relative weight of true-positive term in the objective.")
    parser.add_argument("--optimize_neg_weight", type=float, default=0.5,
                        help="Relative weight of false-positive term in the objective.")
    parser.add_argument("--optimize_reg", type=float, default=0.01,
                        help="L2 regularization strength to keep weights near the starting values.")
    parser.add_argument("--optimize_report", type=str, default=None,
                        help="Optional path to write a TSV report of TP/FP counts and scores for the best weights.")
    parser.add_argument(
        '--breadth_weight',
        metavar="BREADTHSCORE",
        type=float,
        default=0.45,
        help="value of weight for breadth of coverage in final TASS Score",
    )
    parser.add_argument(
        "--minhash_weight",
        metavar="MINHASHSCORE",
        type=float,
        default=0.05,
        help="value of weight for minhash signature reduction in final TASS Score",
    )
    parser.add_argument(
        "--k2_disparity_weight",
        metavar="DISPARITYSCOREWEIGHT",
        type=float,
        default=0.0,
        help="value of weight for disparity of k2 and alignment in final TASS Score",
    )
    parser.add_argument("--disparity_weight", type=float, default=0.0,
        help="Weight applied to disparity_score in tass_score (optimized if --optimize).")
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
        default=0.55,
        help="value of weight for gini coefficient in final TASS Score",
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
    parser.add_argument(
        "-y",
        "--microbert",
        required=False,
        metavar="MICROBERT",
        help="OPTIONAL: Microbert Predictions on downsampled dataset. Contains columns for each rank, e.g. superkingdom, phylum, class, order, family, genus, species and avg, median, and std probablity columns",
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
    parser.add_argument('--minmapq', required=False, type=int, default=0, help="Minimum mapping quality")
    parser.add_argument(
        "--filtered_bam", default=False,  help="Create a filtered bam file of a certain name post sourmash sigfile matching..", type=str
    )
    parser.add_argument(
        "--report_confusion_xlsx",
        action="store_true",
        help="Write an XLSX report containing confusion-matrix stats and remaining false-positive reads."
    )
    parser.add_argument(
        "--confusion_xlsx",
        required=False,
        default=None,
        help="Output XLSX path. If not set, will write to <output_dir>/alignment_confusion_report.xlsx"
    )

    return parser.parse_args(argv)

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

def count_reference_hits(bam_file_path,alignments_to_remove=None, reference_lengths={}):
    """
    Count the number of reads aligned to each reference in a BAM file.

    Args:
    bam_file_path (str): Path to the BAM file.

    Returns:
    dict: A dictionary with reference names as keys and counts of aligned reads as values.
    """
    # Initialize a dictionary to hold the count of reads per reference
    reference_lengths = {}
    unaligned = 0
    aligned_reads = 0
    total_reads = 0

    with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
        # get total reads

        reference_stats = defaultdict(dict)
        for ref in bam_file.header.references:
            # get average read length
            reference_lengths[ref] = bam_file.get_reference_length(ref)
            reference_stats[ref] = dict(
                length = reference_lengths.get(ref, 0),
                depths =  defaultdict(int),
                key = ref,
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
                covered_regions = [],  # Store regions covered by reads
            )
        # Open BAM file
        total_reads = 0
        # check if the alignment was paired end or single end
        start_time = time.time()
        for read in bam_file.fetch():

            if read.is_unmapped:
                unaligned += 1
                continue
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
        for ref, _ in reference_lengths.items():
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
    bam_file.close()
    return reference_stats, aligned_reads, total_reads
def main():
    args = parse_args()
    inputfile = args.input
    pathogenfile = args.pathogens
    # Write to file in a newline format
    weights = {
        'mapq_score': args.mapq_weight,
        'disparity_score': args.disparity_score_weight,
        'hmp_weight': args.hmp_weight,
        'gini_weight': args.gini_weight,
        "breadth_weight": args.breadth_weight,
        "minhash_weight": args.minhash_weight,
        'siblings_score': 0,
        'diamond_identity': args.diamond_identity_weight,
        "k2_disparity_weight": args.k2_disparity_weight,
    }
    total_weight = sum(weights.values())
    if total_weight != 1:
        for key in weights:
            if total_weight != 0:
                weights[key] = weights[key] / total_weight
            else:
                weights[key] = 0
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

    matcher = args.match
    k2_mapping = dict()
    i =0
    if args.k2:
        print("Importing k2 file")
        k2_mapping = import_k2_file(args.k2)
        print("K2 mapping imported")
    comparison_df =  pd.DataFrame()
    alignments_to_remove = defaultdict(set)
    if args.fasta:
        fastas = set(args.fasta)
    if args.match and os.path.exists(matcher):
        # if args.output_dir/organism_ani_matrix.csv doesnt exist then perform it:
        if not args.matrix or not os.path.exists(args.matrix):
            print("generating ani matrix")
            generate_ani_matrix(
                fasta_files = fastas if args.fasta else [],
                matchfile = args.match,
                matchfile_accession_col = args.accessioncol,
                matchfile_taxid_col = args.taxcol,
                matchfile_desc_col = args.namecol,
                output_dir = args.output_dir if args.output_dir else os.path.dirname(args.output),
            )
        else:
            print("Ani matrix already exists across references")
    if args.minhash_weight > 0:
        if args.comparisons and os.path.exists(args.comparisons):
            # if ends with csv
            if args.comparisons.endswith('.csv'):
                comparison_df = pd.read_csv(args.comparisons, sep=',')
            elif args.comparisons.endswith('.tsv'):
                comparison_df = pd.read_csv(args.comparisons, sep='\t')
            else:
                comparison_df = pd.read_excel(args.comparisons)
            # remove the "Total" row if it exists
            comparison_df = comparison_df[comparison_df['Reference'] != 'Total']
        else:
            print("Generating conflict regions info")
            if not args.output_dir:
                args.output_dir = os.path.dirname(args.output)
            alignments_to_remove, comparison_df = determine_conflicts(
                output_dir = args.output_dir,
                input_bam = args.input,
                min_threshold = args.min_threshold,
                fasta_files = fastas if args.fasta else [],
                use_variance = args.use_variance,
                sigfile = args.sigfile,
                bedfile = args.bedgraph,
                scaled = args.scaled,
                kmer_size = args.kmer_size,
                FAST_MODE=args.fast,
                filtered_bam_create=args.filtered_bam,
                sensitive=args.sensitive,
                cpu_count=args.cpu_count,
                jump_threshold = args.jump_threshold,
                gap_allowance=args.gap_allowance,
                compare_to_reference_windows=args.compare_references
            )
            # import the file args.output_dir/region_comparisons.csv
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
    # set reflengths from the Reference Length file if provided, otherwise get the accession lengths from the bam file
    orgn_lengths = dict()
    if comparison_df is not None and not comparison_df.empty and "Reference Length" in comparison_df.columns:
        # set index to Reference column
        print("Using comparison file to get reference lengths")
        comparison_df.set_index('Reference', inplace=True)
        for ref in comparison_df.index.unique():
            # if Reference Length column exists set length else length is 0
            length = 0
            if 'Reference Length' in comparison_df.columns:
                length = comparison_df.loc[ref, 'Reference Length']
            orgn_lengths[ref] = length
    elif args.fasta:
        print("Using fasta files to get reference lengths")
        from Bio import SeqIO
        for fasta in fastas:
            for record in SeqIO.parse(fasta, "fasta"):
                orgn_lengths[record.id] = len(record.seq)
    else:
        print("No fasta file provided, getting reference lengths from BAM file")
        # get reference lengths from bam file
        with pysam.AlignmentFile(inputfile, "rb") as bam_file:
            for ref in bam_file.header.references:
                orgn_lengths[ref] = bam_file.get_reference_length(ref)
        bam_file.close()
    print(len(orgn_lengths), ": total references found")
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
                            contigs= line[1],
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

    mmbert_dict = dict()
    if args.microbert:
        # import the tsv as a dictionary
        mmbert = pd.read_csv(args.microbert, sep='\t', header=0)
        # set the taxid col to str
        mmbert['taxid'] = mmbert['taxid'].astype(str)
        mmbert_dict = mmbert.set_index('taxid').T.to_dict()

    pathogens = import_pathogens(pathogenfile)

    if args.match and os.path.exists(matcher):
        accindex = args.accessioncol
        nameindex = args.namecol
        taxcol = args.taxcol
        orgindex = args.orgcol

        with open(matcher, "r") as f:
            for i, line in enumerate(f):
                line = line.rstrip("\n")
                if i == 0:  # header
                    continue

                splitline = line.split("\t")
                if not splitline or len(splitline) <= accindex:
                    continue

                accession = splitline[accindex].strip() or None
                taxid     = splitline[taxcol].strip() if len(splitline) > taxcol else None
                name      = splitline[nameindex].strip() if len(splitline) > nameindex else None
                organism = splitline[orgindex].strip() if len(splitline) > orgindex else None

                if not accession:
                    continue

                # Optional: only fill taxid if it exists and if this accession is already in reference_hits
                if accession in reference_hits:
                    if taxid and not reference_hits[accession].get("taxid"):
                        reference_hits[accession]["taxid"] = taxid

                    if organism:
                        reference_hits[accession]["name"] = organism
                    elif name:
                        reference_hits[accession]["name"] = name
        f.close()
    taxdump, taxdump_names = {}, {}
    if args.taxdump and os.path.exists(os.path.join(args.taxdump, "nodes.dmp")):
        taxdump = load_taxdump(os.path.join(args.taxdump, "nodes.dmp"))
    if args.taxdump and os.path.exists(os.path.join(args.taxdump, "names.dmp")):
        taxdump_names = load_names(os.path.join(args.taxdump, "names.dmp"))
    # args.rank = False
    # get the ranks for taxid: 198214
    acc_to_parent = dict()
    if args.rank:
        wanted_rank = args.rank  # for example:  species
        for acc, hit in reference_hits.items():
            taxid = hit.get("taxid")
            if not taxid:
                hit['key'] = acc
                hit["toplevelkey"] = acc
                hit["toplevelname"] = hit.get('name', acc)
                hit['strainname'] = hit.get('name', acc)
                continue
            top = taxid_to_rank(taxid, taxdump, wanted_rank )
            acc_to_parent[acc] = taxid
            hit['key'] = taxid
            hit["toplevelkey"] = top if top else taxid  # fallback to itself if rank not found
            hit["rank"] = wanted_rank  # optional: record what you tried
            hit['strainname'] = hit.get('name', '')
            hit['name'] = taxdump_names.get(taxid, hit.get('name', ''))
            hit["toplevelname"] = taxdump_names.get(top, hit.get("name", ""))
    else:
        for acc, hit in reference_hits.items():
            taxid = hit.get("taxid")
            if taxid:
                hit["toplevelkey"] = taxid
                acc_to_parent[acc] = taxid
                hit['key'] = taxid
            else:
                hit["toplevelkey"] = acc
                acc_to_parent[acc] = acc
                hit['key'] = acc
            hit['strainname'] = hit.get('name', '')
            hit["toplevelname"] = taxdump_names.get(hit['key'], hit.get("name", ""))
    species_to_all_accs = defaultdict(set)
    all_readscounts = [x['numreads'] for x in reference_hits.values()]
    total_reads = sum(all_readscounts)
    print(f"Total aligned reads: {total_reads}")
    variance_reads = calculate_var(all_readscounts)
    print(f"\n\tVariance of reads: {variance_reads}")
    if args.sampletype:
        sampletype = body_site_map(args.sampletype.lower())
    else:
        sampletype = "Unknown"
    if sampletype == "sterile":
        # set ALL pathogens in pathogens dict to pathogenic
        for k, v in pathogens.items():
            # if v['callclass'] == "commensal":
            v['callclass'] = "primary (sterile)"

    if args.hmp:
        if sampletype == "Unknown":
            body_sites = []
        else:
            body_sites = [sampletype]
        dists, _ = import_distributions(
            args.hmp,
            "tax_id",
            body_sites
        )
    else:
        dists = {}

    for acc, data in reference_hits.items():
        if not data.get('organism'):
            # try to get organism from taxdump names
            taxid = data.get('taxid')
            if taxid and taxid in taxdump_names:
                data['organism'] = taxdump_names[taxid]
            else:
                data['organism'] = data.get('name', acc)
        top = data.get('key')
        species_to_all_accs[top].add(acc)
        data = compute_scores_per(
            data = data,
            reward_factor = args.reward_factor,
            dispersion_factor = args.dispersion_factor,
            alpha = args.alpha,
            comparison_df = comparison_df,
            fallback_top = top,
            total_reads = total_reads,
        )
    strain_summary = calculate_normalized_groups(
        hits=reference_hits,
        group_field="key",
        reads_key="numreads",
    )

    if args.optimize:
        report_weights = optimize_weights(
            input_bam = args.input,
            final_json = strain_summary,
            accession_to_taxid = acc_to_parent,
            breadth_weight = args.breadth_weight,
            minhash_weight = args.minhash_weight,
            gini_weight = args.gini_weight,
            # disparity_weight = args.disparity_score_weight,
            alpha = args.alpha,
            sampletype = sampletype,
            optimize_pos_weight = args.optimize_pos_weight,
            optimize_neg_weight = args.optimize_neg_weight,
            optimize_reg = args.optimize_reg,
            optimize_seed = args.optimize_seed,
            optimize_maxiter = args.optimize_maxiter,
            optimize_local_maxiter = args.optimize_local_maxiter,
            optimize_report = args.optimize_report,
        )
        weights = report_weights.get("best_weights", weights)
        print("Optimized weights changed: ")
        for k, v in weights.items():
            print(f"\t{k}: {v}")
    # Add sample_name and pathogen annotations to each strain
    for k, data in strain_summary.items():
        data['sample_name'] = args.samplename

        # Annotate strain with pathogen info
        taxid = data.get('key') or data.get('taxid') or k
        rest = calculate_classes(
            rec = data,
            ref = taxid,
            pathogens =pathogens,
            sample_type = sampletype,
            taxdump = taxdump,
        )
        data.update(rest)
        # try to match taxid to names_dict else name is k
        group_reads = [
            dict(reads=x['numreads'], key=x.get('key'))
            for _, x in strain_summary.items()
            if x.get('toplevelkey') == data.get('toplevelkey')
        ]
        calculate_aggregate_scores(
            data = data,
            hmp_dists = dists,
            body_sites = [sampletype],
            k2_mapping = k2_mapping,
            sampletype = sampletype,
            mmbert_dict = mmbert_dict,
            group_reads = group_reads,
            dmnd = dmnd,
        )
        data['tass_score'] = compute_tass_score(
            data = data,
            weights = weights,
        )
    # : Define a function to calculate disparity for each organism
    i=0

    aggregate_dict = calculate_normalized_groups(
        hits=strain_summary,
        group_field="toplevelkey",
        reads_key="numreads",
    )
    # iterate through aggregate_dict, make the accession_to_key dict
    for k, v in aggregate_dict.items():
        for acc in species_to_all_accs.get(k, []):
            acc_to_parent[acc] = k

    # for k, v in strain_summary.items():
    #     print(f"{i}\t{k}\t{v.get('name', '')}\t{100*v.get('tass_score')}")
    # print("\n________________________________\n")
    # exit()
    for _, data in aggregate_dict.items():
        # add all the strains to the strains list from strain_summary if data['toplevelkey'] matches
        data['members'] = [
            x for _, x in strain_summary.items()
            if x.get('toplevelkey') == data.get('toplevelkey')
        ]


        data['name'] = data.get('toplevelname', None)
        data['key'] = data.get('toplevelkey', None)
        data['sample_name'] = args.samplename

        group_reads = [
            dict(reads=x['numreads'], key=x.get('key'))
            for _, x in strain_summary.items()
            if x.get('toplevelkey') == data.get('toplevelkey')
        ]
        calculate_aggregate_scores(
            data = data,
            hmp_dists = dists,
            body_sites = [sampletype],
            k2_mapping = k2_mapping,
            sampletype = sampletype,
            mmbert_dict = mmbert_dict,
            group_reads = group_reads,
        )
        data['tass_score'] = compute_tass_score(
            data = data,
            weights = weights,
        )
    # for k, v in aggregate_dict.items():
    #     print(f"{k}\t{v.get('name', '')}\t{v.get('key', '')}\t{v.get('toplevelkey', '')}\t{v.get('toplevelname', '')}\t{v.get('microbial_category', 0)}\t{v.get('accessions')}\t{v.get('tass_score', 0)}")
    # exit()
    # for values of pathogens, klust the ones with high_cons != ''
    # Next go through the BAM file (inputfile) and see what pathogens match to the reference, use biopython
    # to do this
    if args.min_reads_align:
        # filter the reference_hits based on the minimum number of reads aligned
        print(f"Filtering for minimum reads aligned: {args.min_reads_align}")
        aggregate_dict = {k: v for k, v in aggregate_dict.items() if (v['numreads']) >= int(args.min_reads_align)}
    # if sum of vals in weights isnt 1 then normalize to 1
    if args.readcount:
        total_reads = float(args.readcount)

    print(f"Total Read Count in Entire Sample pre-filter: {total_reads}")

    import json
    def write_to_json(output_path, obj):
        with open(output_path, "w") as f:
            json.dump(obj, f, indent=2)
    final_json = annotate_aggregate_dict(
        aggregate_dict=aggregate_dict,
        pathogens=pathogens,
        sample_type=sampletype,
        taxdump = taxdump,
    )
    for v in final_json:
        v['sampletype'] = sampletype
        v['total_reads'] = total_reads
    # for data in final_json:
    #     print(f"{data.get('name', 'N/A')} ({data.get('key', 'N/A')}):")
    #     print(f"\tGini Coefficient: {data.get('gini_coefficient', 'N/A')},"
    #         f"\n\tMAPQ Score: {data.get('meanmapq', 'N/A')},",
    #         f"\n\t# Reads: {data.get('numreads', 'N/A')},",
    #         f"\n\tBreadth: {data.get('coverage', 'N/A')},",
    #         f"\n\tMinHash Reduction: {data.get('minhash_reduction', 'N/A')},",
    #         f"\n\tBreadth Score: {data.get('breadth_log_score', 'N/A')},",
    #         f"\n\tReads K2: {data.get('k2_reads', 'N/A')},",
    #         f"\n\tK2 Disparity Score: {data.get('k2_disparity_score', 'N/A')},"
    #         f"\n\tDiamond Identity: {data.get('diamond_identity', 'N/A')},"
    #         f"\n\tDisparity Score: {data.get('disparity', 'N/A')},"
    #         f"\n\tHMP Percentile: {data.get('hmp_percentile', 'N/A')},"
    #         f"\n\tAccessions: {data.get('accessions', 0)},"
    #         f"\n\tMicrobert Prob: {data.get('mmbert', 'N/A')},"
    #         f"\n\tMembers : {len(data.get('members', []))},"
    #         f"\n\tHigh Cons: {data.get('high_cons', False)},"
    #         f"\n\tPathogenic Strains: {[x.get('name') for x in data.get('members', []) if x.get('is_pathogen') in ['Primary', 'Potential', 'Opportunistic']]},"
    #         f"\n\tIs Pathogen: {data.get('is_pathogen', 'N/A')},"
    #         f"\n\tAnnotation Class: {data.get('annClass', 'N/A')},"
    #         f"\n\tMicrobial Category: {data.get('microbial_category', 'N/A')},"
    #         f"\n\tFinal Score: {data.get('tass_score', 'N/A')}\n\n+________________________________\n"
    #     )
    write_to_json(args.output.replace(".tsv", ".json"), final_json)


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
            breadth_of_coverage = entry.get('coverage', 0)
            aligned_total = entry.get('total_reads', 0)
            pathogenic_reads = entry.get('pathogenic_reads', 0)
            percent_total = entry.get('percent_total', 0)
            mmbert_proportion = entry.get('mmbert', None)
            k2_disparity_score = entry.get('k2_disparity', 0)
            hmp_percentile = entry.get('hmp_percentile', 0)
            log_breadth_weight = entry.get('log_breadth_weight', 0)
            high_conse = entry.get('high_cons', False)
            mmbert_model = entry.get('mmbert_model', None)
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
                print(f"\tBreadth: {entry.get('coverage', 0):.2f}")
                print(f"\tMinhash score: {entry.get('minhash_score', 0):.2f}")
                print(f"\tDepth (Mean): {entry.get('meandepth', 0):.2f}")
                print(f"\tGT Coverage: {entry.get('gtcov', 0):.2f}")
                print(f"\tFinal Score: {entry.get('tass_score', 0):.2f}")
                print(f"\tCovered Regions: {entry.get('covered_regions', 0)}")
                print(f"\tK2 Reads: {entry.get('k2_reads', 0)}")
                print(f"\tHMP ZScore: {entry.get('zscore', 0)}")
                print(f"\tHMP Percentile: {entry.get('hmp_percentile', 0)}")
                print(f"\tLogWeightBreath: {entry.get('log_breadth_weight', 0)}")
                print(f"\tMicrobeRT Proportion: {entry.get('mmbert', 'N/A')}")
                print(f"\tMicrobeRT Model: {entry.get('mmbert_model', 'N/A')}")
                print()
                total+=1
            file.write(
                f"{formatname}\t{sample_name}\t{sample_type}\t{percent_total}\t{countreads}\t{percent_aligned}\t{breadth_of_coverage:.2f}\t{hmp_percentile:.2f}\t"
                f"{is_annotated}\t{annClass}\t{is_pathogen}\t{high_conse}\t{ref}\t{status}\t{gini_coefficient:.2f}\t"
                f"{meanbaseq:.2f}\t{meanmapq:.2f}\t{meancoverage:.2f}\t{meandepth:.2f}\t{callfamclass}\t"
                f"{k2_reads}\t{k2_parent_reads}\t{mapq_score:.2f}\t{disparity_score:.2f}\t{minhash_score:.2f}\t"
                f"{diamond_identity:.2f}\t{k2_disparity_score:.2f}\t{siblings_score:.2f}\t{log_breadth_weight}\t"
                f"{tass_score:.2f}\t{mmbert_proportion}\t{mmbert_model}\n"
            )
            fulltotal+=1
    print(f"Total pathogenic orgs: {total}, Total entire: {fulltotal}")
if __name__ == "__main__":
    sys.exit(main())
