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
from body_site_normalization import normalize_body_site
import argparse
import json as _json
import re
import csv
import math
import os
from conflict_regions import determine_conflicts, generate_ani_matrix
import pysam
import random
from ground_truth import optimize_weights, compute_tp_fp_counts_by_taxid
from optimize_weights import annotate_aggregate_dict, compute_scores_per, calculate_aggregate_scores, calculate_classes, calculate_normalized_groups, compute_tass_score, pathogen_label, normalize_category, breadth_score_sigmoid, getGiniCoeff, load_control_data, compute_control_comparison, find_missing_positive_controls
from map_taxid import load_taxdump, load_names
from utils import taxid_to_rank, calculate_var, load_matchfile


def _find_sample_only_organisms(final_json, control_index):
    """Identify organisms present in the sample but absent from the control index.

    Returns a set of toplevelkey IDs that exist in the sample but not in the
    control (e.g. insilico) data.
    """
    sample_tlks = set()
    for grp in final_json:
        tlk = str(grp.get('toplevelkey', grp.get('key', '')))
        if tlk:
            sample_tlks.add(tlk)

    ctrl_tlks = set(control_index.get('by_toplevelkey', {}).keys())
    return sample_tlks - ctrl_tlks


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
        type=int,
        help="Index of the column in mapfile (if specified) to match to the reference accession. 0 index start",
    )
    parser.add_argument(
        "--orgcol",
        metavar="ORGCOL",
        type=int,
        default=2, required =False,
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
        default=3,
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
        "-k",  "--rank",  help="Specify the taxonomic rank to group all entries on (e.g., species, genus, family)",  default=None, type=str
    )
    parser.add_argument(
        "--subrank",  help="Specify the taxonomic sub-rank to group all entries on. Default is species. Entries in the final report will have 1 species per row. If set to None or a non-standard ranking, then ignored. If the subrank and rank are equal, subrank is ignore as well.", choices=["phylum", "order", "class", "genus", "species", "strain", "none"], default="species", type=str
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
        "--platform",
        type=str,
        default="unknown",
        help="Sequencing platform (e.g. illumina, nanopore, pacbio, ion_torrent). "
             "Stored in output JSON metadata. Default: unknown.",
    )
    parser.add_argument(
        "--workflow_revision",
        type=str,
        default=None,
        help="Nextflow workflow revision string stored in output JSON metadata.",
    )
    parser.add_argument(
        "--commit_id",
        type=str,
        default=None,
        help="Pipeline commit/version ID stored in output JSON metadata.",
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
        default=0.00,
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
    parser.add_argument("--optimize_entropy", type=float, default=0.1,
                        help="Entropy regularization strength to prevent a single weight from dominating. "
                             "Higher values produce more balanced weights. 0 = disabled. Recommended: 0.05-0.3.")
    parser.add_argument("--optimize_min_weight", type=float, default=0.00,
                        help="Minimum weight floor for each scoring component during optimization. "
                             "Prevents any weight from going to 0. Range: 0.0-0.15. Default: 0.05.")
    parser.add_argument("--optimize_tp_target", type=float, default=None,
                    help="Optional target for TP mean score (0-1). If set, adds a penalty when TP mean is below target.")
    parser.add_argument("--optimize_tp_target_weight", type=float, default=0.0,
                    help="Penalty weight for TP target shortfall. Higher values push TP mean upward.")
    parser.add_argument("--optimize_tp_target_scope", type=str, default="taxa",
                    choices=["taxa", "reads", "hybrid"],
                    help="Which TP mean to target: taxa, reads, or hybrid.")
    parser.add_argument("--optimize_youden_weight", type=float, default=0.0,
                    help="Penalty weight to maximize Youden J (TP retention - FP retention). Higher values push FP down." )
    parser.add_argument("--optimize_fp_cutoff", type=float, default=None,
                    help="Optional TASS cutoff to penalize FP retention at/above this threshold (0-1).")
    parser.add_argument("--optimize_fp_cutoff_weight", type=float, default=0.0,
                    help="Penalty weight for FP retention at the fixed cutoff.")
    parser.add_argument("--optimize_curve_scope", type=str, default="taxa",
                    choices=["taxa", "reads", "hybrid"],
                    help="Which retention curves to use for Youden/cutoff penalties.")
    parser.add_argument("--optimize_tp_floor", type=float, default=0.7,
                    help="TP score floor: all TP organisms should score above this. Default: 0.7")
    parser.add_argument("--optimize_tp_floor_weight", type=float, default=0.0,
                    help="Penalty weight for TP score floor hinge loss. Pushes individual TP organisms "
                         "above the floor. Recommended: 1.0-5.0. 0 = disabled.")
    parser.add_argument("--optimize_fp_ceiling", type=float, default=0.15,
                    help="FP score ceiling: all FP organisms should score below this. Default: 0.15")
    parser.add_argument("--optimize_fp_ceiling_weight", type=float, default=0.0,
                    help="Penalty weight for FP score ceiling hinge loss. Pushes individual FP organisms "
                         "below the ceiling. Recommended: 1.0-5.0. 0 = disabled.")
    parser.add_argument("--optimize_separation_weight", type=float, default=0.0,
                    help="Multi-threshold retention shape penalty. Penalizes TP dropout and FP leakage "
                         "across thresholds [0.1, 0.2, 0.3, 0.5, 0.7]. Recommended: 0.5-2.0. 0 = disabled.")
    parser.add_argument("--optimize_weight_prior", type=str, default=None,
                    help="JSON string or file path specifying target weight values the optimizer "
                         "should bias toward, e.g. '{\"breadth_weight\": 0.4}'. "
                         "Used with --optimize_weight_prior_lambda to control pull strength.")
    parser.add_argument("--optimize_weight_prior_lambda", type=float, default=0.0,
                    help="Strength of the weight-prior pull. Higher values more aggressively "
                         "bias the optimizer toward the target weights. Recommended: 0.5-5.0. 0 = disabled.")
    parser.add_argument("--optimize_report", type=str, default=None,
                        help="Optional path to write a TSV report of TP/FP counts and scores for the best weights.")
    parser.add_argument("--optimize_granularity", type=str, default="subkey",
                        choices=["key", "subkey", "toplevelkey"],
                        help="Preferred granularity level for selecting optimized weights. "
                             "The optimizer runs at all levels but prefers this one unless another "
                             "level has a significantly lower loss. Also controls which "
                             "best_threshold from --thresholds_json is used for the TASS cutoff "
                             "in the report. Default: subkey.")
    parser.add_argument(
        '--breadth_weight',
        metavar="BREADTHSCORE",
        type=float,
        default=0.26,
        help="value of weight for breadth of coverage in final TASS Score",
    )
    parser.add_argument(
        "--minhash_weight",
        metavar="MINHASHSCORE",
        type=float,
        default=0.29,
        help="value of weight for minhash signature reduction in final TASS Score",
    )
    parser.add_argument(
        "--gini_weight",
        metavar="GINIWEIGHT",
        type=float,
        default=0.45,
        help="value of weight for gini coefficient in final TASS Score",
    )
    parser.add_argument(
        "--k2_disparity_weight",
        metavar="DISPARITYSCOREWEIGHT",
        type=float,
        default=0.0,
        help="value of weight for disparity of k2 and alignment in final TASS Score",
    )
    parser.add_argument("--disparity_weight", type=float, default=0.01,
        help="Weight applied to disparity_score in tass_score (optimized if --optimize).")
    parser.add_argument(
        "--diamond_identity_weight",
        metavar="DIAMONDIDENTITYWEIGHT",
        type=float,
        default=0.0,
        help="value of weight for disparity of diamond_identity in final TASS Score",
    )
    parser.add_argument(
        "--hmp_weight",
        metavar="HMPWEIGHT",
        type=float,
        default=0.00,
        help="value of weight for hmp abundance in final TASS Score",
    )
    parser.add_argument(
        "--plasmid_bonus_weight",
        type=float,
        default=0.19,
        help="Additive TASS bonus for strains with strong plasmid coverage "
             "relative to sibling strains in the same species. Applied outside "
             "the normalized weight pool. 0 = disabled. Default: 0.05",
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
    parser.add_argument("--enable_matrix", required=False, action='store_true', help="Enable loading the ANI matrix file")
    parser.add_argument(
        "--ani_threshold",
        type=float,
        default=0.95,
        help="ANI threshold above which two organisms are considered highly similar (default: 0.95). "
             "Used to annotate each member in the output JSON with a 'high_ani_matches' list.",
    )
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
    parser.add_argument('--minmapq', required=False, type=int, default=7,
                    help="MAPQ threshold for high-confidence reads. Reads with MAPQ >= this value "
                         "are counted as high-quality. The fraction of such reads scales breadth score.")
    parser.add_argument('--mapq_breadth_power', required=False, type=float, default=2.0,
                    help="Power exponent for MAPQ-adjusted breadth. breadth *= highmapq_fraction^power. "
                         "Higher values penalize low-MAPQ organisms more aggressively. "
                         "Default: 2.0 (e.g. 10%% high-MAPQ reads → breadth scaled by 0.01).")
    parser.add_argument(
        "--filtered_bam", default=False,  help="Create a filtered bam file of a certain name post sourmash sigfile matching..", type=str
    )
    parser.add_argument(
        "--report_metrics",
        action="store_true",
        help="Compute and output TP/FP metrics per organism using ground-truth read IDs. "
             "Produces a JSON file and stdout table with per-organism TASS scores, TP/FP status, "
             "overall F1/precision/recall, and threshold-based percentile analysis showing "
             "what %% of TPs and FPs pass at each TASS cutoff.",
    )
    parser.add_argument(
        "--thresholds_json",
        type=str,
        default=None,
        help="Path to a JSON file containing per-sample-type best weights and thresholds. "
             "If provided, the weights (breadth_weight, gini_weight, minhash_weight, hmp_weight, "
             "disparity_weight) are automatically set based on the normalized sample type. "
             "Falls back to the 'all' category if the sample type is not found in the JSON.",
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
    parser.add_argument(
        "--taxid_removal_stats",
        action="store_true",
        help="In addition to the standard removal_stats.xlsx, output a taxid-aggregated "
             "removal_stats_by_taxid.xlsx where accessions are grouped by their taxid. "
             "Requires --match to provide accession-to-taxid mappings."
    )

    # ── Control sample arguments ────────────────────────────────────────────
    parser.add_argument(
        "--negative_controls",
        nargs="+",
        default=None,
        help="One or more JSON files (match_paths output format) from negative control "
             "samples.  Organisms whose TASS / reads are within fold-change of these "
             "controls are flagged as 'within_negative' bounds.",
    )
    parser.add_argument(
        "--positive_controls",
        nargs="+",
        default=None,
        help="One or more JSON files (match_paths output format) from positive control "
             "samples.  Used for spark-bar visualisation in the report.",
    )
    parser.add_argument(
        "--control_type",
        type=str,
        default=None,
        choices=["positive", "negative"],
        help="Mark THIS sample as a control.  Stored in output JSON metadata as "
             "'control_type'.  Omit for normal (non-control) samples.",
    )
    parser.add_argument(
        "--control_fold_threshold",
        type=float,
        default=2.0,
        help="Fold-change threshold for negative-control comparison.  If sample "
             "TASS / max(neg TASS) < this value, organism is flagged "
             "'within_negative'.  Default: 2.0.",
    )
    parser.add_argument(
        "--missing_pos_levels",
        nargs="+",
        default=["toplevelkey"],
        choices=["toplevelkey", "key", "subkey"],
        help="Hierarchy level(s) at which to detect missing positive controls. "
             "Default: toplevelkey.  Specify one or more of: toplevelkey, key, subkey.",
    )
    parser.add_argument(
        "--hide_missing_pos_controls",
        action="store_true",
        default=False,
        help="Suppress detection and output of missing positive control organisms.",
    )

    # ── In-silico control arguments ──────────────────────────────────────────
    parser.add_argument(
        "--insilico_controls",
        nargs="+",
        default=None,
        help="One or more JSON files (match_paths output format) from in-silico "
             "simulated control samples.  Used to compare simulated vs real "
             "TASS scores and read counts per organism.",
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
def normalize_weights(weights):
    total = sum(weights.values())

    if total == 0:
        # Avoid division by zero
        return {k: 0 for k in weights}

    return {k: v / total for k, v in weights.items()}
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

def count_reference_hits(bam_file_path,alignments_to_remove=None, reference_lengths={}, args={}):
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
    primary_counts = defaultdict(int)
    secondary_counts = defaultdict(int)

    # Pre-pass: count primary/secondary alignments per read to gate MAPQ=0 reads
    with pysam.AlignmentFile(bam_file_path, "rb") as bam_count:
        for read in bam_count.fetch(until_eof=True):
            if read.is_unmapped:
                continue
            if read.is_secondary or read.is_supplementary:
                secondary_counts[read.query_name] += 1
                continue
            primary_counts[read.query_name] += 1
    # for each of the reads, check which ones have more than 1 count
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
                count_highmapq = 0,  # reads with MAPQ >= threshold
                sum_mapq_filtered = 0,   # MAPQ sum for reads that pass the filter
                count_mapq_filtered = 0, # count of reads that pass the filter
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
            if read.is_secondary or read.is_supplementary:
                continue
            ref = read.reference_name
            total_reads += 1



            # Skip reads that are in alignments_to_remove if provided
            if alignments_to_remove and read.query_name in alignments_to_remove and ref in alignments_to_remove[read.query_name]:
                continue

            # Track MAPQ for the high-quality fraction (computed on ALL mapped reads)
            reference_stats[ref]["sum_mapq"] += read.mapping_quality
            reference_stats[ref]["count_mapq"] += 1
            if read.mapping_quality >= args.minmapq:
                reference_stats[ref]["count_highmapq"] += 1

            # ── Filter: skip reads below --minmapq for all downstream metrics ──
            # These reads don't count toward read totals, coverage, depth, or
            # base quality.  They ARE still counted for highmapq_fraction above
            # so the fraction reflects the full alignment picture.
            if read.mapping_quality < args.minmapq:
                # if read.reference_name == "NC_002695.2":
                #     print(read.mapping_quality, read)
                allow_low_mapq = (
                    read.mapping_quality == 0
                    and secondary_counts.get(read.query_name, 0) > 0
                )
                # if"NC_002695.2" in read.query_name:
                #     print(read.query_qualities)
                #     print(help(read), read.is_mapped)
                #     exit()
                #     print(read.mapping_quality, read.query_name, read.reference_name, allow_low_mapq)
                if not allow_low_mapq:
                    continue

            # Accumulate MAPQ only for reads that passed the filter
            # (so meanmapq reflects the actual reads used for coverage/depth)
            reference_stats[ref]["sum_mapq_filtered"] += read.mapping_quality
            reference_stats[ref]["count_mapq_filtered"] += 1

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
                reference_stats[ref]['highmapq_fraction'] = 0.0
            else:

                # Calculate average read length
                avg_read_length = math.ceil(stats["total_length"] / stats["total_reads"]) if stats["total_reads"] > 0 else 0
                # Calculate average base quality
                avg_baseq = stats["sum_baseq"] / stats["count_baseq"] if stats["count_baseq"] > 0 else 0
                # Calculate average mapping quality from reads that PASSED the
                # minmapq filter (or the allow_low_mapq exception).  Using all
                # reads would drag the mean down with skipped MAPQ=0 alignments.
                _smf = stats.get("sum_mapq_filtered", 0)
                _cmf = stats.get("count_mapq_filtered", 0)
                avg_mapq = _smf / _cmf if _cmf > 0 else 0



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
                # Fraction of reads with MAPQ >= minmapq threshold
                _n_mapped = stats.get("count_mapq", 0)
                _n_hiq = stats.get("count_highmapq", 0)
                reference_stats[ref]['highmapq_fraction'] = (
                    _n_hiq / _n_mapped if _n_mapped > 0 else 0.0)

                # Clean up intermediate fields
                del reference_stats[ref]["sum_baseq"]
                del reference_stats[ref]["count_baseq"]
                del reference_stats[ref]["sum_mapq"]
                del reference_stats[ref]["count_mapq"]
                del reference_stats[ref]["count_highmapq"]
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
    if args.disparity_weight not in (None, 0.0):
        disparity_w = args.disparity_weight
    else:
        disparity_w = args.disparity_score_weight
    if args.platform:
        args.platform = args.platform.lower()
    # ── Load per-sample-type thresholds JSON if provided ─────────────────────
    # Keys in the JSON are "{sampletype}|{platform}", e.g. "blood|illumina".
    # Lookup order (first match wins):
    #   1. {sampletype}|{platform}
    #   2. {sampletype}|all
    #   3. all|{platform}
    #   4. all|all
    # Platform normalisation: pacbio → ont; default platform is "illumina".

    thresholds_config = None
    if args.thresholds_json:
        with open(args.thresholds_json, 'r') as _tf:
            thresholds_config = _json.load(_tf)

        # ─ normalise sample type ─
        _norm_st = normalize_body_site(args.sampletype.lower()) if args.sampletype else "unknown"
        if _norm_st == "sterile":
            _norm_st = "blood"

        # ─ normalise platform ─
        _raw_plat = (args.platform or "unknown").strip().lower()
        _raw_plat = _raw_plat.replace("_", " ").replace("-", " ")
        if _raw_plat in ("pacbio", "pac bio", "pb"):
            _norm_plat = "ont"
        elif _raw_plat in (
            "illumina", "ilumina", "ill", "illum", "miseq", "hi seq", "hiseq",
            "nextseq", "nova seq", "novaseq"
        ):
            _norm_plat = "illumina"
        elif _raw_plat in (
            "ont", "nano", "nanopore", "oxford", "oxford nanopore",
            "oxford nan", "minion", "promethion", "gridion"
        ):
            _norm_plat = "ont"
        elif _raw_plat in ("unknown", ""):
            _norm_plat = "illumina"          # default platform
        else:
            _norm_plat = _raw_plat            # pass through as-is

        # ─ lookup with fallback chain ─
        _candidates = [
            f"{_norm_st}|{_norm_plat}",       # exact match
            f"{_norm_st}|all",                 # any platform for this sampletype
            f"all|{_norm_plat}",               # any sampletype for this platform
            "all|all",                         # universal fallback
        ]
        _st_key = None
        for _cand in _candidates:
            if _cand in thresholds_config:
                _st_key = _cand
                break
        if _st_key is None:
            # Last resort: grab the first key in the JSON
            _st_key = next(iter(thresholds_config))
            print(f"WARNING: No matching thresholds key found; falling back to '{_st_key}'")

        print(f"Thresholds JSON loaded. Sample type '{args.sampletype}' → '{_norm_st}', "
              f"platform '{args.platform}' → '{_norm_plat}', "
              f"using thresholds key: '{_st_key}'")

        _best_w = thresholds_config[_st_key].get("best_weights", {})
        # Override CLI weight defaults with JSON best weights
        if "breadth_weight" in _best_w:
            args.breadth_weight = _best_w["breadth_weight"]
        if "gini_weight" in _best_w:
            args.gini_weight = _best_w["gini_weight"]
        if "minhash_weight" in _best_w:
            args.minhash_weight = _best_w["minhash_weight"]
        if "hmp_weight" in _best_w:
            args.hmp_weight = _best_w["hmp_weight"]
        if "disparity_weight" in _best_w:
            disparity_w = _best_w["disparity_weight"]
            args.disparity_weight = disparity_w
        if "plasmid_bonus_weight" in _best_w:
            args.plasmid_bonus_weight = _best_w["plasmid_bonus_weight"]
        print(f"  Applied weights from JSON: breadth={args.breadth_weight:.6g}, "
              f"gini={args.gini_weight:.6g}, minhash={args.minhash_weight:.6g}, "
              f"hmp={args.hmp_weight:.6g}, disparity={disparity_w:.6g}, "
              f"plasmid_bonus={args.plasmid_bonus_weight:.6g}")

    weights = {
        'mapq_score': args.mapq_weight,
        'disparity_weight': disparity_w,
        'hmp_weight': args.hmp_weight,
        'gini_weight': args.gini_weight,
        "breadth_weight": args.breadth_weight,
        "minhash_weight": args.minhash_weight,
        'siblings_score': 0,
        'diamond_identity': args.diamond_identity_weight,
        "k2_disparity_score_weight": args.k2_disparity_weight,
    }
    # total_weight = sum(weights.values())
    weights = normalize_weights(weights)
    # Plasmid bonus is additive — added AFTER normalization so it doesn't
    # dilute the core weights.  Set to 0 to disable.
    weights['plasmid_bonus_weight'] = args.plasmid_bonus_weight
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
    # ani_df holds the taxid-level ANI DataFrame (index & cols are taxid strings).
    # It is populated either by generating the matrix from FASTA signatures or by
    # loading a pre-existing CSV.  ani_data is the nested dict form used for lookup.
    ani_df = None
    ani_data = {}  # ani_data[taxid1][taxid2] = float ANI value (0..1)

    if args.match and os.path.exists(matcher) and args.enable_matrix:
        # Determine where the default auto-generated matrix would be saved
        default_matrix_path = os.path.join(
            args.output_dir if args.output_dir else os.path.dirname(args.output),
            "organism_ani_matrix.csv",
        )
        matrix_path = args.matrix or default_matrix_path

        if args.matrix and os.path.exists(args.matrix):
            # Use the explicitly-provided pre-computed matrix CSV
            print(f"Loading pre-computed ANI matrix from {args.matrix}")
            ani_df = pd.read_csv(args.matrix, index_col=0)
        elif os.path.exists(default_matrix_path):
            # Auto-generated matrix from a previous run
            print(f"Loading cached ANI matrix from {default_matrix_path}")
            ani_df = pd.read_csv(default_matrix_path, index_col=0)
        else:
            print("Generating ANI matrix from FASTA signatures…")
            ani_df = generate_ani_matrix(
                fasta_files=fastas if args.fasta else [],
                matchfile=args.match,
                matchfile_accession_col=args.accessioncol,
                matchfile_taxid_col=args.taxcol,
                matchfile_desc_col=args.namecol,
                output_dir=args.output_dir if args.output_dir else os.path.dirname(args.output),
            )

        # Convert DataFrame to nested dict for O(1) lookup
        if ani_df is not None:
            for taxid1 in ani_df.index:
                t1 = str(taxid1)
                if t1 not in ani_data:
                    ani_data[t1] = {}
                for taxid2 in ani_df.columns:
                    val = ani_df.loc[taxid1, taxid2]
                    if pd.notna(val):
                        ani_data[t1][str(taxid2)] = float(val)
            print(f"ANI data loaded for {len(ani_data)} taxa (threshold {args.ani_threshold})")

    # Build an early accession→taxid map from the matchfile for use by
    # determine_conflicts (the full acc_to_key map is built later after
    # reference_hits are populated and taxdump is loaded).
    early_acc_to_taxid = {}
    early_taxid_to_desc = {}
    if args.match and os.path.exists(args.match):
        early_acc_to_taxid, early_taxid_to_desc, _ = load_matchfile(
            args.match,
            accession_col=args.accessioncol,
            taxid_col=args.taxcol,
            desc_col=args.namecol,
        )

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
                compare_to_reference_windows=args.compare_references,
                accession_to_taxid=early_acc_to_taxid if early_acc_to_taxid else None,
                taxid_to_name=early_taxid_to_desc if early_taxid_to_desc else None,
                taxid_removal_stats=args.taxid_removal_stats,
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
    reference_hits, _, total_reads = count_reference_hits(
        inputfile,
        alignments_to_remove=alignments_to_remove,
        args=args
    )
    # exit()

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
                splitline = line.split("\t")

                if i == 0 and len(splitline) > 0 and splitline[0] == "Acc":  # header
                    continue

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

                    # Store the raw description from the mapfile BEFORE
                    # overwriting 'name' with the cleaner organism column.
                    # The description often contains "plasmid pXXX" which
                    # we need for plasmid tagging downstream.
                    if name:
                        reference_hits[accession]["description"] = name
                    if organism:
                        reference_hits[accession]["name"] = organism
                    elif name:
                        reference_hits[accession]["name"] = name
        f.close()
    # ── Tag plasmid accessions from the description field ──────────────────
    # Use 'description' (raw mapfile name column) which preserves strings
    # like "E. coli ETEC H10407 plasmid p666, complete sequence".
    # Fall back to 'name' if description wasn't populated.
    _plasmid_count = 0
    for acc, hit in reference_hits.items():
        _desc = (hit.get('description', '') or hit.get('name', '') or '').lower()
        hit['is_plasmid'] = 'plasmid' in _desc
        if hit['is_plasmid']:
            _plasmid_count += 1
    if _plasmid_count:
        print(f"Tagged {_plasmid_count} accessions as plasmid")

    taxdump, taxdump_names = {}, {}
    if args.taxdump and os.path.exists(os.path.join(args.taxdump, "nodes.dmp")):
        taxdump = load_taxdump(os.path.join(args.taxdump, "nodes.dmp"))
    if args.taxdump and os.path.exists(os.path.join(args.taxdump, "names.dmp")):
        taxdump_names = load_names(os.path.join(args.taxdump, "names.dmp"))

    subrank = args.subrank

    if not subrank or subrank.lower() == "none" or subrank.lower() == "strain":
        subrank = None

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
            subrank_value = taxid_to_rank(taxid, taxdump, subrank) if subrank else taxid
            acc_to_parent[acc] = taxid
            hit['subkey'] = subrank_value
            hit['subkeyname'] =taxdump_names.get(subrank_value, hit.get('name', ''))
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
                hit['subkey'] = taxid
                hit['subkeyname'] = hit.get('name', '')
            else:
                hit["toplevelkey"] = acc
                acc_to_parent[acc] = acc
                hit['key'] = acc
                hit['subkey'] = acc
                hit['subkeyname'] = hit.get('name', '')
            hit['strainname'] = hit.get('name', '')
            hit["toplevelname"] = taxdump_names.get(hit['key'], hit.get("name", ""))

    # Build per-accession mappings for optimization and metrics at different granularities.
    acc_to_key = {}
    acc_to_subkey = {}
    acc_to_toplevelkey = {}
    for acc, hit in reference_hits.items():
        acc_to_key[acc] = str(hit.get("key", acc))
        acc_to_subkey[acc] = str(hit.get("subkey", acc))
        acc_to_toplevelkey[acc] = str(hit.get("toplevelkey", acc))
        acc_sub = re.sub(r"\.\d+$", "", acc)
        acc_to_key[acc_sub] = str(hit.get("key", acc))
        acc_to_subkey[acc_sub] = str(hit.get("subkey", acc))
        acc_to_toplevelkey[acc_sub] = str(hit.get("toplevelkey", acc))
    species_to_all_accs = defaultdict(set)
    all_readscounts = [x['numreads'] for x in reference_hits.values()]
    aligned_reads_total = sum(all_readscounts)
    total_reads = aligned_reads_total
    print(f"Total aligned reads: {aligned_reads_total}")
    variance_reads = calculate_var(all_readscounts)
    print(f"\n\tVariance of reads: {variance_reads}")
    if args.sampletype:
        sampletype = body_site_map(args.sampletype.lower())
        # Also check via body_site_normalization for richer sterile detection
        normalized_sampletype = normalize_body_site(args.sampletype.lower())
        if normalized_sampletype == "sterile" and sampletype != "sterile":
            sampletype = "sterile"
            print(f"  Body site normalization detected sterile site from '{args.sampletype}'")
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

    # ── Aggregate comparison_df to subkey (species) level for minhash scoring ──
    # The raw comparison_df is per-accession (contig), but minhash reduction
    # should reflect the species-level conflict picture: sum reads, weighted-avg
    # breadth ratios, so that a species with many contigs gets one composite score.
    subkey_comparison_df = pd.DataFrame()
    if comparison_df is not None and not comparison_df.empty:
        _cdf = comparison_df.copy()
        _cdf['subkey'] = _cdf.index.map(lambda a: acc_to_subkey.get(a, a))
        # Weighted aggregation: weight by Total Reads per accession
        _numeric_cols = [c for c in ['Total Reads', 'Pass Filtered Reads', 'Δ All',
                         'Reference Length'] if c in _cdf.columns]
        _wavg_cols = ['Δ All%', 'Δ^-1 Breadth', 'Breadth Original', 'Breadth New']
        # Force all numeric columns to numeric dtype (some may arrive as strings)
        for _nc in _numeric_cols + _wavg_cols:
            if _nc in _cdf.columns:
                _cdf[_nc] = pd.to_numeric(_cdf[_nc], errors='coerce').fillna(0)
        _grouped = _cdf.groupby('subkey')
        _sums = _grouped[_numeric_cols].sum()
        # Read-weighted averages for ratio/percentage columns
        _wavg_parts = {}
        for col in _wavg_cols:
            if col in _cdf.columns:
                _cdf[f'_w_{col}'] = _cdf[col] * _cdf['Total Reads']
                _wavg_parts[col] = _grouped[f'_w_{col}'].sum() / _sums['Total Reads'].replace(0, 1)
        _wavg_df = pd.DataFrame(_wavg_parts)
        subkey_comparison_df = pd.concat([_sums, _wavg_df], axis=1)
        print(f"Aggregated comparison_df: {len(comparison_df)} accessions → {len(subkey_comparison_df)} subkeys")

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
            comparison_df = subkey_comparison_df,
            fallback_top = top,
            total_reads = total_reads,
            mapq_breadth_power = args.mapq_breadth_power,
        )
    strain_summary = calculate_normalized_groups(
        hits=reference_hits,
        group_field="key",
        reads_key="numreads",
        mapq_breadth_power=args.mapq_breadth_power,
    )
    # for k, v in strain_summary.items():
    #     if v.get('toplevelname') == "Toxoplasma":
    #         print(v.get('toplevelname'), v.get('name'), v.get('minhash_reduction_raw'), v.get('minhash_reduction'), v.get('minhash_score'))
    # exit()

    # ── Plasmid disparity: score strains by relative plasmid presence ────────
    # Within each subkey (species), compare plasmid coverage across strains.
    # A strain whose plasmid(s) have notably better coverage/reads than
    # sibling strains' plasmids gets a higher plasmid_score (0–1).
    # The score also factors in the plasmid's own alignment quality
    # (breadth_log_score and gini_coefficient) so that a plasmid with
    # garbage coverage doesn't inflate the bonus.
    # Strains with no plasmid accessions get plasmid_score = 0 (neutral).
    #
    # 1) Collect plasmid stats per strain, grouped by subkey
    _subkey_plasmid_stats = defaultdict(dict)  # subkey → {strain_key → {reads, covered_bases, length, covered_regions}}
    for acc, hit in reference_hits.items():
        if not hit.get('is_plasmid', False):
            continue
        _sk = str(hit.get('subkey', ''))
        _strain = str(hit.get('key', acc))
        if _strain not in _subkey_plasmid_stats[_sk]:
            _subkey_plasmid_stats[_sk][_strain] = {
                'reads': 0, 'covered_bases': 0, 'length': 0, 'coverage': 0.0,
                'covered_regions': [], '_region_offset': 0,
            }
        _ps = _subkey_plasmid_stats[_sk][_strain]
        _ps['reads'] += hit.get('numreads', 0)
        _ps['covered_bases'] += hit.get('covered_bases', 0)
        # Accumulate covered_regions with offset for multi-plasmid strains
        for _start, _end, _depth in hit.get('covered_regions', []):
            _ps['covered_regions'].append((_start + _ps['_region_offset'], _end + _ps['_region_offset'], _depth))
        _ps['_region_offset'] += int(hit.get('length', 0) or 0)
        _ps['length'] += int(hit.get('length', 0) or 0)

    # 2) Compute coverage, breadth, gini, and rank within each subkey
    _strain_plasmid_score = {}  # strain_key → plasmid_score (0–1)
    for _sk, strains_map in _subkey_plasmid_stats.items():
        # Per-strain plasmid quality metrics
        for _strain, _ps in strains_map.items():
            _plen = max(1, _ps['length'])
            _ps['coverage'] = _ps['covered_bases'] / _plen
            # Compute breadth_log_score for this strain's plasmid(s)
            _ps['breadth'] = breadth_score_sigmoid(_ps['coverage'])
            # Compute gini from combined covered_regions
            if _ps['covered_regions']:
                _ps['gini'] = getGiniCoeff(
                    _ps['covered_regions'], _plen,
                    alpha=1.8, reward_factor=2, beta=0.5)
            else:
                _ps['gini'] = 0.0

        # Max values across siblings for relative comparison
        _max_cov = max((_ps['coverage'] for _ps in strains_map.values()), default=0.0)
        _max_reads = max((_ps['reads'] for _ps in strains_map.values()), default=0)
        _n_strains_with_plasmid = len(strains_map)
        if _max_cov <= 0 and _max_reads <= 0:
            for _strain in strains_map:
                _strain_plasmid_score[_strain] = 0.0
            continue

        for _strain, _ps in strains_map.items():
            # ── Absolute quality: does this plasmid have real coverage? ──
            # breadth (0–1): sigmoid of coverage fraction
            # gini (0–1): evenness of that coverage
            # Both must be decent — geometric mean ensures this.
            _breadth = float(_ps.get('breadth', 0.0))
            _gini = float(_ps.get('gini', 0.0))
            _quality = (_breadth * _gini) ** 0.5  # 0–1

            # ── Relative disparity: how does this strain compare to siblings? ──
            # Only meaningful when multiple strains have plasmids.
            if _n_strains_with_plasmid > 1:
                _cov_ratio = (_ps['coverage'] / _max_cov) if _max_cov > 0 else 0.0
                _read_ratio = (_ps['reads'] / _max_reads) if _max_reads > 0 else 0.0
                _disparity = min(1.0, 0.7 * _cov_ratio + 0.3 * _read_ratio)
            else:
                # Single strain with plasmid: disparity is neutral (1.0).
                # Score is driven entirely by absolute quality — a garbage
                # plasmid with poor breadth/gini will still score low.
                _disparity = 1.0

            # ── Final plasmid score = quality * disparity ──
            # quality gates the score: poor breadth or gini crushes it
            # disparity separates strains when multiple have plasmids
            _strain_plasmid_score[_strain] = min(1.0, _quality * _disparity)

    # 3) Apply plasmid_score to strain_summary
    _plasmid_scored = 0
    for k, data in strain_summary.items():
        _skey = str(data.get('key', k))
        if _skey in _strain_plasmid_score:
            data['plasmid_score'] = _strain_plasmid_score[_skey]
            data['has_plasmid'] = True
            _plasmid_scored += 1
        else:
            data['plasmid_score'] = 0.0
            data['has_plasmid'] = False
    if _plasmid_scored:
        print(f"Plasmid disparity scored for {_plasmid_scored} strains across {len(_subkey_plasmid_stats)} subkeys")

    # ── Pre-annotate strain_summary with microbial_category BEFORE optimization ──
    # This ensures optimize_weights() / build_metrics_df_from_final_json() can
    # read microbial_category for the debug JSON output.
    for k, data in strain_summary.items():
        taxid = data.get('key') or data.get('taxid') or k
        pre_ann = calculate_classes(
            rec=data,
            ref=taxid,
            pathogens=pathogens,
            sample_type=sampletype,
            taxdump=taxdump,
        )
        data.update(pre_ann)

    if args.optimize:
        # Parse weight prior: accepts JSON string or file path
        _weight_prior = None
        if args.optimize_weight_prior:
            _wp_raw = args.optimize_weight_prior.strip()
            if _wp_raw.startswith('{'):
                _weight_prior = _json.loads(_wp_raw)
            else:
                with open(_wp_raw, 'r') as _wpf:
                    _weight_prior = _json.load(_wpf)
            print(f"Weight prior targets: {_weight_prior} "
                  f"(lambda={args.optimize_weight_prior_lambda})")

        report_weights = optimize_weights(
            input_bam = args.input,
            final_json = strain_summary,
            accession_to_taxid = acc_to_key,
            accession_to_subkey = acc_to_subkey,
            accession_to_toplevelkey = acc_to_toplevelkey,
            breadth_weight = args.breadth_weight,
            minhash_weight = args.minhash_weight,
            gini_weight = args.gini_weight,
            disparity_weight = disparity_w,
            hmp_weight = args.hmp_weight,
            alpha = args.alpha,
            sampletype = sampletype,
            optimize_pos_weight = args.optimize_pos_weight,
            optimize_neg_weight = args.optimize_neg_weight,
            optimize_reg = args.optimize_reg,
            optimize_seed = args.optimize_seed,
            optimize_maxiter = args.optimize_maxiter,
            optimize_local_maxiter = args.optimize_local_maxiter,
            optimize_report = args.optimize_report,
            entropy_lambda = args.optimize_entropy,
            min_weight = args.optimize_min_weight,
            optimize_tp_target = args.optimize_tp_target,
            optimize_tp_target_weight = args.optimize_tp_target_weight,
            optimize_tp_target_scope = args.optimize_tp_target_scope,
            optimize_youden_weight = args.optimize_youden_weight,
            optimize_fp_cutoff = args.optimize_fp_cutoff,
            optimize_fp_cutoff_weight = args.optimize_fp_cutoff_weight,
            optimize_curve_scope = args.optimize_curve_scope,
            tp_score_floor = args.optimize_tp_floor,
            tp_floor_weight = args.optimize_tp_floor_weight,
            fp_score_ceiling = args.optimize_fp_ceiling,
            fp_ceiling_weight = args.optimize_fp_ceiling_weight,
            separation_weight = args.optimize_separation_weight,
            weight_prior = _weight_prior,
            weight_prior_lambda = args.optimize_weight_prior_lambda,
            plasmid_bonus_weight = args.plasmid_bonus_weight,
            prefer_granularity = args.optimize_granularity,
            platform=args.platform
        )
        best_weights = report_weights.get("best_weights") or {}
        weights.update(best_weights)
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
            total_reads = total_reads,
        )
        data['tass_score'] = compute_tass_score(
            data = data,
            weights = weights,
        )
    # : Define a function to calculate disparity for each organism
    i=0

    # ── Filter strains by min_reads_align BEFORE species-level aggregation ──
    # This ensures low-read strains don't inflate species aggregate stats
    # (coverage, breadth, Gini, etc.) and are excluded from members.
    _min_ra = int(args.min_reads_align) if args.min_reads_align else 0
    if _min_ra > 0:
        _before = len(strain_summary)
        strain_summary = {
            k: v for k, v in strain_summary.items()
            if v.get('numreads', 0) >= _min_ra
        }
        print(f"Filtered strains by min_reads_align={_min_ra}: {_before} → {len(strain_summary)}")

    aggregate_dict = calculate_normalized_groups(
        hits=strain_summary,
        group_field="toplevelkey",
        reads_key="numreads",
        mapq_breadth_power=args.mapq_breadth_power,
    )
    # iterate through aggregate_dict, make the accession_to_key dict
    for k, v in aggregate_dict.items():
        for acc in species_to_all_accs.get(k, []):
            acc_to_parent[acc] = k

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
            total_reads = total_reads,
        )
        data['tass_score'] = compute_tass_score(
            data = data,
            weights = weights,
        )
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
    # print final weights used:
    print("Final weights used for scoring: ")
    for k, v in weights.items():
        print(f"\t{k}: {v}")
    # ── ANI annotation ────────────────────────────────────────────────────────
    # For each member in every species group, attach a 'high_ani_matches' list
    # containing dicts {key, name, ani_pct} for all other taxa whose ANI with
    # this member exceeds args.ani_threshold.  create_report.py reads this list
    # directly so it does not need to load or parse the ANI matrix itself.
    if ani_data:
        for group in final_json:
            for member in group.get('members', []):
                member_key = str(member.get('key', ''))
                matches = []
                if member_key in ani_data:
                    for other_key, ani_val in ani_data[member_key].items():
                        if other_key == member_key:
                            continue
                        if ani_val >= args.ani_threshold:
                            matches.append({
                                'key': other_key,
                                'ani_pct': round(float(ani_val) * 100, 2),
                            })
                    matches.sort(key=lambda x: x['ani_pct'], reverse=True)
                member['high_ani_matches'] = matches

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
    # ── Control comparison annotation ────────────────────────────────────────
    # If negative and/or positive control JSONs were provided, load them and
    # annotate each group + member with fold-change metrics.  These are purely
    # informational annotations — the TASS score is NOT modified.
    _neg_index = load_control_data(args.negative_controls) if args.negative_controls else None
    _pos_index = load_control_data(args.positive_controls) if args.positive_controls else None

    if _neg_index or _pos_index:
        _neg = _neg_index or {"by_toplevelkey": {}, "by_key": {}, "by_subkey": {}}
        _pos = _pos_index or {"by_toplevelkey": {}, "by_key": {}, "by_subkey": {}}
        _ctrl_threshold = args.control_fold_threshold
        _ctrl_annotated = 0

        for group in final_json:
            # Group-level comparison (toplevelkey)
            grp_ctrl = compute_control_comparison(
                data=group, neg_index=_neg, pos_index=_pos,
                fold_threshold=_ctrl_threshold, level="toplevelkey",
            )
            if grp_ctrl:
                group['control_comparison'] = grp_ctrl
                _ctrl_annotated += 1

            # Member-level comparison (key and subkey)
            for member in group.get('members', []):
                # Compare at key level (strain)
                m_ctrl_key = compute_control_comparison(
                    data=member, neg_index=_neg, pos_index=_pos,
                    fold_threshold=_ctrl_threshold, level="key",
                )
                # Compare at subkey level (species)
                m_ctrl_subkey = compute_control_comparison(
                    data=member, neg_index=_neg, pos_index=_pos,
                    fold_threshold=_ctrl_threshold, level="subkey",
                )
                member['control_comparison'] = {}
                if m_ctrl_key:
                    member['control_comparison']['by_key'] = m_ctrl_key
                if m_ctrl_subkey:
                    member['control_comparison']['by_subkey'] = m_ctrl_subkey
                # Also store a top-level summary using the key-level result
                # (or subkey fallback) for easy access in the report
                _best_ctrl = m_ctrl_key or m_ctrl_subkey
                if _best_ctrl:
                    for _f in ('neg_max_tass', 'neg_max_reads', 'tass_fold_over_neg',
                               'reads_fold_over_neg', 'tass_fold_over_pos',
                               'reads_fold_over_pos', 'control_flag',
                               'neg_control_values', 'pos_control_values',
                               'pos_min_tass', 'pos_min_reads'):
                        member['control_comparison'][_f] = _best_ctrl.get(_f)

        if _ctrl_annotated:
            print(f"Control comparison: annotated {_ctrl_annotated} species groups "
                  f"(neg={len(args.negative_controls or [])} files, "
                  f"pos={len(args.positive_controls or [])} files, "
                  f"fold_threshold={_ctrl_threshold})")

        # Detect positive control organisms missing from the sample
        if _pos and not args.hide_missing_pos_controls:
            _missing_pos = find_missing_positive_controls(
                final_json, _pos, levels=args.missing_pos_levels)
            if _missing_pos:
                print(f"Missing positive controls: {len(_missing_pos)} organism(s) "
                      f"in positive control not found in sample "
                      f"(levels={args.missing_pos_levels})")
        else:
            _missing_pos = []
    else:
        _missing_pos = []

    # ── In-silico control comparison annotation ──────────────────────────────
    # Same pattern as lab controls: load the insilico JSON(s), index by key
    # levels, compute fold-change metrics, and annotate each organism.
    _insilico_index = load_control_data(args.insilico_controls) if args.insilico_controls else None
    _missing_insilico = []

    if _insilico_index:
        # Treat insilico as a positive-like control: we compare sample vs
        # what the simulation produced.  Use positive-control slots in the
        # metrics so fold-over-pos semantics apply (sample / insilico).
        _empty_neg = {"by_toplevelkey": {}, "by_key": {}, "by_subkey": {}}
        _insilico_annotated = 0

        for group in final_json:
            grp_isil = compute_control_comparison(
                data=group, neg_index=_empty_neg, pos_index=_insilico_index,
                fold_threshold=args.control_fold_threshold, level="toplevelkey",
            )
            if grp_isil:
                group['insilico_comparison'] = grp_isil
                _insilico_annotated += 1

            for member in group.get('members', []):
                m_isil_key = compute_control_comparison(
                    data=member, neg_index=_empty_neg, pos_index=_insilico_index,
                    fold_threshold=args.control_fold_threshold, level="key",
                )
                m_isil_subkey = compute_control_comparison(
                    data=member, neg_index=_empty_neg, pos_index=_insilico_index,
                    fold_threshold=args.control_fold_threshold, level="subkey",
                )
                member['insilico_comparison'] = {}
                if m_isil_key:
                    member['insilico_comparison']['by_key'] = m_isil_key
                if m_isil_subkey:
                    member['insilico_comparison']['by_subkey'] = m_isil_subkey
                _best_isil = m_isil_key or m_isil_subkey
                if _best_isil:
                    for _f in ('pos_min_tass', 'pos_min_reads',
                               'tass_fold_over_pos', 'reads_fold_over_pos',
                               'pos_control_values'):
                        member['insilico_comparison'][_f] = _best_isil.get(_f)
                    # Rename for clarity: these are insilico values
                    member['insilico_comparison']['insilico_tass'] = _best_isil.get('pos_min_tass')
                    member['insilico_comparison']['insilico_reads'] = _best_isil.get('pos_min_reads')
                    member['insilico_comparison']['tass_fold_over_insilico'] = _best_isil.get('tass_fold_over_pos')
                    member['insilico_comparison']['reads_fold_over_insilico'] = _best_isil.get('reads_fold_over_pos')

        if _insilico_annotated:
            print(f"In-silico comparison: annotated {_insilico_annotated} species groups "
                  f"(insilico={len(args.insilico_controls or [])} files, "
                  f"fold_threshold={args.control_fold_threshold})")

        # Detect organisms present in insilico but missing from sample
        _missing_insilico = find_missing_positive_controls(
            final_json, _insilico_index, levels=args.missing_pos_levels)
        if _missing_insilico:
            # Tag as insilico-origin
            for entry in _missing_insilico:
                entry['control_source'] = 'insilico'
            print(f"Missing in-silico controls: {len(_missing_insilico)} organism(s) "
                  f"in in-silico control not found in sample")

        # Also detect organisms in sample but missing from insilico
        _sample_only = _find_sample_only_organisms(final_json, _insilico_index)
        if _sample_only:
            print(f"Sample-only organisms: {len(_sample_only)} organism(s) "
                  f"in sample not found in in-silico control")
            # Store on each relevant group/member
            for group in final_json:
                tlk = str(group.get('toplevelkey', group.get('key', '')))
                if tlk and tlk in _sample_only:
                    if 'insilico_comparison' not in group:
                        group['insilico_comparison'] = {}
                    group['insilico_comparison']['missing_from_insilico'] = True

    # ── Build structured output with metadata ────────────────────────────────
    _all_keys = set()
    _all_subkeys = set()
    _all_toplevelkeys = set()
    _total_organism_reads = 0
    for grp in final_json:
        _all_toplevelkeys.add(grp.get('toplevelkey', grp.get('key', '')))
        for m in grp.get('members', []):
            _all_keys.add(m.get('key', ''))
            _all_subkeys.add(m.get('subkey', m.get('key', '')))
            _total_organism_reads += float(m.get('numreads', 0))
    # Strip covered_regions from JSON output — large and only needed internally
    import copy as _copy
    _json_out = _copy.deepcopy(final_json)
    for _grp in _json_out:
        _grp.pop('covered_regions', None)
        for _m in _grp.get('members', []):
            _m.pop('covered_regions', None)
    # ── Resolve best cutoffs from thresholds JSON (if available) ──────────
    # Extract per-group best_threshold values so create_report.py can use
    # them instead of its hard-coded sampletype defaults.
    _best_cutoffs = None
    if thresholds_config and _st_key in thresholds_config:
        _groups_block = thresholds_config[_st_key].get("groups", {})
        if _groups_block:
            _best_cutoffs = {}
            for _gname, _gdata in _groups_block.items():
                _best_cutoffs[_gname] = {
                    "best_threshold": _gdata.get("best_threshold"),
                    "fp_le_0_1pct": _gdata.get("fp_le_0_1pct", {}).get("threshold"),
                    "tp_ge_99_5pct": _gdata.get("tp_ge_99_5pct", {}).get("threshold"),
                }
            print(f"  Best cutoffs from thresholds JSON ({_st_key}):")
            for _gn, _gc in _best_cutoffs.items():
                print(f"    {_gn}: best_threshold={_gc['best_threshold']}")

    output_json = {
        "metadata": {
            "sample_name": args.samplename,
            "sample_type": sampletype,
            "platform": args.platform,
            "workflow_revision": args.workflow_revision,
            "commit_id": args.commit_id,
            "total_reads": total_reads,
            "aligned_reads": aligned_reads_total,
            "total_organism_reads": int(_total_organism_reads),
            "num_species_groups": len(final_json),
            "num_keys": len(_all_keys),
            "num_subkeys": len(_all_subkeys),
            "num_toplevelkeys": len(_all_toplevelkeys),
            "minmapq": args.minmapq,
            "mapq_breadth_power": args.mapq_breadth_power,
            "weights": dict(weights),
            "min_conf_applied": None,  # populated by create_report per-sample
            "best_cutoffs": _best_cutoffs,  # from thresholds JSON; None if not available
            "best_cutoffs_source": _st_key if _best_cutoffs else None,
            "preferred_granularity": args.optimize_granularity,
            "control_type": args.control_type,  # "positive", "negative", or None
            "negative_controls_used": [
                os.path.basename(p) for p in (args.negative_controls or [])
            ] or None,
            "positive_controls_used": [
                os.path.basename(p) for p in (args.positive_controls or [])
            ] or None,
            "control_fold_threshold": args.control_fold_threshold if (
                args.negative_controls or args.positive_controls) else None,
            "missing_positive_controls": _missing_pos if _missing_pos else None,
            "insilico_controls_used": [
                os.path.basename(p) for p in (args.insilico_controls or [])
            ] or None,
            "missing_insilico_controls": _missing_insilico if _missing_insilico else None,
        },
        "organisms": _json_out,
    }
    write_to_json(args.output.replace(".tsv", ".json"), output_json)

    # ── Report Metrics (TP/FP analysis with threshold percentiles) ───────────
    if args.report_metrics:
        print("\n" + "=" * 80)
        print("REPORT METRICS: Computing TP/FP metrics from ground-truth read IDs")
        print("=" * 80)

        tp_fp_counts = compute_tp_fp_counts_by_taxid(
            inputfile,
            accession_to_taxid=acc_to_key,
            debug_n=10,
        )

        # Build per-organism table from strain_summary
        organism_rows = []
        for k, data in strain_summary.items():
            taxid_str = str(data.get('key', k))
            tp_reads, fp_reads = tp_fp_counts.get(taxid_str, (0, 0))
            total_rd = tp_reads + fp_reads
            is_tp = tp_reads > 0
            organism_rows.append({
                "name": data.get('name', ''),
                "key": data.get('key', k),
                "subkey": data.get('subkey', ''),
                "subkeyname": data.get('subkeyname', ''),
                "tass_score": round(float(data.get('tass_score', 0)), 6),
                "tp_reads": int(tp_reads),
                "fp_reads": int(fp_reads),
                "total_reads": int(total_rd),
                "classification": "TP" if is_tp else "FP",
                "microbial_category": data.get('microbial_category', 'Unknown'),
            })

        # Filter out organisms with zero TP and zero FP reads –
        # they carry no ground-truth signal and skew mean scores.
        organism_rows = [r for r in organism_rows if (r['tp_reads'] + r['fp_reads']) >= 1]

        organism_rows.sort(key=lambda r: r['tass_score'], reverse=True)

        # Compute overall metrics at taxid level
        tp_taxa = [r for r in organism_rows if r['classification'] == 'TP']
        fp_taxa = [r for r in organism_rows if r['classification'] == 'FP']
        n_tp = len(tp_taxa)
        n_fp = len(fp_taxa)
        # FN: taxids that have ground-truth reads but aren't in strain_summary
        strain_taxids = {str(data.get('key', k)) for k, data in strain_summary.items()}
        fn_taxids = [t for t, (tp, fp) in tp_fp_counts.items() if tp > 0 and t not in strain_taxids]
        n_fn = len(fn_taxids)

        precision = n_tp / (n_tp + n_fp) if (n_tp + n_fp) > 0 else 0.0
        recall = n_tp / (n_tp + n_fn) if (n_tp + n_fn) > 0 else 0.0
        f1 = (2 * precision * recall / (precision + recall)) if (precision + recall) > 0 else 0.0

        # Threshold-based percentile analysis
        # At each TASS cutoff, what % of TPs pass and what % of FPs pass
        import numpy as np
        thresholds = sorted(set(
            [round(x * 0.01, 2) for x in range(1, 100)] +
            [round(x * 0.1, 1) for x in range(1, 10)]
        ))
        threshold_analysis = []
        for thresh in thresholds:
            tp_pass = sum(1 for r in tp_taxa if r['tass_score'] >= thresh)
            fp_pass = sum(1 for r in fp_taxa if r['tass_score'] >= thresh)
            tp_pct = (tp_pass / n_tp * 100) if n_tp > 0 else 0.0
            fp_pct = (fp_pass / n_fp * 100) if n_fp > 0 else 0.0
            threshold_analysis.append({
                "tass_cutoff": thresh,
                "tp_passing": tp_pass,
                "tp_total": n_tp,
                "tp_percent": round(tp_pct, 2),
                "fp_passing": fp_pass,
                "fp_total": n_fp,
                "fp_percent": round(fp_pct, 2),
            })

        def _build_group_metrics(group_field, acc_map, label):
            counts = compute_tp_fp_counts_by_taxid(
                inputfile,
                accession_to_taxid=acc_map,
                debug_n=10,
            )
            groups = defaultdict(list)
            for _, data in strain_summary.items():
                gval = str(data.get(group_field, data.get("key", "")))
                groups[gval].append(data)

            rows = []
            for gval, members in groups.items():
                tp_reads, fp_reads = counts.get(gval, (0, 0))
                total_rd = tp_reads + fp_reads
                tass_score = max(float(m.get("tass_score", 0)) for m in members) if members else 0.0
                rows.append({
                    "group": gval,
                    "tass_score": round(tass_score, 6),
                    "tp_reads": int(tp_reads),
                    "fp_reads": int(fp_reads),
                    "total_reads": int(total_rd),
                    "classification": "TP" if tp_reads > 0 else "FP",
                    "member_count": len(members),
                })

            # Filter out groups with zero TP and zero FP reads
            rows = [r for r in rows if (r['tp_reads'] + r['fp_reads']) >= 1]

            rows.sort(key=lambda r: r["tass_score"], reverse=True)
            tp_groups = [r for r in rows if r["classification"] == "TP"]
            fp_groups = [r for r in rows if r["classification"] == "FP"]
            n_tp_g = len(tp_groups)
            n_fp_g = len(fp_groups)

            present_groups = set(groups.keys())
            fn_groups = [g for g, (tp, _) in counts.items() if tp > 0 and g not in present_groups]
            n_fn_g = len(fn_groups)

            precision_g = n_tp_g / (n_tp_g + n_fp_g) if (n_tp_g + n_fp_g) > 0 else 0.0
            recall_g = n_tp_g / (n_tp_g + n_fn_g) if (n_tp_g + n_fn_g) > 0 else 0.0
            f1_g = (2 * precision_g * recall_g / (precision_g + recall_g)) if (precision_g + recall_g) > 0 else 0.0

            print(f"\n{label.upper()} METRICS: Precision={precision_g:.4f}  Recall={recall_g:.4f}  F1={f1_g:.4f}")
            print(f"         TP groups={n_tp_g}  FP groups={n_fp_g}  FN groups={n_fn_g}")
            print(f"\n{'GROUP':<28s} {'TASS':>8s} {'CLASS':>5s} {'TP_RD':>8s} {'FP_RD':>8s} {'MEMBERS':>8s}")
            print("-" * 70)
            for r in rows:
                print(f"{r['group'][:27]:<28s} {r['tass_score']:>8.4f} {r['classification']:>5s} "
                      f"{r['tp_reads']:>8d} {r['fp_reads']:>8d} {r['member_count']:>8d}")

            return {
                "summary": {
                    "total_groups": len(rows),
                    "true_positives": n_tp_g,
                    "false_positives": n_fp_g,
                    "false_negatives": n_fn_g,
                    "precision": round(precision_g, 6),
                    "recall": round(recall_g, 6),
                    "f1_score": round(f1_g, 6),
                    "tp_total_reads": sum(r['tp_reads'] for r in rows),
                    "fp_total_reads": sum(r['fp_reads'] for r in rows),
                },
                "groups": rows,
            }

        metrics_report = {
            "summary": {
                "total_organisms": len(organism_rows),
                "true_positives": n_tp,
                "false_positives": n_fp,
                "false_negatives": n_fn,
                "precision": round(precision, 6),
                "recall": round(recall, 6),
                "f1_score": round(f1, 6),
                "tp_total_reads": sum(r['tp_reads'] for r in organism_rows),
                "fp_total_reads": sum(r['fp_reads'] for r in organism_rows),
            },
            "organisms": organism_rows,
            "threshold_analysis": threshold_analysis,
            "subkey_metrics": _build_group_metrics("subkey", acc_to_subkey, "subkey"),
            "toplevelkey_metrics": _build_group_metrics("toplevelkey", acc_to_toplevelkey, "toplevelkey"),
            "weights_used": {k: round(v, 6) for k, v in weights.items()},
        }

        # Write JSON
        metrics_json_path = args.output.replace(".tsv", ".metrics.json")
        write_to_json(metrics_json_path, metrics_report)
        print(f"\nMetrics JSON written to: {metrics_json_path}")

        # Print to stdout
        print(f"\n{'─' * 80}")
        print(f"SUMMARY: Precision={precision:.4f}  Recall={recall:.4f}  F1={f1:.4f}")
        print(f"         TP taxa={n_tp}  FP taxa={n_fp}  FN taxa={n_fn}")
        print(f"         TP reads={metrics_report['summary']['tp_total_reads']}  "
              f"FP reads={metrics_report['summary']['fp_total_reads']}")
        print(f"{'─' * 80}")

        # Organism table
        print(f"\n{'NAME':<40s} {'KEY':<12s} {'SUBKEY':<12s} {'TASS':>8s} {'CLASS':>5s} "
              f"{'TP_RD':>8s} {'FP_RD':>8s} {'CATEGORY':<20s}")
        print("─" * 125)
        for r in organism_rows:
            print(f"{r['name'][:39]:<40s} {str(r['key']):<12s} {str(r['subkey']):<12s} "
                  f"{r['tass_score']:>8.4f} {r['classification']:>5s} "
                  f"{r['tp_reads']:>8d} {r['fp_reads']:>8d} {r['microbial_category']:<20s}")

        # Threshold table
        print(f"\n{'─' * 80}")
        print(f"THRESHOLD ANALYSIS: %% of TP/FP taxa passing at each TASS cutoff")
        print(f"{'─' * 80}")
        print(f"{'TASS_CUTOFF':>12s} {'TP_PASS':>8s} {'TP_%':>8s} {'FP_PASS':>8s} {'FP_%':>8s}")
        print("─" * 50)
        # Print at 0.05 increments for readability on stdout
        for t in threshold_analysis:
            cutoff = t['tass_cutoff']
            if cutoff * 100 % 5 == 0 or cutoff in (0.01, 0.99):
                print(f"{cutoff:>12.2f} {t['tp_passing']:>8d} {t['tp_percent']:>7.1f}% "
                      f"{t['fp_passing']:>8d} {t['fp_percent']:>7.1f}%")

        print(f"\n(Full per-0.01 breakdown available in {metrics_json_path})")
        print("=" * 80)


def random_tweak(weights, scale=0.2):
    """
    Returns a dict, each weight randomly offset by up to ±scale.
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
            # Compute percent_aligned relative to total reads in sample
            _entry_total = entry.get('total_reads', 0)
            _entry_reads = entry.get('reads_aligned', 0)
            percent_aligned = (100 * _entry_reads / _entry_total) if _entry_total > 0 and _entry_reads > 0 else 0
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
