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

# ---------------------------------------------------------------------------
# FIPS-safe environment setup (must run before ANY imports that touch OpenSSL)
# On HPC systems with FIPS mode enforced at the kernel level, the OpenSSL
# library inside containers may fail the FIPS self-test and SIGABRT the
# process.  Setting OPENSSL_CONF to /dev/null skips the FIPS provider
# entirely.  We also ensure matplotlib and fontconfig have writable cache
# directories so they don't emit warnings on read-only home filesystems.
# ---------------------------------------------------------------------------
import os as _os
import tempfile as _tempfile

if not _os.environ.get("OPENSSL_CONF"):
    _os.environ["OPENSSL_CONF"] = "/dev/null"

if not _os.environ.get("MPLCONFIGDIR"):
    _mpl_tmp = _os.path.join(_tempfile.gettempdir(), "matplotlib_cache")
    _os.makedirs(_mpl_tmp, exist_ok=True)
    _os.environ["MPLCONFIGDIR"] = _mpl_tmp

if not _os.environ.get("FONTCONFIG_PATH"):
    _os.environ["FONTCONFIG_PATH"] = _tempfile.gettempdir()

if not _os.environ.get("XDG_CACHE_HOME"):
    _os.environ["XDG_CACHE_HOME"] = _tempfile.gettempdir()
# ---------------------------------------------------------------------------

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
        "--annotate_report",
        default=None,
        type=str,
        help="Annotation report TSV (merged DIAMOND hits + metadata from annotate_report.py). "
             "Used to attach protein annotation data (gene names, classification, resistance info) "
             "to each organism in the output JSON, matched by species_taxid to the subkey.",
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
    parser.add_argument("--optimize_mode", type=str, default="taxids",
                    choices=["taxids", "reads", "hybrid"],
                    help="Optimization objective mode. 'taxids' weights each organism equally, "
                            "'reads' weights by read count, 'hybrid' blends both via --optimize_hybrid_lambda. "
                            "Default: taxids.")
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
    parser.add_argument("--optimize_hybrid_lambda", type=float, default=0.5,
                        help="Blend factor when --optimize_mode hybrid. 1.0 = pure taxid-level, "
                             "0.0 = pure read-level. Default: 0.5.")
    parser.add_argument("--optimize_categories", nargs='+', default=None,
                        choices=["Primary", "Potential", "Opportunistic", "Commensal", "Unknown"],
                        help="Restrict weight optimization to only these microbial categories. "
                             "Accepts one or more of: Primary, Potential, Opportunistic, Commensal, Unknown. "
                             "If not set, all categories are included in the optimization objective.")
    # ── Site-aware breadth sigmoid tuning ────────────────────────────────
    parser.add_argument("--breadth_midpoint", type=float, default=0.0005,
                        help="Breadth sigmoid midpoint (fraction of genome covered at 50%% score). "
                             "Lower values make the sigmoid more sensitive to low-coverage organisms. "
                             "For sterile/blood sites use 0.001 or 0.0001. Default: 0.0005 (0.05%% coverage).")
    parser.add_argument("--breadth_steepness", type=float, default=50000,
                        help="Breadth sigmoid steepness. Higher = sharper transition. "
                             "When lowering --breadth_midpoint, increase steepness proportionally "
                             "(e.g. midpoint=0.001 → steepness=120000). Default: 50000.")
    # ── Low-abundance confidence (sterile-site boost) ────────────────────
    parser.add_argument("--abundance_confidence_weight", type=float, default=None,
                        help="Weight for the low-abundance confidence component in TASS score. "
                             "This component uses log-RPM to boost organisms that are meaningful "
                             "at low read counts (e.g. pathogens in sterile/blood sites). "
                             "Recommended: 0.15-0.30 for sterile sites. 0 = disabled. Default: 0.20.")
    parser.add_argument("--abundance_rpm_midpoint", type=float, default=5.0,
                        help="Expected RPM for a meaningful detection in this site type. "
                             "For sterile/blood: 1.0-5.0. For gut/skin: 50-200. Default: 5.0.")
    parser.add_argument("--abundance_rpm_steepness", type=float, default=2.0,
                        help="Steepness of the log-RPM sigmoid for abundance confidence. Default: 2.0.")
    parser.add_argument("--abundance_gate", action="store_true", default=False,
                        help="Use abundance_confidence as a multiplicative gate on the entire "
                             "TASS score.  Organisms with trivially low RPM (e.g. 3 reads in a "
                             "deep sample) get their score crushed toward 0, preventing noise "
                             "from accumulating small metric contributions into inflated scores. "
                             "The gate uses the same log-RPM sigmoid controlled by "
                             "--abundance_rpm_midpoint and --abundance_rpm_steepness. "
                             "Default: disabled.")
    parser.add_argument("--score_power", type=float, default=None,
                        help="Power transform (gamma) applied to TASS scores. "
                             "Values < 1 lift compressed scores: 0.09^0.5=0.30, 0.09^0.3=0.52. "
                             "Preserves monotonic ordering so thresholds still separate TP/FP. "
                             "Default: 1.0 (no transform).")
    parser.add_argument("--auto_score_power", action="store_true", default=False,
                        help="Automatically compute score_power after optimization so that "
                             "the mean TP TASS score hits --optimize_tp_target. "
                             "Requires --optimize and --optimize_tp_target. "
                             "Overrides --score_power if both are set.")
    # ── Youden J minimum threshold floor ─────────────────────────────────
    parser.add_argument("--youden_min_threshold", type=float, default=None,
                        help="Minimum allowed TASS threshold for Youden J cutoff. "
                             "Prevents the optimizer from selecting unreasonably low cutoffs "
                             "for sterile sites. E.g. 0.05 means Youden J can never return "
                             "a threshold below 0.05. Default: None (no floor).")
    parser.add_argument(
        '--breadth_weight',
        metavar="BREADTHSCORE",
        type=float,
        default=None,
        help="value of weight for breadth of coverage in final TASS Score",
    )
    parser.add_argument(
        "--minhash_weight",
        metavar="MINHASHSCORE",
        type=float,
        default=None,
        help="value of weight for minhash signature reduction in final TASS Score",
    )
    parser.add_argument(
        "--gini_weight",
        metavar="GINIWEIGHT",
        type=float,
        default=None,
        help="value of weight for gini coefficient in final TASS Score",
    )
    parser.add_argument(
        "--k2_disparity_weight",
        metavar="DISPARITYSCOREWEIGHT",
        type=float,
        default=0.0,
        help="value of weight for disparity of k2 and alignment in final TASS Score",
    )
    parser.add_argument("--disparity_weight", type=float, default=None,
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
        default=None,
        help="value of weight for hmp abundance in final TASS Score",
    )
    parser.add_argument(
        "--plasmid_bonus_weight",
        type=float,
        default=None,
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
    parser.add_argument("--matrix", required=False, help="A Matrix file for ANI in long format from fastANI")
    parser.add_argument("--enable_matrix", required=False, action='store_true', help="Enable loading the ANI matrix file")
    parser.add_argument(
        "--ani_threshold",
        type=float,
        default=0.97,
        help="ANI threshold above which two organisms are considered highly similar (default=0.97). "
             "Used to annotate each member in the output JSON with a 'high_ani_matches' list.",
    )
    parser.add_argument("--comparisons", required=False, help="Skip comparison metrics if present, can be either csv, tsv, or xlsx")
    parser.add_argument("--failed_reads", required=False, help="Load a 2 col tsv of reference   read_id that is to be the passed reads. Remove all others and update bedgraph and cov file(s)")
    parser.add_argument("--scaled", type=int, default=5000, help="scaled factor for MinHash.")
    parser.add_argument("--kmer_size", type=int, default=51, help="k-mer size for MinHash.")
    parser.add_argument("--window_size", type=int, default=900_000,
                        help="Window size (bp) for shared-region sketching in conflict detection. "
                             "Smaller values give finer positional resolution. Default: 500,000.")
    parser.add_argument("--step_size", type=int, default=900_000,
                        help="Step size (bp) between windows for shared-region sketching. "
                             "Use window_size/2 for 50%% overlap. Default: 500,000.")
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
    parser.add_argument('--mapq_breadth_power', required=False, type=float, default=0.3,
                    help="Power exponent for MAPQ-adjusted breadth. breadth *= highmapq_fraction^power. "
                         "Higher values penalize low-MAPQ organisms more aggressively. "
                         "Default: 2 (e.g. 10%% high-MAPQ reads → breadth scaled by 0.1).")
    parser.add_argument('--mapq_gini_power', required=False, type=float, default=0.3,
                    help="Power exponent for MAPQ-adjusted Gini. gini *= highmapq_fraction^power. "
                         "Penalizes Gini score for organisms dominated by low-MAPQ reads. "
                         "Default: 0.85. Set to 0 to disable.")
    parser.add_argument('--contig_penalty_power', required=False, type=float, default=0.3,
                    help="Power exponent for contig utilization penalty on Gini. "
                         "gini *= (covered_contigs/total_contigs)^power. "
                         "Penalizes organisms where reads concentrate on few contigs "
                         "(e.g. 2/2000 contigs covered → Gini drops to ~6%% of original). "
                         "Default: 0.3. Set to 0 to disable.")
    parser.add_argument('--depth_concentration_power', required=False, type=float, default=0.3,
                    help="Power exponent for depth-concentration penalty on Gini. "
                         "Detects the conserved-human-reads pattern: many reads map to a "
                         "tiny region of a large genome (high depth, negligible coverage). "
                         "Measures actual_coverage / expected_coverage from read count. "
                         "e.g. 80K reads on 65Mbp genome with 0.15%% coverage: "
                         "efficiency=0.008, penalty^0.3=0.22 -> Gini crushed to 22%%. "
                         "Default: 0.1. Set to 0 to disable.")
    parser.add_argument(
        '--dominance_protect_ratio',
        type=float,
        default=2.0,
        help="Read-count ratio (dominant / runner-up) within a conflict group above which the "
             "dominant reference is fully exempt from the base (min_cov) removal.  Dominance "
             "is determined by raw read count in the shared region.  Lower values protect more "
             "aggressively (e.g. 2.0 = full protection when winner has ≥2× reads).  Higher "
             "values are more conservative (e.g. 5.0 requires a very lopsided group).  "
             "Default: 3.0.  Set to 0 or a very large number (e.g. 999) to disable protection "
             "and restore the original unconditional base-removal behaviour.",
    )
    parser.add_argument(
        '--skip_read_removal_scoring',
        action='store_true',
        default=False,
        help="When set, conflict detection still runs (and breadth-change stats are computed "
             "and written to the removal_stats report), but the removed read IDs are NOT "
             "applied when counting reads for scoring.  Coverage, Gini, numreads, and all "
             "TASS-score inputs are derived from the full original read set, exactly as if "
             "no reads had been removed.  Useful for comparing Gini/breadth behaviour with "
             "and without removal applied to the scoring pipeline.",
    )
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

    # ── Description-based exclusion filters ─────────────────────────────────
    parser.add_argument(
        "--exclude_descriptions",
        nargs="*",
        default=None,
        help="One or more regex patterns to match against reference descriptions "
             "from the matchfile.  Any accession whose description matches ANY of "
             "these patterns will be removed from the analysis (reads discarded).  "
             "Useful for removing unplaced genomic scaffolds, mitochondrial sequences, "
             "etc.  Example: --exclude_descriptions 'unplaced genomic scaffold' "
             "'mitochondrion'.  If the flag is provided with no patterns, nothing is "
             "filtered.  If the flag is omitted entirely, nothing is filtered.",
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
    # ── Three-layer param priority: defaults < thresholds_json < explicit CLI ─
    # Params that are still None were not passed on the CLI; record which ones
    # the user DID explicitly supply so the JSON block never overwrites them.
    _WEIGHT_DEFAULTS = {
        "breadth_weight":              0.27,
        "minhash_weight":              0.31,
        "gini_weight":                 0.42,
        "disparity_weight":            0.00,
        "hmp_weight":                  0.00,
        "plasmid_bonus_weight":        0.0,
        "abundance_confidence_weight": 0.0,
        "score_power":                 0.7,
    }
    # Which of these did the user pass explicitly on the CLI?
    _cli_provided = {p for p, _ in _WEIGHT_DEFAULTS.items()
                     if getattr(args, p, None) is not None}
    # Fill in defaults for everything not explicitly set
    for _p, _v in _WEIGHT_DEFAULTS.items():
        if getattr(args, _p, None) is None:
            setattr(args, _p, _v)

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
        # Override defaults with JSON best weights — but never overwrite a param
        # the user explicitly passed on the CLI (tracked in _cli_provided).
        if "breadth_weight" in _best_w and "breadth_weight" not in _cli_provided:
            args.breadth_weight = _best_w["breadth_weight"]
        if "gini_weight" in _best_w and "gini_weight" not in _cli_provided:
            args.gini_weight = _best_w["gini_weight"]
        if "minhash_weight" in _best_w and "minhash_weight" not in _cli_provided:
            args.minhash_weight = _best_w["minhash_weight"]
        if "hmp_weight" in _best_w and "hmp_weight" not in _cli_provided:
            args.hmp_weight = _best_w["hmp_weight"]
        if "disparity_weight" in _best_w and "disparity_weight" not in _cli_provided:
            disparity_w = _best_w["disparity_weight"]
            args.disparity_weight = disparity_w
        if "plasmid_bonus_weight" in _best_w and "plasmid_bonus_weight" not in _cli_provided:
            args.plasmid_bonus_weight = _best_w["plasmid_bonus_weight"]
        if "abundance_confidence_weight" in _best_w and "abundance_confidence_weight" not in _cli_provided:
            args.abundance_confidence_weight = _best_w["abundance_confidence_weight"]
        if "score_power" in _best_w and "score_power" not in _cli_provided:
            args.score_power = _best_w["score_power"]
        print(f"  Applied weights from JSON: breadth={args.breadth_weight:.6g}, "
              f"gini={args.gini_weight:.6g}, minhash={args.minhash_weight:.6g}, "
              f"hmp={args.hmp_weight:.6g}, disparity={disparity_w:.6g}, "
              f"plasmid_bonus={args.plasmid_bonus_weight:.6g}, "
              f"abundance_conf={args.abundance_confidence_weight:.6g}, "
              f"score_power={args.score_power:.6g}")

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
    # Abundance confidence is also additive (outside normalized pool).
    # It boosts organisms meaningful at low read counts (sterile/blood).
    weights['abundance_confidence_weight'] = args.abundance_confidence_weight
    # Multiplicative abundance gate: when enabled, the abundance_confidence
    # sigmoid is used as a multiplier on the entire TASS score, crushing
    # noise organisms with trivially low RPM.
    weights['abundance_gate'] = args.abundance_gate
    # Score power transform: applied to TASS scores to lift compressed values.
    # score_power < 1 expands the low-score range (0.09^0.5 = 0.30).
    weights['score_power'] = args.score_power
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
                window_size=args.window_size,
                step_size=args.step_size,
                accession_to_taxid=early_acc_to_taxid if early_acc_to_taxid else None,
                taxid_to_name=early_taxid_to_desc if early_taxid_to_desc else None,
                taxid_removal_stats=args.taxid_removal_stats,
                dominance_protect_ratio=args.dominance_protect_ratio,
            )
            # ── Skip-removal-for-scoring mode ────────────────────────────────
            # conflict detection + breadth-change stats are always computed and
            # written to disk (removal_stats.xlsx, cluster_dominance.json, etc.).
            # When --skip_read_removal_scoring is set we simply discard the
            # removal dict so that count_reference_hits uses every read for
            # numreads / coverage / gini, giving a clean baseline for comparison.
            if getattr(args, 'skip_read_removal_scoring', False):
                print(
                    "[skip_read_removal_scoring] Conflict removal computed but NOT applied "
                    "to scoring reads.  All original reads will be used for coverage, "
                    "Gini, numreads, and TASS-score inputs."
                )
                alignments_to_remove = defaultdict(set)
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
    print("Counting reference hits from BAM file…")
    reference_hits, _, total_reads = count_reference_hits(
        inputfile,
        alignments_to_remove=alignments_to_remove,
        args=args
    )
    # exit()

    mmbert_dict = dict()
    if args.microbert:
        print("Loading MicroBERT data from", args.microbert)
        # import the tsv as a dictionary
        mmbert = pd.read_csv(args.microbert, sep='\t', header=0)
        # set the taxid col to str
        mmbert['taxid'] = mmbert['taxid'].astype(str)
        mmbert_dict = mmbert.set_index('taxid').T.to_dict()

    taxdump, taxdump_names = {}, {}
    if args.taxdump and os.path.exists(os.path.join(args.taxdump, "nodes.dmp")):
        print("Loading taxdump from", args.taxdump)
        taxdump = load_taxdump(os.path.join(args.taxdump, "nodes.dmp"))
    if args.taxdump and os.path.exists(os.path.join(args.taxdump, "names.dmp")):
        print("Loading taxdump names from", args.taxdump)
        taxdump_names = load_names(os.path.join(args.taxdump, "names.dmp"))

    # When the matchfile lacks a name/description column (e.g. only
    # accession + taxid), taxid_to_desc entries will be empty strings.
    # Back-fill from taxdump/names.dmp so downstream code has real
    # organism names even for minimal 2-column matchfiles.
    if taxdump_names and early_taxid_to_desc:
        print("Backfilling taxid→description from taxdump names.dmp for taxa missing descriptions in the matchfile…")
        for _tid, _desc in early_taxid_to_desc.items():
            if not _desc and _tid in taxdump_names:
                early_taxid_to_desc[_tid] = taxdump_names[_tid]
    # Also fill in any taxids from early_acc_to_taxid that aren't in
    # early_taxid_to_desc at all (shouldn't happen, but defensive).
    if taxdump_names and early_acc_to_taxid:
        print("Backfilling taxid→description for taxids from early accession→taxid map that are missing from the early taxid→desc map…")
        for _acc, _tid in early_acc_to_taxid.items():
            if _tid not in early_taxid_to_desc and _tid in taxdump_names:
                early_taxid_to_desc[_tid] = taxdump_names[_tid]

    pathogens = import_pathogens(pathogenfile)

    if args.match and os.path.exists(matcher):
        print("Processing matchfile to enrich reference hits with taxid and description information…")
        accindex = args.accessioncol
        nameindex = args.namecol
        taxcol = args.taxcol
        orgindex = args.orgcol

        # Map unversioned accession -> canonical BAM reference accession so
        # mapfiles that omit version suffixes (e.g., .1) still match.
        _ref_by_unversioned = {}
        for _ref_acc in reference_hits.keys():
            _ref_unv = re.sub(r"\.\d+$", "", str(_ref_acc))
            _ref_by_unversioned[_ref_unv] = _ref_acc

        with open(matcher, "r") as f:
            _sample = f.read(8192)
            f.seek(0)
            _delim = '\t' if _sample.count('\t') >= _sample.count(',') else ','
            _reader = csv.reader(f, delimiter=_delim)

            # ── Header auto-detection ────────────────────────────────
            # Recognised header synonyms (all lowercased).
            _ACC_HDRS = {"acc", "accession", "accession_version", "refseq",
                         "accession.version", "ref"}
            _TAX_HDRS = {"taxid", "tax_id", "mapped_value", "staxids",
                         "taxon_id", "ncbi_taxid", "species_taxid"}
            _NAME_HDRS = {"name", "description", "desc", "seqname",
                          "sequence_name", "refseq_name"}
            _ORG_HDRS  = {"organism", "org", "species", "scientific_name",
                          "organism_name"}

            _skip_header = False
            _ncols = None

            for i, splitline in enumerate(_reader):
                # First row: try to auto-detect columns from header
                if i == 0 and splitline:
                    _ncols = len(splitline)
                    _hdr_lower = [c.strip().lower() for c in splitline]

                    # Check if this looks like a header row
                    _is_header = any(
                        h in (_ACC_HDRS | _TAX_HDRS | _NAME_HDRS | _ORG_HDRS)
                        for h in _hdr_lower
                    )

                    if _is_header:
                        _skip_header = True
                        # Auto-detect column indices from header names
                        for _ci, _h in enumerate(_hdr_lower):
                            if _h in _ACC_HDRS:
                                accindex = _ci
                            elif _h in _TAX_HDRS:
                                taxcol = _ci
                            elif _h in _NAME_HDRS:
                                nameindex = _ci
                            elif _h in _ORG_HDRS:
                                orgindex = _ci
                        print(f"Matchfile header detected ({_ncols} columns): "
                              f"acc={accindex}, taxid={taxcol}, "
                              f"name={nameindex}, org={orgindex}")
                        continue

                if not splitline or len(splitline) <= accindex or accindex < 0:
                    continue
                # Track column count from first data row if header wasn't present
                if _ncols is None:
                    _ncols = len(splitline)

                accession = splitline[accindex].strip() or None
                taxid     = splitline[taxcol].strip() if len(splitline) > taxcol else None
                name      = splitline[nameindex].strip() if len(splitline) > nameindex else None
                organism = splitline[orgindex].strip() if len(splitline) > orgindex else None
                if not accession:
                    continue

                accession = str(accession)
                accession_unv = re.sub(r"\.\d+$", "", accession)
                target_accession = accession
                if target_accession not in reference_hits:
                    target_accession = _ref_by_unversioned.get(accession_unv)
                if not target_accession:
                    continue

                # Optional: only fill taxid if it exists and if this accession is already in reference_hits
                if target_accession in reference_hits:
                    if taxid and not reference_hits[target_accession].get("taxid"):
                        reference_hits[target_accession]["taxid"] = taxid

                    # If org/name columns are empty or missing, prefer a taxdump
                    # scientific name over accession fallback when taxid is known.
                    _resolved_taxid = taxid or reference_hits[target_accession].get("taxid")
                    _taxdump_name = (
                        taxdump_names.get(str(_resolved_taxid))
                        if _resolved_taxid and taxdump_names
                        else None
                    )

                    # Store the raw description from the mapfile BEFORE
                    # overwriting 'name' with the cleaner organism column.
                    # The description often contains "plasmid pXXX" which
                    # we need for plasmid tagging downstream.
                    if name:
                        reference_hits[target_accession]["description"] = name
                    if organism:
                        reference_hits[target_accession]["name"] = organism
                    elif name:
                        reference_hits[target_accession]["name"] = name
                    elif _taxdump_name:
                        reference_hits[target_accession]["name"] = _taxdump_name
        f.close()

    # ── Resolve missing names from taxdump/names.dmp ─────────────────────
    # When the matchfile has only accession + taxid (no name/org columns),
    # entries may still have name == accession.  Resolve them from
    # taxdump_names so downstream JSON has real organism names.
    if taxdump_names:
        print("Resolving missing organism names from taxdump names.dmp for entries with taxids but no names in the matchfile…")
        _resolved_count = 0
        for _acc, _hit in reference_hits.items():
            _cur_name = _hit.get('name', _acc)
            _tid = _hit.get('taxid')
            # If the name is still just the accession (never overwritten),
            # look it up from names.dmp via the taxid.
            if _tid and (_cur_name == _acc or not _cur_name):
                _sci_name = taxdump_names.get(str(_tid))
                if _sci_name:
                    _hit['name'] = _sci_name
                    _resolved_count += 1
        if _resolved_count:
            print(f"Resolved {_resolved_count} organism names from names.dmp (matchfile had no name/org columns)")

    # ── Exclude accessions whose description matches --exclude_descriptions ──
    # Patterns are compiled as regex and matched case-insensitively against
    # the 'description' (or fallback 'name') field populated from the matchfile.
    if args.exclude_descriptions:
        _exclude_patterns = [re.compile(p, re.IGNORECASE) for p in args.exclude_descriptions]
        _excluded_accs = []
        _excluded_reads = 0
        for acc, hit in list(reference_hits.items()):
            _desc = hit.get('description', '') or hit.get('name', '') or ''
            for _pat in _exclude_patterns:
                if _pat.search(_desc):
                    _excluded_accs.append(acc)
                    _excluded_reads += hit.get('numreads', 0)
                    break
        for acc in _excluded_accs:
            del reference_hits[acc]
        if _excluded_accs:
            print(f"--exclude_descriptions: removed {len(_excluded_accs)} accessions "
                  f"({_excluded_reads} reads) matching {len(_exclude_patterns)} pattern(s)")

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

    subrank = args.subrank

    if not subrank or subrank.lower() == "none" or subrank.lower() == "strain":
        subrank = None

    # get the ranks for taxid: 198214
    acc_to_parent = dict()
    if args.rank:
        print(f"Resolving taxonomic ranks from taxdump for rank '{args.rank}' (subrank '{subrank}')…")
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
    # variance_reads = calculate_var(all_readscounts)
    # print(f"\n\tVariance of reads: {variance_reads}")
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
        print("Aggregating comparison_df to subkey level for minhash scoring…")
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

    # ── Build cross-genus conflict clusters from shared_windows_report.csv ──
    # The standard minhash winner/loser logic groups organisms by toplevelkey
    # (genus).  Highly similar organisms from different genera — e.g. E. coli
    # (Escherichia) and Shigella — compete for the same reads but land in
    # separate genus buckets and are scored as independent solos.  This causes
    # the true winner to receive less boost than it deserves and the true loser
    # to go unpenalised at the cross-genus level.
    #
    # If shared_windows_report.csv exists in the output dir, we use its pairwise
    # window-match data to build a union-find conflict graph.  Every accession
    # that shares ≥1 window with another (across any genus) is assigned a shared
    # conflict_cluster_id.  This ID supersedes toplevelkey in the parent_buckets
    # computation inside calculate_normalized_groups, so cross-genus competitors
    # are grouped together for the winner/loser minhash adjustment.
    #
    # When FASTA files are not provided the shared_windows_report.csv will be
    # absent or empty; in that case the fallback toplevelkey grouping applies.
    _acc_to_conflict_cluster: dict = {}
    _sw_cluster_dominance: dict = {}
    _output_dir_for_sw = args.output_dir if args.output_dir else os.path.dirname(args.output)
    _sw_report_path = os.path.join(_output_dir_for_sw, "shared_windows_report.csv")
    if os.path.exists(_sw_report_path):
        try:
            _sw_df = pd.read_csv(_sw_report_path)
            if not _sw_df.empty and "query_contig" in _sw_df.columns and "match_contig" in _sw_df.columns:
                # Build union-find across (query_contig, match_contig) pairs
                _uf_parent: dict = {}

                def _uf_find(x):
                    while _uf_parent.get(x, x) != x:
                        _uf_parent[x] = _uf_parent.get(_uf_parent[x], _uf_parent[x])
                        x = _uf_parent[x]
                    return x

                def _uf_union(a, b):
                    ra, rb = _uf_find(a), _uf_find(b)
                    if ra != rb:
                        _uf_parent[ra] = rb

                for _, row in _sw_df.iterrows():
                    qa = str(row.get("query_contig", "")).strip()
                    mb = str(row.get("match_contig", "")).strip()
                    if qa and mb and qa != mb:
                        _uf_union(qa, mb)

                # Map each accession to its cluster root
                _all_accs_in_sw = set(_sw_df["query_contig"].astype(str)) | set(_sw_df["match_contig"].astype(str))
                for _acc in _all_accs_in_sw:
                    _acc_to_conflict_cluster[_acc] = f"conflict_cluster_{_uf_find(_acc)}"

                _n_clusters = len(set(_acc_to_conflict_cluster.values()))
                print(f"Loaded shared_windows_report.csv: {len(_all_accs_in_sw)} accessions "
                      f"→ {_n_clusters} cross-genus conflict clusters")

                # Derive per-accession cluster_dominance from shared windows
                # so the penalty works regardless of --compare_references mode.
                _sw_partner_count: dict = defaultdict(set)
                _sw_region_count: dict = defaultdict(int)
                for _, _swr in _sw_df.iterrows():
                    _qa = str(_swr.get("query_contig", "")).strip()
                    _mb = str(_swr.get("match_contig", "")).strip()
                    if _qa and _mb and _qa != _mb:
                        _sw_partner_count[_qa].add(_mb)
                        _sw_partner_count[_mb].add(_qa)
                        _sw_region_count[_qa] += 1
                        _sw_region_count[_mb] += 1
                _sw_cluster_dominance: dict = {}
                for _acc in _all_accs_in_sw:
                    _sw_cluster_dominance[_acc] = {
                        "n_shared_regions": _sw_region_count.get(_acc, 0),
                        "n_cluster_partners": len(_sw_partner_count.get(_acc, set())),
                        "cluster_removal_frac": 0.0,
                    }
        except Exception as _sw_err:
            print(f"Warning: could not load shared_windows_report.csv for conflict clustering: {_sw_err}")

    # Annotate each reference_hit with its conflict_cluster_id (if any)
    _clustered_count = 0
    for _acc, _data in reference_hits.items():
        _cid = _acc_to_conflict_cluster.get(_acc)
        if _cid:
            _data["conflict_cluster_id"] = _cid
            _clustered_count += 1
    if _clustered_count:
        print(f"Annotated {_clustered_count} accessions with cross-genus conflict_cluster_id")

    # ── Load cluster_dominance.json for minhash penalty scaling ────────────
    # Written by finalize_proportional_removal when compare_to_reference_windows
    # is False.  Contains per-accession stats: n_shared_regions,
    # n_cluster_partners, cluster_removal_frac.
    _cluster_dominance: dict = {}
    _cd_path = os.path.join(_output_dir_for_sw, "cluster_dominance.json")
    if os.path.exists(_cd_path):
        try:
            import json as _json_cd
            with open(_cd_path) as _cd_fh:
                _cluster_dominance = _json_cd.load(_cd_fh)
            print(f"Loaded cluster_dominance.json: {len(_cluster_dominance)} accessions")
        except Exception as _cd_err:
            print(f"Warning: could not load cluster_dominance.json: {_cd_err}")

    for acc, data in reference_hits.items():
        if not data.get('organism'):
            # try to get organism from taxdump names
            taxid = data.get('taxid')
            if taxid and taxid in taxdump_names:
                data['organism'] = taxdump_names[taxid]
            else:
                data['organism'] = data.get('name', acc)
        # Annotate with cluster dominance stats.
        # Prefer cluster_dominance.json (from False path — has removal_frac),
        # fall back to shared_windows-derived counts (from True path).
        if acc in _cluster_dominance:
            data["cluster_dominance"] = _cluster_dominance[acc]
        elif acc in _sw_cluster_dominance:
            data["cluster_dominance"] = _sw_cluster_dominance[acc]
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
            breadth_midpoint = args.breadth_midpoint,
            breadth_steepness = args.breadth_steepness,
            abundance_rpm_midpoint = args.abundance_rpm_midpoint,
            abundance_rpm_steepness = args.abundance_rpm_steepness,
        )
    strain_summary = calculate_normalized_groups(
        hits=reference_hits,
        group_field="key",
        reads_key="numreads",
        mapq_breadth_power=args.mapq_breadth_power,
        mapq_gini_power=args.mapq_gini_power,
        contig_penalty_power=args.contig_penalty_power,
        depth_concentration_power=args.depth_concentration_power,
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
            abundance_confidence_weight = args.abundance_confidence_weight,
            youden_min_threshold = args.youden_min_threshold,
            prefer_granularity = args.optimize_granularity,
            platform=args.platform,
            abundance_gate = args.abundance_gate,
            optimize_mode = args.optimize_mode,
            hybrid_lambda = args.optimize_hybrid_lambda,
            optimize_categories = args.optimize_categories,
            score_power = args.score_power,
            auto_score_power = args.auto_score_power,
        )
        best_weights = report_weights.get("best_weights") or {}
        weights.update(best_weights)
        print("Optimized weights changed: ")
        for k, v in weights.items():
            print(f"\t{k}: {v}")
    # ── Pre-compute score_power_scale for each strain ────────────────────
    # score_power_scale ∈ [0, 1] controls how much the power transform
    # applies to each organism.  Dominant organisms in their group get the
    # full boost (scale≈1); minor siblings sharing reads with many close
    # relatives get almost none (scale→0).
    #
    # Grouping strategy:
    #   1. If ani_data is available, group organisms that share ANI ≥
    #      ani_threshold (transitive closure of pairwise edges).
    #   2. Otherwise fall back to toplevelkey grouping.
    #
    # Within each group, score_power_scale = reads / group_total_reads,
    # so the effective exponent becomes:
    #   effective_power = 1.0 - (1.0 - score_power) * score_power_scale
    # A sole organism (proportion=1) gets full score_power; one of five
    # equal Shigella (proportion≈0.2) gets almost no boost.

    def _build_ani_groups(strain_dict, ani_data, threshold):
        """Build groups of taxa connected by ANI ≥ threshold (union-find)."""
        parent = {}

        def find(x):
            while parent.get(x, x) != x:
                parent[x] = parent.get(parent[x], parent[x])
                x = parent[x]
            return x

        def union(a, b):
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb

        taxids = set()
        for _k, _d in strain_dict.items():
            t = str(_d.get('key', _k))
            taxids.add(t)

        for t1 in taxids:
            if t1 not in ani_data:
                continue
            for t2, ani_val in ani_data[t1].items():
                if t2 in taxids and ani_val >= threshold:
                    union(t1, t2)

        # Collect groups: root -> list of taxids
        groups = defaultdict(list)
        for t in taxids:
            groups[find(t)].append(t)
        return groups

    if ani_data:
        _ani_groups = _build_ani_groups(strain_summary, ani_data, args.ani_threshold)
        # Invert: taxid -> group root
        _taxid_to_ani_group = {}
        for root, members in _ani_groups.items():
            for m in members:
                _taxid_to_ani_group[m] = root

    for _k, _d in strain_summary.items():
        _taxid = str(_d.get('key', _k))
        _reads = float(_d.get('numreads', 0))

        if ani_data and _taxid in _taxid_to_ani_group:
            # ANI-based grouping
            _group_root = _taxid_to_ani_group[_taxid]
            _group_total = sum(
                float(_d2.get('numreads', 0))
                for _d2 in strain_summary.values()
                if _taxid_to_ani_group.get(str(_d2.get('key', ''))) == _group_root
            )
        else:
            # Fallback: toplevelkey grouping
            _tlk = _d.get('toplevelkey', '')
            _group_total = sum(
                float(_d2.get('numreads', 0))
                for _d2 in strain_summary.values()
                if _d2.get('toplevelkey', '') == _tlk
            )

        _d['score_power_scale'] = (_reads / _group_total) if _group_total > 0 else 1.0

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

    # ── Subkey (species) level aggregation ─────────────────────────────────
    # Intermediate aggregation: strain → subkey (species).  Each subkey group
    # aggregates all strains sharing the same species-level taxid and will
    # later be nested *inside* the toplevelkey (genus) group.
    subkey_summary = calculate_normalized_groups(
        hits=strain_summary,
        group_field="subkey",
        reads_key="numreads",
        mapq_breadth_power=args.mapq_breadth_power,
        mapq_gini_power=args.mapq_gini_power,
        contig_penalty_power=args.contig_penalty_power,
        depth_concentration_power=args.depth_concentration_power,
    )
    # Attach strain-level members to each subkey group
    for _, sk_data in subkey_summary.items():
        sk_data['members'] = [
            x for _, x in strain_summary.items()
            if str(x.get('subkey', x.get('key', ''))) == str(sk_data.get('subkey', sk_data.get('key', '')))
        ]
        # Ensure subkey groups carry correct identity fields
        sk_data['name'] = sk_data.get('subkeyname', sk_data.get('name', None))
        sk_data['key'] = sk_data.get('subkey', sk_data.get('key', None))
    print(f"Subkey (species) aggregation: {len(strain_summary)} strains → {len(subkey_summary)} subkeys")

    aggregate_dict = calculate_normalized_groups(
        hits=strain_summary,
        group_field="toplevelkey",
        reads_key="numreads",
        mapq_breadth_power=args.mapq_breadth_power,
        mapq_gini_power=args.mapq_gini_power,
        contig_penalty_power=args.contig_penalty_power,
        depth_concentration_power=args.depth_concentration_power,
    )
    # iterate through aggregate_dict, make the accession_to_key dict
    for k, v in aggregate_dict.items():
        for acc in species_to_all_accs.get(k, []):
            acc_to_parent[acc] = k

    # Pre-compute score_power_scale for aggregate (species-level) entries.
    # At this level each entry IS a toplevelkey group, so default scale=1.0.
    # But if ANI data links multiple toplevelkeys, compute proportion within
    # the ANI super-group.
    if ani_data:
        _agg_ani_groups = _build_ani_groups(aggregate_dict, ani_data, args.ani_threshold)
        _agg_taxid_to_ani_group = {}
        for root, members in _agg_ani_groups.items():
            for m in members:
                _agg_taxid_to_ani_group[m] = root

        for _k, _d in aggregate_dict.items():
            _taxid = str(_d.get('toplevelkey', _k))
            _reads = float(_d.get('numreads', 0))
            if _taxid in _agg_taxid_to_ani_group:
                _group_root = _agg_taxid_to_ani_group[_taxid]
                _group_total = sum(
                    float(_d2.get('numreads', 0))
                    for _d2 in aggregate_dict.values()
                    if _agg_taxid_to_ani_group.get(str(_d2.get('toplevelkey', ''))) == _group_root
                )
                _d['score_power_scale'] = (_reads / _group_total) if _group_total > 0 else 1.0
            else:
                _d['score_power_scale'] = 1.0
    else:
        for _k, _d in aggregate_dict.items():
            _d['score_power_scale'] = 1.0  # no ANI info at species level, no penalty

    # ── Compute TASS scores for subkey (species) groups ────────────────────
    for _, sk_data in subkey_summary.items():
        sk_data['sample_name'] = args.samplename
        sk_group_reads = [
            dict(reads=x['numreads'], key=x.get('key'))
            for x in sk_data.get('members', [])
        ]
        calculate_aggregate_scores(
            data=sk_data,
            hmp_dists=dists,
            body_sites=[sampletype],
            k2_mapping=k2_mapping,
            sampletype=sampletype,
            mmbert_dict=mmbert_dict,
            group_reads=sk_group_reads,
            total_reads=total_reads,
        )
        sk_data['tass_score'] = compute_tass_score(
            data=sk_data,
            weights=weights,
        )

    for _, data in aggregate_dict.items():
        # Nest subkey groups as members of the toplevelkey group.
        # Each subkey group already contains its own 'members' list of strains.
        tlk = str(data.get('toplevelkey', ''))
        data['members'] = [
            sk_data for _, sk_data in subkey_summary.items()
            if str(sk_data.get('toplevelkey', '')) == tlk
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
            for subkey_member in group.get('members', []):
                # ANI at the subkey (species) level
                sk_key = str(subkey_member.get('key', subkey_member.get('subkey', '')))
                sk_matches = []
                if sk_key in ani_data:
                    for other_key, ani_val in ani_data[sk_key].items():
                        if other_key == sk_key:
                            continue
                        if ani_val >= args.ani_threshold:
                            sk_matches.append({
                                'key': other_key,
                                'ani_pct': round(float(ani_val) * 100, 2),
                            })
                    sk_matches.sort(key=lambda x: x['ani_pct'], reverse=True)
                subkey_member['high_ani_matches'] = sk_matches
                # ANI at the strain level (nested members)
                for strain in subkey_member.get('members', []):
                    strain_key = str(strain.get('key', ''))
                    strain_matches = []
                    if strain_key in ani_data:
                        for other_key, ani_val in ani_data[strain_key].items():
                            if other_key == strain_key:
                                continue
                            if ani_val >= args.ani_threshold:
                                strain_matches.append({
                                    'key': other_key,
                                    'ani_pct': round(float(ani_val) * 100, 2),
                                })
                        strain_matches.sort(key=lambda x: x['ani_pct'], reverse=True)
                    strain['high_ani_matches'] = strain_matches

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

            # Subkey-level comparison (species) and strain-level (nested)
            for subkey_member in group.get('members', []):
                # Compare at subkey level (species)
                sk_ctrl_subkey = compute_control_comparison(
                    data=subkey_member, neg_index=_neg, pos_index=_pos,
                    fold_threshold=_ctrl_threshold, level="subkey",
                )
                subkey_member['control_comparison'] = {}
                if sk_ctrl_subkey:
                    subkey_member['control_comparison']['by_subkey'] = sk_ctrl_subkey
                    for _f in ('neg_max_tass', 'neg_max_reads', 'tass_fold_over_neg',
                               'reads_fold_over_neg', 'tass_fold_over_pos',
                               'reads_fold_over_pos', 'control_flag',
                               'neg_control_values', 'pos_control_values',
                               'pos_min_tass', 'pos_min_reads'):
                        subkey_member['control_comparison'][_f] = sk_ctrl_subkey.get(_f)

                # Strain-level comparison (nested members of the subkey)
                for strain in subkey_member.get('members', []):
                    s_ctrl_key = compute_control_comparison(
                        data=strain, neg_index=_neg, pos_index=_pos,
                        fold_threshold=_ctrl_threshold, level="key",
                    )
                    s_ctrl_subkey = compute_control_comparison(
                        data=strain, neg_index=_neg, pos_index=_pos,
                        fold_threshold=_ctrl_threshold, level="subkey",
                    )
                    strain['control_comparison'] = {}
                    if s_ctrl_key:
                        strain['control_comparison']['by_key'] = s_ctrl_key
                    if s_ctrl_subkey:
                        strain['control_comparison']['by_subkey'] = s_ctrl_subkey
                    _best_ctrl = s_ctrl_key or s_ctrl_subkey
                    if _best_ctrl:
                        for _f in ('neg_max_tass', 'neg_max_reads', 'tass_fold_over_neg',
                                   'reads_fold_over_neg', 'tass_fold_over_pos',
                                   'reads_fold_over_pos', 'control_flag',
                                   'neg_control_values', 'pos_control_values',
                                   'pos_min_tass', 'pos_min_reads'):
                            strain['control_comparison'][_f] = _best_ctrl.get(_f)

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
    # We split by simulator type (ISS vs NanoSim) for per-type metrics tables.
    _missing_insilico = []

    def _classify_insilico_files(file_list):
        """Split insilico control files into simulator-type buckets."""
        buckets = {}
        for fpath in (file_list or []):
            bname = os.path.basename(fpath).lower()
            if '_insilico_iss' in bname or '_iss.' in bname:
                buckets.setdefault('iss', []).append(fpath)
            elif '_insilico_nanosim' in bname or '_nanosim.' in bname:
                buckets.setdefault('nanosim', []).append(fpath)
            else:
                buckets.setdefault('unknown', []).append(fpath)
        return buckets

    _isil_buckets = _classify_insilico_files(args.insilico_controls)
    # Also build a combined index for the merged Ctrl spark bar annotation
    _insilico_index = load_control_data(args.insilico_controls) if args.insilico_controls else None
    _per_type_missing = {}  # sim_type -> list of missing entries

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

            for subkey_member in group.get('members', []):
                # Subkey-level insilico comparison
                sk_isil = compute_control_comparison(
                    data=subkey_member, neg_index=_empty_neg, pos_index=_insilico_index,
                    fold_threshold=args.control_fold_threshold, level="subkey",
                )
                subkey_member['insilico_comparison'] = {}
                if sk_isil:
                    subkey_member['insilico_comparison']['by_subkey'] = sk_isil
                    for _f in ('pos_min_tass', 'pos_min_reads',
                               'tass_fold_over_pos', 'reads_fold_over_pos',
                               'pos_control_values'):
                        subkey_member['insilico_comparison'][_f] = sk_isil.get(_f)
                    subkey_member['insilico_comparison']['insilico_tass'] = sk_isil.get('pos_min_tass')
                    subkey_member['insilico_comparison']['insilico_reads'] = sk_isil.get('pos_min_reads')
                    subkey_member['insilico_comparison']['tass_fold_over_insilico'] = sk_isil.get('tass_fold_over_pos')
                    subkey_member['insilico_comparison']['reads_fold_over_insilico'] = sk_isil.get('reads_fold_over_pos')

                # Strain-level insilico comparison (nested members)
                for strain in subkey_member.get('members', []):
                    s_isil_key = compute_control_comparison(
                        data=strain, neg_index=_empty_neg, pos_index=_insilico_index,
                        fold_threshold=args.control_fold_threshold, level="key",
                    )
                    s_isil_subkey = compute_control_comparison(
                        data=strain, neg_index=_empty_neg, pos_index=_insilico_index,
                        fold_threshold=args.control_fold_threshold, level="subkey",
                    )
                    strain['insilico_comparison'] = {}
                    if s_isil_key:
                        strain['insilico_comparison']['by_key'] = s_isil_key
                    if s_isil_subkey:
                        strain['insilico_comparison']['by_subkey'] = s_isil_subkey
                    _best_isil = s_isil_key or s_isil_subkey
                    if _best_isil:
                        for _f in ('pos_min_tass', 'pos_min_reads',
                                   'tass_fold_over_pos', 'reads_fold_over_pos',
                                   'pos_control_values'):
                            strain['insilico_comparison'][_f] = _best_isil.get(_f)
                        strain['insilico_comparison']['insilico_tass'] = _best_isil.get('pos_min_tass')
                        strain['insilico_comparison']['insilico_reads'] = _best_isil.get('pos_min_reads')
                        strain['insilico_comparison']['tass_fold_over_insilico'] = _best_isil.get('tass_fold_over_pos')
                        strain['insilico_comparison']['reads_fold_over_insilico'] = _best_isil.get('reads_fold_over_pos')

        if _insilico_annotated:
            print(f"In-silico comparison: annotated {_insilico_annotated} species groups "
                  f"(insilico={len(args.insilico_controls or [])} files, "
                  f"fold_threshold={args.control_fold_threshold})")

        # ── Per-simulator-type annotation (for separate metrics tables) ────
        # Load each simulator type's files separately and compute per-type
        # comparison data stored under insilico_comparison_<type>
        for _sim_type, _sim_files in _isil_buckets.items():
            _type_index = load_control_data(_sim_files)
            if not _type_index:
                continue
            _type_key = f"insilico_comparison_{_sim_type}"
            for group in final_json:
                grp_isil_t = compute_control_comparison(
                    data=group, neg_index=_empty_neg, pos_index=_type_index,
                    fold_threshold=args.control_fold_threshold, level="toplevelkey",
                )
                if grp_isil_t:
                    group[_type_key] = grp_isil_t
                for subkey_member in group.get('members', []):
                    sk_isil_t = compute_control_comparison(
                        data=subkey_member, neg_index=_empty_neg, pos_index=_type_index,
                        fold_threshold=args.control_fold_threshold, level="subkey",
                    )
                    if sk_isil_t:
                        subkey_member[_type_key] = {
                            'insilico_tass': sk_isil_t.get('pos_min_tass'),
                            'insilico_reads': sk_isil_t.get('pos_min_reads'),
                            'tass_fold_over_insilico': sk_isil_t.get('tass_fold_over_pos'),
                            'reads_fold_over_insilico': sk_isil_t.get('reads_fold_over_pos'),
                            'pos_control_values': sk_isil_t.get('pos_control_values'),
                        }
                    for strain in subkey_member.get('members', []):
                        s_isil_t = compute_control_comparison(
                            data=strain, neg_index=_empty_neg, pos_index=_type_index,
                            fold_threshold=args.control_fold_threshold, level="key",
                        ) or compute_control_comparison(
                            data=strain, neg_index=_empty_neg, pos_index=_type_index,
                            fold_threshold=args.control_fold_threshold, level="subkey",
                        )
                        if s_isil_t:
                            strain[_type_key] = {
                                'insilico_tass': s_isil_t.get('pos_min_tass'),
                                'insilico_reads': s_isil_t.get('pos_min_reads'),
                                'tass_fold_over_insilico': s_isil_t.get('tass_fold_over_pos'),
                                'reads_fold_over_insilico': s_isil_t.get('reads_fold_over_pos'),
                                'pos_control_values': s_isil_t.get('pos_control_values'),
                            }

            # Per-type missing organisms
            _type_missing = find_missing_positive_controls(
                final_json, _type_index, levels=["key", "subkey"])
            if _type_missing:
                for entry in _type_missing:
                    entry['control_source'] = f'insilico_{_sim_type}'
                    entry['simulator_type'] = _sim_type
                    _cat = 'Unknown'
                    _entry_id = str(entry.get('id', ''))
                    _entry_name = entry.get('name', '')
                    _pinfo = pathogens.get(_entry_id) or pathogens.get(_entry_name)
                    if _pinfo:
                        _raw_cat = _pinfo.get('callclass', 'Unknown')
                        if _raw_cat:
                            _raw_lower = str(_raw_cat).strip().lower()
                            if 'primary' in _raw_lower:
                                _cat = 'Primary'
                            elif 'opportunistic' in _raw_lower:
                                _cat = 'Opportunistic'
                            elif 'potential' in _raw_lower:
                                _cat = 'Potential'
                            elif 'commensal' in _raw_lower:
                                _cat = 'Commensal'
                    entry['microbial_category'] = _cat
                _per_type_missing[_sim_type] = _type_missing

            # Per-type sample-only organisms
            _type_sample_only = _find_sample_only_organisms(final_json, _type_index)
            if _type_sample_only:
                for group in final_json:
                    tlk = str(group.get('toplevelkey', group.get('key', '')))
                    if tlk and tlk in _type_sample_only:
                        if _type_key not in group:
                            group[_type_key] = {}
                        if isinstance(group[_type_key], dict):
                            group[_type_key]['missing_from_insilico'] = True

        # Detect organisms present in insilico but missing from sample (combined)
        # Use key/subkey levels for species-level resolution (not toplevelkey/genus)
        _missing_insilico = find_missing_positive_controls(
            final_json, _insilico_index, levels=["key", "subkey"])
        if _missing_insilico:
            # Tag as insilico-origin and enrich with microbial category from pathogens dict
            for entry in _missing_insilico:
                entry['control_source'] = 'insilico'
                # Look up microbial category using taxid or name from pathogens dict
                _cat = 'Unknown'
                _entry_id = str(entry.get('id', ''))
                _entry_name = entry.get('name', '')
                _pinfo = pathogens.get(_entry_id) or pathogens.get(_entry_name)
                if _pinfo:
                    _raw_cat = _pinfo.get('callclass', 'Unknown')
                    # Normalize to standard category names
                    if _raw_cat:
                        _raw_lower = str(_raw_cat).strip().lower()
                        if 'primary' in _raw_lower:
                            _cat = 'Primary'
                        elif 'opportunistic' in _raw_lower:
                            _cat = 'Opportunistic'
                        elif 'potential' in _raw_lower:
                            _cat = 'Potential'
                        elif 'commensal' in _raw_lower:
                            _cat = 'Commensal'
                entry['microbial_category'] = _cat
            print(f"Missing in-silico controls: {len(_missing_insilico)} organism(s) "
                  f"in in-silico control not found in sample")

        # Also detect organisms in sample but missing from insilico (combined)
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

    # ── Protein annotation enrichment ─────────────────────────────────────────
    # If an annotate_report TSV is provided (from annotate_report.py), load it
    # and attach matching annotations to each organism group/member by
    # species_taxid ↔ subkey.  Also try genus_taxid ↔ toplevelkey fallback.
    if args.annotate_report and os.path.exists(args.annotate_report):
        print(f"Loading protein annotation report from {args.annotate_report}")
        if args.annotate_report.endswith('.xlsx') or args.annotate_report.endswith('.xls'):
            _annot_df = pd.read_excel(args.annotate_report, dtype=str)
        elif args.annotate_report.endswith('.csv'):
            _annot_df = pd.read_csv(args.annotate_report, dtype=str)
        else:
            _annot_df = pd.read_csv(args.annotate_report, sep='\t', dtype=str)
        _annot_df.columns = _annot_df.columns.str.strip()
        # Build lookup: species_taxid → list of annotation dicts
        _annot_by_species = defaultdict(list)
        _annot_by_genus = defaultdict(list)
        _annot_cols = [
            'sseqid', 'pident', 'evalue', 'bitscore', 'source_id',
            'gene_name', 'product', 'classification', 'antibiotics_class',
            'antibiotics', 'organism', 'genus', 'species', 'property',
            'source', 'level', 'host_name',
        ]
        _annot_numeric_cols = {'pident', 'evalue', 'bitscore'}
        def _coerce(col, val):
            if pd.isna(val):
                return None
            if col in _annot_numeric_cols:
                try:
                    return float(val)
                except (ValueError, TypeError):
                    return val
            return val
        _available_cols = [c for c in _annot_cols if c in _annot_df.columns]
        for _, row in _annot_df.iterrows():
            entry = {c: _coerce(c, row.get(c)) for c in _available_cols}
            sp_taxid = str(row.get('species_taxid', '')).strip()
            gn_taxid = str(row.get('genus_taxid', '')).strip()
            if sp_taxid and sp_taxid != 'nan':
                _annot_by_species[sp_taxid].append(entry)
            if gn_taxid and gn_taxid != 'nan':
                _annot_by_genus[gn_taxid].append(entry)

        _annot_matched = 0
        for group in final_json:
            tlk = str(group.get('toplevelkey', group.get('key', '')))
            # Attach at genus (toplevel) level
            if tlk in _annot_by_genus:
                # Deduplicate by gene_name
                _seen_genes = set()
                _unique = []
                for a in _annot_by_genus[tlk]:
                    gn = a.get('gene_name', '')
                    if gn not in _seen_genes:
                        _seen_genes.add(gn)
                        _unique.append(a)
                group['protein_annotations_genus'] = _unique

            for subkey_member in group.get('members', []):
                sk = str(subkey_member.get('subkey', subkey_member.get('key', '')))
                if sk in _annot_by_species:
                    _seen_genes = set()
                    _unique = []
                    for a in _annot_by_species[sk]:
                        gn = a.get('gene_name', '')
                        if gn not in _seen_genes:
                            _seen_genes.add(gn)
                            _unique.append(a)
                    subkey_member['protein_annotations'] = _unique
                    _annot_matched += 1

                # Also propagate to individual strains via taxdump hierarchy
                for strain in subkey_member.get('members', []):
                    strain_key = str(strain.get('key', ''))
                    # Try direct match first, then walk up to species via taxdump
                    matched_annots = []
                    if strain_key in _annot_by_species:
                        matched_annots = _annot_by_species[strain_key]
                    elif taxdump and strain_key in taxdump:
                        # Walk up the taxonomy to find species-level match
                        _current = strain_key
                        _visited = set()
                        while _current and _current not in _visited:
                            _visited.add(_current)
                            if _current in _annot_by_species:
                                matched_annots = _annot_by_species[_current]
                                break
                            _parent = taxdump.get(_current, {}).get('parent', '')
                            if _parent == _current:
                                break
                            _current = _parent
                    if matched_annots:
                        _seen_genes = set()
                        _unique = []
                        for a in matched_annots:
                            gn = a.get('gene_name', '')
                            if gn not in _seen_genes:
                                _seen_genes.add(gn)
                                _unique.append(a)
                        strain['protein_annotations'] = _unique

        print(f"Protein annotations: matched {_annot_matched} species/subkey groups, "
              f"{len(_annot_by_species)} species taxids in annotation file, "
              f"{len(_annot_by_genus)} genus taxids")
    elif args.annotate_report:
        print(f"WARNING: annotate_report file not found: {args.annotate_report}")

    # ── Build structured output with metadata ────────────────────────────────
    _all_keys = set()
    _all_subkeys = set()
    _all_toplevelkeys = set()
    _total_organism_reads = 0
    for grp in final_json:
        _all_toplevelkeys.add(grp.get('toplevelkey', grp.get('key', '')))
        for sk_m in grp.get('members', []):
            _all_subkeys.add(sk_m.get('subkey', sk_m.get('key', '')))
            for strain in sk_m.get('members', []):
                _all_keys.add(strain.get('key', ''))
                _total_organism_reads += float(strain.get('numreads', 0))
    # Strip covered_regions from JSON output — large and only needed internally
    import copy as _copy
    _json_out = _copy.deepcopy(final_json)
    for _grp in _json_out:
        _grp.pop('covered_regions', None)
        for _sk_m in _grp.get('members', []):
            _sk_m.pop('covered_regions', None)
            for _strain in _sk_m.get('members', []):
                _strain.pop('covered_regions', None)
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
            "insilico_simulator_types": list(_isil_buckets.keys()) if _isil_buckets else None,
            "missing_insilico_controls": _missing_insilico if _missing_insilico else None,
            "missing_insilico_by_type": _per_type_missing if _per_type_missing else None,
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
            minhash_score = entry.get('minhash_reduction', entry.get('minhash_score', 0))
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
                print(f"\tMinhash score: {entry.get('minhash_reduction', entry.get('minhash_score', 0)):.2f}")
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
