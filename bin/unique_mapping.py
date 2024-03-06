import pysam
import sys
import numpy as np

import argparse

def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and transform a tabular samplesheet.",
        epilog="Example: python check_samplesheet.py samplesheet.csv samplesheet.valid.csv",
    )
    parser.add_argument(
        "-i",
        "--input",
        metavar="BAM",
        help="BAM File to process",
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="FILE",
        help="Output file",
    )
    return parser.parse_args(argv)

def coverage_stats_per_reference(bam_path):
    # Open the BAM file without using an index
    bamfile = pysam.AlignmentFile(bam_path, "rb", index_filename=None)

    # Initialize variables
    reference_coverage_stats = {}

    # Function to process reads
    def process_reads(reads, ref_name):
        coverage = [0] * (max(read.reference_end for read in reads if not read.is_unmapped) + 1)
        for read in reads:
            if not read.is_unmapped:
                for pos in range(read.reference_start, read.reference_end):
                    coverage[pos] += 1
        mean_coverage = np.mean(coverage) if coverage else 0
        std_dev_coverage = np.std(coverage) if coverage else 0
        reference_coverage_stats[ref_name] = (mean_coverage, std_dev_coverage)

    # Group reads by reference
    current_ref = None
    grouped_reads = []
    for read in bamfile:
        if read.reference_name != current_ref:
            if grouped_reads:
                process_reads(grouped_reads, current_ref)
            grouped_reads = []
            current_ref = read.reference_name
        grouped_reads.append(read)

    # Process the last group
    if grouped_reads:
        process_reads(grouped_reads, current_ref)

    bamfile.close()

    return reference_coverage_stats


def evaluate_alignment_improved(mean_coverage, std_dev_coverage, target_coverage=30, max_std_dev=20, scale_factor=1000):
    """
    Evaluate alignment quality with adjusted sensitivity for low coverage values.

    :param mean_coverage: Mean coverage across the reference, expected to be very low in this context.
    :param std_dev_coverage: Standard deviation of the coverage.
    :param target_coverage: Expected target mean coverage, not directly used here due to the scaling for low coverage.
    :param max_std_dev: Maximum acceptable standard deviation for optimal alignment.
    :param scale_factor: Factor to scale up low coverage values to make the function more responsive.
    :return: Alignment quality score between 0 and 1.
    """

    # Scale and normalize the mean coverage to exaggerate the differences at low coverage
    scaled_coverage = mean_coverage * scale_factor
    coverage_score = np.exp(-((scaled_coverage - target_coverage / scale_factor) ** 2) / (2 * ((target_coverage / scale_factor) / 10) ** 2))

    # Inverse linear normalization for standard deviation
    std_dev_score = max(0, 1 - (std_dev_coverage / max_std_dev))

    # Combine metrics
    alignment_quality = (coverage_score + std_dev_score) / 2

    return alignment_quality

# Example usage
if __name__ == "__main__":
    args = parse_args()
    bam_file = args.input
    # uniqueness = calculate_uniqueness(bam_file)
    coverage_stats = coverage_stats_per_reference(bam_file)

    for ref, stats in coverage_stats.items():
        align_quality = evaluate_alignment_improved(stats[0], stats[1], 30, 20, 2000)
        print(f"{ref} - MeanCov: {stats[0]:.2f}, CovStdDev: {stats[1]:.2f}, QualConf: {align_quality:.2f}")

