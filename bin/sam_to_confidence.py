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

import sys
import os
import pysam
import math
import argparse


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Merge a paf to confidence output mqc.tsv with the kraken report it corresponds to include tax info",
        epilog="Example: ",
    )
    parser.add_argument(
        "-i",
        "--bamfile",
        metavar="BAMFILE",
        type=str,
        required=True,
        help="BAM alignment file",
    )
    parser.add_argument(
        "-d",
        "--depth",
        metavar="DEPTHFILE",
        type=str,
        default=None,
        required=True,
        help="OPTIONAL, samplename analyzed",
    )
    parser.add_argument(
        "-m",
        "--mapping",
        metavar="MAPPING",
        default=None,
        required=False,
        help="Output accession to name tsv. 2 columns minimum. must be in 1 col for accession second is name",
    )
    parser.add_argument(
        "-s",
        "--splitidx",
        metavar="SPLITIDX",
        default=None,
        required=False,
        type=int,
        help="If using -m mapping, split the index of the chr and get element. Split on | by default. 0 based index.",
    )
    parser.add_argument(
        "-o",
        "--file_out",
        metavar="FILE_OUT",
        type=str,
        required=True,
        help="Name of the output tsv file containing a mapping of top n organisms at individual taxa levels",
    )

    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)
# Load the depth file into a dictionary


def load_depth_file(depth_file):
    depth = {}
    total_positions = {}


    with open(depth_file, 'r') as f:
        for line in f:
            cols = line.strip().split()
            ref, pos, dp = cols[0], int(cols[1]), int(cols[2])

            if ref not in depth:
                depth[ref] = {}
                total_positions[ref] = 0

            depth[ref][pos] = dp
            total_positions[ref] += 1

    return depth, total_positions

# Process BAM file


def process_bam_file(bam_file, depth, total_positions):
    thresholds = [1, 10, 50, 100, 300]
    coverage = {}
    count = {}
    sumsq = {}
    ireads = {}
    countreads = {}
    ilen = {}

    total_reads_aligned = 0

    with pysam.AlignmentFile(bam_file, 'rb') as f:
        for read in f.fetch():
            if read.is_unmapped:
                continue

            ref_name = read.reference_name
            start = read.reference_start
            end = read.reference_end

            # Update coverage and counts
            if ref_name not in coverage:
                coverage[ref_name] = 0
                count[ref_name] = 0
                ireads[ref_name] = 0
                countreads[ref_name] = 0
            ireads[ref_name] += 1
            total_reads_aligned += 1

            for pos in range(start, end):
                coverage[ref_name] = coverage.get(ref_name, 0) + 1
                count[ref_name] = count.get(ref_name, 0) + 1

                if ref_name in depth and pos in depth[ref_name]:
                    dp = depth[ref_name][pos]
                    sumsq[ref_name] = sumsq.get(ref_name, 0) + (dp ** 2)

        for ref in f.references:
            ilen[ref] = f.lengths[f.references.index(ref)]

    return coverage, count, sumsq, ireads, ilen, total_reads_aligned


def fm(x):
    # if the number is 0, then return 0 else return the number as is
    if x == 0:
        return "0"
    # elif x > 0 and x < 0.01:
    #     return "<0.01"
    else:
        return f"{x:.2f}"
# Calculate statistics and output


def output_statistics(depth, total_positions, coverage, count, sumsq, ireads, ilen, total_reads_aligned, FILE_OUT, mapping, splitidx):
    headers = ["Accession",
        "Organism",
        "Full Name",
        "Mean Depth",
        "Average Coverage",
        "Ref Size",
        "Total Reads Aligned",
        "% Reads Aligned",
        "Stdev",
        "Abundance Aligned",
        "1:10:50:100:300X Cov."
    ]
    # print the statistics to a file
    data = []
    with open(FILE_OUT, "w") as f:
        f.write("\t".join(headers)+"\n")
        for ref in ilen:
            if ref in coverage and ref in ilen:
                avg_coverage = coverage[ref] / ilen[ref]

                if ref in depth:
                    total_depth = sum(depth[ref].values())
                    avgDepth = total_depth / total_positions[ref]

                    newThresholds = {t: sum(1 for v in depth[ref].values() if v > t) for t in [
                        1, 10, 50, 100, 300]}
                    xCov = [100 * newThresholds[t] / ilen[ref]
                            for t in [1, 10, 50, 100, 300]]
                    stdev = math.sqrt(
                        abs(sumsq[ref] / count[ref] - (total_depth / count[ref]) ** 2))
                    abu_total_aligned = ireads[ref] / total_reads_aligned
                    check_idx = ref
                    if splitidx or splitidx == 0:
                        check_idx = check_idx.split("|")[splitidx]
                    if check_idx in mapping:
                        organism = mapping[check_idx]['organism']
                        fullname = mapping[check_idx]['fullname']
                    else:
                        organism = "N/A"
                        fullname = "N/A"
                    outstring = f"{check_idx}\t{organism}\t{fullname}\t{fm(avgDepth)}\t{fm(avg_coverage)}\t{ilen[ref]}\t{ireads[ref]}\t{fm(100*(ireads[ref]/total_reads_aligned))}\t{fm(stdev)}\t{fm(abu_total_aligned)}\t{':'.join([fm(x) for x in xCov])}"
                    f.write(outstring+"\n")


def main(argv=None):

    args = parse_args(argv)
    DEPTHFILE = args.depth
    BAMFILE = args.bamfile
    mapping = {}

    if args.mapping:
        # import mapping file. Split each line on "\t" and store in a dictionary
        # where the key is the accession and the value is the name
        with open(args.mapping, 'r') as f:
            for line in f:
                cols = line.strip().split("\t")
                idx = cols[0]
                organism = "N/A"
                fullname = "N/A"
                obj = dict(organism=organism, fullname=fullname)
                if len(cols) > 2:
                    obj['organism'] = cols[2]
                if len(cols) > 3:
                    obj['fullname'] = cols[3]
                mapping[idx] = obj
        f.close()
    depth, total_positions = load_depth_file(DEPTHFILE)

    def ensure_bam_index(bam_path):
        if not os.path.exists(bam_path + ".csi"):
            pysam.index(bam_path, csi = True, force = True)


    # Usage
    # ensure_bam_index(BAMFILE)
    coverage, count, sumsq, ireads, ilen, total_reads_aligned = process_bam_file(
        BAMFILE, depth, total_positions)
    output_statistics(depth, total_positions, coverage, count,
        sumsq, ireads, ilen, total_reads_aligned, args.file_out, mapping, args.splitidx)


if __name__ == "__main__":
    sys.exit(main())
