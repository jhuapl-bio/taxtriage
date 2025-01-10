#!/usr/bin/env python3

import argparse
import os
import subprocess
import csv
import pysam
from sourmash import MinHash, SourmashSignature, save_signatures, load_file_as_signatures

from collections import defaultdict

def apply_proportional_scaling(conflict_sets, region_coverage, all_reads, fetch_reads_in_region, passed_regions):
    reads_to_keep = set()

    if not conflict_sets:
        # No conflict sets at all, just keep all reads from passed_regions.
        for reg in passed_regions:
            r, s, e = reg
            reads = fetch_reads_in_region(all_reads, r, s, e)
            for (read_id, refname, seq) in reads:
                reads_to_keep.add(read_id)
        return reads_to_keep

    # If we do have conflict sets:
    for conflict_set in conflict_sets:
        coverages = [region_coverage.get(reg,0) for reg in conflict_set]
        if not coverages:
            # No coverage info, keep all reads from these regions in this set
            for reg in conflict_set:
                r, s, e = reg
                reads = fetch_reads_in_region(all_reads, r, s, e)
                for (read_id, refname, seq) in reads:
                    reads_to_keep.add(read_id)
            continue

        min_cov = min(coverages)
        # Compute fractions
        fractions = {}
        for reg in conflict_set:
            c = region_coverage.get(reg,0)
            if c > 0:
                fractions[reg] = min_cov / c
                ref, s, end = reg
                if ref == "NC_066642.1":
                    print(conflict_set)
                    print(f"Region: {reg}, Coverage: {c}, Fraction: {fractions[reg]}")
                    print(coverages,"<")
            else:
                fractions[reg] = 1.0
        # Apply scaling
        for reg in conflict_set:
            (r, s, e) = reg
            reads = fetch_reads_in_region(all_reads, r, s, e)
            fraction = fractions[reg]
            # Ensure at least 1 if fraction > 0, or skip if fraction=0
            if fraction > 0 and len(reads) > 0:
                num_to_keep = max(1, int(len(reads)*fraction + 0.9999))
                sampled_reads = reads[:num_to_keep]
                for (read_id, refname, seq) in sampled_reads:
                    reads_to_keep.add(read_id)

    # Now, for regions that are passed and not in any conflict sets, we still need to keep their reads.
    # Identify all regions that are not part of any conflict set:
    conflict_region_set = set().union(*conflict_sets) if conflict_sets else set()
    non_conflict_regions = [reg for reg in passed_regions if reg not in conflict_region_set]
    for reg in non_conflict_regions:
        r, s, e = reg
        reads = fetch_reads_in_region(all_reads, r, s, e)
        for (read_id, refname, seq) in reads:
            reads_to_keep.add(read_id)

    return reads_to_keep
def transform_sum_comparisons(sum_comparisons):
    """
    Transform from:
    {
      'NC_066642.1': [
        {'jaccard': ..., 'to': 'NC_003310.1', 's1':..., 'e1':..., 's2':..., 'e2':...},
        ...
      ],
      ...
    }
    to:
    {
      (ref, s1, e1): [
        {'jaccard': ..., 'to': to_ref, 's2':..., 'e2':..., 's1':..., 'e1':...}
      ],
      ...
    }
    """
    new_sum = defaultdict(list)
    for ref, conflicts in sum_comparisons.items():
        for c in conflicts:
            # Construct the key region for 'ref' using s1, e1 from the conflict dict
            region = (ref, c['s1'], c['e1'])
            new_sum[region].append(c)
    return dict(new_sum)

def build_conflict_graph(passed_regions, sum_comparisons, min_threshold=0.1):
    """
    Build a graph where each node is a region (ref, start, end) and edges represent conflicts.
    Only consider edges if both regions are in passed_regions and jaccard >= min_threshold.
    """
    graph = defaultdict(set)
    for region, conflicts in sum_comparisons.items():

        if region not in passed_regions:
            continue
        for c in conflicts:
            jaccard = c['jaccard']
            if jaccard >= min_threshold:
                to_ref = c['to']
                s2 = c['s2']
                e2 = c['e2']
                other_region = (to_ref, s2, e2)
                if other_region in passed_regions:
                    graph[region].add(other_region)
                    graph[other_region].add(region)
    return graph

def identify_conflict_sets(graph):
    visited = set()
    conflict_sets = []

    def dfs(node, comp):
        comp.add(node)
        visited.add(node)
        for neigh in graph[node]:
            if neigh not in visited:
                dfs(neigh, comp)

    for node in graph:
        if node not in visited:
            comp = set()
            dfs(node, comp)
            conflict_sets.append(comp)

    return conflict_sets

def compute_region_coverage(passed_regions, all_reads, fetch_reads_in_region):
    region_coverage = {}
    for (ref, start, end) in passed_regions:
        reads = fetch_reads_in_region(all_reads, ref, start, end)
        coverage_approx = len(reads)
        region_coverage[(ref, start, end)] = coverage_approx
    return region_coverage

def compute_coverage_with_bedtools(bam_path, outpath=None, bedtools_path='bedtools'):
    """
    Use bedtools genomecov to compute per-base coverage.
    Returns the output as a string.
    """

    # remove the etensions only using os path
    cmd = [bedtools_path, 'genomecov', '-ibam', bam_path, '-d', ]
    coverage_output = run_subprocess(cmd, "bedtools genomecov",
                        stdout_file = outpath,
                        read_stdout=True if not outpath else False
                    )
    return coverage_output

def coverage_to_bed(coverage_output, min_coverage=1):
    """
    Convert bedtools genomecov output to BED format for positions with coverage >= min_coverage.
    Returns a list of BED intervals as tuples: (chrom, start, end)
    """
    bed_intervals = []
    for line in coverage_output.strip().split('\n'):
        parts = line.strip().split()
        if len(parts) < 3:
            continue
        chrom, pos, depth = parts[0], int(parts[1])-1, int(parts[2])  # bedtools genomecov is 1-based
        if depth >= min_coverage:
            bed_intervals.append((chrom, pos, pos+1))  # Each covered base as a separate interval
    return bed_intervals

def write_bed_file(bed_intervals, bed_file_path):
    """
    Write BED intervals to a BED file.
    """
    with open(bed_file_path, 'w') as f:
        for interval in bed_intervals:
            chrom, start, end = interval
            f.write(f"{chrom}\t{start}\t{end}\n")

def merge_bed_with_gap(bed_file_path, merged_bed_path, coverage_spacing, bedtools_path='bedtools'):
    """
    Merge BED intervals allowing gaps up to coverage_spacing.
    """
    cmd = [bedtools_path, 'merge', '-i', bed_file_path, '-d', str(coverage_spacing)]
    merged_output = run_subprocess(cmd, "bedtools merge", read_stdout=False, stdout_file=merged_bed_path)
    return
def split_regions_into_windows(merged_bed_path, windows_bed_path, coverage_length, bedtools_path='bedtools'):
    """
    Split merged BED regions into fixed-length windows.
    """
    cmd = [bedtools_path, 'makewindows', '-b', merged_bed_path, '-w', str(coverage_length)]
    run_subprocess(cmd, "bedtools makewindows", read_stdout=False, stdout_file=windows_bed_path)

def filter_windows_with_coverage(windows_bed_path, bed_file_path, filtered_bed_path, bedtools_path='bedtools'):
    """
    Use bedtools intersect to keep only windows that overlap with coverage regions.
    """
    cmd = [bedtools_path, 'intersect', '-a', windows_bed_path, '-b', bed_file_path, '-wa']
    run_subprocess(cmd, "bedtools intersect", stdout_file=filtered_bed_path, read_stdout=False)

def parse_bed_file(bed_file_path):
    """
    Parse a BED file and return a list of regions as tuples: (chrom, start, end)
    """
    regions = []
    with open(bed_file_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            chrom, start, end, depth = parts[0], int(parts[1]), int(parts[2]), int(parts[3])
            regions.append((chrom, start, end, depth))
    return regions

def load_reads_from_bam(bam_path):
    """
    Load all eligible reads from a BAM file into memory.
    Returns a dictionary keyed by reference name:
    {
      refname: [
        {
           'id': read_name,
           'seq': read_sequence,
           'start': ref_start,
           'aligned_ref': ....,
           'end': ref_end
        }, ...
      ]
    }
    """
    bam = pysam.AlignmentFile(bam_path, "rb")
    all_reads = {}
    for read in bam:
        if read.is_unmapped:
            continue
        if read.is_secondary or read.is_supplementary:
            continue

        refname = bam.get_reference_name(read.reference_id)
        seq = read.query_sequence
        if seq:
            if refname not in all_reads:
                all_reads[refname] = []
            all_reads[refname].append({
                'id': read.query_name,
                'seq': seq,
                'start': read.reference_start,
                'end': read.reference_end
            })
    bam.close()
    return all_reads

def fetch_reads_in_region(all_reads, refname, start, end):
    """
    Given the pre-loaded reads dictionary, return all reads overlapping [start, end)
    of the given refname. Overlap means the read's alignment region intersects [start, end).

    start and end are 0-based coordinates.
    """
    results = []
    if refname not in all_reads:
        return results
    for read in all_reads[refname]:
        # Check overlap
        # A read overlaps the region if read.end > start and read.start < end
        if read['end'] > start and read['start'] < end:
            results.append((read['id'], refname, read['seq']))
    return results

def create_signatures_for_regions(regions, all_reads, bam_path, kmer_size, scaled):
    """
    For each region, fetch reads and create a signature.
    Returns a dictionary of the form:
    {
        (refname, start, end): SourmashSignature,
        ...
    }
    """
    signatures = {}
    seenrefs = set()
    sets = dict()
    for (refname, start, end, depth, val) in regions:
        combined = f"{refname}:{start}-{end}"
        if combined not in seenrefs:
            # print(f"Processing {combined}...")
            seenrefs.add(combined)
        else:
            continue
        reads = fetch_reads_in_region(all_reads, refname, start, end)
        if len(reads) == 0:
            # No reads in this region
            continue
        # Create a single MinHash for this region
        mh = MinHash(n=0, ksize=kmer_size, scaled=scaled)
        for read_id, refname, read_seq in reads:
            mh.add_sequence(read_seq, force=True)
        # Create signature
        region_name = f"{refname}:{start+1}-{end}"  # 1-based coords for display
        sig = SourmashSignature(mh, name=region_name)
        name=f"{refname}:{start+1}-{end}"
        signatures[name] = sig
    return signatures

def compare_signatures(signatures, output_csv, min_threshold=0.1, clusters=None):
    """
    signatures: dict of region_string -> SourmashSignature
                where region_string looks like "NC_003310.1:1234-5678"
    output_csv: path to CSV file for saving results
    min_threshold: float, Jaccard threshold for reporting comparisons
    clusters: list of sets, e.g. [
                 {'NC_008291.1', 'NC_006998.1', 'NC_066642.1', 'NC_001611.1'},
                 {'NC_003310.1'},
                 {'NC_055230.1'}
              ]
    """

    # ------------------------------------------------
    # 1) Build a reference->cluster map if clusters given
    # ------------------------------------------------
    # Example: ref2cluster['NC_008291.1'] = 0  (meaning cluster 0)
    #          ref2cluster['NC_003310.1'] = 1
    #          ref2cluster['NC_055230.1'] = 2
    ref2cluster = {}
    if clusters is not None:
        for cluster_idx, cluster_set in enumerate(clusters):
            for ref in cluster_set:
                ref2cluster[ref] = cluster_idx

    # ------------------------------------------------
    # 2) Prepare data structures
    # ------------------------------------------------
    siglist = list(signatures.items())  # [ (region_string, signature), ... ]

    # Initialize passed_regions with all region keys.
    all_regions = set(signatures.keys())
    passed_regions = set(all_regions)

    sum_comparisons = defaultdict(list)

    # Utility to parse "ref:start-end" into (ref, start, end)
    def parse_split(region_str):
        # region_str example: "NC_003310.1:1234-5678"
        # Split on the last colon, then split on '-'
        ref_part, coords = region_str.rsplit(":", 1)
        s, e = coords.split("-")
        return ref_part, int(s), int(e)

    # ------------------------------------------------
    # 3) Compare signatures pairwise
    # ------------------------------------------------
    with open(output_csv, "w", newline="") as outfh:
        writer = csv.writer(outfh)
        writer.writerow([
            "reference1","start1","end1",
            "reference2","start2","end2",
            "jaccard","containment_1_in_2","containment_2_in_1"
        ])

        for i in range(len(siglist)):
            region1, sig1 = siglist[i]
            r1, s1, e1 = parse_split(region1)
            mh1 = sig1.minhash

            for j in range(i+1, len(siglist)):
                region2, sig2 = siglist[j]
                r2, s2, e2 = parse_split(region2)

                # 3A) Skip if same reference (as in your original code)
                if r2 == r1:
                    continue

                # 3B) Skip if not in the same cluster
                if clusters is not None:
                    # If either reference isn't in the map, or they differ in cluster, skip
                    if (r1 not in ref2cluster or
                        r2 not in ref2cluster or
                        ref2cluster[r1] != ref2cluster[r2]):
                        continue

                # 3C) Compute similarities
                mh2 = sig2.minhash
                jaccard = mh1.jaccard(mh2)
                if jaccard >= min_threshold:
                    c1_in_2 = mh1.avg_containment(mh2)
                    c2_in_1 = mh2.avg_containment(mh1)

                    # Write to CSV
                    writer.writerow([r1, s1, e1, r2, s2, e2, jaccard, c1_in_2, c2_in_1])

                    # Record comparison
                    sum_comparisons[r1].append(dict(
                        jaccard=jaccard, to=r2,
                        s2=s2, e2=e2, s1=s1, e1=e1
                    ))
                    sum_comparisons[r2].append(dict(
                        jaccard=jaccard, to=r1,
                        s2=s1, e2=e1, s1=s2, e1=e2
                    ))

                    # (Optional) If you still want to remove from passed_regions
                    # if (region1 in passed_regions):
                    #     passed_regions.remove(region1)
                    # if (region2 in passed_regions):
                    #     passed_regions.remove(region2)

    # ------------------------------------------------
    # 4) Print aggregated comparisons
    # ------------------------------------------------
    # for ref, jaccards in sum_comparisons.items():
    #     print(f"{ref}:")
    #     for entry in jaccards:
    #         print(
    #             f"\t{entry['s1']}-{entry['e1']}\t"
    #             f"Jaccard: {entry['jaccard']:.4f}\n"
    #             f"\t\tto {entry['to']}:{entry['s2']}-{entry['e2']}"
    #         )

    return passed_regions, sum_comparisons
def create_filtered_bam(original_bam_path, output_bam_path, read_ids_to_keep):
    """
    Create a new BAM file containing only reads with IDs in read_ids_to_keep.
    """
    bam_in = pysam.AlignmentFile(original_bam_path, "rb")
    bam_out = pysam.AlignmentFile(output_bam_path, "wb", template=bam_in)

    count_total = 0
    count_written = 0
    for read in bam_in:
        count_total += 1
        if read.query_name not in read_ids_to_keep:
            bam_out.write(read)
            count_written += 1

    bam_in.close()
    bam_out.close()
    print(f"Total reads processed: {count_total}")
    print(f"Reads written to filtered BAM: {count_written}")

def get_references_from_bam(bam_path):
    """
    Return a set of references that have at least one read in the bam.
    """
    bam = pysam.AlignmentFile(bam_path, "rb")
    ref_presence = {r:False for r in bam.references}
    # A quick check: if a single read maps to the reference, mark it True
    for read in bam:
        if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:
            refname = bam.get_reference_name(read.reference_id)
            ref_presence[refname] = True
    bam.close()
    # Return only those references that are present
    present_refs = {r for r, present in ref_presence.items() if present}
    return present_refs

def compute_performance(ground_truth, filtered_bam_path):
    """
    Compute precision, recall, F1-score using ground truth and predictions.

    Ground truth: If coverage > 0 => True presence, else False.
    Prediction: If filtered bam has reads for that ref => True presence, else False.

    Metrics:
    - TP: ground truth = True, predicted = True
    - FP: ground truth = False, predicted = True
    - FN: ground truth = True, predicted = False
    """
    ground_truth_presence = {ref:(cov>0) for ref,cov in ground_truth.items()}
    predicted_refs = get_references_from_bam(filtered_bam_path)
    # Predicted presence is True if ref is in predicted_refs, else False

    TP = 0
    FP = 0
    FN = 0

    # Consider only references we have ground truth for
    for ref, truth_present in ground_truth_presence.items():
        pred_present = (ref in predicted_refs)
        if truth_present and pred_present:
            TP += 1
        elif not truth_present and pred_present:
            FP += 1
        elif truth_present and not pred_present:
            FN += 1
    # Precision = TP/(TP+FP) if TP+FP>0
    # Recall = TP/(TP+FN) if TP+FN>0
    # F1 = 2*(Precision*Recall)/(Precision+Recall)

    precision = TP/(TP+FP) if (TP+FP)>0 else 0
    recall = TP/(TP+FN) if (TP+FN)>0 else 0
    f1 = 2*(precision*recall)/(precision+recall) if (precision+recall)>0 else 0

    return TP, FP, FN, precision, recall, f1

def compute_breadth_and_depth(bam_path, outpath=None, bedtools_path='bedtools'):
    """
    Compute breadth of coverage and mean depth for each reference in the BAM file.

    Returns a dictionary:
    {
        refname: {
            'breadth': float,  # proportion of genome covered
            'mean_depth': float
        },
        ...
    }
    """
    coverage_output = compute_coverage_with_bedtools(bam_path, outpath=outpath, bedtools_path=bedtools_path)

    coverage_stats = {}
    # Get reference lengths from BAM header
    bam = pysam.AlignmentFile(bam_path, "rb")
    ref_lengths = {ref: length for ref, length in zip(bam.references, bam.lengths)}
    bam.close()

    # Initialize stats
    for ref in ref_lengths:
        coverage_stats[ref] = {'covered_positions': 0, 'total_coverage': 0}

    for line in coverage_output.strip().split('\n'):
        parts = line.strip().split()
        if len(parts) < 3:
            continue
        chrom, pos, depth = parts[0], int(parts[1])-1, int(parts[2])
        if chrom not in coverage_stats:
            continue
        coverage_stats[chrom]['total_coverage'] += depth
        if depth >=1:
            coverage_stats[chrom]['covered_positions'] +=1

    # Now, compute breadth and mean depth
    final_stats = {}
    for ref, stats in coverage_stats.items():
        length = ref_lengths[ref]
        if length ==0:
            final_stats[ref] = {'breadth':0, 'mean_depth':0}
            continue
        breadth = stats['covered_positions'] / length
        mean_depth = stats['total_coverage'] / length
        final_stats[ref] = {'breadth': breadth, 'mean_depth': mean_depth}

    return final_stats

def parse_ground_truth(abu_file):
    """
    Parse the ground truth coverage file.
    File format:
    refname <tab> coverage
    coverage is a float, representing mean depth per position across the genome.
    """
    gt = {}
    with open(abu_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            ref, cov = parts[0], float(parts[1])
            gt[ref] = cov
    return gt

import random


def finalize_removed_reads(remove_counts, all_reads, fetch_reads_in_region,
                          remove_mode='first',
                          random_seed=None):
    """
    Convert the 'remove_counts' constraints into a global set of read IDs to remove,
    and also gather stats on how many reads are removed per reference.

    Args:
        remove_counts: dict keyed by (ref, start, end) -> int
                       indicating how many reads to remove from that region.
        all_reads: dict from refname -> list of read dicts, as loaded from the BAM.
        fetch_reads_in_region: function that returns a list of (read_id, ref, seq)
                               for the given region.
        remove_mode: 'first' or 'random' determines how we pick which reads to remove.
        random_seed: optional seed for reproducible random removal.

    Returns:
        removed_read_ids: a set of read IDs that are removed globally.
        removal_stats: a dict keyed by refname -> number_of_reads_removed
    """

    # Optionally set a random seed if we want reproducible sampling
    if random_seed is not None:
        random.seed(random_seed)

    removed_read_ids = set()
    removal_stats = defaultdict(int)  # refname -> total removed

    for (ref, start, end), remove_n in remove_counts.items():
        # If remove_n is None or 0, skip
        if not remove_n:
            continue

        reads = fetch_reads_in_region(all_reads, ref, start, end)
        total_reads = len(reads)

        if remove_n >= total_reads:
            # We remove them all
            for (r_id, r_ref, r_seq) in reads:
                if r_id not in removed_read_ids:
                    removed_read_ids.add(r_id)
                    removal_stats[ref] += 1
        else:
            # We remove exactly remove_n reads from this region
            if remove_mode == 'random':
                # Shuffle or sample
                sampled = random.sample(reads, remove_n)
            else:
                # By default, just remove the first remove_n
                sampled = reads[:remove_n]

            for (r_id, r_ref, r_seq) in sampled:
                if r_id not in removed_read_ids:
                    removed_read_ids.add(r_id)
                    removal_stats[ref] += 1

    return removed_read_ids, dict(removal_stats)

def report_coverage_stats(ground_truth, original_stats, filtered_stats, output_dir):
    """
    Report coverage stats compared to ground truth.
    Writes a CSV file with columns:
    reference, ground_truth_coverage, original_breadth, original_mean_depth, filtered_breadth, filtered_mean_depth
    """
    report_file = os.path.join(output_dir, "coverage_stats.csv")
    with open(report_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["reference", "ground_truth_coverage", "original_breadth", "original_mean_depth", "filtered_breadth", "filtered_mean_depth"])
        for ref in ground_truth:
            gt_cov = ground_truth[ref]/100
            orig_breadth = original_stats.get(ref, {}).get('breadth', 0)
            orig_mean_depth = original_stats.get(ref, {}).get('mean_depth', 0)
            filt_breadth = filtered_stats.get(ref, {}).get('breadth', 0)
            filt_mean_depth = filtered_stats.get(ref, {}).get('mean_depth', 0)
            mean_Depth_delta = orig_breadth - filt_breadth
            if orig_mean_depth != 0:
                perc_delta = (orig_breadth - filt_breadth)  * 100
            else:
                perc_delta = 0
            breadth_change = orig_breadth - filt_breadth
            writer.writerow([ref, gt_cov, orig_breadth, orig_mean_depth, filt_breadth, filt_mean_depth])
            print(f"{ref}:\n\tGT Cov: {gt_cov:.4f}\n\tOrig Breadth: {orig_breadth:.4f}\n\tOrig Depth: {orig_mean_depth:.4f}\n\tFilt Breadth: {filt_breadth:.4f}\n\tFilt Depth: {filt_mean_depth:.4f}\n\tBreadth Delta: {mean_Depth_delta:.4f}\n\t% Delta: {perc_delta:.2f}%")
    print(f"Coverage statistics written to {report_file}")


def pairwise_remove_reads(sum_comparisons, all_reads):
    """
    Go through each conflict pair and figure out how many reads we need to remove
    from each region in order to match (scale down to) the smaller coverage.

    We'll store the removal counts in `remove_counts`.
    If `remove_counts[(ref, start, end)]` is None => no constraint yet (default keep all).
    Once we set it to an integer, we interpret that as "remove at least this many reads".

    If a given region appears in multiple conflicts, we take the maximum removal
    so that it's guaranteed to remove enough reads to satisfy all conflicts.
    """
    remove_counts = defaultdict(lambda: None)  # None => default remove zero reads

    for ref, conflict_list in sum_comparisons.items():
        for conflict in conflict_list:
            jaccard = conflict['jaccard']
            if jaccard > 0:
                # Region 1
                region1 = (ref, conflict['s1'], conflict['e1'])
                # Region 2
                region2 = (conflict['to'], conflict['s2'], conflict['e2'])

                # Fetch reads (lists of (read_id, refname, seq), presumably)
                reads_r1 = fetch_reads_in_region(all_reads, region1[0], region1[1], region1[2])
                reads_r2 = fetch_reads_in_region(all_reads, region2[0], region2[1], region2[2])

                cov_r1 = len(reads_r1)
                cov_r2 = len(reads_r2)

                # We'll scale them both to min_cov
                min_cov = min(cov_r1, cov_r2)

                # fraction for each (to keep)
                frac_r1 = (min_cov / cov_r1) if cov_r1 else 1.0
                frac_r2 = (min_cov / cov_r2) if cov_r2 else 1.0

                # number to keep for each
                keep_r1 = max(1, int(frac_r1 * cov_r1 + 0.9999)) if cov_r1 else 0
                keep_r2 = max(1, int(frac_r2 * cov_r2 + 0.9999)) if cov_r2 else 0

                # number to remove for each
                remove_r1 = cov_r1 - keep_r1
                remove_r2 = cov_r2 - keep_r2

                # Combine constraints by taking the maximum removal so that
                # if another conflict says we need to remove more, we abide by that.
                if remove_counts[region1] is None:
                    remove_counts[region1] = remove_r1
                else:
                    remove_counts[region1] = min(remove_counts[region1], remove_r1)

                if remove_counts[region2] is None:
                    remove_counts[region2] = remove_r2
                else:
                    remove_counts[region2] = min(remove_counts[region2], remove_r2)

                # Debug prints
                print(f"Conflict between {region1} and {region2} (jaccard={jaccard:.4f})")
                print(f"  cov_r1={cov_r1}, keep={keep_r1}, remove={remove_r1}")
                print(f"  cov_r2={cov_r2}, keep={keep_r2}, remove={remove_r2}\n")
    return remove_counts

def build_conflict_groups(sum_comparisons, min_jaccard=0.0):
    """
    sum_comparisons is in the original format:
      {
        'NC_066642.1': [
          {
            'jaccard': 0.0137,
            'to': 'NC_003310.1',
            's1': 133813, 'e1': 134113,
            's2': 122124, 'e2': 123123
          },
          ...
        ],
        ...
      }

    min_jaccard: only treat as conflicts if jaccard >= min_jaccard.

    Returns: A list of conflict sets, each set is { (ref, start, end), (ref2, start2, end2), ... }
    """
    from collections import defaultdict, deque

    # Step 1: Build adjacency
    graph = defaultdict(set)

    for ref, conflict_list in sum_comparisons.items():
        for conflict in conflict_list:
            jacc = conflict['jaccard']
            if jacc >= min_jaccard:
                region1 = (ref, conflict['s1'], conflict['e1'])
                region2 = (conflict['to'], conflict['s2'], conflict['e2'])
                graph[region1].add(region2)
                graph[region2].add(region1)

    # Step 2: Find connected components
    visited = set()
    conflict_groups = []

    def bfs(start):
        queue = deque([start])
        component = set([start])
        visited.add(start)
        while queue:
            node = queue.popleft()
            for neighbor in graph[node]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    component.add(neighbor)
                    queue.append(neighbor)
        return component

    for node in graph:
        if node not in visited:
            comp = bfs(node)
            conflict_groups.append(comp)

    return conflict_groups



def remove_symmetric(conflict_group, all_reads, fetch_reads_in_region, remove_n=1):
    """
    conflict_group: a set of regions, e.g. { (r1, s1, e1), (r2, s2, e2), ... }
    remove_n: the number of reads to remove from each region in this group.

    Returns: dict { region -> set_of_removed_read_ids }

    If a region has fewer than remove_n reads, we'll remove them all (or remove as many as possible).
    """
    removed = {}
    for region in conflict_group:
        (ref, start, end) = region
        reads = fetch_reads_in_region(all_reads, ref, start, end)
        total = len(reads)
        if total == 0:
            removed[region] = set()  # nothing to remove
            continue

        # If total < remove_n, remove them all
        actual_remove = min(remove_n, total)
        # For simplicity, remove the first N reads (or random.sample if you prefer)
        to_remove = set(r_id for (r_id, _, _) in reads[:actual_remove])
        removed[region] = to_remove
    return removed

def finalize_proportional_removal(conflict_groups, all_reads, fetch_reads_in_region, remove_mode='random', random_seed=None):
    """
    For each conflict group, we:
      1) Gather all sub-regions belonging to each reference in that group.
      2) Sum coverage for each reference.
      3) Let 'min_cov' be the smallest coverage among references in that group.
      4) Remove exactly 'min_cov' reads from each reference's union of reads in that group
         (or remove them all if coverage < min_cov).

    remove_mode = 'random' or 'first'
    random_seed can be set for reproducible removals.

    Returns:
      removed_read_ids: a global set of read IDs that are removed from all conflict groups.

    Example:
      If conflict group has:
        - NC_003310.1 sub-regions total 5 reads
        - NC_006998.1 sub-regions total 70 reads
      => min_cov = 5
      => remove 5 from NC_003310.1 (all) and 5 from NC_006998.1.
    """
    if random_seed is not None:
        random.seed(random_seed)

    removed_read_ids = defaultdict(list)

    for group in conflict_groups:
        # 1) Group sub-regions by reference
        #    Example: ref_subregions[ref] = [(ref, s, e), (ref, s2, e2), ...]
        ref_subregions = defaultdict(list)
        for (ref, s, e) in group:
            ref_subregions[ref].append((s, e))

        # 2) For each reference, gather all reads from those sub-regions
        #    Then sum coverage => coverage_by_ref[ref] = total #reads
        coverage_by_ref = {}
        reads_by_ref = {}
        for ref, subregs in ref_subregions.items():
            # unify all reads for these sub-regions
            unified_reads = []
            for (start, end) in subregs:
                region_reads = fetch_reads_in_region(all_reads, ref, start, end)
                unified_reads.extend(region_reads)
            # remove duplicates if the same read appears in multiple sub-regions
            # using a dict or set keyed by read_id
            unique_ids = {}
            for (r_id, r_ref, r_seq) in unified_reads:
                if r_id not in unique_ids:
                    unique_ids[r_id] = (r_id, r_ref, r_seq)
            final_reads = list(unique_ids.values())

            coverage_by_ref[ref] = len(final_reads)
            reads_by_ref[ref] = final_reads

        if not coverage_by_ref:
            # no coverage in this group => skip
            continue

        # 3) min_cov = smallest coverage among references
        min_cov = min(coverage_by_ref.values())

        # 4) For each reference, remove min_cov reads (or all if coverage < min_cov)
        for ref, cov_count in coverage_by_ref.items():
            # reads for this ref in the group
            rlist = reads_by_ref[ref]
            if cov_count <= min_cov:
                # remove them all
                for (read_id, r_ref, r_seq) in rlist:
                    removed_read_ids[read_id].append(r_ref)
            else:
                # remove exactly min_cov reads
                if remove_mode == 'random':
                    sampled = random.sample(rlist, min_cov)
                else:
                    # remove the first min_cov
                    sampled = rlist[:min_cov]
                for (read_id, r_ref, r_seq) in sampled:
                    removed_read_ids[read_id].append(r_ref)

    return removed_read_ids

def import_bedfile(coverage_file):
    """
    Import coverage file in bedtools genomecov format.
    """
    with open(coverage_file, 'r') as f:
        return f.read()
def merge_bedgraph_regions(
    regions,
    statistic_func,
    max_stat_threshold,
    max_group_size=1000
):
    """
    Merges consecutive regions if the statistic (variance/Gini) is below the given threshold.

    :param regions: List of (chrom, start, end, coverage)
    :param statistic_func: A function that takes a list of coverage values and returns a float (e.g., variance or Gini)
    :param max_stat_threshold: The maximum allowed statistic to keep merging
    :param max_group_size: (Optional) max number of intervals to merge at once
                           before forcing a cutoff to avoid merging everything into one region.
    :return: A list of merged regions with format (chrom, merged_start, merged_end, average_coverage, <statistic>)
    """
    # Sort the input by (chrom, start)
    regions = sorted(regions, key=lambda x: (x[0], x[1]))

    merged = []

    # Accumulators for the current merge group
    current_chrom = None
    current_start = None
    current_end = None
    coverage_values = []
    region_count = 0

    def close_and_save_region():
        # Called when we need to finalize the current group.
        if not coverage_values:
            return

        avg_coverage = sum(coverage_values) / len(coverage_values)
        stat_value = statistic_func(coverage_values)
        merged.append((current_chrom, current_start, current_end, avg_coverage, stat_value))

    for i, (chrom, start, end, cov) in enumerate(regions):
        if current_chrom is None:
            # Initialize first region
            current_chrom = chrom
            current_start = start
            current_end = end
            coverage_values = [cov]
            region_count = 1
            continue

        # Check if new interval is contiguous and on same chrom
        # For "merging", we typically want (chrom matches) AND (start == current_end)
        # but you could be flexible if you want to merge even non-contiguous intervals.
        if (chrom == current_chrom) and (start == current_end):
            # Tentatively add to the current group
            extended_coverage_values = coverage_values + [cov]
            stat_value = statistic_func(extended_coverage_values)

            # Check if adding this interval exceeds the threshold
            if stat_value <= max_stat_threshold and region_count < max_group_size:
                # OK to merge
                coverage_values.append(cov)
                current_end = end
                region_count += 1
            else:
                # The statistic went over the threshold, so close the current group
                close_and_save_region()

                # Start a new group
                current_chrom = chrom
                current_start = start
                current_end = end
                coverage_values = [cov]
                region_count = 1
        else:
            # Different chromosome or non-contiguous => close previous region, start fresh
            close_and_save_region()
            current_chrom = chrom
            current_start = start
            current_end = end
            coverage_values = [cov]
            region_count = 1

    # Close the final region if needed
    close_and_save_region()

    return merged
import statistics

def compute_variance(values):
    return statistics.pvariance(values)  # or statistics.variance(values)


def compute_gini(values):
    """
    Computes the Gini index for a list of numeric values.
    Returns a value between 0 and 1 (0 = perfect equality).
    """
    if len(values) == 0:
        return 0

    sorted_values = sorted(values)
    n = len(values)
    cumulative = 0
    partial_sum = 0

    for i, val in enumerate(sorted_values, start=1):
        partial_sum += val
        cumulative += partial_sum

    gini = 1 - (2 * cumulative) / (n * partial_sum)
    return gini


def compute_f1_scores(all_reads):
    """
    Computes precision, recall, and F1 scores by comparing ground truth (gt) to predicted references.

    Parameters:
        all_reads (dict): A dictionary where each key is a reference and each value is a list of read dictionaries.
                          Each read dictionary must contain at least 'id' and 'gt' keys.

    Returns:
        per_ref_results (dict): A dictionary mapping each reference to its (precision, recall, F1) tuple.
        overall_scores (tuple): A tuple containing (micro_precision, micro_recall, micro_f1).
    """
    # Initialize dictionaries
    per_ref_predictions = defaultdict(set)
    per_ref_gt = defaultdict(set)
    per_ref_counts = defaultdict(lambda: defaultdict(int))

    # Populate per_ref_predictions and per_ref_gt
    for ref, reads in all_reads.items():
        for read in reads:
            read_id = read['id']
            gt_ref = read['gt']

            # Add read_id to the predicted reference
            per_ref_predictions[ref].add(read_id)

            # Add read_id to the ground truth reference
            per_ref_gt[gt_ref].add(read_id)

            # Optional: Populate per_ref_counts if needed
            per_ref_counts[gt_ref][ref] += 1

    # Initialize dictionary to store per-reference metrics
    per_ref_results = {}

    # Compute TP, FP, FN for each reference using set operations
    for ref in per_ref_predictions.keys():
        predicted = per_ref_predictions[ref]
        gt_set = per_ref_gt.get(ref, set())

        # True Positives: Reads correctly predicted to ref
        tp = len(predicted & gt_set)

        # False Positives: Reads predicted to ref but ground truth is different
        fp = len(predicted - gt_set)

        # False Negatives: Reads that are actually ref but not predicted as ref
        fn = len(gt_set - predicted)
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        f1 = (2 * precision * recall) / (precision + recall) if (precision + recall) > 0 else 0.0
        per_ref_results[ref] = {
            'precision': precision,
            'recall': recall,
            'f1_score': f1,
            'read_count': len(predicted),
            "FP": fp,
            "TP": tp,
            "FN": fn,
        }


    # Calculate micro-averaged metrics
    tp_total = sum(len(per_ref_predictions[ref] & per_ref_gt.get(ref, set())) for ref in per_ref_predictions)
    fp_total = sum(len(per_ref_predictions[ref] - per_ref_gt.get(ref, set())) for ref in per_ref_predictions)
    fn_total = sum(len(per_ref_gt.get(ref, set()) - per_ref_predictions[ref]) for ref in per_ref_predictions)

    micro_precision = tp_total / (tp_total + fp_total) if (tp_total + fp_total) > 0 else 0.0
    micro_recall = tp_total / (tp_total + fn_total) if (tp_total + fn_total) > 0 else 0.0
    micro_f1 = (2 * micro_precision * micro_recall) / (micro_precision + micro_recall) if (micro_precision + micro_recall) > 0 else 0.0

    overall_scores = {
        'micro_precision': micro_precision,
        'micro_recall': micro_recall,
        'micro_f1': micro_f1
    }

    return per_ref_results, overall_scores

def compare_metrics(original, new):
    """
    Compares original and new metrics and returns a pandas DataFrame.

    Parameters:
    - original (dict): Original metrics.
    - new (dict): New metrics.

    Returns:
    - comparison_df (DataFrame): DataFrame containing original, new, and delta metrics.
    """
    # Get all unique references
    all_refs = set(original.keys()).union(new.keys())

    # Prepare data for DataFrame
    data = []
    for ref in sorted(all_refs):
        orig = original.get(ref, {'TP': 0, 'FP': 0, 'FN': 0})
        nw = new.get(ref, {'TP': 0, 'FP': 0, 'FN': 0})
        delta_tp = nw['TP'] - orig['TP']
        delta_fp = nw['FP'] - orig['FP']
        delta_fn = nw['FN'] - orig['FN']
        perc_tp_change = (delta_tp / orig['TP']) * 100 if orig['TP'] > 0 else 0.0
        perc_fp_change = (delta_fp / orig['FP']) * 100 if orig['FP'] > 0 else 0.0
        total_perc_removed = (delta_fp + delta_fn + delta_tp) / (orig['TP'] + orig['FP'] + orig['FN']) * 100
        data.append({
            'Reference': ref,
            'TP Original': orig['TP'],
            'FP Original': orig['FP'],
            'FN Original': orig['FN'],
            'TP New': nw['TP'],
            'FP New': nw['FP'],
            'FN New': nw['FN'],
            'Δ TP': delta_tp,
            'Δ TP%': f"{perc_tp_change:.2f}",
            'Δ FP': delta_fp,
            'Δ FP%': f"{perc_fp_change:.2f}",
            'Δ All%': f"{total_perc_removed:.2f}",
            'Δ FN': delta_fn
        })

    # Create DataFrame
    import pandas as pd
    comparison_df = pd.DataFrame(data)
    # get the total number of FP original and new
    total_fp_original = comparison_df['FP Original'].sum()
    total_fp_new = comparison_df['FP New'].sum()
    print(f"Total FP Original: {total_fp_original}")
    print(f"Total FP New: {total_fp_new}")
    if total_fp_original > 0:
        print(f"Drop in FPs: {total_fp_original - total_fp_new} ({((total_fp_original - total_fp_new) / total_fp_original) * 100:.2f}%)")
    return comparison_df
def rebuild_sig_dict(loaded_sigs):
    """
    loaded_sigs: list of (SourmashSignature, name_string)
    Returns a dict of (ref, start, end) -> signature
    """
    sig_dict = {}
    for sig in loaded_sigs:
        # e.g. name_str = 'NC_003310.1:1000-2000'
        name = sig.name
        sig_dict[name] = sig
    return sig_dict

def save_signatures_sourmash(signatures, output_sigfile):
    """
    signatures: dict of (ref, start, end) -> SourmashSignature
    Saves all signatures into one `.sig` file.
    """
    os.makedirs(os.path.dirname(output_sigfile), exist_ok=True)
    # Extract just the signature objects
    siglist = []
    for ref, sig in signatures.items():
        # Optionally embed region info in the signature's name or metadata
        # so you can retrieve it later
        sig._name = ref
        # or set custom metadata: sig._md['region_info'] = (ref, start, end)
        siglist.append(sig)

    with open(output_sigfile, "wt") as fp:
        save_signatures(siglist, fp)
    fp.close()
    print(f"Signatures saved to {output_sigfile}")

def load_signatures_sourmash(input_sigfile):
    """
    Loads multiple signatures from a single `.sig` file.
    Returns a list of (signature, original_name).
    If you stored region info in metadata, you can reconstruct your dict.
    """
    loaded_sigs = load_file_as_signatures(input_sigfile)
    return loaded_sigs
def parse_sourmash_compare_csv(csv_file):
    """
    Parses the CSV output of `sourmash compare --csv`.
    Returns:
      references: list of reference names in row/column order
      dist_dict: dict of dicts, dist_dict[ref1][ref2] = similarity (float)
    """
    with open(csv_file, mode='r', newline='', encoding='utf-8') as file:
        # Use csv.reader to parse the file
        parsed_data = list(csv.reader(file, quotechar='"', delimiter=','))
        # assigne references to header
        references = parsed_data[0][0:]
        references = [ref.split(" ")[0] for ref in references]
        datalines = []
        if len(parsed_data ) > 1:
            datalines = parsed_data[1:]
        matrix = []
        for line in datalines:
            # Print the result
            matrix.append([float(x) for x in line])
    # Build dist_dict
    dist_dict = {}
    for i, ref1 in enumerate(references):
        dist_dict[ref1] = {}
        for j, ref2 in enumerate(references):
            dist_dict[ref1][ref2] = matrix[i][j]
    return references, dist_dict


def find_reference_clusters(dist_dict, min_sim=0.7):
    """
    Returns a list of clusters (each cluster is a set of reference names).
    A cluster = references that are connected by >= min_sim edges.

    Example:
      If A & C are >= 0.7, B & D are >= 0.8, E is alone =>
      we get [{A, C}, {B, D}, {E}]
    """
    references = list(dist_dict.keys())

    # 1) Build adjacency list: adjacency[ref] = set of neighbors >= min_sim
    adjacency = {}
    for i, ref1 in enumerate(references):
        adjacency[ref1] = set()
        for j, ref2 in enumerate(references):
            if ref1 == ref2:
                continue
            if dist_dict[ref1][ref2] >= min_sim:
                adjacency[ref1].add(ref2)

    visited = set()
    clusters = []

    # 2) BFS or DFS to find connected components
    def bfs(start):
        queue = [start]
        visited.add(start)
        cluster = {start}
        while queue:
            current = queue.pop(0)
            for neighbor in adjacency[current]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append(neighbor)
                    cluster.add(neighbor)
        return cluster

    for ref in references:
        if ref not in visited:
            cluster = bfs(ref)
            clusters.append(cluster)

    return clusters
def filter_singleton_clusters(clusters):
    """
    Returns only the clusters that have size >= 2.
    """
    return [c for c in clusters if len(c) >= 2]
def filter_regions_for_cluster(merged_regions, cluster):
    """
    Keep only regions whose ref is in the given cluster set.
    """
    return [(ref, s, e) for (ref, s, e) in merged_regions if ref in cluster]

def create_signatures_for_cluster(merged_regions, cluster, all_reads, bam_path, kmer_size, scaled):
    """
    1) Filter regions to the cluster's references.
    2) Create signatures for those filtered regions.
    Returns a dict of (ref, start, end) -> SourmashSignature.
    """
    cluster_regions = filter_regions_for_cluster(merged_regions, cluster)
    if not cluster_regions:
        return {}
    # then call your existing function that does minhashing
    signatures = create_signatures_for_regions(cluster_regions, all_reads, bam_path, kmer_size, scaled)
    return signatures


def import_matrix_ani(matrix_file):
    import pandas as pd
    # Load the ANI table
    ani_df = pd.read_csv(matrix_file, sep='\t', header=None, names=['Genome1', 'Genome2', 'ANI', 'AF', 'BF'])

    ani_df['Genome1'] = ani_df['Genome1'].apply(lambda x: os.path.basename(x).replace('.fasta', ''))
    ani_df['Genome2'] = ani_df['Genome2'].apply(lambda x: os.path.basename(x).replace('.fasta', ''))

    # Pivot the table to create a matrix
    ani_matrix = ani_df.pivot(index='Genome1', columns='Genome2', values='ANI')
    # cahnge genome1 and genome2 to be basename and remove .fasta ext using os path formatting


    # Since ANI is symmetric, mirror the upper triangle to the lower triangle
    ani_matrix = ani_matrix.combine_first(ani_matrix.T)
    # set values of cells to 0-1 rather than 0-100
    ani_matrix = ani_matrix / 100
    # Fill diagonal with 100% ANI (if not already)
    for genome in ani_matrix.columns:
        ani_matrix.at[genome, genome] = 1
    # convert to dictionary
    ani_matrix = ani_matrix.to_dict()
    return ani_matrix



def determine_conflicts(
        output_dir = None,
        input_bam = None,
        matrix = None,
        min_threshold = 0.2,
        abu_file = None,
        min_similarity_comparable = 0.0,
        use_variance = False,
        apply_ground_truth = False,
        use_reads_gt = False,
        sigfile = None,
        bedfile = None,
        reference_signatures = None,
        scaled = 100,
        kmer_size=31,
):
    # parser = argparse.ArgumentParser(description="Use bedtools to define coverage-based regions and compare read signatures using sourmash MinHash.")
    # parser.add_argument("--input_bam", required=True, help="Input BAM file.")
    # parser.add_argument("--matrix", required=False, help="A Matrix file for ANI in long format from fastANI")
    # parser.add_argument("--output_dir", required=True, help="Output directory for results.")
    # parser.add_argument("--kmer_size", type=int, default=51, help="k-mer size for MinHash.")
    # parser.add_argument("--scaled", type=int, default=2000, help="scaled factor for MinHash.")
    # parser.add_argument("--bedtools_path", default='bedtools', help="Path to bedtools executable.")
    # parser.add_argument("--bedfile", default=None, help="Coverage file from bedtools")
    # parser.add_argument("--min_threshold", type=float, default=0.2, help="Min Jaccard similarity to report.")
    # parser.add_argument("--abu_file", required=False, help="Path to ground truth coverage file (abu.txt).")
    # parser.add_argument("--use_variance", required=False, action='store_true', help="Use variance instead of Gini index.")
    # parser.add_argument("--apply_ground_truth", required=False, action='store_true', help="Overwrites abu file for F1 score metrics. Parses the read name for the actual reference. Assigns --use_reads_gt as true")
    # parser.add_argument("--use_reads_gt", required=False, action='store_true', help="Calculates the f1 score based on the reads rather than coverage i.e. abu.txt file from --abu_file")
    # parser.add_argument("--sigfile", required=False, type=str, help="Skip signatures comparison if provided ad load it instead")
    # parser.add_argument(
    #     "--reference_signatures",
    #     help="Path to the CSV file generated by `sourmash compare --csv ...`",
    #     required=False
    # )
    # parser.add_argument(
    #     "--min_similarity_comparable",
    #     type=float,
    #     default=0.0,
    #     help="Minimum ANI threshold to consider references comparable (default=0.7)"
    # )
    # args = parser.parse_args()

    os.makedirs(output_dir, exist_ok=True)
    import time
    start = time.time()
    basename = os.path.basename(input_bam).replace(".bam", "")

    # Step 1: Compute coverage using bedtools genomecov, if not present then do it with depthfile
    if not bedfile:
        bed_path = os.path.join(output_dir, f"{basename}.bed")
        bed_path = None
        if not bed_path:
            coverage_output = compute_coverage_with_bedtools(
                input_bam,
                None,
            )
    else:
        coverage_output = import_bedfile(bedfile)
    end = time.time()
    print(f"Time to complete import of bed information: {end-start} seconds")
    # Step 8: Load all reads from BAM once
    print("Loading reads from BAM...")
    all_reads = load_reads_from_bam(input_bam)
    print(f"Total references with reads: {len(all_reads)}")
    gt_performances = dict()
    # 3. If user provided a reference comparison CSV


    if matrix:
        dist_dict = import_matrix_ani(matrix)
        clusters = find_reference_clusters(dist_dict, min_sim=min_similarity_comparable)
    elif reference_signatures:  # the CSV from `sourmash compare --csv`
        _, dist_dict = parse_sourmash_compare_csv(reference_signatures)
        # 3A. Build clusters (connected components)
        clusters = find_reference_clusters(dist_dict, min_sim=min_similarity_comparable)
        # 3B. Filter out singletons
        # clusters = filter_singleton_clusters(clusters)
        # print(f"Found {len(clusters)} multi-ref clusters above similarity {min_similarity_comparable}")
        # for c in clusters:
        #     print("  Cluster:", c)
    else:
        # If no CSV is provided, each reference stands alone,
        # or you can treat everything as a single big cluster
        # (depending on your logic).
        all_references = set(all_reads.keys())
        clusters = [all_references]

    # Step 7: Parse filtered BED to get regions
    regions = parse_bed_file(bedfile)
    print(f"Total regions defined: {len(regions)}")
    # 2. Choose statistic function
    if use_variance:
        from statistics import pvariance
        statistic_func = lambda vals: pvariance(vals)  # population variance
        stat_name = "variance"
    else:
        statistic_func = compute_gini
        stat_name = "gini"
    # 3. Merge
    threshold=0.9
    merged_regions = merge_bedgraph_regions(
        regions,
        statistic_func=statistic_func,
        max_stat_threshold=threshold,
        max_group_size=400  # Example: limit merges to 200 intervals
    )

    # 4. Print or save results
    print(f"Merged regions (using {stat_name} <= {threshold}):")
    print(f"Length of original regions : {len(regions)})")
    print(f"Length of merged regions: {len(merged_regions)}")

    if not sigfile:
        # Step 9: Create signatures for each region
        print("Creating signatures for regions...")
        signatures = create_signatures_for_regions(merged_regions, all_reads, input_bam, kmer_size, scaled)
        print(f"Total signatures created: {len(signatures)}")
        # Step 10: Save signatures to files
        sig_dir = os.path.join(output_dir, "signatures")
        single_sigfile = os.path.join(sig_dir, "merged_regions.sig")
        save_signatures_sourmash(signatures, single_sigfile)
    else:
        # read in the signatures
        loaded_sigs = load_signatures_sourmash(sigfile)
        signatures = rebuild_sig_dict(loaded_sigs)

    output_csv = os.path.join(output_dir, "region_comparisons.csv")
    print(f"Comparison results written to: {output_csv}")
    print("Comparing signatures pairwise...")
    print(f"Clusters to compare", )
    for i, c in enumerate(clusters):
        print(f"C{i}:", c)
    passed_regions, sum_comparisons = compare_signatures(
        signatures,
        output_csv,
        min_threshold,
        clusters
    )
    # 1) Build conflict groups from sum_comparisons
    conflict_groups = build_conflict_groups(sum_comparisons, min_jaccard=0.0)
    # for g in conflict_groups:
    #     print(f"Conflict group: {g}")

    # 2) For each group, remove 1 read from each region
    removed_read_ids = finalize_proportional_removal(
        conflict_groups,
        all_reads,
        fetch_reads_in_region,
        remove_mode='random'
    )




    # 3) Write final BAM
    filtered_bam = os.path.join(output_dir, "filtered.hash.bam")
    create_filtered_bam(
        input_bam,
        filtered_bam,
        removed_read_ids.keys()
    )


   # Create filtered BAM
    if use_reads_gt:
        for ref, reads in all_reads.items():
            for i, read in enumerate(reads):
                gt_read = read['id'].split("_")
                gt_name = "_".join(gt_read[:-2])
                all_reads[ref][i]['gt'] = gt_name
        print("Calculating f1 scores of alignment first....")
        per_ref, overall = compute_f1_scores(all_reads)
        # Print per-reference scores
         # Display per-reference results
        # Display micro-averaged results
        gt_performances['original'] = dict(overall=overall, per_ref = per_ref)
        all_reads = load_reads_from_bam(filtered_bam)
        for ref, reads in all_reads.items():
            for i, read in enumerate(reads):
                gt_read = read['id'].split("_")
                gt_name = "_".join(gt_read[:-2])
                all_reads[ref][i]['gt'] = gt_name
        print("Calculating f1 scores of alignment first....")
        per_ref, overall = compute_f1_scores(all_reads)
        gt_performances['new'] = dict(overall=overall, per_ref = per_ref)
        # Compare the metrics
        comparison_df = compare_metrics(gt_performances['original']['per_ref'], per_ref)

        # Display the comparison
        print("\n=== Metrics Comparison ===\n")
        print(comparison_df.to_string(index=False))
        print("Original Metrics Overall: ")
        print(f"\tPrecision: {gt_performances['original']['overall']['micro_precision']:.4f}")
        print(f"\tRecall: {gt_performances['original']['overall']['micro_recall']:.4f}")
        print(f"\tF1: {gt_performances['original']['overall']['micro_f1']:.4f}")
        print("New Filtered Metrics Overall: ")
        print(f"\tPrecision: {gt_performances['new']['overall']['micro_precision']:.4f}")
        print(f"\tRecall: {gt_performances['new']['overall']['micro_recall']:.4f}")
        print(f"\tF1: {gt_performances['new']['overall']['micro_f1']:.4f}")
        # print df to output_dir/removal_stats.xlsx
        comparison_df.to_excel(os.path.join(output_dir, "removal_stats.xlsx"), index=False)
    if abu_file:
        # Parse ground truth coverage
        ground_truth = parse_ground_truth(abu_file)
        # Step 14: Compute coverage stats before and after filtering
        print("Computing coverage statistics before filtering...")
        original_coverage_stats = compute_breadth_and_depth(input_bam)
        print("Computing coverage statistics after filtering...")
        filtered_coverage_stats = compute_breadth_and_depth(filtered_bam)

        # Step 15: Report coverage stats compared to ground truth
        print("Reporting coverage statistics...")
        report_coverage_stats(ground_truth, original_coverage_stats, filtered_coverage_stats, output_dir)

        # Step 16: Compute performance metrics
        print("Computing performance metrics based on ground truth and filtered BAM...")
        TP, FP, FN, precision, recall, f1 = compute_performance(ground_truth, filtered_bam)
        print("Performance Metrics:")
        print(f"TP: {TP}, FP: {FP}, FN: {FN}")
        print(f"Precision: {precision:.4f}")
        print(f"Recall: {recall:.4f}")
        print(f"F1-Score: {f1:.4f}")
    return comparison_df
