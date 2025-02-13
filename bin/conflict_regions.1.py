#!/usr/bin/env python3
import os
import csv
import pysam
import time
from sourmash.sbt import SBT, GraphFactory
from sourmash.sbtmh import SigLeaf
from sourmash import MinHash, SourmashSignature, save_signatures, load_file_as_signatures
from itertools import groupby
import statistics

import subprocess

from tqdm import tqdm  # For progress bar
from collections import defaultdict

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
    reflengths = {}
    bam = pysam.AlignmentFile(bam_path, "rb")
    all_reads = {}
    # get reflengths all references
    for ref in bam.references:
        reflengths[ref] = bam.get_reference_length(ref)
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
    return all_reads, reflengths

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


def create_signatures_for_regions(regions, all_reads, bam_path, kmer_size, scaled, num_workers = 1):
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
    i=0
    for (refname, start, end, val) in tqdm(regions):
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
        # i+=1
        # if i % 10_000 == 0:
        #     print(f"Processed {i} regions...")
    return signatures

#######################################################
# Helper functions
#######################################################
def parse_split(region_name):
    """
    Example: region_name = 'NC_003310.1:1000-2000'
    Returns (ref, start, end).
    Adjust to match your actual naming scheme.
    """
    ref_region = region_name.split(':')
    refname = ref_region[0]
    coords = ref_region[1].split('-')
    start = int(coords[0])
    end = int(coords[1])
    return (refname, start, end)

def build_sbt_index(siglist, ksize=31, sbt_name="my_sbt", clusters={}):
    """
    Build an SBT index in memory (or on disk).
    - siglist: list of (region_name, sourmash_signature)
    - ksize, scaled: specify how your MinHashes were created.

    Returns an SBT object.
    """
    # Create an empty SBT
    sbt_index = dict()
    factory = GraphFactory(ksize=ksize,  n_tables=1, starting_size=len(siglist))
    ref2cluster = {}
    if clusters is not None:
        for cluster_idx, cluster_set in enumerate(clusters):
            for ref in cluster_set:
                ref2cluster[ref] = cluster_idx
    # Add each signature as a leaf

    for i, (region_name, sig) in enumerate(siglist):
        r1, s1, e1 = parse_split(region_name)
        # sig.name = region_name
        leaf = SigLeaf(str(region_name), sig)
        cluster = ref2cluster.get(r1, "Unknown")
        if cluster in sbt_index:
            sbt = sbt_index[cluster]
        else:
            sbt = SBT(factory)
            sbt_index[cluster] = sbt
        sbt.add_node(leaf)


    # If you want to save it to disk, do:
    # sbt.save(sbt_name)
    return sbt_index

def fast_mode_sbt(
    siglist,
    sbt_index,
    output_csv,
    min_threshold=0.05,
    cluster_map=None
):
    """
    Use SBT index to find similar signatures above min_threshold (Jaccard).
    cluster_map: dictionary { refname -> cluster_id } (optional).
    """
    sum_comparisons = {}
    for (region_name, _) in siglist:
        ref, s1, e1 = parse_split(region_name)
        sum_comparisons.setdefault(ref, [])
    print(f"Building a cluster map to each reference input, total indices types: {len(sbt_index)}")
    ref2cluster = {}
    if cluster_map is not None:
        for cluster_idx, cluster_set in enumerate(cluster_map):
            for ref in cluster_set:
                ref2cluster[ref] = cluster_idx
    start_time = time.time()
    seen_comparisons = defaultdict()
    with open(output_csv, "w", newline="") as outfh:
        writer = csv.writer(outfh)
        writer.writerow([
            "reference1", "start1", "end1",
            "reference2", "start2", "end2",
            "jaccard", "containment_1_in_2", "containment_2_in_1"
        ])

        # Query each signature in the SBT
        for i, (region1, sig1) in enumerate(siglist):
            r1, s1, e1 = parse_split(region1)
            if i % 1000 == 0 and i > 0:  # Optional: progress update
                print(f"Processed {i} signatures in {time.time()-start_time} seconds...")
            if region1 in seen_comparisons:
                # print(f"{region1}. repeat comparison, skipping...")
                continue

            # SBT search returns a generator of SearchResult objects
            # at or above min_threshold (Jaccard).
            sbt = sbt_index.get(ref2cluster.get(r1, "Unknown"))
            if sbt is None:
                continue
            results = sbt.search(
                sig1, threshold=min_threshold, best_only=True
            )

            # We'll iterate over these matches
            # if results:
            #     print(region1)
            for sr in results:
                # sr is a SearchResult with .similarity, .location, etc.
                jaccard = sr.score  #  jaccard similarity score is here
                desc = sr.signature.name    # we stored region_name in the leaf constructor
                mh2 = sr.signature.minhash
                r2, s2, e2 = parse_split(desc)
                if r2 == r1:
                    continue  # skip self match to the same accession

                # print(f"\tFound {desc} with jaccard {jaccard:.4f}")
                if cluster_map is not None:
                    # If either reference isn't in the map, or they differ in cluster, skip
                    if (r1 not in ref2cluster or
                        r2 not in ref2cluster or
                        ref2cluster[r1] != ref2cluster[r2]):
                        continue
                ani = mh2.containment_ani(sig1.minhash)
                # print(f"{ani.dist:.2f}", desc, region1, f"{jaccard:.2f}")
                # # If you need actual containments, compute them:
                # mh1 = sig1.minhash
                # # We'll load the matching leaf's signature from the SBT
                # # or you can keep a dictionary of region_name -> signature.
                # # But let's suppose we can find it from siglist in memory:
                # # (In an actual pipeline, store a dict for quick lookup.)
                # sig2 = None
                # # for (rname2, s2sig) in siglist:
                # #     if rname2 == region2:
                # #         sig2 = s2sig
                # #         break

                # if not sig2:
                #     continue  # couldn't find the matching sig in memory
                # mh2 = sig2.minhash

                # c1_in_2 = mh1.avg_containment(mh2)
                # c2_in_1 = mh2.avg_containment(mh1)
                c1_in_2 = c2_in_1 = jaccard

                seen_comparisons[desc] = True
                seen_comparisons[region1] = True

                # Write to CSV
                writer.writerow([r1, s1, e1, r2, s2, e2, jaccard, c1_in_2, c2_in_1])

                # Optional: store in sum_comparisons
                sum_comparisons[r1].append(dict(
                    jaccard=jaccard, to=r2, s2=s2, e2=e2, s1=s1, e1=e1
                ))
                sum_comparisons[r2].append(dict(
                    jaccard=jaccard, to=r1, s2=s1, e2=e1, s1=s2, e1=e2
                ))
    print(f"Processed all {i} signatures in {time.time()- start_time} seconds...")
    return sum_comparisons

def slow_mode_linear(siglist, output_csv, min_threshold=0.1, cluster_map=None):
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
    print("Building a cluster map to each reference input")
    ref2cluster = {}
    if cluster_map is not None:
        for cluster_idx, cluster_set in enumerate(cluster_map):
            for ref in cluster_set:
                ref2cluster[ref] = cluster_idx
    seen_comparisons = defaultdict()

    # ------------------------------------------------
    # 2) Prepare data structures
    # ------------------------------------------------
    print("\tPrepping the data structures")

    # Initialize passed_regions with all region keys.
    all_regions = set([region for (region, _) in siglist])
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
    print("\tWriting comparison to outfile")
    start_time = time.time()
    with open(output_csv, "w", newline="") as outfh:
        writer = csv.writer(outfh)
        writer.writerow([
            "reference1","start1","end1",
            "reference2","start2","end2",
            "jaccard","containment_1_in_2","containment_2_in_1"
        ])

        for i in range(len(siglist)):
            region1, sig1 = siglist[i]
            if i > 0 and i % 1_000 == 0:
                print(f"Processed {i} signatures in {time.time()- start_time} seconds...")
            # if region1 in seen_comparisons:
                # print(f"{region1}. repeat comparison, skipping...")
                # continue

            r1, s1, e1 = parse_split(region1)
            mh1 = sig1.minhash
            for j in range(i+1, len(siglist)):
                region2, sig2 = siglist[j]
                r2, s2, e2 = parse_split(region2)

                # 3A) Skip if same reference (as in your original code)
                if r2 == r1:
                    continue

                # 3B) Skip if not in the same cluster
                if cluster_map is not None:
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
                    seen_comparisons[region1] = True
                    seen_comparisons[region2] = True
                    # # Write to CSV
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
    print(f"Processed all {i} signatures in {time.time()- start_time} seconds...")
    # ------------------------------------------------
    # 4) Print aggregated comparisons
    # ------------------------------------------------
    # for ref, jaccards in sum_comparisons.items():
    #     if ref == "NC_003310.1":
    #         print(f"{ref}:")

    #         for entry in jaccards:
    #             if entry.get('to') == "NC_006998.1":
    #                 print(
    #                     f"\t{entry['s1']}-{entry['e1']}\t"
    #                     f"Jaccard: {entry['jaccard']:.4f}\n"
    #                     f"\t\tto {entry['to']}:{entry['s2']}-{entry['e2']}"
    #                 )

    return sum_comparisons

def create_filtered_bam(original_bam_path, output_bam_path, alignments_to_keep):
    """
    Create a new BAM file containing only reads with IDs in read_ids_to_keep.
    """
    bam_in = pysam.AlignmentFile(original_bam_path, "rb")
    bam_out = pysam.AlignmentFile(output_bam_path, "wb", template=bam_in)
    set_reads_refs = set()
    count_total = 0
    count_written = 0
    for read in bam_in:
        count_total += 1
        if (read.reference_name, read.query_name) in alignments_to_keep:
            bam_out.write(read)
            count_written += 1
        set_reads_refs.add(read.query_name)
    bam_in.close()
    bam_out.close()
    # subprocess.run(["samtools", "index", output_bam_path])
    pysam.index(output_bam_path)
    print(f"Total reads processed: {count_total}")
    print(f"Reads written to filtered BAM: {count_written}")
    print(f"Unique Query Read IDS: {len(set_reads_refs)}")

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
        merged.append((current_chrom, current_start, current_end, stat_value))

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

def merge_bedgraph_regions(
    intervals,
    statistic_func,
    max_stat_threshold,
    max_group_size=20000,
    max_length=None,
    value_diff_tolerance=None
):
    """
    intervals: list of tuples (chrom, start, end, value)
               e.g. [('NC_003310.1', 6085, 11112, 86.58), ...]
    Returns a list of merged intervals, each of which is (chrom, start, end, stat_value, <optional>).
    """

    if not intervals:
        return []

    # Sort intervals by (chrom, start) to ensure a sensible merge order
    # (If you already have them sorted by contig, then by start, skip.)
    intervals = sorted(intervals, key=lambda x: (x[0], x[1]))

    merged = []
    current_chrom = None
    current_buffer = []  # Temporary storage for the intervals to be merged

    for interval in intervals:
        chrom, start, end, depth_val = interval

        # If we are on a new chromosome, flush the old buffer
        if chrom != current_chrom:
            # If there's anything in the current_buffer, finalize it
            if current_buffer:
                merged.extend(_finalize_buffer(current_buffer, statistic_func))
            current_chrom = chrom
            current_buffer = [interval]
            continue

        # Otherwise, we are on the same chromosome.
        # Check if we can merge this interval into the current_buffer
        if should_merge(
            current_buffer,
            interval,
            statistic_func,
            max_stat_threshold,
            max_group_size,
            max_length,
            value_diff_tolerance
        ):
            # Merge by appending
            current_buffer.append(interval)
        else:
            # Finalize the current_buffer (turn it into one merged region) and start anew
            merged.extend(_finalize_buffer(current_buffer, statistic_func))
            current_buffer = [interval]

    # End of loop – finalize whatever is left
    if current_buffer:
        merged.extend(_finalize_buffer(current_buffer, statistic_func))

    return merged


def _finalize_buffer(buffer_intervals, statistic_func):
    """
    Takes all intervals currently in buffer_intervals and merges them into one region.
    Returns a list with a single “merged” entry, or you could break it into multiple if you prefer.
    """
    if not buffer_intervals:
        return []

    chrom = buffer_intervals[0][0]
    start = min(iv[1] for iv in buffer_intervals)
    end   = max(iv[2] for iv in buffer_intervals)
    values = [iv[3] for iv in buffer_intervals]
    stat_val = statistic_func(values)

    # Return a single merged interval (you can store additional info if desired)
    return [(chrom, start, end, stat_val)]


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

def compare_metrics(original, new, breadth_old, breadth_new):
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
        total_p =  (orig['TP'] + orig['FP'] + orig['FN'])
        total_perc_removed =  ( (delta_fp + delta_fn + delta_tp) / total_p * 100 ) if total_p > 0 else 0.0
        data.append({
            'Reference': ref,
            'TP Original': orig['TP'],
            'FP Original': orig['FP'],
            'FN Original': orig['FN'],
            'Breadth Original': breadth_old.get(ref, 0),
            'Breadth New': breadth_new.get(ref, 0),
            'TP New': nw['TP'],
            'FP New': nw['FP'],
            'FN New': nw['FN'],
            # 'Δ TP': delta_tp,
            'Δ TP%': f"{perc_tp_change:.2f}",
            # 'Δ FP': delta_fp,
            'Δ FP%': f"{perc_fp_change:.2f}",
            'Δ All%': f"{total_perc_removed:.2f}",
            'Δ Breadth': (breadth_new.get(ref, 0) - breadth_old.get(ref, 0)),
            'Δ Breadth%': (breadth_new.get(ref, 0)  / breadth_old.get(ref, 0)) if breadth_old.get(ref, 0) > 0 else 0,
            # 'Δ Regions %': f"{100*regions_stat.get('regions_remain', 0):.2f}",
            # 'Sum Rs %': f"{ ( ( 100+total_perc_removed ) + ( 100*regions_stat.get('regions_remain', 0) ) ) /2  :.2f}",
            # 'Δ FN': delta_fn
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

def create_breadth_coverage_pysam(bamfile=None, covfile_output=None, coverage_threshold=1):
    """
    Generate breadth of coverage statistics using pysam.

    Parameters:
    - bamfile (str): Path to the input BAM file.
    - covfile_output (str, optional): Path to the output coverage file.
    - coverage_threshold (int, optional): Minimum coverage depth to consider a base as covered. Default is 1.

    Returns:
    - None: If coverage stats are printed to stdout.
    - dict: Dictionary of coverage stats {reference: breadth_of_coverage_percentage} if the output file is empty or the command fails.
    """
    try:
        if not bamfile:
            raise ValueError("Parameter 'bamfile' is required.")

        if not os.path.exists(bamfile):
            raise FileNotFoundError(f"BAM file '{bamfile}' does not exist.")

        # Open the BAM file
        with pysam.AlignmentFile(bamfile, "rb") as bam:
            references = bam.references
            lengths = bam.lengths

            breadth_coverage = {}

            for ref, length in zip(references, lengths):
                # Initialize covered bases count
                covered_bases = 0

                # Using count_coverage to get coverage per base
                # Returns a tuple of four lists (A, C, G, T) coverage
                coverage = bam.count_coverage(ref, quality_threshold=0)

                # Sum coverage across all four bases to get total coverage per base
                total_coverage = [sum(base) for base in zip(*coverage)]

                # Calculate the number of bases with coverage >= threshold
                for cov in total_coverage:
                    if cov >= coverage_threshold:
                        covered_bases += 1

                # Calculate breadth of coverage percentage
                breadth = (covered_bases / length) * 100 if length > 0 else 0
                breadth_coverage[ref] = breadth

        if covfile_output:
            # Write coverage statistics to the specified file
            with open(covfile_output, 'w') as outfile:
                outfile.write("Reference\tBreadth_of_Coverage(%)\n")
                for ref, breadth in breadth_coverage.items():
                    outfile.write(f"{ref}\t{breadth:.2f}\n")

            # Check if the output file is non-empty
            if os.path.getsize(covfile_output) > 0:
                # Print the coverage statistics to stdout
                with open(covfile_output, 'r') as infile:
                    coverage_data = infile.read()
                    print(coverage_data)
                return
            else:
                # Output file is empty; return the coverage dictionary
                print("Coverage output file is empty. Returning coverage data as a dictionary.")
                return breadth_coverage
        else:
            # No output file specified; return the coverage dictionary
            return breadth_coverage

    except FileNotFoundError as fnf_error:
        print(f"Error: {fnf_error}")
        raise
    except ValueError as ve:
        print(f"Error: {ve}")
        raise
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        raise

def calculate_breadth_of_coverage_dict(all_reads_dict, reference_lengths=None, skip_ids={}):
    """
    Calculate the breadth of coverage for each reference in a dictionary of reads.

    Parameters:
    - all_reads_dict (dict): Keys are reference IDs (str), and values are lists of read dictionaries.
                            Each read dictionary must have 'start' and 'end' keys.
    - reference_lengths (dict, optional): Keys are reference IDs (str), and values are the total length (int) of the reference.
                                        If not provided, the reference length is inferred from the reads.

    Returns:
    - dict: Keys are reference IDs, and values are breadth of coverage percentages (float).
    """
    def calculate_breadth_of_coverage(reads, ref_length=None):
        """
        Helper function to calculate coverage for a single reference.

        Parameters:
        - reads (list of dict): Each dict has 'start' and 'end' keys.
        - ref_length (int, optional): Total length of the reference. If not provided, inferred from reads.

        Returns:
        - float: Breadth of coverage as a percentage.
        """
        if not reads:
            return 0.0  # No reads provided

        # Extract intervals
        intervals = [(read['start'], read['end']) for read in reads]

        # Sort intervals by start position
        sorted_intervals = sorted(intervals, key=lambda x: x[0])

        # Merge overlapping or adjacent intervals
        merged_intervals = []
        current_start, current_end = sorted_intervals[0]

        for start, end in sorted_intervals[1:]:
            if start <= current_end + 1:  # Overlapping or adjacent
                current_end = max(current_end, end)
            else:
                merged_intervals.append((current_start, current_end))
                current_start, current_end = start, end
        merged_intervals.append((current_start, current_end))  # Add the last interval

        # Calculate total covered bases
        total_covered = sum(end  - start   for start, end in merged_intervals)

        # Define reference length
        if ref_length is not None:
            reference_length = ref_length
            # Optionally, you might want to adjust merged_intervals to fit within reference_length
            # Here, we assume reads are already within the reference
        else:
            min_start = min(start for start, end in intervals)
            max_end = max(end for start, end in intervals)
            reference_length = max_end - min_start + 1

        # Calculate breadth of coverage
        breadth = (total_covered / reference_length) * 100  # Percentage

        return breadth

    coverage_dict = {}
    for ref_id, reads in all_reads_dict.items():
        ref_length = reference_lengths.get(ref_id) if reference_lengths else None
        # filter reads where != skip_ids
        reads_to_skip = []
        for r in reads:
            if r.get('id') in skip_ids and ref_id in skip_ids.get(r.get('id'), []):
                reads_to_skip.append(r.get('id'))
        reads = [r for r in reads if r.get('id') not in reads_to_skip]
        coverage = calculate_breadth_of_coverage(reads, ref_length)
        coverage_dict[ref_id] = coverage

    return coverage_dict

def calculate_breadth_of_coverage_dict_regions(all_reads_dict, reference_lengths=None, removed_read_ids=None, regions_dict=None):
    """
    Calculate the breadth of coverage for each reference over specified regions,
    excluding specified reads per reference.

    Parameters:
    - all_reads_dict (dict): Keys are reference IDs (str), and values are lists of read dictionaries.
    - reference_lengths (dict, optional): Keys are reference IDs (str), and values are the total length (int) of the reference.
    - removed_read_ids (defaultdict(list), optional): Keys are read IDs (str), and values are lists of reference IDs (str) to exclude.
    - regions_dict (dict, optional): Keys are reference IDs (str), and values are lists of [start, end] regions.

    Returns:
    - dict: Keys are reference IDs, and values are breadth of coverage percentages (float).
    """
    from collections import defaultdict

    def merge_intervals(intervals):
        """
        Merge overlapping or adjacent intervals.
        """
        if not intervals:
            return []
        sorted_intervals = sorted(intervals, key=lambda x: x[0])
        merged = [sorted_intervals[0]]
        for current in sorted_intervals[1:]:
            last = merged[-1]
            if current[0] <= last[1] + 1:
                merged[-1] = (last[0], max(last[1], current[1]))
            else:
                merged.append(current)
        return merged

    def calculate_coverage(reads, regions):
        """
        Calculate the total covered bases within specified regions.
        """
        total_covered = 0
        total_length = 0
        for region_start, region_end in regions:
            total_length += (region_end - region_start + 1)
            # Find reads that overlap with this region
            overlapping_reads = [
                (max(read['start'], region_start), min(read['end'], region_end))
                for read in reads
                if not (read['end'] < region_start or read['start'] > region_end)
            ]
            # Merge overlapping reads within this region
            merged_reads = merge_intervals(overlapping_reads)
            # Sum covered bases
            region_covered = sum(end - start  for start, end in merged_reads)
            total_covered += region_covered
        if total_length == 0:
            return 0.0
        return (total_covered / total_length) * 100

    coverage_dict = {}
    for ref_id, reads in all_reads_dict.items():
        # Filter reads based on removed_read_ids
        if removed_read_ids:
            filtered_reads = [
                read for read in reads
                if not (
                    read['id'] in removed_read_ids and
                    ref_id in removed_read_ids[read['id']]
                )
            ]
        else:
            filtered_reads = reads  # No reads to remove

        # Retrieve regions for the current reference
        if regions_dict and ref_id in regions_dict:
            regions = regions_dict[ref_id]
        else:
            # If no regions provided, consider the entire reference
            if reference_lengths and ref_id in reference_lengths:
                regions = [[1, reference_lengths[ref_id]]]  # 1-based coordinates
            else:
                # Infer from reads
                if not filtered_reads:
                    coverage_dict[ref_id] = 0.0
                    continue
                min_start = min(read['start'] for read in filtered_reads)
                max_end = max(read['end'] for read in filtered_reads)
                regions = [[min_start, max_end]]

        # Calculate coverage over the specified regions
        coverage = calculate_coverage(filtered_reads, regions)
        coverage_dict[ref_id] = coverage

    return coverage_dict

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
    num_workers = os.cpu_count()
    signatures = create_signatures_for_regions(
        cluster_regions,
        all_reads,
        bam_path,
        kmer_size,
        scaled,
        num_workers=num_workers
    )
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


def should_merge(
    current_intervals,
    new_interval,
    statistic_func,
    max_stat_threshold,
    max_group_size,
    max_length=None,
    value_diff_tolerance=None
):
    """
    Decide whether a new_interval can be merged into the current_intervals buffer
    under these constraints:
        - The combined statistic must be < max_stat_threshold
        - The number of intervals in the merged chunk must be <= max_group_size
        - (Optional) The total basepair span <= max_length
        - (Optional) The difference in average (4th-col) values must be <= value_diff_tolerance
    """
    if len(current_intervals) >= max_group_size:
        # Already at capacity for how many intervals we want to merge at once
        return False

    # Combine current_intervals with new_interval
    candidate = current_intervals + [new_interval]

    # 1) Check total length if max_length is defined
    if max_length is not None:
        min_start = min(interval[1] for interval in candidate)
        max_end   = max(interval[2] for interval in candidate)
        if (max_end - min_start) > max_length:
            return False

    # 2) Check combined statistic
    values = [interval[3] for interval in candidate]
    combined_stat = statistic_func(values)
    if combined_stat > max_stat_threshold:
        return False

    # 3) (Optional) Check difference in average value between "neighbors"
    #    In simpler terms, you might just check the difference between
    #    new_interval's value vs. the average of current_intervals' values
    #    if value_diff_tolerance is set.
    if value_diff_tolerance is not None:
        current_avg = sum(iv[3] for iv in current_intervals) / len(current_intervals)
        new_val = new_interval[3]
        if abs(new_val - current_avg) > value_diff_tolerance:
            return False

    # If all checks pass, we can merge
    return True

def merge_bedgraph_regions_by_ref(
    regions,
    statistic_func,
    max_stat_threshold,
    max_group_size=1000
):
    """
    Merges consecutive bedgraph-like intervals per reference (chrom).

    :param regions: List of (chrom, start, end, coverage)
    :param statistic_func: Function to compute a statistic from coverage values (e.g., Gini or variance)
    :param max_stat_threshold: Maximum allowed statistic to keep merging
    :param max_group_size: Max intervals to merge at once before forcing a cutoff
    :return: List of (chrom, merged_start, merged_end, avg_coverage, statistic)
    """
    # Sort first by chrom, then by start
    regions = sorted(regions, key=lambda x: (x[0], x[1]))

    merged_all_refs = []

    # Group the intervals by their reference (chrom)
    for chrom, group_iter in groupby(regions, key=lambda x: x[0]):
        group_list = list(group_iter)  # All intervals for this chrom
        # Merge intervals for this chrom
        merged_for_this_chrom = merge_intervals_for_single_ref(
            group_list,
            statistic_func,
            max_stat_threshold,
            max_group_size
        )
        merged_all_refs.extend(merged_for_this_chrom)

    return merged_all_refs


def merge_intervals_for_single_ref(regions, statistic_func, max_stat_threshold, max_group_size=1000):
    """
    Merges consecutive regions for a single reference, using statistic_func for coverage
    (e.g., Gini or variance).
    Returns a list of (chrom, start, end, avg_coverage, statistic).
    """
    # Sort by start, though we assume they are all same chrom
    regions = sorted(regions, key=lambda x: x[1])

    merged = []

    # Accumulators for the current merge group
    current_chrom = None
    current_start = None
    current_end = None
    coverage_values = []
    region_count = 0

    def finalize_group():
        """Finalize the current merge group and append to merged list."""
        if not coverage_values:
            return
        avg_coverage = sum(coverage_values) / len(coverage_values)
        stat_value = statistic_func(coverage_values)
        merged.append((current_chrom, current_start, current_end, avg_coverage))

    for (chrom, start, end, cov) in regions:
        if current_chrom is None:
            # Initialize the first group
            current_chrom = chrom
            current_start = start
            current_end = end
            coverage_values = [cov]
            region_count = 1
            continue

        # Check if the next interval is contiguous with the current
        # (you could relax this condition if you want to merge even if not contiguous)
        if start == current_end:
            # Tentatively add coverage
            extended_values = coverage_values + [cov]
            test_stat = statistic_func(extended_values)

            # Check if adding this interval keeps us under threshold and max_group_size
            if test_stat <= max_stat_threshold and region_count < max_group_size:
                coverage_values.append(cov)
                current_end = end
                region_count += 1
            else:
                # Statistic or group size exceeded
                finalize_group()

                # Start a new group
                current_chrom = chrom
                current_start = start
                current_end = end
                coverage_values = [cov]
                region_count = 1
        else:
            # Different segment or gap => finalize the current group and start fresh
            finalize_group()
            current_chrom = chrom
            current_start = start
            current_end = end
            coverage_values = [cov]
            region_count = 1

    # Final group
    finalize_group()

    return merged

def determine_conflicts(
        output_dir = None,
        input_bam = None,
        matrix = None,
        min_threshold = 0.2,
        abu_file = None,
        min_similarity_comparable = 0.0,
        use_variance = False,
        apply_ground_truth = False,
        sigfile = None,
        bedfile = None,
        reference_signatures = None,
        scaled = 100,
        kmer_size=31,
        FAST_MODE=True
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

    # Step 8: Load all reads from BAM once
    print("Loading reads from BAM...")
    all_reads, reflengths = load_reads_from_bam(input_bam)
    end = time.time()
    print(f"Time to complete import of bed information: {end-start} seconds")

    # filter all _reads to be only NC_006998.1 or NC_003310.1
    # all_reads = {k: v for k, v in all_reads.items() if k in ['NC_006998.1', 'NC_003310.1']}
    # get the total number of reads in all_reads (len of values in all_reads)
    # total_reads = sum([len(v) for k, v in all_reads.items()])
    start_stime = time.time()
    try:

        breadth_old = calculate_breadth_of_coverage_dict(all_reads, reference_lengths=reflengths)
    except Exception as e:
        print(f"Error calculating coverage with samtools, attempting another way internally...: {e}")
        try:
            breadth_old = create_breadth_coverage_pysam(input_bam, None)
        except Exception as e:
            print(f"Error calculating coverage with internal method: {e}")
            breadth_old = {}
    print(f"Total references with reads: {len(all_reads)}. Done in {time.time()-start_stime} seconds")
    gt_performances = dict()



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
        clusters = [set(all_reads.keys())]
    # Step 7: Parse filtered BED to get regions
    regions = parse_bed_file(bedfile)
    # filter regions only present in all_reads
    regions = [r for r in regions if r[0] in all_reads.keys()]
    print(f"Total regions defined: {len(regions)}")
    # get unique count of read_ids in all_reads
    unique_read_ids = defaultdict(int)
    for ref, reads in all_reads.items():
        for read in reads:
            unique_read_ids[read.get('id')]+=1
    total_gt1 = sum([1 for k, v in unique_read_ids.items() if v > 1])
    total_sum = sum([v for k, v in unique_read_ids.items()])
    print(f"Total unique reads: {len(unique_read_ids)}, sum {total_sum}, sumgt1: {total_gt1}")
    # 2. Choose statistic function
    use_variance=False
    if use_variance:
        from statistics import pvariance
        statistic_func = lambda vals: pvariance(vals)  # population variance
        stat_name = "variance"
    else:
        statistic_func = compute_gini
        stat_name = "gini"
    # 3. Merge

    threshold=0.8

    merged_regions = merge_bedgraph_regions(
        regions,
        statistic_func=statistic_func,
        max_stat_threshold=threshold,
        max_group_size=600  #Limit merges to x intervals
    )

    # merged_regions = merge_bedgraph_regions_by_ref(
    #     regions,
    #     statistic_func=statistic_func,
    #     max_stat_threshold=threshold,
    #     max_group_size=40  #Limit merges to x intervals
    # )



    # 4. Print or save results
    print(f"Merged regions (using {stat_name} <= {threshold}):")
    print(f"Length of original regions : {len(regions)})")
    print(f"Length of merged regions: {len(merged_regions)}")



    if not sigfile or not os.path.exists(sigfile):
        # Step 9: Create signatures for each region
        print("Creating signatures for regions...")
        num_workers = os.cpu_count()
        signatures = create_signatures_for_regions(
            merged_regions,
            all_reads,
            input_bam,
            kmer_size,
            scaled,
            num_workers=num_workers
        )
        print(f"Total signatures created: {len(signatures)}")
        # Step 10: Save signatures to files
        sig_dir = os.path.join(output_dir, "signatures")
        single_sigfile = os.path.join(sig_dir, "merged_regions.sig")
        save_signatures_sourmash(signatures, single_sigfile)
    else:
        # read in the signatures
        print(f"loading signatures from sigfile: {sigfile}")
        loaded_sigs = load_signatures_sourmash(sigfile)
        signatures = rebuild_sig_dict(loaded_sigs)

    output_csv = os.path.join(output_dir, "region_comparisons.csv")
    print(f"Comparison results written to: {output_csv}")
    print("Comparing signatures pairwise...")
    print(f"Clusters to compare", )
    for i, c in enumerate(clusters):
        print(f"C{i}:", c)

    # # get count signautres by reference
    # sig_count = defaultdict(int)
    # for region1, sigs in signatures.items():
    #     ref, start_end = region1.split(":")
    #     sig_count[ref] +=1
    # print(f"Total signatures by reference")
    # for k, v in sig_count.items():
    #     print(f"\t{k}: {v}")
    signatures = list(signatures.items())
    # get index of all_reads where read.id is NC_003310.1:178327-181685
    # sigs_to_compare = []
    # for sig in signatures:
    #     if sig[0] == 'NC_003310.1:178327-181685' or sig[0] == 'NC_006998.1:9528-183447':
    #         print(sig)
    #         sigs_to_compare.append(sig)
    # # compare signatures in sigs_to_compare
    # print(f"Comparing {len(sigs_to_compare)} signatures")
    # for i, sig in enumerate(sigs_to_compare):
    #     print(f"Comparing {i+1} of {len(sigs_to_compare)}")
    #     mh1 = sig[1].minhash
    #     for j, sig2 in enumerate(sigs_to_compare):
    #         if i == j:
    #             continue
    #         mh2 = sig2[1].minhash
    #         jaccard = mh1.jaccard(mh2)
    #         print(jaccard)
    #         if jaccard >= min_threshold:
    #             c1_in_2 = mh1.avg_containment(mh2)
    #             c2_in_1 = mh2.avg_containment(mh1)



    # exit()
    if FAST_MODE:
        print("Building SBT index for fast mode...")
        sbt_index = build_sbt_index(
            signatures,
            ksize=51,
            clusters=clusters,
        )
        print("Searching SBT...")
        sum_comparisons = fast_mode_sbt(
            signatures,
            sbt_index,
            output_csv,
            min_threshold,
            clusters
        )
    else:
        print("Using linear pairwise comparison (slow mode)...")
        sum_comparisons = slow_mode_linear(
            signatures,
            output_csv,
            min_threshold,
            cluster_map=clusters
        )

    print("Building Conflict Groups")
    # 1) Build conflict groups from sum_comparisons
    start_time = time.time()
    conflict_groups = build_conflict_groups(
        sum_comparisons,
        min_jaccard=0.0
    )
    print(f"Conflict groups {len(conflict_groups)} built in {time.time() - start_time:.2f} seconds. Next up is proportion removal")
    start_time = time.time()
    # 2) For each group, remove 1 read from each region
    removed_read_ids = finalize_proportional_removal(
        conflict_groups,
        all_reads,
        fetch_reads_in_region,
        remove_mode='random'
    )
    print(removed_read_ids)
    exit()


    start_time = time.time()
    filtered_bam = os.path.join(output_dir, "filtered.hash.bam")
    comparison_df = None
    includable_read_ids = dict()
    for ref, reads in all_reads.items():
        for read in reads:
            if read.get('id') not in removed_read_ids:
                if ref not in includable_read_ids:
                    includable_read_ids[ref] = []
                includable_read_ids[(ref, read.get('id'))] = True
    try:

        print(f"Completed Bread of coverage setting in {time.time() - start_time}. Filtering to a new bamfile")
        create_filtered_bam(
            input_bam,
            filtered_bam,
            includable_read_ids.keys()
        )
    except Exception as e:
        print(f"Error while filtering BAM: {e}")

    finally:

        print(f"Done with proportional removal in built in {time.time() - start_time:.2f} seconds. Calculating Bread of coverage for the new set of reads")
        try:
            breadth = calculate_breadth_of_coverage_dict(all_reads, reference_lengths=reflengths, skip_ids=removed_read_ids)
        except Exception as e:
            print(f"Error calculating coverage with samtools, attempting another way internally...: {e}")
            try:
                breadth = create_breadth_coverage_pysam(filtered_bam, None)
            except Exception as e:
                print(f"Error calculating coverage with internal method: {e}")
                breadth = {}
        # # 3. If user provided a reference comparison CSV
        # for k, v in breadth_old.items():
        #     print(k, v)
        # print("-----------")
        # for k, v in breadth.items():
        #     print(k, v)
   # Create filtered BAM
    # if use_reads_gt:
    for ref, reads in all_reads.items():
        for i, read in enumerate(reads):
            gt_read = read['id'].split("_")
            gt_name = "_".join(gt_read[:-2])
            all_reads[ref][i]['gt'] = gt_name

    per_ref, overall = compute_f1_scores(all_reads)
    # Print per-reference scores
        # Display per-reference results
    # Display micro-averaged results
    gt_performances['original'] = dict(overall=overall, per_ref = per_ref)
    all_reads, _ = load_reads_from_bam(filtered_bam)
    for ref, reads in all_reads.items():
        for i, read in enumerate(reads):
            gt_read = read['id'].split("_")
            gt_name = "_".join(gt_read[:-2])
            all_reads[ref][i]['gt'] = gt_name
    print("Calculating f1 scores of alignment first....")
    per_ref, overall = compute_f1_scores(
        all_reads
    )
    gt_performances['new'] = dict(overall=overall, per_ref = per_ref)
    # Compare the metrics
    comparison_df = compare_metrics(
        gt_performances['original']['per_ref'],
        per_ref,
        breadth_old,
        breadth
    )
    try:
        # print df to output_dir/removal_stats.xlsx
        comparison_df.to_excel(os.path.join(output_dir, "removal_stats.xlsx"), index=False)
            # Display the comparison
        print("\n=== Metrics Comparison ===\n")
        # convert 'Δ All%' to float
        comparison_df['Δ All%'] = comparison_df['Δ All%'].astype(float)

        print(comparison_df[comparison_df['Δ All%'] != 0 ][['Reference', 'Δ All%', 'Δ Breadth%', 'Breadth New', 'Breadth Original', 'TP New', 'TP Original', 'FP Original', 'FP New']].to_string(index=False))
    except Exception as ex:
        print(f"Error while saving comparison to Excel and printing it to console: {ex}")
    if apply_ground_truth:
        print("Original Metrics Overall: ")
        print(f"\tPrecision: {gt_performances['original']['overall']['micro_precision']:.4f}")
        print(f"\tRecall: {gt_performances['original']['overall']['micro_recall']:.4f}")
        print(f"\tF1: {gt_performances['original']['overall']['micro_f1']:.4f}")
        print("New Filtered Metrics Overall: ")
        print(f"\tPrecision: {gt_performances['new']['overall']['micro_precision']:.4f}")
        print(f"\tRecall: {gt_performances['new']['overall']['micro_recall']:.4f}")
        print(f"\tF1: {gt_performances['new']['overall']['micro_f1']:.4f}")
    # if abu_file:
    #     # Parse ground truth coverage
    #     ground_truth = parse_ground_truth(abu_file)
    #     # Step 14: Compute coverage stats before and after filtering
    #     print("Computing coverage statistics before filtering...")
    #     original_coverage_stats = compute_breadth_and_depth(input_bam)
    #     print("Computing coverage statistics after filtering...")
    #     filtered_coverage_stats = compute_breadth_and_depth(filtered_bam)

    #     # Step 15: Report coverage stats compared to ground truth
    #     print("Reporting coverage statistics...")
    #     report_coverage_stats(ground_truth, original_coverage_stats, filtered_coverage_stats, output_dir)

    #     # Step 16: Compute performance metrics
    #     print("Computing performance metrics based on ground truth and filtered BAM...")
    #     TP, FP, FN, precision, recall, f1 = compute_performance(ground_truth, filtered_bam)
    #     print("Performance Metrics:")
    #     print(f"TP: {TP}, FP: {FP}, FN: {FN}")
    #     print(f"Precision: {precision:.4f}")
    #     print(f"Recall: {recall:.4f}")
    #     print(f"F1-Score: {f1:.4f}")
        return filtered_bam, comparison_df
    else:
        return filtered_bam, comparison_df
