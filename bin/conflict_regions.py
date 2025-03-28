#!/usr/bin/env python3
import os
import csv
import pysam
import time
from sourmash.sbt import SBT, GraphFactory
from sourmash.sbtmh import SigLeaf
from sourmash import MinHash, SourmashSignature, save_signatures, load_file_as_signatures
import statistics
from typing import List, Tuple, Optional, Dict
import numpy as np
import math
import concurrent.futures
from collections import defaultdict
from functools import partial
import pandas as pd
import pysam
from tqdm import tqdm
import time




def parse_bed_file(bed_file_path):
    """
    Parse a BED file using pandas and return a DataFrame with optimized dtypes.
    """
    # Define column names
    column_names = ['chrom', 'start', 'end', 'depth']

    # Read the BED file with specified dtypes
    df = pd.read_csv(
        bed_file_path,
        sep='\t',
        header=None,
        names=column_names,
        dtype={
            'chrom': 'category',       # Use 'category' dtype for repeated strings
            'start': 'int32',          # Use smallest sufficient integer type
            'end': 'int32',
            'depth': 'float32'         # Assuming 'depth' can be a float; adjust as needed
        }
    )
    return df

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
    reads_map = {}
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
            if refname not in reads_map:
                reads_map[refname] = []
            reads_map[refname].append({
                'id': read.query_name,
                # 'seq': seq,
                'start': read.reference_start,
                'end': read.reference_end,
                'ref': refname
            })
    bam.close()
    return reads_map, reflengths

def fetch_reads_in_region(reads_map, refname, start, end):
    """
    Given the pre-loaded reads dictionary, return all reads overlapping [start, end)
    of the given refname. Overlap means the read's alignment region intersects [start, end).

    start and end are 0-based coordinates.
    """
    results = []
    if refname not in reads_map:
        return results
    for read in reads_map[refname]:
        # Check overlap
        # A read overlaps the region if read.end > start and read.start < end
        if read['end'] > start and read['start'] < end:
            results.append((read['id'], refname))
    return results


def create_signature_for_single_region(args) -> Tuple[str, 'SourmashSignature']:
    """
    Worker function to create a signature for a single region.
    This function is designed to be used with parallel processing.

    :param args: Tuple containing (refname, start, end, kmer_size, scaled, reads_map)
    :return: Tuple of (region_name, SourmashSignature)
    """
    refname, start, end, kmer_size, scaled, reads_map = args
    # Fetch reads in the specified region
    if not reads_map:
        return None  # No reads in this region

    # Create a single MinHash for this region
    mh = MinHash(n=0, ksize=kmer_size, scaled=scaled)
    for read in reads_map:
        mh.add_sequence(read.get('seq'), force=True)

    # Create signature
    region_name = f"{refname}:{start + 1}-{end}"  # 1-based coords for display
    sig = SourmashSignature(mh, name=region_name)
    return (region_name, sig)


global_bam = None
global_fastas  = []

def init_worker(bam_path, fasta_paths):
    """Initializer for each worker process; open the BAM file once."""
    global global_bam
    global global_fastas
    global_bam = pysam.AlignmentFile(bam_path, "rb")
    # append all fasta to global_fasta
    for fasta_path in fasta_paths:
        # if ends with .gz treat as compressed
        if fasta_path.endswith('.gz'):
            # index fasta file
            fasta = pysam.FastxFile(fasta_path )
        else:
            fasta = pysam.FastxFile(fasta_path)
        global_fastas.append(fasta)

def process_region(region, kmer_size, scaled):
    """
    Process a single region using the globally opened BAM file.
    :param region: Tuple (chrom, start, end, mean_depth)
    :param kmer_size: K-mer size for signature.
    :param scaled: Scaling factor for signature.
    :return: (region_name, signature, chrom, region_reads) or None if no signature.
    """
    global global_bam
    global global_fastas
    chrom, start, end, mean_depth = region
    reads_in_region = []
    # Use the globally opened BAM file.

    result=None
    if len(global_fastas) > 0 :
        # fetch the sequence from fasta
        seen_references = False
        for fasta in global_fastas:
            # check if chrom is in fasta first
            # get references from pysam.fastxfile
            for record in fasta:
                # Check if the record corresponds to the desired chromosome
                if record.name != chrom:
                    continue

                # Slice the sequence from 'start' to 'end'
                seq = record.sequence[start:end]
                if seq:
                    reads_in_region.append({
                        "id": f"{chrom}:{start}-{end}",
                        "seq": seq,
                        "start": start,
                        "end": end
                    })
                args = [chrom, start, end, kmer_size, scaled, reads_in_region]
                result = create_signature_for_single_region(args)
                seen_references = True
                # Once the desired record is found and processed, you can break out of the loop.
                break
    else:
        try:
            for read in global_bam.fetch(chrom, start, end):
                reads_in_region.append({
                    "id": read.query_name,
                    "seq": read.seq,
                    "start": read.reference_start,
                    "end": read.reference_end
                })
            args = [chrom, start, end, kmer_size, scaled, reads_in_region]
            result = create_signature_for_single_region(args)
        except Exception as e:
            print(f"Error processing region {region}: {e}")
    if result is None:
        return None
    region_name, sig = result
    return (region_name, sig, chrom, reads_in_region)


def create_signatures_for_regions(
    regions_df: pd.DataFrame,
    bam_path: str,
    fasta_paths: list,
    kmer_size: int,
    scaled: float,
    num_workers: int = 8
) -> (dict, dict, dict): # type: ignore
    """
    For each region in the DataFrame, fetch reads and create a Sourmash signature.
    Parallelized over multiple processes.

    :param regions_df: DataFrame with columns ['chrom', 'start', 'end', 'mean_depth']
    :param bam_path: Path to the BAM file.
    :param kmer_size: Size of k-mers for MinHash.
    :param scaled: Scaling factor for MinHash.
    :param num_workers: Number of parallel processes to use.
    :return: Tuple (signatures, reads_map, reflengths)
             signatures: dict mapping region names to signatures.
             reads_map: dict mapping chrom to list of read dicts.
             reflengths: dict mapping reference names to lengths.
    """
    # Validate columns.
    required_columns = {'chrom', 'start', 'end'}
    if not required_columns.issubset(regions_df.columns):
        raise ValueError(f"Input DataFrame must contain columns: {required_columns}")

    # Ensure proper types.
    regions_df = regions_df.copy()
    regions_df['chrom'] = regions_df['chrom'].astype(str)
    regions_df['start'] = regions_df['start'].astype(int)
    regions_df['end'] = regions_df['end'].astype(int)

    # Convert DataFrame rows into a list of tuples.
    regions = regions_df.itertuples(index=False, name=None)  # (chrom, start, end, mean_depth)

    signatures = {}
    reads_map = defaultdict(list)


    num_workers = 1
    init_worker(bam_path, fasta_paths)
    for region in tqdm(regions, total=regions_df.shape[0], miniters=1_000, desc="Processing regions"):
        result = process_region(region, kmer_size, scaled)
        if result is not None:
            region_name, sig, _, _ = result
            signatures[region_name] = sig
    print(f"\tStart processing pool now... with: {num_workers}")

    # with concurrent.futures.ProcessPoolExecutor(
    #         max_workers=num_workers,
    #         initializer=init_worker,
    #         initargs=(bam_path,fasta_paths, )) as executor:
    #     #  bind kmer_size and scaled to process_region.
    #     process_region_partial = partial(process_region, kmer_size=kmer_size, scaled=scaled)
    #     # map over regions without using a lambda.
    #     results = executor.map(process_region_partial, regions, chunksize=1)
    #     for result in tqdm(results, total=regions_df.shape[0], desc="Processing regions"):
    #         if result is not None:
    #             region_name, sig, _, _ = result
    #             signatures[region_name] = sig
    # close the global BAM file.
    global global_bam
    if global_bam is not None:
        global_bam.close()
    return signatures

from collections import Counter

def consensus_reference_from_bam(bam_path, region):
    """
    Construct a consensus "reference" sequence for a given region from a BAM file.

    Parameters:
      bam_path (str): Path to the BAM file.
      region (tuple): (chrom, start, end, mean_depth) defining the region.

    Returns:
      str: The consensus sequence across the region.
    """
    chrom, start, end, mean_depth = region
    # Dictionary to collect bases observed at each reference position
    bases_at_pos = defaultdict(list)

    # Open the BAM file and iterate over all reads overlapping the region
    with pysam.AlignmentFile(bam_path, "rb") as bam_fs:
        for read in bam_fs.fetch(chrom, start, end):
            # print the reference_sequence
            # get_aligned_pairs(with_seq=True) returns tuples: (query_pos, ref_pos, base)
            for qpos, rpos, base in read.get_aligned_pairs(with_seq=True):
                # Skip positions not aligned to the reference
                if rpos is None:
                    continue
                # Only consider positions within the region of interest
                if rpos < start or rpos >= end:
                    continue
                if base is not None:
                    bases_at_pos[rpos].append(base)

    # Build the consensus sequence for each position in the region
    consensus_bases = []
    for pos in range(start, end):
        if pos in bases_at_pos:
            # Majority vote at this position
            most_common_base, count = Counter(bases_at_pos[pos]).most_common(1)[0]
            consensus_bases.append(most_common_base)
        else:
            # If no read covers this position, insert an 'N'
            consensus_bases.append('N')

    return ''.join(consensus_bases)
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

    for region_name, sig in tqdm(siglist, miniters=1_000, desc="Adding signatures to SBT"):
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
    print(f"Processing all signatures in SBT index")
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
                sig1, threshold=min_threshold, best_only=False
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

def create_filtered_bam(bam_in, output_bam_path, alignments_to_skip):
    """
    Create a new BAM file containing only reads with IDs in read_ids_to_keep.
    """
    bam_out = pysam.AlignmentFile(output_bam_path, "wb", template=bam_in)
    # referesh bam iteration
    bam_in.reset()
    set_reads_refs = set()
    count_total = 0
    count_written = 0
    for read in bam_in:
        count_total += 1
        if (read.query_name not in alignments_to_skip) and (read.reference_name not in alignments_to_skip.get(read.query_name, [])):
            bam_out.write(read)
            count_written += 1
        set_reads_refs.add(read.query_name)
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


def finalize_removed_reads(remove_counts, reads_map, fetch_reads_in_region,
                          remove_mode='first',
                          random_seed=None):
    """
    Convert the 'remove_counts' constraints into a global set of read IDs to remove,
    and also gather stats on how many reads are removed per reference.

    Args:
        remove_counts: dict keyed by (ref, start, end) -> int
                       indicating how many reads to remove from that region.
        reads_map: dict from refname -> list of read dicts, as loaded from the BAM.
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

        reads = fetch_reads_in_region(reads_map, ref, start, end)
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


def pairwise_remove_reads(sum_comparisons, reads_map):
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
                reads_r1 = fetch_reads_in_region(reads_map, region1[0], region1[1], region1[2])
                reads_r2 = fetch_reads_in_region(reads_map, region2[0], region2[1], region2[2])

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



def remove_symmetric(conflict_group, reads_map, fetch_reads_in_region, remove_n=1):
    """
    conflict_group: a set of regions, e.g. { (r1, s1, e1), (r2, s2, e2), ... }
    remove_n: the number of reads to remove from each region in this group.

    Returns: dict { region -> set_of_removed_read_ids }

    If a region has fewer than remove_n reads, we'll remove them all (or remove as many as possible).
    """
    removed = {}
    for region in conflict_group:
        (ref, start, end) = region
        reads = fetch_reads_in_region(reads_map, ref, start, end)
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

def finalize_proportional_removal(conflict_groups, bam_fs, fetch_reads_in_region, remove_mode='random', random_seed=None):
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
    global global_bam
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
                # region_reads = fetch_reads_in_region(reads_map, ref, start, end)
                reads = bam_fs.fetch(ref, start, end)
                region_reads = [(read.query_name, ref) for read in reads]
                unified_reads.extend(region_reads)
            # remove duplicates if the same read appears in multiple sub-regions
            # using a dict or set keyed by read_id
            unique_ids = {}
            for (r_id, r_ref) in unified_reads:
                if r_id not in unique_ids:
                    unique_ids[r_id] = (r_id, r_ref)
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
                for (read_id, r_ref) in rlist:
                    removed_read_ids[read_id].append(r_ref)
            else:
                # remove exactly min_cov reads
                if remove_mode == 'random':
                    sampled = random.sample(rlist, min_cov)
                else:
                    # remove the first min_cov
                    sampled = rlist[:min_cov]
                for (read_id, r_ref) in sampled:
                    removed_read_ids[read_id].append(r_ref)
    return removed_read_ids

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


def compute_f1_scores(reads_map):
    """
    Computes precision, recall, and F1 scores by comparing ground truth (gt) to predicted references.

    Parameters:
        reads_map (dict): A dictionary where each key is a reference and each value is a list of read dictionaries.
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
    for ref, reads in reads_map.items():
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

def compare_metrics(reads_map):
    """
    Create a DataFrame comparing metrics for each reference from the reads_map dictionary.

    The expected keys in each reference's dictionary are:
      - TP: True Positives
      - FP: False Positives
      - FN: False Negatives
      - total_reads
      - proportion_aligned
      - precision
      - recall
      - f1
      - breadth_old: Original breadth (from an external computation)
      - breadth: New breadth (from an external computation)

    This function also calculates the delta in breadth and the relative percentage change.

    Returns:
      pd.DataFrame: A DataFrame with one row per reference.
    """
    data = []
    for ref, stats in sorted(reads_map.items()):
        TP = stats.get('TP Original', 0)
        FP = stats.get('FP Original', 0)
        FN = stats.get('FN Original', 0)
        FN_NEW = stats.get('FN New', 0)
        TP_NEW = stats.get('TP New', 0)
        total_reads = stats.get('total_reads', 0)
        pass_filtered_reads = stats.get('pass_filtered_reads', 0)
        prop_aligned = stats.get('proportion_aligned', 0)
        delta_reads = pass_filtered_reads - total_reads
        precision = stats.get('precision', 0)
        recall = stats.get('recall', 0)
        f1 = stats.get('f1', 0)
        breadth_old_val = stats.get('breadth_old', 0)
        breadth_new_val = stats.get('breadth', 0)
        delta_breadth = breadth_new_val - breadth_old_val
        delta_breadth_perc = (breadth_new_val / breadth_old_val) if breadth_old_val > 0 else 0
        delta_all_readspercent = 100*(delta_reads / total_reads) if total_reads > 0 else 0
        # if not whole number set to 2 decimal places else set to 0 decimal places
        # delta_all_readspercent = round(delta_all_readspercent, 2) if delta_all_readspercent % 1 == 0 else int(delta_all_readspercent)

        data.append({
            'Reference': ref,
            'TP Original': TP,
            'FP Original': FP,
            'FN Original': FN,
            'TP New': TP_NEW,
            'FP New': FN_NEW,
            'Total Reads': total_reads,
            'Pass Filtered Reads': pass_filtered_reads,
            'Proportion Aligned': prop_aligned,
            'Precision': precision,
            'Recall': recall,
            'F1': f1,
            'Δ All': delta_reads,
            'Δ All%': delta_all_readspercent,
            'Breadth Original': breadth_old_val,
            'Breadth New': breadth_new_val,
            'Δ Breadth': delta_breadth,
            'Δ^-1 Breadth': delta_breadth_perc
        })

    df = pd.DataFrame(data)

    # # Optionally, print summary stats
    # total_fp = df['FP'].sum()
    # print(f"Total FP: {total_fp}")

    return df


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
        bam.close()
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

def calculate_breadth_coverage_from_bam(bam_fs, reflengths, removed_ids):
    """
    Calculate the breadth of coverage for each reference in the BAM file,
    iterating over each reference in reflengths. For each reference,
    it fetches reads across the entire reference, skipping any read whose ID
    appears in removed_ids for that reference.

    Parameters:
      bam_fs (pysam.AlignmentFile): An open pysam AlignmentFile (BAM) object.
      reflengths (dict): Dictionary with keys as reference names (str) and values
                         as the total length (int) of that reference.
      removed_ids (dict): Dictionary where keys are read IDs (str) and values are lists
                          of reference names (str) for which the read should be skipped.

    Returns:
      dict: A dictionary with keys as reference names and values as the breadth
            of coverage (percentage) computed over the full reference.
    """
    coverage_dict = {}
    bam_fs.reset()
    # Iterate over each reference in reflengths.
    for ref, ref_length in reflengths.items():
        total_covered = 0
        current_start = None
        current_end = None

        # Fetch all reads spanning the entire reference.
        for read in bam_fs.fetch(ref, 0, ref_length):
            read_id = read.query_name
            # Skip this read if it is in the removed_ids for this reference.
            if read_id in removed_ids and ref in removed_ids[read_id]:
                continue

            # Clamp the read's aligned region to the reference boundaries.
            read_start = max(read.reference_start, 0)
            read_end = min(read.reference_end, ref_length)
            if read_start >= read_end:
                continue  # Ignore reads that don't contribute coverage

            # Merge intervals on the fly.
            if current_start is None:
                current_start, current_end = read_start, read_end
            else:
                if read_start <= current_end + 1:
                    # Overlapping or adjacent; extend current interval.
                    current_end = max(current_end, read_end)
                else:
                    # No overlap; add the current interval length and start a new one.
                    total_covered += current_end - current_start
                    current_start, current_end = read_start, read_end

        # Add the final interval if any.
        if current_start is not None:
            total_covered += current_end - current_start

        # Calculate breadth as a percentage of the reference length.
        breadth = (total_covered / ref_length) * 100 if ref_length > 0 else 0.0
        coverage_dict[ref] = breadth

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

def compute_gini(values: List[float]) -> float:
    """Compute the Gini coefficient of a list of values."""
    sorted_values = sorted(values)
    n = len(values)
    cumulative = 0
    cumulative_values = 0
    for i, value in enumerate(sorted_values, 1):
        cumulative += value
        cumulative_values += cumulative
    if cumulative == 0:
        return 0.0
    return (2 * cumulative_values) / (n * cumulative) - (n + 1) / n
def merge_bedgraph_regions(
    intervals: pd.DataFrame,
    merging_method: str = 'jump',
    max_stat_threshold: Optional[float] = None,
    max_group_size: int = 20000,
    max_length: Optional[int] = None,
    value_diff_tolerance: Optional[float] = None,
    breadth_allowance: Optional[int] = 1000,
    gap_allowance: Optional[float] = 0.1,
    reflengths: Optional[Dict[str, int]] = None,
) -> pd.DataFrame:
    """
    Merges consecutive regions in a BEDGRAPH-like DataFrame based on the specified merging method.

    :param intervals: pandas DataFrame with columns ['chrom', 'start', 'end', 'depth']
    :param merging_method: 'variance', 'gini', or 'jump'
    :param max_stat_threshold: Threshold for the statistic (or for the jump method, it can be used as a jump threshold)
    :param max_group_size: Maximum number of intervals to merge at once
    :param max_length: (Optional) Maximum length of merged interval
    :param value_diff_tolerance: (Optional) Maximum allowed difference in average values
    :param breadth_allowance: Allowable gap (in bp) between intervals to consider merging.
                              (Used only when reflengths are provided.)
    :param gap_allowance: If reflengths is provided, the allowed gap is computed as ref_length * gap_allowance.
                          Otherwise, if reflengths is None, gap condition will always pass.
    :param reflengths: (Optional) Dictionary mapping chromosome names to their reference lengths.
                      If not provided, the gap check will not restrict merging.
    :return: pandas DataFrame with merged intervals and mean depth.
    """
    if intervals.empty:
        return pd.DataFrame(columns=['chrom', 'start', 'end', 'mean_depth'])

    # Ensure the DataFrame is sorted by chromosome and start position.
    intervals = intervals.sort_values(['chrom', 'start']).reset_index(drop=True)

    # If reference lengths are provided, use them. Otherwise, set to infinity so gap check does not interfere.
    if reflengths is not None:
        intervals['ref_length'] = intervals['chrom'].map(reflengths)
    else:
        intervals['ref_length'] = float('inf')

    # Compute allowed gap per chromosome.
    # (If reflengths is provided, allowed gap = ref_length * gap_allowance.
    #  Otherwise, since ref_length is infinity, gap condition is never triggered.)
    percent_dict = {}
    for chrom in intervals['chrom'].unique():
        ref_length = intervals[intervals['chrom'] == chrom]['ref_length'].iloc[0]
        percent_dict[chrom] = ref_length * gap_allowance

    # Precompute jump threshold for the 'jump' method.
    if merging_method == 'jump':
        if max_stat_threshold is not None:
            effective_jump_threshold = max_stat_threshold
        else:
            # Fill NA depths with 0 and compute the 97.5th quantile of depth differences per chromosome.
            intervals['depth'] = intervals['depth'].fillna(0)
            try:
                q = intervals.groupby('chrom', observed=True)['depth'].quantile(0.975).mean()
                effective_jump_threshold = math.ceil(q) + 1
            except Exception as ex:
                print(ex)
                effective_jump_threshold = 200
    else:
        effective_jump_threshold = None

    merged_regions = []
    # Group intervals by chromosome.
    grouped = intervals.groupby('chrom', observed=True)
    for chrom, group in grouped:
        group = group.reset_index(drop=True)
        if group.empty:
            continue

        # Initialize a buffer with the first interval.
        buffer = {
            'chrom': chrom,
            'start': group.at[0, 'start'],
            'end': group.at[0, 'end'],
            'depths': [group.at[0, 'depth']],
        }
        per_ref = percent_dict.get(chrom, breadth_allowance)

        for i, row in enumerate(group.itertuples(index=False)):
            if i == 0:
                continue  # first row is already in the buffer

            new_start, new_end, new_depth = row.start, row.end, row.depth
            prev_end = buffer['end']
            gap = new_start - prev_end  # gap between intervals

            # Start with merging allowed.
            can_merge = True

            # Check maximum group size.
            if len(buffer['depths']) >= max_group_size:
                can_merge = False

            # Check maximum merged region length.
            if max_length is not None:
                total_length = new_end - buffer['start']
                if total_length > max_length:
                    can_merge = False

            # Apply method-specific merging criteria.
            if can_merge:
                if merging_method == 'jump':
                    # Calculate the jump based on the last depth in the current buffer.
                    last_depth = buffer['depths'][-1]
                    jump = abs(new_depth - last_depth)
                    # Merge only if the gap is within the allowed per_ref and the jump is within effective_jump_threshold.
                    if gap > per_ref or jump > effective_jump_threshold:
                        can_merge = False
                elif merging_method in ['variance', 'gini']:
                    current_depths = buffer['depths'] + [new_depth]
                    if merging_method == 'variance':
                        stat_val = statistics.pvariance([float(x) for x in current_depths])
                    elif merging_method == 'gini':
                        stat_val = compute_gini(current_depths)
                    if max_stat_threshold is not None and stat_val > max_stat_threshold:
                        can_merge = False
                    # Optional: Check average value difference tolerance.
                    if value_diff_tolerance is not None:
                        current_avg = np.mean(buffer['depths'])
                        if abs(new_depth - current_avg) > value_diff_tolerance:
                            can_merge = False
                else:
                    raise ValueError(f"Unknown merging method: {merging_method}")

            if can_merge:
                # Merge the current interval into the buffer.
                buffer['end'] = new_end
                buffer['depths'].append(new_depth)
            else:
                # Finalize the current buffer.
                mean_depth = np.mean(buffer['depths'])
                merged_regions.append({
                    'chrom': buffer['chrom'],
                    'start': buffer['start'],
                    'end': buffer['end'],
                    'mean_depth': mean_depth
                })
                # Start a new buffer with the current interval.
                buffer = {
                    'chrom': chrom,
                    'start': new_start,
                    'end': new_end,
                    'depths': [new_depth],
                }

        # Finalize the last buffer for this chromosome.
        if buffer['depths']:
            mean_depth = np.mean(buffer['depths'])
            merged_regions.append({
                'chrom': buffer['chrom'],
                'start': buffer['start'],
                'end': buffer['end'],
                'mean_depth': mean_depth
            })

    merged_df = pd.DataFrame(merged_regions)
    return merged_df



def determine_conflicts(
        output_dir = None,
        input_bam = None,
        matrix = None,
        min_threshold = 0.2,
        fasta_files = [],
        abu_file = None,
        min_similarity_comparable = 0.0,
        use_variance = False,
        apply_ground_truth = False,
        sigfile = None,
        bedfile = None,
        reference_signatures = None,
        scaled = 100,
        kmer_size=31,
        filtered_bam_create=False,
        FAST_MODE=True,
        sensitive=False,
        cpu_count=None,
        jump_threshold=None,
        gap_allowance = 0.1,
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
    print(f"Starting conflict region detection at {time.ctime()}")
    if len(fasta_files) == 0:
        print("No fasta files provided, using sensitive mode")
    elif sensitive:
        print("Fasta files provided but sensitive mode enabled, using BAM alignments for signature generations")
        fasta_files = []
    else:
        print(f"Using {len(fasta_files)} fasta files for regional comparisons")
    # Step 7: Parse filtered BED to get regions
    regions = parse_bed_file(bedfile)
    # filter regions only present in reads_map
    # regions = [r for r in regions if r[0] in reads_map.keys()]
    # convert regions to a list, each row is an element of the list
    # pause for 20 seconds
    print(f"Total regions defined: {len(regions)} in default mode for signature generation")

    # 2. Choose statistic function
    use_jump = True
    use_variance=False
    if use_variance:
        from statistics import pvariance
        statistic_func = lambda vals: pvariance(vals)  # population variance
        stat_name = "variance"
        threshold=0.8
    elif use_jump:
        statistic_func = lambda vals: max(vals) - min(vals)
        stat_name="jump"
        threshold = jump_threshold
    else:
        statistic_func = compute_gini
        stat_name = "gini"
        threshold=0.8
    # 3. Merge
    bam_fs = pysam.AlignmentFile(input_bam, "rb")
    reflengths = {ref: bam_fs.get_reference_length(ref) for ref in bam_fs.references}

    print(f"Merging regions (using {stat_name} <= {threshold}, Gap: {gap_allowance})...")
    start_time = time.time()

    merged_regions = merge_bedgraph_regions(
        regions,
        max_stat_threshold=threshold,
        max_group_size=4_000_000, # Limit merges to x intervals
        merging_method=stat_name,
        reflengths=reflengths,
        gap_allowance = gap_allowance
    )
    print(f"Regions merged in {time.time() - start_time:.2f} seconds.")


    # 4. Print or save results
    print(f"Merged regions (using {stat_name} <= {threshold}):")
    print(f"Length of original regions : {len(regions)})")
    print(f"Length of merged regions: {len(merged_regions)}")
    start = time.time()
    gt_performances = dict()

    reads_map = defaultdict(list)
    if not sigfile or not os.path.exists(sigfile):
        # Step 9: Create signatures for each region
        # get cpu count / 2, only if it is > 1, round down
        if cpu_count:
            num_workers = cpu_count
        else:
            num_workers = max(1, int(os.cpu_count() / 2))
        # print("Creating signatures for regions...", f"parallelized across {num_workers} workers")

        start_time = time.time()
        signatures = create_signatures_for_regions(
            regions_df=merged_regions,
            bam_path = input_bam,
            fasta_paths = fasta_files,
            kmer_size=kmer_size,
            scaled=scaled,
            num_workers=num_workers
        )

        print(f"Signatures created in {time.time() - start_time:.2f} seconds.")
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
        # set to the names of all signatures
        clusters = [set([x.split(":")[0] for x in (signatures)])]
        # clusters = None

    output_csv = os.path.join(output_dir, "region_comparisons.csv")
    print(f"Comparison results written to: {output_csv}")
    print("Comparing signatures pairwise...")
    print(f"Clusters to compare", len(clusters))
    # for i, c in enumerate(clusters):
    #     print(f"C{i}:", c)
    signatures = list(signatures.items())
    start_stime = time.time()

    print(f"Total references with reads: {len(reads_map)}. Done in {time.time()-start_stime} seconds")


    if FAST_MODE:
        print("Building SBT index for fast mode...")
        start_time = time.time()
        sbt_index = build_sbt_index(
            signatures,
            ksize=51,
            clusters=clusters,
        )
        print(f"SBT index built in {time.time() - start_time:.2f} seconds.")
        start_time = time.time()
        print("Searching SBT...")
        sum_comparisons = fast_mode_sbt(
            signatures,
            sbt_index,
            output_csv,
            min_threshold,
            clusters
        )
        print(f"Comparisons done in {time.time() - start_time:.2f} seconds.")
    else:
        print("Using linear pairwise comparison (slow mode)...")
        start_time = time.time()
        sum_comparisons = slow_mode_linear(
            signatures,
            output_csv,
            min_threshold,
            cluster_map=clusters
        )
        print(f"Comparisons done in {time.time() - start_time:.2f} seconds.")
    print("Building Conflict Groups")

    # 1) Build conflict groups from sum_comparisons
    start_time = time.time()
    conflict_groups = build_conflict_groups(
        sum_comparisons,
        min_jaccard=0.0
    )
    print(f"Conflict groups length: {len(conflict_groups)} built in {time.time() - start_time:.2f} seconds. Next up is proportion removal")
    start_time = time.time()
    # 2) For each group, remove 1 read from each region
    removed_read_ids = finalize_proportional_removal(
        conflict_groups,
        bam_fs,
        fetch_reads_in_region,
        remove_mode='random'
    )
    start_time = time.time()
    comparison_df = None
    includable_read_ids = dict()
    sum_of_ref_aligned = defaultdict(int)
    output_file = os.path.join(output_dir, "failed_reads.txt")
    with open (output_file, 'w') as f:
        for read, refs in removed_read_ids.items():
            for ref in refs:
                f.write(f"{ref}\t{read}\n")
        f.close()

    try:
        # Create new variabels of only the filtered reads that passed based on the dict
        # new_filtered_reads = []
        # size_old_reads_map = sum([len(v) for k, v in reads_map.items()])
        # print(f"Size of old all reads: {size_old_reads_map}")
        # for k, v in sum_of_ref_aligned.items():
        #     print(f"Reference {k} has {v} reads aligned")
        # print(f"Size of new filtered reads: {len(new_filtered_reads)}")
        # print(f"Completed Bread of coverage setting in {time.time() - start_time}")
        if filtered_bam_create:
            create_filtered_bam(
                bam_fs,
                filtered_bam_create,
                removed_read_ids
            )
    except Exception as e:
        print(f"Error while filtering BAM: {e}")

    finally:
        try:
            # breadth_old = calculate_breadth_of_coverage_dict(reads_map, reference_lengths=reflengths)
            breadth_old = calculate_breadth_coverage_from_bam(bam_fs=bam_fs,  reflengths=reflengths, removed_ids=dict())
        except Exception as e:
            print(f"Error calculating coverage with samtools, attempting another way internally...: {e}")
            try:
                breadth_old = create_breadth_coverage_pysam(input_bam, None)
            except Exception as e:
                print(f"Error calculating coverage with internal method: {e}")
                breadth_old = {}
        print(f"Done with proportional removal in built in {time.time() - start_time:.2f} seconds. Calculating Bread of coverage for the new set of reads")
        try:
            breadth = calculate_breadth_coverage_from_bam(bam_fs=bam_fs,  reflengths=reflengths, removed_ids=removed_read_ids)
        except Exception as e:
            print(f"Error calculating coverage with samtools, attempting another way internally...: {e}")
            try:
                breadth = create_breadth_coverage_pysam(filtered_bam_create, None)
            except Exception as e:
                print(f"Error calculating coverage with internal method: {e}")
                breadth = {}
        print(f"Done with the breadth creation in {time.time()-start_time} seconds")

    reads_map = defaultdict(lambda: defaultdict(int))
    print(f"Calculating reads map information from original and new filters")
    for ref, ref_length in reflengths.items():
        reads = bam_fs.fetch(ref, 0, ref_length)
        total_reads = 0
        pass_filtered_reads = 0

        for read in reads:
            total_reads += 1

            # Parse the read query name and get the 'gt_name'
            gt_read = read.query_name.split("_")
            gt_name = "_".join(gt_read[:-2])

            # Original mapping stats (each alignment is counted)
            if gt_name == read.reference_name:
                reads_map[ref]['TP Original'] += 1
            else:
                reads_map[ref]['FP Original'] += 1

            # Determine if this alignment should be counted for the "new" stats.
            # We consider an alignment as "passing" if either:
            #   1. The read is not in removed_read_ids, OR
            #   2. The read is in removed_read_ids, but the current reference is not one of those removed.
            if (read.query_name not in removed_read_ids or
                read.reference_name not in removed_read_ids.get(read.query_name, [])):
                pass_filtered_reads += 1
                if gt_name == read.reference_name:
                    reads_map[ref]['TP New'] += 1
                else:
                    reads_map[ref]['FP New'] += 1
            else:
                # The read was removed for this reference
                if gt_name == read.reference_name:
                    reads_map[ref]['FN New'] += 1
                else:
                    reads_map[ref]['TN New'] += 1

        # Save overall counts for this reference
        reads_map[ref]['total_reads'] = total_reads
        reads_map[ref]['pass_filtered_reads'] = pass_filtered_reads

        # Calculate additional stats
        reads_map[ref]['proportion_aligned'] = (
            reads_map[ref]['TP Original'] / total_reads if total_reads > 0 else 0
        )
        fft = reads_map[ref]['TP Original'] + reads_map[ref]['FP Original']
        reads_map[ref]['precision'] = (reads_map[ref]['TP Original'] / fft if fft > 0 else 0)
        reads_map[ref]['recall'] = (reads_map[ref]['TP Original'] / total_reads if total_reads > 0 else 0)
        prrecallplus = reads_map[ref]['precision'] + reads_map[ref]['recall']
        reads_map[ref]['f1'] = (
            2 * (reads_map[ref]['precision'] * reads_map[ref]['recall']) / prrecallplus
            if prrecallplus > 0 else 0
        )
        reads_map[ref]['breadth_old'] = breadth_old.get(ref, 0)
        reads_map[ref]['breadth'] = breadth.get(ref, 0)
    comparison_df = compare_metrics(
        reads_map
    )
    try:
        # print df to output_dir/removal_stats.xlsx
        comparison_df.to_excel(os.path.join(output_dir, "removal_stats.xlsx"), index=False)
            # Display the comparison
        print("\n=== Metrics Comparison ===\n")
        # convert 'Δ All%' to float
        # comparison_df['Δ All%'] = comparison_df['Δ All%'].astype(float)

        print(comparison_df[comparison_df['Δ All'] != 0 ][['Reference', 'Δ All%', 'Δ^-1 Breadth', 'Breadth New', 'Breadth Original', 'TP New', 'TP Original', 'FP Original', 'FP New']].to_string(index=False))
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

        return removed_read_ids, comparison_df
    else:
        return removed_read_ids, comparison_df
