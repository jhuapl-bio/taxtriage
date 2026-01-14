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
import bisect
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
    global global_bam, global_fastas
    global_bam = pysam.AlignmentFile(bam_path, "rb")
    global_fastas = []

    for fasta_path in fasta_paths:
        # pysam.FastaFile expects an *indexed* fasta (has .fai).
        # It can work with bgzipped fasta too (bgzip) + .gzi typically.
        global_fastas.append(pysam.FastaFile(fasta_path))
def process_region(region, kmer_size, scaled):
    global global_bam, global_fastas
    chrom, start, end, mean_depth = region

    reads_in_region = []
    result = None

    if global_fastas and len(global_fastas) > 0:
        seq = None
        for fasta in global_fastas:
            if chrom in fasta.references:
                seq = fasta.fetch(chrom, start, end)
                break

        if not seq:
            return None

        reads_in_region.append({
            "id": f"{chrom}:{start}-{end}",
            "seq": seq,
            "start": start,
            "end": end
        })

        args = [chrom, start, end, kmer_size, scaled, reads_in_region]
        result = create_signature_for_single_region(args)

    else:
        try:
            for read in global_bam.fetch(chrom, start, end):
                reads_in_region.append({
                    "id": read.query_name,
                    "seq": read.query_sequence,   # slightly safer than read.seq
                    "start": read.reference_start,
                    "end": read.reference_end
                })
            args = [chrom, start, end, kmer_size, scaled, reads_in_region]
            result = create_signature_for_single_region(args)
        except Exception as e:
            print(f"Error processing region {region}: {e}")
            return None

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

def build_reference_window_index(
    fasta_files,
    *,
    ksize=31,
    scaled=200,
    window=10_000,
    step=1_000,
    max_n_frac=0.05,
):
    """
    Returns:
      ref_sigs: dict[window_id] = signature
      sbt_by_contig: dict[contig] = SBT containing ONLY windows from that contig
    """
    ref_sigs = {}
    sigs_by_contig = defaultdict(dict)
    # only do unique fasta files
    for fp in fasta_files:
        label = os.path.basename(fp)
        sigs = window_sigs_from_fasta(
            fp,
            fasta_label=label,
            ksize=ksize,
            scaled=scaled,
            window=window,
            step=step,
            max_n_frac=max_n_frac,
        )
        ref_sigs.update(sigs)

        # group by contig (chrom)
        for wid, sig in sigs.items():
            _, contig, _, _ = parse_window_id(wid)
            sigs_by_contig[contig][wid] = sig

    # build SBT per contig
    sbt_by_contig = {}
    for contig, contig_sigs in sigs_by_contig.items():
        factory = GraphFactory(ksize=ksize, n_tables=1, starting_size=max(1, len(contig_sigs)))
        sbt = SBT(factory)
        for wid, sig in contig_sigs.items():
            sbt.add_node(SigLeaf(wid, sig))
        sbt_by_contig[contig] = sbt

    return ref_sigs, sbt_by_contig

def compare_bed_regions_to_reference_windows(
    bed_siglist,                 # list[(region_name, sig)]
    ref_sbt_by_contig,           # dict[contig] -> SBT
    output_csv,
    *,
    min_threshold=0.10,
    top_k_per_contig=1,          # keep only best hit per contig to avoid spam
    max_total_hits=10,           # cap total hits per query region
):
    """
    Writes CSV like your region_comparisons.csv BUT region2 is a reference window id.
    Also returns sum_comparisons in your usual format:
      sum_comparisons[r1].append({jaccard,to,s2,e2,s1,e1})
    Where 'to' is the matched contig and s2/e2 are window coords.
    """
    sum_comparisons = defaultdict(list)

    with open(output_csv, "w", newline="") as outfh:
        w = csv.writer(outfh)
        w.writerow([
            "reference1","start1","end1",
            "reference2","start2","end2",
            "jaccard","containment_1_in_2","containment_2_in_1",
            "match_window_id"
        ])

        for region1, sig1 in bed_siglist:
            r1, s1, e1 = parse_split(region1)  # your existing parser (ref, start, end)
            mh1 = sig1.minhash

            all_hits = []

            # search all contigs except the query contig
            for contig, sbt in ref_sbt_by_contig.items():
                if contig == r1:
                    continue

                # collect best hits for this contig
                contig_hits = []
                for sr in sbt.search(sig1, threshold=min_threshold, best_only=False):
                    wid = sr.signature.name  # window_id
                    _, m_contig, m_s, m_e = parse_window_id(wid)
                    j = float(sr.score)
                    mh2 = sr.signature.minhash
                    c12 = mh1.avg_containment(mh2)
                    c21 = mh2.avg_containment(mh1)
                    contig_hits.append((j, m_contig, m_s, m_e, c12, c21, wid))

                if contig_hits:
                    contig_hits.sort(key=lambda x: x[0], reverse=True)
                    all_hits.extend(contig_hits[:top_k_per_contig])

            if not all_hits:
                continue

            all_hits.sort(key=lambda x: x[0], reverse=True)
            all_hits = all_hits[:max_total_hits]

            for (j, r2, s2, e2, c12, c21, wid) in all_hits:
                w.writerow([r1, s1, e1, r2, s2, e2, j, c12, c21, wid])

                # your downstream expects this shape
                sum_comparisons[r1].append(dict(
                    jaccard=j,
                    to=r2,
                    s1=s1, e1=e1,
                    s2=s2, e2=e2
                ))
                # optional: make symmetric entry
                sum_comparisons[r2].append(dict(
                    jaccard=j,
                    to=r1,
                    s1=s2, e1=e2,
                    s2=s1, e2=e1
                ))

    return sum_comparisons

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
def filter_sum_comparisons_by_mat(sum_comparisons, mat, min_mat_sim=0.0):
    if not mat or min_mat_sim <= 0:
        return sum_comparisons

    filtered = defaultdict(list)
    for ref, conflicts in sum_comparisons.items():
        for c in conflicts:
            other = c["to"]
            if mat_sim(mat, ref, other) >= min_mat_sim:
                filtered[ref].append(c)
    return filtered



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

def build_region_sigs_from_fasta(regions_df, fasta_path, ksize=31, scaled=200):
    """
    regions_df: columns chrom,start,end (0-based half-open)
    returns: dict region_id -> signature
    region_id example: "chr1:1000-2000"
    """
    fa = pysam.FastaFile(fasta_path)
    sigs = {}

    for row in regions_df.itertuples(index=False):
        chrom, start, end = str(row.chrom), int(row.start), int(row.end)
        if chrom not in fa.references:
            continue
        seq = fa.fetch(chrom, start, end)
        if not seq or len(seq) < ksize:
            continue

        mh = MinHash(n=0, ksize=ksize, scaled=scaled)
        mh.add_sequence(seq, force=True)
        region_id = f"{chrom}:{start}-{end}"
        sigs[region_id] = SourmashSignature(mh, name=region_id)

    fa.close()
    return sigs

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
        # for ref, cov_count in coverage_by_ref.items():
        #     rlist = reads_by_ref[ref]

        #     # how many we need to remove to bring this ref down to min_cov
        #     remove_n = max(0, cov_count - min_cov)
        #     if remove_n == 0:
        #         continue

        #     if remove_mode == 'random':
        #         sampled = random.sample(rlist, remove_n)
        #     else:
        #         sampled = rlist[:remove_n]

        #     for (read_id, r_ref) in sampled:
        #         removed_read_ids[read_id].append(r_ref)
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
            try:
                read_id = read.query_name
                if read.is_unmapped or read.reference_start is None or read.reference_end is None:
                    continue
                # Skip this read if it is in the removed_ids for this reference.
                if read_id in removed_ids and ref in removed_ids[read_id]:
                    continue

                # Clamp the read's aligned region to the reference boundaries.
                read_start = max(read.reference_start, 0)
                read_end = min(read.reference_end, ref_length)
                if read_start >= read_end:
                    continue

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
            except Exception as tr:
                print(f"Error processing read {read.query_name} on reference {ref}: {tr}")
                continue
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


def iter_windows(length: int, window: int, step: int):
    for s in range(0, max(0, length - window + 1), step):
        yield s, s + window


def make_window_id(fasta_label: str, contig: str, start: int, end: int) -> str:
    # Unique, parseable ID for each window
    return f"{fasta_label}::{contig}:{start}-{end}"


def parse_window_id(wid: str) -> Tuple[str, str, int, int]:
    # wid format: fasta_label::contig:start-end
    fasta_label, rest = wid.split("::", 1)
    contig, coords = rest.split(":", 1)
    s, e = coords.split("-", 1)
    return fasta_label, contig, int(s), int(e)


def window_sigs_from_fasta(
    fasta_path: str,
    *,
    fasta_label: Optional[str] = None,
    ksize: int = 31,
    scaled: int = 200,
    window: int = 2000,
    step: int = 500,
    max_n_frac: float = 0.05,
    contigs: Optional[List[str]] = None,
) -> Dict[str, SourmashSignature]:
    """
    Build windowed signatures for all contigs in a FASTA (or a specified subset).
    Returns dict window_id -> signature.
    """
    if fasta_label is None:
        fasta_label = os.path.basename(fasta_path)

    fa = pysam.FastaFile(fasta_path)
    sigs: Dict[str, SourmashSignature] = {}

    contig_list = contigs if contigs is not None else list(fa.references)

    for contig in contig_list:
        if contig not in fa.references:
            continue
        seq = fa.fetch(contig).upper()
        L = len(seq)
        if L < window:
            continue

        for s, e in iter_windows(L, window, step):
            subseq = seq[s:e]
            if len(subseq) < ksize:
                continue
            if subseq.count("N") / len(subseq) > max_n_frac:
                continue

            mh = MinHash(n=0, ksize=ksize, scaled=scaled)
            mh.add_sequence(subseq, force=True)

            wid = make_window_id(fasta_label, contig, s, e)
            sigs[wid] = SourmashSignature(mh, name=wid)

    fa.close()
    return sigs


def build_sbt(sigs: Dict[str, SourmashSignature], ksize: int) -> SBT:
    factory = GraphFactory(ksize=ksize, n_tables=1, starting_size=max(1, len(sigs)))
    sbt = SBT(factory)
    for wid, sig in sigs.items():
        sbt.add_node(SigLeaf(wid, sig))
    return sbt


def report_shared_windows_across_fastas(
    fasta_files: List[str],
    output_csv: str,
    *,
    ksize: int = 31,
    scaled: int = 200,
    window: int = 2000,
    step: int = 500,
    jaccard_threshold: float = 0.10,
    max_hits_per_query: int = 3,
    skip_self_same_fasta: bool = True,
    skip_self_same_contig: bool = True,
    max_windows_per_fasta: Optional[int] = None,
) -> None:
    """
    Compare all window signatures across a list of FASTA files and write a report of shared regions.

    Output CSV columns:
      query_fasta, query_contig, q_start, q_end,
      match_fasta, match_contig, m_start, m_end,
      jaccard, containment_q_in_m, containment_m_in_q
    """
    # 1) Sketch all windows across all FASTAs
    all_sigs: Dict[str, SourmashSignature] = {}
    print(f"Sketching windows across {len(fasta_files)} FASTA(s)...")
    for fp in fasta_files:
        label = os.path.basename(fp)
        sigs = window_sigs_from_fasta(
            fp,
            fasta_label=label,
            ksize=ksize,
            scaled=scaled,
            window=window,
            step=step,
        )
        if max_windows_per_fasta is not None and len(sigs) > max_windows_per_fasta:
            # keep only first N (quick throttling for huge genomes/tons of windows)
            kept = dict(list(sigs.items())[:max_windows_per_fasta])
            sigs = kept
        print(f"  {label}: {len(sigs)} windows")
        all_sigs.update(sigs)

    if not all_sigs:
        raise ValueError("No signatures were generated. Check window/step/ksize or FASTA content.")

    # 2) Build one global SBT index
    print(f"Building SBT over {len(all_sigs)} total windows...")
    sbt = build_sbt(all_sigs, ksize=ksize)

    # 3) Search each window vs the index, report cross-fasta / cross-contig hits
    print(f"Searching (threshold={jaccard_threshold}) and writing report to {output_csv} ...")
    with open(output_csv, "w", newline="") as out:
        w = csv.writer(out)
        w.writerow([
            "query_fasta", "query_contig", "q_start", "q_end",
            "match_fasta", "match_contig", "m_start", "m_end",
            "jaccard", "containment_q_in_m", "containment_m_in_q"
        ])

        for q_wid, q_sig in all_sigs.items():
            q_fa, q_contig, q_s, q_e = parse_window_id(q_wid)

            results = sbt.search(q_sig, threshold=jaccard_threshold, best_only=False)

            hits = []
            mh_q = q_sig.minhash

            for sr in results:
                m_wid = sr.signature.name
                if m_wid == q_wid:
                    continue
                m_fa, m_contig, m_s, m_e = parse_window_id(m_wid)

                if skip_self_same_fasta and (m_fa == q_fa):
                    continue
                if skip_self_same_contig and (m_fa == q_fa and m_contig == q_contig):
                    continue

                j = float(sr.score)
                mh_m = sr.signature.minhash
                c_q_in_m = mh_q.avg_containment(mh_m)
                c_m_in_q = mh_m.avg_containment(mh_q)
                hits.append((j, m_fa, m_contig, m_s, m_e, c_q_in_m, c_m_in_q))

            if not hits:
                continue

            hits.sort(key=lambda x: x[0], reverse=True)
            hits = hits[:max_hits_per_query]

            for (j, m_fa, m_contig, m_s, m_e, c1, c2) in hits:
                w.writerow([
                    q_fa, q_contig, q_s, q_e,
                    m_fa, m_contig, m_s, m_e,
                    f"{j:.6f}", f"{c1:.6f}", f"{c2:.6f}"
                ])

    print("Done.")



def build_similarity_matrix_from_shared_windows(
    report_csv: str,
    output_matrix_csv: str | None = None,
    accession_level: str = "contig",   # "contig" | "fasta" | "fasta_contig"
    score_col: str = "jaccard",        # "jaccard" or "containment_q_in_m" etc.
    agg: str = "max",                  # "max" | "mean" | "median"
    fill_diagonal: float = 1.0,
) -> pd.DataFrame:
    """
    Build a symmetric similarity matrix across accessions that appear in the report.

    - accession_level controls what becomes a node in the matrix:
        * "contig": uses query_contig / match_contig
        * "fasta": uses query_fasta / match_fasta
        * "fasta_contig": uses "fasta::contig"
    - agg controls how multiple window hits between the same pair are reduced to one number.
    """
    df = pd.read_csv(report_csv)

    # Ensure score is numeric (your report writes it as formatted string)
    df[score_col] = pd.to_numeric(df[score_col], errors="coerce")
    df = df.dropna(subset=[score_col])

    def node(row, side: str) -> str:
        fa = row[f"{side}_fasta"]
        ct = row[f"{side}_contig"]
        if accession_level == "fasta":
            return str(fa)
        if accession_level == "fasta_contig":
            return f"{fa}::{ct}"
        # default: contig
        return str(ct)

    df["A"] = df.apply(lambda r: node(r, "query"), axis=1)
    df["B"] = df.apply(lambda r: node(r, "match"), axis=1)

    # Keep only accessions that have at least one hit (already true by construction)
    accessions = sorted(set(df["A"]).union(set(df["B"])))

    # Aggregate per undirected pair (A,B)
    # Normalize so (A,B) == (B,A)
    df["X"] = df[["A", "B"]].min(axis=1)
    df["Y"] = df[["A", "B"]].max(axis=1)

    if agg == "mean":
        pair_scores = df.groupby(["X", "Y"], as_index=False)[score_col].mean()
    elif agg == "median":
        pair_scores = df.groupby(["X", "Y"], as_index=False)[score_col].median()
    else:
        pair_scores = df.groupby(["X", "Y"], as_index=False)[score_col].max()

    # Initialize matrix
    mat = pd.DataFrame(0.0, index=accessions, columns=accessions)

    # Fill scores symmetrically
    for _, row in pair_scores.iterrows():
        x, y, s = row["X"], row["Y"], float(row[score_col])
        mat.at[x, y] = max(mat.at[x, y], s)
        mat.at[y, x] = max(mat.at[y, x], s)

    # Fill diagonal
    np.fill_diagonal(mat.values, fill_diagonal)

    if output_matrix_csv:
        mat.to_csv(output_matrix_csv)

    return mat
import pandas as pd
import numpy as np

def build_support_weighted_matrix(
    report_csv: str,
    *,
    accession_level: str = "contig",   # "contig" | "fasta_contig"
    jaccard_col: str = "jaccard",
    hit_threshold: float = 0.10,
    min_hit_windows: int = 5,          # require at least this many windows supporting the pair
    window_size: int = 2000,
    step: int = 500,
    output_csv: str | None = None,
) -> pd.DataFrame:
    """
    Similarity(i,j) = fraction of i's windows that have >=1 hit to j with jaccard>=hit_threshold.
    Then symmetrize by averaging i->j and j->i.
    """
    df = pd.read_csv(report_csv)
    df[jaccard_col] = pd.to_numeric(df[jaccard_col], errors="coerce")
    df = df.dropna(subset=[jaccard_col])

    # nodes
    if accession_level == "fasta_contig":
        df["A"] = df["query_fasta"].astype(str) + "::" + df["query_contig"].astype(str)
        df["B"] = df["match_fasta"].astype(str) + "::" + df["match_contig"].astype(str)
    else:
        df["A"] = df["query_contig"].astype(str)
        df["B"] = df["match_contig"].astype(str)

    # window id for query side (so we can count unique windows that matched)
    df["A_win"] = df["A"] + ":" + df["q_start"].astype(str) + "-" + df["q_end"].astype(str)

    # keep only sufficiently good hits
    df = df[df[jaccard_col] >= hit_threshold].copy()
    if df.empty:
        raise ValueError("No hits left after applying hit_threshold.")

    # Count how many unique query windows matched each (A,B)
    win_hits = (
        df.groupby(["A", "B"])["A_win"]
          .nunique()
          .reset_index(name="hit_windows")
    )

    # Estimate number of windows per accession A from observed max coordinate range in the report
    # (better: compute from FASTA lengths; but this works for “matrix from report” only)
    # Here we approximate n_windows_A as count of unique query windows appearing for A in the report.
    nwin_A = (
        df.groupby("A")["A_win"]
          .nunique()
          .to_dict()
    )

    win_hits["frac_A_to_B"] = win_hits.apply(
        lambda r: r["hit_windows"] / max(1, nwin_A.get(r["A"], 1)),
        axis=1
    )

    # Apply minimum support requirement (prevents “one perfect window”)
    win_hits = win_hits[win_hits["hit_windows"] >= min_hit_windows].copy()

    accs = sorted(set(win_hits["A"]).union(set(win_hits["B"])))
    mat = pd.DataFrame(0.0, index=accs, columns=accs)

    # Fill directed fractions
    for _, r in win_hits.iterrows():
        mat.at[r["A"], r["B"]] = max(mat.at[r["A"], r["B"]], float(r["frac_A_to_B"]))

    # Symmetrize by averaging A->B and B->A
    mat_sym = mat.copy()
    for i in accs:
        for j in accs:
            if i == j:
                mat_sym.at[i, j] = 1.0
            else:
                mat_sym.at[i, j] = 0.5 * (mat.at[i, j] + mat.at[j, i])

    if output_csv:
        mat_sym.to_csv(output_csv)

    return mat_sym
def mat_sim(mat: dict | None, ref1: str, ref2: str) -> float:
    if not mat:
        return 0.0
    # mat could be dict-of-dict or pd.DataFrame.to_dict(); try both directions safely
    if ref1 in mat and ref2 in mat[ref1]:
        return float(mat[ref1][ref2])
    if ref2 in mat and ref1 in mat[ref2]:
        return float(mat[ref2][ref1])
    return 0.0
def infer_gt_from_readname(readname: str) -> str | None:
    parts = readname.split("_")
    if len(parts) < 3:
        return None
    return "_".join(parts[:-2])  # your existing convention

def resolve_conflict_groups_gt_biased(
    conflict_groups,
    bam_fs,
    reflengths,
    mat: dict | None = None,
    min_mat_sim: float = 0.0,
    remove_mode: str = "all_wrong",  # or "random_wrong"
    random_seed: int | None = None,
):
    """
    For each conflict group (regions across refs), pull reads from those regions.
    If GT exists, mark read as removable ONLY on refs != GT (and preferably only if mat says they’re similar).
    Returns removed_read_ids: dict read_id -> [refs_to_remove_from]
    """
    if random_seed is not None:
        random.seed(random_seed)

    removed_read_ids = defaultdict(list)

    for group in conflict_groups:
        # group: {(ref,s,e), ...}
        by_ref = defaultdict(list)
        for (ref, s, e) in group:
            by_ref[ref].append((s, e))

        # Gather reads per ref from the group’s regions
        reads_per_ref = {}
        for ref, segs in by_ref.items():
            uniq = {}
            for (s, e) in segs:
                for read in bam_fs.fetch(ref, s, e):
                    if read.is_unmapped:
                        continue
                    if read.is_secondary or read.is_supplementary:
                        continue
                    uniq[read.query_name] = True
            reads_per_ref[ref] = list(uniq.keys())

        # For each ref, decide which reads to remove from that ref
        refs = list(reads_per_ref.keys())
        for ref in refs:
            for rid in reads_per_ref[ref]:
                gt = infer_gt_from_readname(rid)
                if gt is None:
                    continue  # no GT => skip here; you can add MAPQ/AS fallback later

                # Only remove from the "wrong" ref, and only if ref is similar to gt (optional)
                if gt != ref:
                    if (not mat) or (mat_sim(mat, gt, ref) >= min_mat_sim):
                        removed_read_ids[rid].append(ref)

    return removed_read_ids
def overlap_len(a0, a1, b0, b1) -> int:
    return max(0, min(a1, b1) - max(a0, b0))

def region_ambiguity_hits(
    chrom: str,
    start: int,
    end: int,
    shared_idx: dict[str, list[tuple[int,int,str,int,int,float]]],
) -> list[dict]:
    """
    Returns a list of hits, each hit includes:
      other_chrom, jaccard, overlap_bp, overlap_frac_of_region
    """
    hits = []
    if chrom not in shared_idx:
        return hits

    region_len = max(1, end - start)

    # shared_idx[chrom] sorted by window start; cheap scan.
    # (If huge, swap to interval tree later.)
    for (ws, we, other_chr, os, oe, j) in shared_idx[chrom]:
        if we <= start:
            continue
        if ws >= end:
            break
        ol = overlap_len(start, end, ws, we)
        if ol > 0:
            hits.append({
                "other_chrom": other_chr,
                "jaccard": j,
                "overlap_bp": ol,
                "overlap_frac": ol / region_len,
                "win_start": ws,
                "win_end": we,
                "other_start": os,
                "other_end": oe,
            })
    return hits

def ambiguity_score_from_hits(hits: list[dict]) -> float:
    """
    Turn overlap hits into a single [0..1+] score.
    This version: sum(overlap_frac * jaccard), capped at 1.0.
    """
    if not hits:
        return 0.0
    s = sum(h["overlap_frac"] * h["jaccard"] for h in hits)
    return min(1.0, float(s))

from dataclasses import dataclass

@dataclass(frozen=True)
class SharedWindow:
    # window on "this" contig
    contig: str
    start: int
    end: int
    # corresponding window on the alternative contig
    alt_contig: str
    alt_start: int
    alt_end: int
    jaccard: float


def load_shared_windows_csv(
    csv_path: str,
    *,
    min_jaccard: float = 1.0,
    skip_same_contig: bool = True
) -> Dict[str, List[SharedWindow]]:
    """
    Load shared window pairs and index them by contig.

    Returns:
      index[contig] = list of SharedWindow objects sorted by start
    Includes BOTH directions:
      query->match and match->query
    """
    idx: Dict[str, List[SharedWindow]] = {}

    with open(csv_path, newline="") as fp:
        r = csv.DictReader(fp)
        required = {"query_contig", "q_start", "q_end", "match_contig", "m_start", "m_end", "jaccard"}
        missing = required - set(r.fieldnames or [])
        if missing:
            raise ValueError(f"CSV missing required columns: {missing}")

        for row in r:
            try:
                j = float(row["jaccard"])
            except Exception:
                continue
            if j < min_jaccard:
                continue

            qc = row["query_contig"]
            mc = row["match_contig"]
            if skip_same_contig and qc == mc:
                continue

            qs, qe = int(row["q_start"]), int(row["q_end"])
            ms, me = int(row["m_start"]), int(row["m_end"])

            # store forward
            idx.setdefault(qc, []).append(SharedWindow(qc, qs, qe, mc, ms, me, j))
            # store reverse too (so we can query by either contig)
            idx.setdefault(mc, []).append(SharedWindow(mc, ms, me, qc, qs, qe, j))

    # sort each contig list by start for interval querying
    for contig in idx:
        idx[contig].sort(key=lambda w: (w.start, w.end, w.alt_contig))
    return idx


def _windows_overlapping(
    windows_sorted: List[SharedWindow],
    pos_start: int,
    pos_end: int
) -> List[SharedWindow]:
    """
    Return windows that overlap [pos_start, pos_end) on a contig.
    windows_sorted must be sorted by start.
    """
    if not windows_sorted:
        return []

    starts = [w.start for w in windows_sorted]
    # Candidate windows start before pos_end
    i = bisect.bisect_left(starts, pos_end)

    hits = []
    # scan backwards a bit (typical window sizes make this cheap)
    # stop when window.end <= pos_start and since we're scanning backward by start, it will only get worse
    for w in reversed(windows_sorted[:i]):
        if w.end <= pos_start:
            break
        if w.start < pos_end and w.end > pos_start:
            hits.append(w)
    return hits


def mark_reads_in_shared_regions(
    bam_path: str,
    shared_idx: Dict[str, List[SharedWindow]],
    *,
    only_primary: bool = True,
) -> List[dict]:
    """
    Scan BAM alignments. If an alignment overlaps any SharedWindow on its contig,
    emit a record with read_id and alternative contig info.

    Returns: list of dicts with:
      read_id, aligned_contig, aln_start, aln_end,
      shared_contig, shared_start, shared_end,
      alt_contig, alt_start, alt_end, jaccard
    """
    out: List[dict] = []
    bam = pysam.AlignmentFile(bam_path, "rb")

    for read in bam:
        if read.is_unmapped:
            continue
        if only_primary and (read.is_secondary or read.is_supplementary):
            continue
        if read.reference_name is None or read.reference_start is None or read.reference_end is None:
            continue

        contig = read.reference_name
        if contig not in shared_idx:
            continue

        aln_s = int(read.reference_start)
        aln_e = int(read.reference_end)

        overlaps = _windows_overlapping(shared_idx[contig], aln_s, aln_e)
        if not overlaps:
            continue

        for w in overlaps:
            out.append({
                "read_id": read.query_name,
                "aligned_contig": contig,
                "aln_start": aln_s,
                "aln_end": aln_e,
                "shared_contig": w.contig,
                "shared_start": w.start,
                "shared_end": w.end,
                "alt_contig": w.alt_contig,
                "alt_start": w.alt_start,
                "alt_end": w.alt_end,
                "jaccard": w.jaccard,
            })

    bam.close()
    return out


def write_marked_reads_tsv(marked: List[dict], out_tsv: str) -> None:
    if not marked:
        # still write header
        header = [
            "read_id","aligned_contig","aln_start","aln_end",
            "shared_contig","shared_start","shared_end",
            "alt_contig","alt_start","alt_end","jaccard"
        ]
        with open(out_tsv, "w", newline="") as fp:
            fp.write("\t".join(header) + "\n")
        return

    header = list(marked[0].keys())
    with open(out_tsv, "w", newline="") as fp:
        fp.write("\t".join(header) + "\n")
        for row in marked:
            fp.write("\t".join(str(row[h]) for h in header) + "\n")

def annotate_regions_with_ambiguity(merged_regions_df, shared_idx):
    merged = merged_regions_df.copy()
    scores = []
    partners = []
    for r in merged.itertuples(index=False):
        hits = region_ambiguity_hits(str(r.chrom), int(r.start), int(r.end), shared_idx)
        scores.append(ambiguity_score_from_hits(hits))
        partners.append(",".join(sorted(set(h["other_chrom"] for h in hits))) if hits else "")
    merged["ambiguity_score"] = scores
    merged["ambiguous_partners"] = partners
    return merged


def compute_remove_counts_from_ambiguity(
    merged_regions_df,
    reads_map,
    *,
    alpha: float = 1.0,           # 1.0 = remove up to N*A reads
    min_remove: int = 0,
    max_remove_frac: float = 0.95, # never remove >95% from a region
):
    """
    Returns remove_counts: dict[(chrom,start,end)] -> remove_n
    """
    remove_counts = {}
    for r in merged_regions_df.itertuples(index=False):
        chrom = str(r.chrom); start = int(r.start); end = int(r.end)
        A = float(getattr(r, "ambiguity_score", 0.0))
        if A <= 0:
            continue

        region_reads = fetch_reads_in_region(reads_map, chrom, start, end)
        N = len(region_reads)
        if N == 0:
            continue

        remove_n = int(math.floor(N * A * alpha))
        remove_n = max(min_remove, remove_n)
        remove_n = min(remove_n, int(math.floor(N * max_remove_frac)))

        if remove_n > 0:
            remove_counts[(chrom, start, end)] = remove_n

    return remove_counts

def load_ambiguous_windows(report_csv: str, jaccard_min: float = 0.75):
    df = pd.read_csv(report_csv)

    # force numeric
    df["jaccard"] = pd.to_numeric(df["jaccard"], errors="coerce")
    df = df.dropna(subset=["jaccard"])
    df = df[df["jaccard"] >= jaccard_min].copy()

    ambig = defaultdict(list)

    # Store BOTH directions so lookup works no matter which chrom is the region chrom
    for r in df.itertuples(index=False):
        qchrom, qs, qe = str(r.query_contig), int(r.q_start), int(r.q_end)
        mchrm, ms, me = str(r.match_contig), int(r.m_start), int(r.m_end)
        j = float(r.jaccard)

        ambig[qchrom].append((qs, qe, mchrm, ms, me, j))
        ambig[mchrm].append((ms, me, qchrom, qs, qe, j))

    # sort by start for each contig
    for chrom in list(ambig.keys()):
        ambig[chrom].sort(key=lambda x: x[0])

    return ambig

def overlaps(a_start, a_end, b_start, b_end):
    return (a_end > b_start) and (a_start < b_end)

def find_overlapping_ambig(ambig_list, start, end):
    """
    ambig_list is sorted by interval start:
      [(s,e, pchr, ps, pe, j), ...]
    Return subset that overlaps [start,end)
    """
    if not ambig_list:
        return []

    starts = [x[0] for x in ambig_list]
    # find candidate insertion point near region start
    i = bisect_left(starts, start)
    i = max(0, i - 1)

    hits = []
    for k in range(i, len(ambig_list)):
        s, e, pchr, ps, pe, j = ambig_list[k]
        if s >= end:
            break
        if overlaps(s, e, start, end):
            hits.append((s, e, pchr, ps, pe, j))
    return hits


def flag_reads_in_ambiguous_ranges(
    bam_path: str,
    merged_regions_df: pd.DataFrame,   # columns: chrom,start,end
    ambig_map: dict,                   # output of load_ambiguous_windows
    *,
    require_primary_only: bool = True,
):
    bam = pysam.AlignmentFile(bam_path, "rb")
    flagged = defaultdict(list)

    for row in merged_regions_df.itertuples(index=False):
        chrom = str(row.chrom)
        start = int(row.start)
        end   = int(row.end)

        ambig_list = ambig_map.get(chrom)
        if not ambig_list:
            continue

        hits = find_overlapping_ambig(ambig_list, start, end)
        if not hits:
            continue

        for (as_, ae, pchr, ps, pe, j) in hits:
            # intersection sub-range on this chrom
            ovl_s = max(start, as_)
            ovl_e = min(end, ae)
            if ovl_s >= ovl_e:
                continue

            for read in bam.fetch(chrom, ovl_s, ovl_e):
                if read.is_unmapped:
                    continue
                if require_primary_only and (read.is_secondary or read.is_supplementary):
                    continue

                flagged[read.query_name].append({
                    "on_chrom": chrom,
                    "read_ref_start": read.reference_start,
                    "read_ref_end": read.reference_end,
                    "ambig_span_start": ovl_s,
                    "ambig_span_end": ovl_e,
                    "partner_chrom": pchr,
                    "partner_span_start": ps,
                    "partner_span_end": pe,
                    "jaccard": j,
                })

    bam.close()
    return flagged
from bisect import bisect_left
def alignment_alt_contigs(shared_idx, contig: str, aln_s: int, aln_e: int):
    """
    Returns (alt_set, overlap_bp_sum) for this alignment on contig.
    alt_set: set of alternative contigs that share any window overlap with [aln_s, aln_e)
    overlap_bp_sum: total overlapped bp across windows (double-counting possible if windows overlap)
    """
    if contig not in shared_idx:
        return set(), 0

    overlaps = _windows_overlapping(shared_idx[contig], aln_s, aln_e)
    if not overlaps:
        return set(), 0

    alt_set = set()
    overlap_bp_sum = 0
    for w in overlaps:
        # overlap on the contig coordinates
        ovl_s = max(aln_s, w.start)
        ovl_e = min(aln_e, w.end)
        if ovl_s < ovl_e:
            alt_set.add(w.alt_contig)
            overlap_bp_sum += (ovl_e - ovl_s)

    # Don't count self as alternative even if CSV had it
    alt_set.discard(contig)
    return alt_set, overlap_bp_sum
from collections import defaultdict

def build_removed_ids_best_alignment(
    bam_path: str,
    shared_idx: Dict[str, List[SharedWindow]],
    *,
    penalize_weight: float = 50.0,    # how harsh ambiguity is (MAPQ scale ~0-60)
    as_weight: float = 0.0,           # set >0 if AS is present and useful
    drop_contigs: set[str] | None = None,   # contigs you *know* are wrong
    drop_if_ambiguous: bool = True,         # for drop_contigs: drop if ambiguous
    min_alt_count: int = 1,                 # ambiguous if >= this many alt contigs
    only_primary: bool = False,             # set True if you truly only want primary
) -> Dict[str, List[str]]:
    """
    Returns removed_read_ids[read_id] = [refs_to_remove_from]

    Logic:
      - For each read, score each alignment.
      - Keep best scoring alignment(s).
      - Mark all others for removal.
      - Additionally, if drop_contigs is provided, and alignment is on those contigs AND ambiguous,
        mark for removal regardless of score (aggressive targeting).
    """
    if drop_contigs is None:
        drop_contigs = set()

    # Collect alignments per read
    per_read = defaultdict(list)
    bam = pysam.AlignmentFile(bam_path, "rb")

    for r in bam:
        if r.is_unmapped:
            continue
        if only_primary and (r.is_secondary or r.is_supplementary):
            continue
        if r.reference_name is None or r.reference_start is None or r.reference_end is None:
            continue

        contig = r.reference_name
        s = int(r.reference_start)
        e = int(r.reference_end)

        alt_set, _ovl_bp = alignment_alt_contigs(shared_idx, contig, s, e)
        alt_n = len(alt_set)

        # Pull AS if present
        AS = 0
        try:
            AS = int(r.get_tag("AS"))
        except Exception:
            AS = 0

        mapq = int(r.mapping_quality)

        # higher is better
        score = mapq + as_weight * AS - penalize_weight * alt_n

        per_read[r.query_name].append({
            "contig": contig,
            "score": score,
            "mapq": mapq,
            "AS": AS,
            "alt_n": alt_n,
            "is_secondary": bool(r.is_secondary or r.is_supplementary),
        })

    bam.close()

    removed = defaultdict(list)

    # Decide best alignment(s) per read
    for rid, alns in per_read.items():
        if not alns:
            continue

        best = max(a["score"] for a in alns)

        # Keep everything tied for best (rare, but possible)
        keep_contigs = {a["contig"] for a in alns if a["score"] == best}

        for a in alns:
            c = a["contig"]
            ambiguous = (a["alt_n"] >= min_alt_count)

            # Aggressive rule: if this contig is known-wrong and ambiguous, drop it no matter what
            if (c in drop_contigs) and drop_if_ambiguous and ambiguous:
                removed[rid].append(c)
                continue

            # Otherwise keep only best contig(s)
            if c not in keep_contigs:
                removed[rid].append(c)

    # de-dupe
    removed = {rid: sorted(set(refs)) for rid, refs in removed.items() if refs}
    return removed
def report_removed_read_stats(
    bam_path: str,
    removed_read_ids: dict[str, list[str]],
):
    """
    Prints:
      - % of unique reads affected
      - % of alignments removed
      - per-contig removal percentages
    """
    bam = pysam.AlignmentFile(bam_path, "rb")

    total_reads = set()
    total_alignments = 0

    removed_reads = set(removed_read_ids.keys())
    removed_alignments = 0

    per_contig_total = defaultdict(int)
    per_contig_removed = defaultdict(int)

    for r in bam:
        if r.is_unmapped:
            continue
        total_alignments += 1
        total_reads.add(r.query_name)

        contig = r.reference_name
        per_contig_total[contig] += 1

        if r.query_name in removed_read_ids:
            if contig in removed_read_ids[r.query_name]:
                removed_alignments += 1
                per_contig_removed[contig] += 1

    bam.close()

    n_reads = len(total_reads)
    n_removed_reads = len(removed_reads)

    print("\n=== Removal Summary ===")
    print(f"Total unique reads:        {n_reads}")
    print(f"Reads with ≥1 removal:     {n_removed_reads} "
          f"({100.0 * n_removed_reads / max(1, n_reads):.2f}%)")

    print(f"\nTotal alignments:          {total_alignments}")
    print(f"Alignments removed:        {removed_alignments} "
          f"({100.0 * removed_alignments / max(1, total_alignments):.2f}%)")

    print("\nPer-contig alignment removal:")
    for contig in sorted(per_contig_total):
        tot = per_contig_total[contig]
        rem = per_contig_removed.get(contig, 0)
        if tot == 0:
            continue
        pct = 100.0 * rem / tot
        print(f"  {contig:20s}  {rem:8d}/{tot:8d}  ({pct:6.2f}%)")

    return {
        "total_reads": n_reads,
        "removed_reads": n_removed_reads,
        "pct_reads_removed": 100.0 * n_removed_reads / max(1, n_reads),
        "total_alignments": total_alignments,
        "removed_alignments": removed_alignments,
        "pct_alignments_removed": 100.0 * removed_alignments / max(1, total_alignments),
    }

def build_read_to_refs(bam_path: str, include_secondary=True):
    bam = pysam.AlignmentFile(bam_path, "rb")
    read2refs = defaultdict(set)
    for r in bam:
        if r.is_unmapped:
            continue
        if (not include_secondary) and (r.is_secondary or r.is_supplementary):
            continue
        read2refs[r.query_name].add(r.reference_name)
    bam.close()
    return read2refs

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
        shared_windows_report_csv = None,
        compare_to_reference_windows = False

):
    os.makedirs(output_dir, exist_ok=True)
    import time
    print(f"Starting conflict region detection at {time.ctime()}")
    print(fasta_files)


    if len(fasta_files) == 0:
        print("No fasta files provided, using sensitive mode")
    elif sensitive:
        print("Fasta files provided but sensitive mode enabled, using BAM alignments for signature generations")
        fasta_files = []
    else:
        print(f"Using {len(fasta_files)} fasta files for regional comparisons")
    # Step 7: Parse filtered BED to get regions
    regions = parse_bed_file(bedfile)
    print(f"Total regions defined: {len(regions)} in default mode for signature generation")

    # 2. Choose statistic function
    use_jump = True
    use_variance=False
    jump_threshold = 1
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
        print("Creating signatures for regions")

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

    compare_to_reference_windows = True
    if compare_to_reference_windows:
        # Build FASTA window index
        # if not os.path.exists("shared_windows_report.csv"):
        report_shared_windows_across_fastas(
            fasta_files=fasta_files,
            output_csv="shared_windows_report.csv",
            ksize=21,
            scaled=8000,
            window=2000,
            step=2000,
            jaccard_threshold=0.80,
            max_hits_per_query=3,
            skip_self_same_fasta=False,
        )
        # import shared_windows_report.csv as a dictionary
        shared_idx = load_shared_windows_csv(
            "shared_windows_report.csv",
            min_jaccard=1,          # or 0.99 if you want near-identical too
            skip_same_contig=True
        )


        print("Reference SBT by contig:")



        marked = mark_reads_in_shared_regions(
            bam_path=input_bam,
            shared_idx=shared_idx,
            only_primary=True
        )

        write_marked_reads_tsv(marked, "reads_in_shared_regions.tsv")
        print("Flagged alignments:", len(marked))

        # creeate merged_regions from shared_idx
        merged_regions_list = []
        seen=dict()
        for chr, windows in shared_idx.items():
            for w in windows:
                key = (w.contig, w.start, w.end)
                if key in seen:
                    continue
                seen[key]=True
                merged_regions_list.append({
                    "chrom": w.contig,
                    "start": w.start,
                    "end": w.end,
                    "depth": 1
                })
                key_alt = (w.alt_contig, w.alt_start, w.alt_end)
                if key_alt in seen:
                    continue
                seen[key_alt]=True
                merged_regions_list.append({
                    "chrom": w.alt_contig,
                    "start": w.alt_start,
                    "end": w.alt_end,
                    "depth": 1
                })


        # convert to dataframe
        merged_regions_df = pd.DataFrame(merged_regions_list)

        signatures_df = create_signatures_for_regions(
            regions_df=merged_regions_df,
            bam_path = input_bam,
            fasta_paths = fasta_files,
            kmer_size=kmer_size,
            scaled=scaled,
            num_workers=num_workers
        )
        # merge the two signature dicts
        signatures.update(signatures_df)

    if matrix:
        dist_dict = import_matrix_ani(matrix)
        clusters = find_reference_clusters(dist_dict, min_sim=min_similarity_comparable)
    elif reference_signatures:  # the CSV from `sourmash compare --csv`
        _, dist_dict = parse_sourmash_compare_csv(reference_signatures)
        # 3A. Build clusters (connected components)
        clusters = find_reference_clusters(dist_dict, min_sim=min_similarity_comparable)
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
    # removed_read_ids = finalize_proportional_removal(
    #     conflict_groups,
    #     bam_fs,
    #     fetch_reads_in_region,
    #     remove_mode='random'
    # )
    bad = {}

    removed_read_ids = build_removed_ids_best_alignment(
        bam_path=input_bam,
        shared_idx=shared_idx,
        penalize_weight=1.0,     # start here; raise to be more aggressive
        as_weight=0.0,            # set 1.0 if you trust AS and it exists
        drop_contigs=bad,
        drop_if_ambiguous=True,
        min_alt_count=1,          # ambiguous if ANY alternative exists
        only_primary=False        # set True if you only have primaries anyway
    )
    stats = report_removed_read_stats(
        bam_path=input_bam,
        removed_read_ids=removed_read_ids
    )

    start_time = time.time()
    comparison_df = None
    includable_read_ids = dict()
    sum_of_ref_aligned = defaultdict(int)
    output_file = os.path.join(output_dir, "failed_reads.txt")

    print("Writing failed reads to outputfile:", output_file)
    with open (output_file, 'w') as f:
        for read, refs in removed_read_ids.items():
            for ref in refs:
                f.write(f"{ref}\t{read}\n")
        f.close()

    try:
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

        return removed_read_ids, comparison_df
    else:
        return removed_read_ids, comparison_df
