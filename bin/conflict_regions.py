#!/usr/bin/env python3
"""
Conflict-region detection and read-removal planning.

This is a cleaned-up version of your script with:
- duplicate imports removed
- unused functions removed
- comments rewritten to be short and neutral
- formatting normalized
- only functions required by determine_conflicts() and its called paths kept

Notes:
- FASTA paths are optional; if provided (and not in sensitive mode), signatures are built from FASTA slices.
- compare_to_reference_windows path generates shared-window ambiguity via sourmash window sketches and uses an
  alignment-scoring heuristic to decide which alignments to drop per read.
"""

from __future__ import annotations

import bisect
import csv
import itertools
import re
import math
import os
import random
import shutil
import statistics
import tempfile
import time
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import pysam
from tqdm import tqdm

from sourmash import MinHash, SourmashSignature, save_signatures, load_file_as_signatures
from sourmash.sbt import SBT, GraphFactory
from sourmash.sbtmh import SigLeaf
from utils import load_matchfile

# -----------------------------
# IO / parsing helpers
# -----------------------------

def parse_bed_file(bed_file_path: str) -> pd.DataFrame:
    """Read bedgraph-like file: chrom, start, end, depth."""
    cols = ["chrom", "start", "end", "depth"]
    return pd.read_csv(
        bed_file_path,
        sep="\t",
        header=None,
        names=cols,
        dtype={
            "chrom": "category",
            "start": "int32",
            "end": "int32",
            "depth": "float32",
        },
    )


def parse_split(region_name: str) -> Tuple[str, int, int]:
    """
    Region string format: ref:start-end
    start/end may be 0-based or 1-based depending on your upstream; this preserves numeric values as-is.
    """
    ref, coords = region_name.split(":", 1)
    s, e = coords.split("-", 1)
    return ref, int(s), int(e)


def save_signatures_sourmash(signatures: Dict[str, SourmashSignature], output_sigfile: str) -> None:
    """Save a dict of signatures to a single .sig file."""
    os.makedirs(os.path.dirname(output_sigfile), exist_ok=True)
    sigs = []
    for name, sig in signatures.items():
        sig._name = name
        sigs.append(sig)
    with open(output_sigfile, "wt") as fp:
        save_signatures(sigs, fp)


def load_signatures_sourmash(input_sigfile: str):
    """Load multiple signatures from a single .sig file."""
    return load_file_as_signatures(input_sigfile)


def rebuild_sig_dict(loaded_sigs) -> Dict[str, SourmashSignature]:
    """Convert loaded signatures iterator/list into {sig.name: sig}."""
    out: Dict[str, SourmashSignature] = {}
    for sig in loaded_sigs:
        out[sig.name] = sig
    return out


# -----------------------------
# Region merging (bedgraph -> merged regions)
# -----------------------------

def compute_gini(values: List[float]) -> float:
    """Gini coefficient; 0 = equal, higher = more unequal."""
    if not values:
        return 0.0
    vals = sorted(values)
    n = len(vals)
    cum = 0.0
    cumvals = 0.0
    for i, v in enumerate(vals, 1):
        cum += v
        cumvals += cum
    if cum == 0:
        return 0.0
    return (2 * cumvals) / (n * cum) - (n + 1) / n


def merge_bedgraph_regions(
    intervals: pd.DataFrame,
    merging_method: str = "jump",
    max_stat_threshold: Optional[float] = None,
    max_group_size: int = 20000,
    max_length: Optional[int] = None,
    value_diff_tolerance: Optional[float] = None,
    breadth_allowance: int = 1000,
    gap_allowance: float = 0.1,
    reflengths: Optional[Dict[str, int]] = None,
) -> pd.DataFrame:
    """
    Merge adjacent intervals within each chrom using a rule:
      - jump: merge if abs(depth - last_depth) <= threshold and gap <= allowed_gap
      - variance/gini: compute stat on buffered depths and merge if <= threshold

    Optimised: uses numpy arrays per-chrom instead of itertuples, and
    tracks running variance via Welford's algorithm (avoids O(n^2) pvariance calls).
    """
    if intervals.empty:
        return pd.DataFrame(columns=["chrom", "start", "end", "mean_depth"])

    df = intervals.sort_values(["chrom", "start"]).reset_index(drop=True).copy()
    df["depth"] = df["depth"].fillna(0)

    # allowed gap per chrom
    if reflengths is None:
        allowed_gap = {c: float("inf") for c in df["chrom"].astype(str).unique()}
    else:
        allowed_gap = {}
        for chrom in df["chrom"].astype(str).unique():
            L = reflengths.get(chrom)
            if L is None:
                allowed_gap[chrom] = breadth_allowance
            else:
                allowed_gap[chrom] = max(breadth_allowance, int(L * gap_allowance))

    # pick default jump threshold if none passed
    if merging_method == "jump":
        if max_stat_threshold is not None:
            jump_thr = float(max_stat_threshold)
        else:
            try:
                q = df.groupby("chrom", observed=True)["depth"].quantile(0.975).mean()
                jump_thr = float(math.ceil(q) + 1)
            except Exception:
                jump_thr = 200.0
    else:
        jump_thr = None

    merged = []
    use_jump = merging_method == "jump"
    use_variance = merging_method == "variance"
    use_gini = merging_method == "gini"

    for chrom, g in df.groupby("chrom", observed=True):
        if g.empty:
            continue

        chrom_str = str(chrom)
        gap_thr = allowed_gap.get(chrom_str, float("inf"))

        # Extract numpy arrays for fast iteration (avoids itertuples overhead)
        starts = g["start"].values   # int32 ndarray
        ends = g["end"].values
        depths = g["depth"].values.astype(np.float64)
        n = len(starts)

        buf_start = int(starts[0])
        buf_end = int(ends[0])
        last_depth = depths[0]
        buf_count = 1
        buf_sum = last_depth

        # Welford's running variance state (for variance method)
        _wf_mean = last_depth
        _wf_m2 = 0.0

        # For gini we still need the full list (rare path, usually not perf-critical)
        buf_depths_list = [last_depth] if use_gini else None

        for i in range(1, n):
            new_start = int(starts[i])
            new_end = int(ends[i])
            new_depth = depths[i]
            gap = new_start - buf_end

            can_merge = True
            if buf_count >= max_group_size:
                can_merge = False
            if can_merge and max_length is not None and (new_end - buf_start) > max_length:
                can_merge = False
            if can_merge and gap > gap_thr:
                can_merge = False

            if can_merge:
                if use_jump:
                    if abs(new_depth - last_depth) > jump_thr:
                        can_merge = False
                elif use_variance:
                    # Welford's online variance: O(1) per step instead of O(n)
                    new_count = buf_count + 1
                    delta = new_depth - _wf_mean
                    new_mean = _wf_mean + delta / new_count
                    delta2 = new_depth - new_mean
                    new_m2 = _wf_m2 + delta * delta2
                    pvar = new_m2 / new_count
                    if max_stat_threshold is not None and pvar > max_stat_threshold:
                        can_merge = False
                    if can_merge and value_diff_tolerance is not None:
                        if abs(new_depth - (buf_sum / buf_count)) > value_diff_tolerance:
                            can_merge = False
                elif use_gini:
                    vals = buf_depths_list + [new_depth]
                    stat = compute_gini(vals)
                    if max_stat_threshold is not None and stat > max_stat_threshold:
                        can_merge = False
                    if can_merge and value_diff_tolerance is not None:
                        if abs(new_depth - (buf_sum / buf_count)) > value_diff_tolerance:
                            can_merge = False
                else:
                    raise ValueError(f"Unknown merging method: {merging_method}")

            if can_merge:
                buf_end = new_end
                last_depth = new_depth
                buf_count += 1
                buf_sum += new_depth
                # Update Welford state
                if use_variance:
                    _wf_mean = new_mean
                    _wf_m2 = new_m2
                if use_gini:
                    buf_depths_list.append(new_depth)
            else:
                merged.append(
                    {
                        "chrom": chrom_str,
                        "start": buf_start,
                        "end": buf_end,
                        "mean_depth": buf_sum / buf_count,
                    }
                )
                buf_start = new_start
                buf_end = new_end
                last_depth = new_depth
                buf_count = 1
                buf_sum = new_depth
                _wf_mean = new_depth
                _wf_m2 = 0.0
                if use_gini:
                    buf_depths_list = [new_depth]

        merged.append(
            {
                "chrom": chrom_str,
                "start": buf_start,
                "end": buf_end,
                "mean_depth": buf_sum / buf_count,
            }
        )

    return pd.DataFrame(merged)


# -----------------------------
# Signature generation for merged regions
# -----------------------------

def create_signature_for_single_region(
    refname: str,
    start: int,
    end: int,
    kmer_size: int,
    scaled: int,
    seq: str,
) -> Tuple[str, SourmashSignature]:
    mh = MinHash(n=0, ksize=kmer_size, scaled=scaled)
    mh.add_sequence(seq, force=True)
    region_name = f"{refname}:{start}-{end}"
    sig = SourmashSignature(mh, name=region_name)
    return region_name, sig



# -----------------------------
# BAM-only region sketching (no FASTA) — single-pass per chromosome
# -----------------------------

def _sketch_bam_single_pass_chrom(args: Tuple) -> List[Tuple]:
    """Worker: stream all reads for *one chromosome* and assign to regions.

    Instead of calling bam.fetch(chrom, start, end) per region (N random I/O
    lookups), we call bam.fetch(chrom) once and use binary search on the
    sorted, non-overlapping region list to find which region(s) each read
    overlaps.  This turns millions of index lookups into one sequential scan.

    Returns list of (region_name, SourmashSignature) pairs.
    """
    chrom, sorted_regions, bam_path, kmer_size, scaled = args
    # sorted_regions: list of (start, end) sorted by start — non-overlapping
    n_regions = len(sorted_regions)
    region_mhs = [MinHash(n=0, ksize=kmer_size, scaled=scaled) for _ in range(n_regions)]

    # Pre-extract ends for bisect lookups
    ends = [r[1] for r in sorted_regions]

    bam = pysam.AlignmentFile(bam_path, "rb")
    try:
        for read in bam.fetch(chrom):
            seq = read.query_sequence
            if not seq or len(seq) < kmer_size:
                continue
            r_start = read.reference_start
            r_end = read.reference_end
            if r_end is None:
                r_end = r_start + len(seq)

            # First region whose end > r_start (could overlap the read)
            lo = bisect.bisect_right(ends, r_start)
            for i in range(lo, n_regions):
                s, e = sorted_regions[i]
                if s >= r_end:
                    break  # all remaining regions are past the read
                region_mhs[i].add_sequence(seq, force=True)
    except ValueError:
        pass
    finally:
        bam.close()

    results: List[Tuple] = []
    for i, (s, e) in enumerate(sorted_regions):
        mh = region_mhs[i]
        if mh.hashes:
            name = f"{chrom}:{s}-{e}"
            results.append((name, SourmashSignature(mh, name=name)))
    return results


def create_signatures_from_bam(
    regions_df: pd.DataFrame,
    bam_path: str,
    kmer_size: int,
    scaled: int,
    num_workers: int = 1,
) -> Dict[str, SourmashSignature]:
    """Build signatures for merged regions using BAM read sequences only (no FASTA).

    Optimised: one sequential BAM scan per chromosome with binary-search
    assignment to regions, parallelised across chromosomes.
    """
    required = {"chrom", "start", "end"}
    if not required.issubset(set(regions_df.columns)):
        raise ValueError(f"regions_df must contain columns: {required}")

    df = regions_df[["chrom", "start", "end"]].copy()
    df["chrom"] = df["chrom"].astype(str)
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)

    # Group by chromosome; sort regions within each group
    grouped = df.groupby("chrom", observed=True)
    chrom_args: List[Tuple] = []
    for chrom, grp in grouped:
        sorted_regions = sorted(zip(grp["start"], grp["end"]))
        chrom_args.append((str(chrom), sorted_regions, bam_path, kmer_size, scaled))

    n_chroms = len(chrom_args)
    n_rows = len(df)
    sigs: Dict[str, SourmashSignature] = {}

    if num_workers <= 1 or n_chroms <= 1:
        for ca in tqdm(chrom_args, desc="Sketching chroms (BAM single-pass)", unit="chrom"):
            for name, sig in _sketch_bam_single_pass_chrom(ca):
                sigs[name] = sig
    else:
        actual_workers = min(num_workers, n_chroms)
        print(f"  Sketching {n_rows} regions across {n_chroms} chrom(s) from BAM "
              f"(single-pass, {actual_workers} workers)")
        with ProcessPoolExecutor(max_workers=actual_workers) as pool:
            futures = [pool.submit(_sketch_bam_single_pass_chrom, ca) for ca in chrom_args]
            for fut in tqdm(as_completed(futures), total=len(futures),
                            desc="Sketching chroms (BAM single-pass)", unit="chrom"):
                for name, sig in fut.result():
                    sigs[name] = sig

    print(f"  Built {len(sigs)} region signature(s) from BAM reads "
          f"({n_rows} input rows, {n_chroms} chrom(s))")
    return sigs


# -----------------------------
# FASTA-based region sketching
# -----------------------------

def _sketch_region_batch(args: Tuple) -> List[Tuple]:
    """Worker: fetch sequences from FASTA files and sketch a batch of regions.

    Opens each FASTA file only once per batch.  Returns a list of
    (region_name, hashes_tuple) so that the main process can reconstruct
    SourmashSignature objects without sending unpicklable file handles.
    """
    batch, chrom_to_fasta_path, kmer_size, scaled = args
    open_fastas: Dict[str, pysam.FastaFile] = {}
    results: List[Tuple] = []
    try:
        for chrom, start, end in batch:
            fp = chrom_to_fasta_path.get(chrom)
            if fp is None:
                continue
            if fp not in open_fastas:
                open_fastas[fp] = pysam.FastaFile(fp)
            fa = open_fastas[fp]
            seq = fa.fetch(chrom, start, end)
            if not seq or len(seq) < kmer_size:
                continue
            mh = MinHash(n=0, ksize=kmer_size, scaled=scaled)
            mh.add_sequence(seq, force=True)
            if mh.hashes:
                results.append((f"{chrom}:{start}-{end}", tuple(mh.hashes)))
    finally:
        for fa in open_fastas.values():
            fa.close()
    return results


def create_signatures_for_regions(
    regions_df: pd.DataFrame,
    bam_path: str,
    fasta_paths: List[str],
    kmer_size: int,
    scaled: int,
    num_workers: int = 1,
) -> Dict[str, SourmashSignature]:
    """Build signatures for each merged region.

    Uses FASTA files for sequence retrieval.  Workers are parallelised across
    batches of rows so that ``num_workers`` is actually honoured (previously
    the function was entirely serial regardless of this parameter).

    A pre-built ``chrom → fasta_path`` dict replaces the O(N) linear scan of
    ``fa.references`` that was previously done per-region, which caused
    progressive slowdown as the reference list grew.
    """
    required = {"chrom", "start", "end"}
    if not required.issubset(set(regions_df.columns)):
        raise ValueError(f"regions_df must contain columns: {required}")

    if not fasta_paths:
        print("WARNING: No FASTA files provided; create_signatures_for_regions "
              "cannot sketch regions without sequence. Returning empty dict.")
        return {}

    # ── pre-build O(1) chrom→fasta_path lookup ───────────────────────────
    # Previously: `chrom in fa.references` was an O(len(references)) linear
    # scan done for every region.  With 2.8M regions and 200+ references that
    # adds up to hundreds of millions of string comparisons and causes the
    # progressive slowdown seen toward the end of a large run.
    chrom_to_fasta_path: Dict[str, str] = {}
    for fp in fasta_paths:
        fa = pysam.FastaFile(fp)
        for ref in fa.references:
            if ref not in chrom_to_fasta_path:
                chrom_to_fasta_path[ref] = fp
        fa.close()
    print(f"  Region sketching: {len(chrom_to_fasta_path)} contig(s) across "
          f"{len(fasta_paths)} FASTA file(s)")

    # ── build row list once ───────────────────────────────────────────────
    df = regions_df[["chrom", "start", "end"]].copy()
    df["chrom"] = df["chrom"].astype(str)
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    rows: List[Tuple] = list(df.itertuples(index=False, name=None))
    n_rows = len(rows)

    # ── choose batch size: aim for ~4× more chunks than workers ──────────
    # Small batches → good load balance; large batches → less IPC overhead.
    # 2000 rows/batch is a good default; scale up for very large datasets.
    batch_size = max(500, min(5000, n_rows // max(1, num_workers * 4)))
    batches = [rows[i:i + batch_size] for i in range(0, n_rows, batch_size)]
    args_list = [(b, chrom_to_fasta_path, kmer_size, scaled) for b in batches]

    sigs: Dict[str, SourmashSignature] = {}

    if num_workers <= 1:
        # Serial path — simple loop, no IPC overhead
        for batch_args in tqdm(args_list, desc="Sketching regions", unit="batch"):
            for region_name, hashes in _sketch_region_batch(batch_args):
                mh = MinHash(n=0, ksize=kmer_size, scaled=scaled)
                mh.add_many(hashes)
                sigs[region_name] = SourmashSignature(mh, name=region_name)
    else:
        # Parallel path: workers do FASTA fetch + MinHash, main process
        # reconstructs SourmashSignature from hashes tuples (which are
        # picklable, unlike FastaFile or MinHash objects).
        print(f"  Sketching {n_rows} regions in {len(batches)} batch(es) "
              f"across {num_workers} worker(s) (batch_size={batch_size})")
        with ProcessPoolExecutor(max_workers=num_workers) as pool:
            futures = [pool.submit(_sketch_region_batch, a) for a in args_list]
            for fut in tqdm(as_completed(futures), total=len(futures),
                            desc="Sketching regions", unit="batch"):
                for region_name, hashes in fut.result():
                    mh = MinHash(n=0, ksize=kmer_size, scaled=scaled)
                    mh.add_many(hashes)
                    sigs[region_name] = SourmashSignature(mh, name=region_name)

    print(f"  Built {len(sigs)} region signature(s) from {n_rows} input rows")
    return sigs


# -----------------------------
# SBT clustering + search
# -----------------------------

def build_sbt_index(
    siglist: List[Tuple[str, SourmashSignature]],
    ksize: int = 31,
    clusters: Optional[List[set]] = None,
) -> Dict[str, SBT]:
    """
    Build an SBT per cluster label. Cluster label is derived from clusters list membership.
    """
    factory = GraphFactory(ksize=ksize, n_tables=1, starting_size=max(1, len(siglist)))
    ref2cluster: Dict[str, str] = {}

    if clusters:
        for idx, cset in enumerate(clusters):
            for ref in cset:
                ref2cluster[ref] = str(idx)

    sbts: Dict[str, SBT] = {}

    for region_name, sig in tqdm(siglist, miniters=1000, desc="Indexing signatures"):
        ref, _, _ = parse_split(region_name)
        cluster_id = ref2cluster.get(ref, "Unknown")
        if cluster_id not in sbts:
            sbts[cluster_id] = SBT(factory)
        sbts[cluster_id].add_node(SigLeaf(region_name, sig))

    return sbts


# ── Parallel SBT search ──────────────────────────────────────────────────────
# Build SBT once (fast), save to a temp directory, then let N workers each
# load their own copy from disk and search a chunk of queries in parallel.
# This turns the sequential O(Q * search_cost) into O(Q * search_cost / N).

_SBT_WORKER_TREE: Optional[SBT] = None
_SBT_WORKER_KSIZE: int = 31
_SBT_WORKER_SCALED: int = 100


def _init_sbt_parallel_worker(sbt_path: str, ksize: int, scaled: int) -> None:
    """Worker initializer: load SBT from disk once per worker process."""
    global _SBT_WORKER_TREE, _SBT_WORKER_KSIZE, _SBT_WORKER_SCALED
    _SBT_WORKER_TREE = SBT.load(sbt_path, leaf_loader=SigLeaf.load)
    _SBT_WORKER_KSIZE = ksize
    _SBT_WORKER_SCALED = scaled


def _sbt_search_chunk(args: Tuple) -> List[Tuple]:
    """Worker: search a chunk of query signatures against the pre-loaded SBT.

    Reconstructs MinHash objects from serialised hash tuples so that nothing
    unpicklable crosses the process boundary.
    """
    query_items, min_threshold = args
    results: List[Tuple] = []
    for region_name, q_hashes in query_items:
        mh = MinHash(n=0, ksize=_SBT_WORKER_KSIZE, scaled=_SBT_WORKER_SCALED)
        mh.add_many(q_hashes)
        if not mh.hashes:
            continue
        sig = SourmashSignature(mh, name=region_name)
        q_ref, q_s, q_e = parse_split(region_name)

        for sr in _SBT_WORKER_TREE.search(sig, threshold=min_threshold, best_only=False):
            m_name = sr.signature.name
            m_ref, m_s, m_e = parse_split(m_name)
            if m_ref == q_ref:
                continue
            j = float(sr.score)
            # Compute containment from the underlying hashes
            m_mh = sr.signature.minhash
            n_shared = len(set(mh.hashes) & set(m_mh.hashes))
            c12 = (n_shared / len(mh.hashes)) if mh.hashes else 0.0
            c21 = (n_shared / len(m_mh.hashes)) if m_mh.hashes else 0.0
            results.append((q_ref, q_s, q_e, m_ref, m_s, m_e, j, c12, c21))
    return results


def search_sbt_parallel(
    sig_items: List[Tuple[str, SourmashSignature]],
    output_csv: str,
    min_threshold: float = 0.05,
    n_jobs: int = 1,
    kmer_size: int = 31,
    scaled: int = 100,
) -> Dict[str, List[dict]]:
    """Build SBT, save to disk, and search in parallel across worker processes.

    The SBT index is built once (fast) and serialised to a temp directory.
    Each worker process loads its own copy of the tree from disk and searches
    an independent chunk of query signatures.  This gives near-linear speedup
    with the number of workers.
    """
    n_sigs = len(sig_items)
    print(f"Building SBT index for {n_sigs} signatures (ksize={kmer_size}, scaled={scaled})...")
    t0 = time.time()
    # Build a single SBT containing all signatures
    factory = GraphFactory(ksize=kmer_size, n_tables=1, starting_size=max(1, n_sigs))
    sbt = SBT(factory)
    for region_name, sig in tqdm(sig_items, miniters=1000, desc="Indexing SBT"):
        sbt.add_node(SigLeaf(region_name, sig))
    print(f"  SBT built in {time.time() - t0:.2f}s")

    # Save SBT to temp directory for worker loading
    tmp_dir = tempfile.mkdtemp(prefix="sbt_parallel_")
    sbt_path = os.path.join(tmp_dir, "index.sbt.json")
    sbt.save(sbt_path)
    del sbt  # free memory in the main process
    print(f"  SBT saved to {sbt_path}")

    # Serialise queries as (name, hashes_tuple) for pickling
    serialized = [(name, tuple(sig.minhash.hashes)) for name, sig in sig_items]

    sum_comparisons: Dict[str, List[dict]] = defaultdict(list)

    with open(output_csv, "w", newline="") as outfh:
        w = csv.writer(outfh)
        w.writerow(["reference1", "start1", "end1",
                    "reference2", "start2", "end2",
                    "jaccard", "containment_1_in_2", "containment_2_in_1"])

        if n_sigs <= 200 or n_jobs <= 1:
            # Serial path: load SBT once, search sequentially
            _init_sbt_parallel_worker(sbt_path, kmer_size, scaled)
            results = _sbt_search_chunk((serialized, min_threshold))
            for (r1, s1, e1, r2, s2, e2, j, c12, c21) in results:
                w.writerow([r1, s1, e1, r2, s2, e2,
                            f"{j:.6f}", f"{c12:.6f}", f"{c21:.6f}"])
                sum_comparisons[r1].append(dict(jaccard=j, to=r2, s1=s1, e1=e1, s2=s2, e2=e2))
                sum_comparisons[r2].append(dict(jaccard=j, to=r1, s1=s2, e1=e2, s2=s1, e2=e1))
        else:
            # Parallel path: each worker loads its own SBT from disk
            chunk_size = max(1, -(-n_sigs // (n_jobs * 4)))
            chunks = [serialized[i:i + chunk_size]
                      for i in range(0, n_sigs, chunk_size)]
            actual_workers = min(n_jobs, len(chunks))
            print(f"  Searching {n_sigs} signatures: {len(chunks)} chunk(s) "
                  f"across {actual_workers} worker(s) (SBT loaded per worker)")
            with ProcessPoolExecutor(
                max_workers=actual_workers,
                initializer=_init_sbt_parallel_worker,
                initargs=(sbt_path, kmer_size, scaled),
            ) as pool:
                jobs = [pool.submit(_sbt_search_chunk, (chunk, min_threshold))
                        for chunk in chunks]
                for fut in tqdm(as_completed(jobs), total=len(jobs),
                                desc="  SBT search", unit="chunk"):
                    for (r1, s1, e1, r2, s2, e2, j, c12, c21) in fut.result():
                        w.writerow([r1, s1, e1, r2, s2, e2,
                                    f"{j:.6f}", f"{c12:.6f}", f"{c21:.6f}"])
                        sum_comparisons[r1].append(
                            dict(jaccard=j, to=r2, s1=s1, e1=e1, s2=s2, e2=e2))
                        sum_comparisons[r2].append(
                            dict(jaccard=j, to=r1, s1=s2, e1=e2, s2=s1, e2=e1))

    # Clean up temp SBT files
    shutil.rmtree(tmp_dir, ignore_errors=True)

    print(f"  SBT search complete: {sum(len(v) for v in sum_comparisons.values()) // 2} "
          f"cross-reference pairs found")
    return sum_comparisons


# ── Worker state for parallel region–region search ───────────────────────────
# Mirrors the pattern used in report_shared_windows_across_fastas but keyed by
# reference name instead of FASTA label so same-reference pairs are skipped O(1).
_REGION_WORKER_TARGETS: Dict[str, List[Tuple]] = {}


def _init_region_search_worker(serialized: List[Tuple]) -> None:
    """Worker initializer: parse & group all region targets once per process.

    Each entry in serialized is (region_name, hashes_tuple).
    Targets are grouped by reference name (the part before the first ':') so
    _search_region_chunk can skip entire same-reference groups in O(1).
    """
    global _REGION_WORKER_TARGETS
    tgt: Dict[str, List[Tuple]] = defaultdict(list)
    for region_name, hashes in serialized:
        ref, s, e = parse_split(region_name)
        m_set = frozenset(hashes)
        if m_set:
            tgt[ref].append((region_name, m_set, ref, s, e))
    _REGION_WORKER_TARGETS = dict(tgt)


def _search_region_chunk(args: Tuple) -> List[Tuple]:
    """Worker: frozenset Jaccard between a query chunk and all pre-loaded targets.

    Same-reference pairs are skipped in O(1) by iterating target groups keyed
    by reference.  Only pairs sharing ≥1 hash are evaluated for full Jaccard.
    Returns list of (r1,s1,e1, r2,s2,e2, jaccard) tuples.
    """
    query_items, min_threshold = args
    results: List[Tuple] = []
    for region_name, q_hashes in query_items:
        q_ref, q_s, q_e = parse_split(region_name)
        q_set = frozenset(q_hashes)
        if not q_set:
            continue
        q_size = len(q_set)
        hits = []
        for t_ref, t_items in _REGION_WORKER_TARGETS.items():
            if t_ref == q_ref:
                continue  # skip same-reference group in O(1)
            for m_name, m_set, m_ref, m_s, m_e in t_items:
                n_shared = len(q_set & m_set)
                if n_shared == 0:
                    continue
                n_union = q_size + len(m_set) - n_shared
                j = n_shared / n_union
                if j < min_threshold:
                    continue
                hits.append((j, m_ref, m_s, m_e, n_shared / q_size, n_shared / len(m_set)))
        for j, m_ref, m_s, m_e, c12, c21 in hits:
            results.append((q_ref, q_s, q_e, m_ref, m_s, m_e, j, c12, c21))
    return results


def compare_regions_parallel(
    sig_items: List[Tuple[str, SourmashSignature]],
    output_csv: str,
    min_threshold: float = 0.05,
    n_jobs: int = 1,
) -> Dict[str, List[dict]]:
    """Replace the SBT-based search with parallel frozenset Jaccard.

    The sourmash SBT is not designed for 2M+ signatures: search time per query
    grows with tree depth and the build itself exhausts memory.  This function
    replicates the initializer-based pattern from report_shared_windows_across_fastas:

    - Targets are serialised as (region_name, hashes_tuple) and loaded ONCE per
      worker process at startup (not once per job).
    - Same-reference pairs are skipped O(1) by grouping targets by reference name.
    - Most inter-reference pairs share zero hashes → early-exit in one bitwise op.
    - Fully parallelised across n_jobs workers.
    """
    # Serialise: drop full SourmashSignature objects, keep only hashes tuple
    serialized = [(name, tuple(sig.minhash.hashes)) for name, sig in sig_items]
    n_queries = len(serialized)

    sum_comparisons: Dict[str, List[dict]] = defaultdict(list)

    with open(output_csv, "w", newline="") as outfh:
        w = csv.writer(outfh)
        w.writerow(["reference1", "start1", "end1",
                    "reference2", "start2", "end2",
                    "jaccard", "containment_1_in_2", "containment_2_in_1"])

        if n_queries <= 200 or n_jobs <= 1:
            # Serial path
            _init_region_search_worker(serialized)
            results = _search_region_chunk((serialized, min_threshold))
            for (r1, s1, e1, r2, s2, e2, j, c12, c21) in results:
                w.writerow([r1, s1, e1, r2, s2, e2,
                            f"{j:.6f}", f"{c12:.6f}", f"{c21:.6f}"])
                sum_comparisons[r1].append(dict(jaccard=j, to=r2, s1=s1, e1=e1, s2=s2, e2=e2))
                sum_comparisons[r2].append(dict(jaccard=j, to=r1, s1=s2, e1=e2, s2=s1, e2=e1))
        else:
            chunk_size = max(1, -(-n_queries // n_jobs))
            chunks = [serialized[i:i + chunk_size]
                      for i in range(0, n_queries, chunk_size)]
            actual_workers = min(n_jobs, len(chunks))
            print(f"  Searching {n_queries} regions: {len(chunks)} chunk(s) "
                  f"across {actual_workers} worker(s) "
                  f"(targets pre-loaded via initializer; no SBT)")
            with ProcessPoolExecutor(
                max_workers=actual_workers,
                initializer=_init_region_search_worker,
                initargs=(serialized,),
            ) as pool:
                jobs = [pool.submit(_search_region_chunk, (chunk, min_threshold))
                        for chunk in chunks]
                for fut in tqdm(as_completed(jobs), total=len(jobs),
                                desc="  Region search", unit="chunk"):
                    for (r1, s1, e1, r2, s2, e2, j, c12, c21) in fut.result():
                        w.writerow([r1, s1, e1, r2, s2, e2,
                                    f"{j:.6f}", f"{c12:.6f}", f"{c21:.6f}"])
                        sum_comparisons[r1].append(
                            dict(jaccard=j, to=r2, s1=s1, e1=e1, s2=s2, e2=e2))
                        sum_comparisons[r2].append(
                            dict(jaccard=j, to=r1, s1=s2, e1=e2, s2=s1, e2=e1))

    print(f"  Region search complete: {sum(len(v) for v in sum_comparisons.values()) // 2} "
          f"cross-reference pairs found")
    return sum_comparisons


def fast_mode_sbt(
    siglist: List[Tuple[str, SourmashSignature]],
    sbt_index: Dict[str, SBT],
    output_csv: str,
    min_threshold: float = 0.05,
    cluster_map: Optional[List[set]] = None,
) -> Dict[str, List[dict]]:
    """
    Search each region signature against its cluster SBT and record matches.
    Returns sum_comparisons[ref] = [{jaccard,to,s1,e1,s2,e2}, ...]
    """
    ref2cluster: Dict[str, str] = {}
    if cluster_map:
        for idx, cset in enumerate(cluster_map):
            for ref in cset:
                ref2cluster[ref] = str(idx)

    sum_comparisons: Dict[str, List[dict]] = defaultdict(list)

    with open(output_csv, "w", newline="") as outfh:
        w = csv.writer(outfh)
        w.writerow(
            [
                "reference1",
                "start1",
                "end1",
                "reference2",
                "start2",
                "end2",
                "jaccard",
                "containment_1_in_2",
                "containment_2_in_1",
            ]
        )

        for region1, sig1 in siglist:
            r1, s1, e1 = parse_split(region1)
            sbt = sbt_index.get(ref2cluster.get(r1, "Unknown"))
            if sbt is None:
                continue

            for sr in sbt.search(sig1, threshold=min_threshold, best_only=False):
                region2 = sr.signature.name
                r2, s2, e2 = parse_split(region2)
                if r2 == r1:
                    continue

                j = float(sr.score)

                # containment values can be computed, but you were using j as proxy in fast mode
                c12 = j
                c21 = j

                w.writerow([r1, s1, e1, r2, s2, e2, j, c12, c21])

                sum_comparisons[r1].append(dict(jaccard=j, to=r2, s1=s1, e1=e1, s2=s2, e2=e2))
                sum_comparisons[r2].append(dict(jaccard=j, to=r1, s1=s2, e1=e2, s2=s1, e2=e1))

    return sum_comparisons


def load_region_comparisons_csv(csv_path: str) -> Dict[str, List[dict]]:
    """Reload sum_comparisons from a previously-written region_comparisons.csv.

    This mirrors the dict structure produced by fast_mode_sbt so that
    build_conflict_groups can consume it directly.
    """
    sum_comparisons: Dict[str, List[dict]] = defaultdict(list)
    with open(csv_path, "r", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            r1 = row["reference1"]
            s1 = int(row["start1"])
            e1 = int(row["end1"])
            r2 = row["reference2"]
            s2 = int(row["start2"])
            e2 = int(row["end2"])
            j = float(row["jaccard"])
            sum_comparisons[r1].append(dict(jaccard=j, to=r2, s1=s1, e1=e1, s2=s2, e2=e2))
            sum_comparisons[r2].append(dict(jaccard=j, to=r1, s1=s2, e1=e2, s2=s1, e2=e1))
    return sum_comparisons


def build_conflict_groups(sum_comparisons: Dict[str, List[dict]], min_jaccard: float = 0.0) -> List[set]:
    """Connected components on region nodes using jaccard>=min_jaccard edges."""
    from collections import deque

    graph = defaultdict(set)
    for ref, lst in sum_comparisons.items():
        for c in lst:
            j = float(c.get("jaccard", 0.0))
            if j < min_jaccard:
                continue
            a = (ref, int(c["s1"]), int(c["e1"]))
            b = (str(c["to"]), int(c["s2"]), int(c["e2"]))
            graph[a].add(b)
            graph[b].add(a)

    visited = set()
    comps = []

    for node in graph:
        if node in visited:
            continue
        q = deque([node])
        visited.add(node)
        comp = {node}
        while q:
            cur = q.popleft()
            for nb in graph[cur]:
                if nb not in visited:
                    visited.add(nb)
                    comp.add(nb)
                    q.append(nb)
        comps.append(comp)

    return comps


# -----------------------------
# Shared-window ambiguity (reference windowing)
# -----------------------------

def iter_windows(length: int, window: int, step: int):
    for s in range(0, max(0, length - window + 1), step):
        yield s, s + window


def make_window_id(fasta_label: str, contig: str, start: int, end: int) -> str:
    return f"{fasta_label}::{contig}:{start}-{end}"


def parse_window_id(wid: str) -> Tuple[str, str, int, int]:
    fasta_label, rest = wid.split("::", 1)
    contig, coords = rest.split(":", 1)
    s, e = coords.split("-", 1)
    return fasta_label, contig, int(s), int(e)

def _sketch_window(args):
    """Top-level worker: sketch a single window (must be picklable for ProcessPoolExecutor)."""
    wid, subseq, ksize, scaled = args
    mh = MinHash(n=0, ksize=ksize, scaled=scaled)
    mh.add_sequence(subseq, force=True)
    return wid, SourmashSignature(mh, name=wid)


def _sketch_window_batch(args):
    """Top-level worker: sketch a batch of windows to reduce IPC overhead.

    Accepts (batch_of_tasks, ksize, scaled) where batch_of_tasks is a list of
    (wid, subseq) tuples.  Returns a list of (wid, SourmashSignature).
    """
    batch, ksize, scaled = args
    results = []
    for wid, subseq in batch:
        mh = MinHash(n=0, ksize=ksize, scaled=scaled)
        mh.add_sequence(subseq, force=True)
        results.append((wid, SourmashSignature(mh, name=wid)))
    return results


def _iter_window_tasks(fasta_path, fasta_label, ksize, window, step, max_n_frac, contigs):
    """Generator that yields (wid, subseq) one contig at a time.

    This avoids holding all subsequences in memory simultaneously -- only one
    contig's worth of data is live at a time.
    """
    fa = pysam.FastaFile(fasta_path)
    contig_list = contigs if contigs is not None else list(fa.references)
    for contig in contig_list:
        if contig not in fa.references:
            continue

        seq = fa.fetch(contig).upper()
        L = len(seq)
        if L < ksize:
            del seq
            continue

        if L < window:
            win_list = [(0, L)]
        else:
            win_list = list(iter_windows(L, window, step))

        for s, e in win_list:
            subseq = seq[s:e]
            if len(subseq) < ksize:
                continue
            if max_n_frac is not None and max_n_frac >= 0:
                if subseq.count("N") / len(subseq) > max_n_frac:
                    continue
            wid = make_window_id(fasta_label, contig, s, e)
            yield wid, subseq

        del seq  # free contig memory before moving to next

    fa.close()


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
    n_jobs: Optional[int] = None,
    batch_size: int = 500,
) -> Dict[str, SourmashSignature]:
    """Sketch fixed windows across contigs; if contig shorter than window, sketch the full contig once.

    Parameters
    ----------
    n_jobs : int or None
        Number of parallel workers for sketching.  ``None`` (default) uses all CPUs.
        Set to 1 to disable parallelism.
    batch_size : int
        Number of windows per worker batch.  Limits peak memory by streaming
        tasks to workers in fixed-size batches instead of materializing all
        subsequences up front.
    """
    if fasta_label is None:
        fasta_label = os.path.basename(fasta_path)
    if n_jobs is None:
        n_jobs = os.cpu_count() or 1

    sigs: Dict[str, SourmashSignature] = {}

    task_iter = _iter_window_tasks(fasta_path, fasta_label, ksize, window, step, max_n_frac, contigs)

    if n_jobs <= 1:
        # Serial: process one window at a time -- no list accumulation
        for wid, subseq in task_iter:
            mh = MinHash(n=0, ksize=ksize, scaled=scaled)
            mh.add_sequence(subseq, force=True)
            sigs[wid] = SourmashSignature(mh, name=wid)
    else:
        # Parallel: stream batches to workers to cap memory
        def _batch_iter():
            batch = []
            for wid, subseq in task_iter:
                batch.append((wid, subseq))
                if len(batch) >= batch_size:
                    yield (batch, ksize, scaled)
                    batch = []
            if batch:
                yield (batch, ksize, scaled)

        with ProcessPoolExecutor(max_workers=n_jobs) as pool:
            for result_batch in pool.map(_sketch_window_batch, _batch_iter()):
                for wid, sig in result_batch:
                    sigs[wid] = sig

    return sigs


def build_sbt(sigs: Dict[str, SourmashSignature], ksize: int) -> SBT:
    factory = GraphFactory(ksize=ksize, n_tables=1, starting_size=max(1, len(sigs)))
    sbt = SBT(factory)
    for wid, sig in sigs.items():
        sbt.add_node(SigLeaf(wid, sig))
    return sbt


# ---- helpers for parallel report_shared_windows_across_fastas ----

def _generate_sigs_for_fasta(args):
    """Top-level worker for parallel signature generation (must be picklable).

    When called from a ProcessPoolExecutor (multi-FASTA path), n_jobs_inner
    should be 1 to avoid nested pools.  When called from the serial path,
    n_jobs_inner enables intra-file parallelism.

    aligned_refs (optional set/frozenset): when provided, only contigs whose
    ID appears in this set are sketched.  Contigs with no aligned reads are
    silently skipped, dramatically reducing work for large reference databases.
    """
    fp, ksize, scaled, window, step, max_windows_per_fasta, n_jobs_inner, aligned_refs = args
    label = os.path.basename(fp)
    # Build the contig allowlist for this FASTA: only keep contigs present in
    # aligned_refs so we don't spend time sketching zero-read references.
    contigs_filter = None
    if aligned_refs is not None:
        contigs_filter = [c for c in aligned_refs]  # window_sigs_from_fasta expects a list
    sigs = window_sigs_from_fasta(
        fp, fasta_label=label, ksize=ksize, scaled=scaled,
        window=window, step=step, n_jobs=n_jobs_inner,
        contigs=contigs_filter,
    )
    if max_windows_per_fasta is not None and len(sigs) > max_windows_per_fasta:
        sigs = dict(list(sigs.items())[:max_windows_per_fasta])
    return sigs


def _search_chunk(args):
    """Top-level worker: compare a chunk of query windows against a chunk of target sigs.

    Rebuilds MinHash objects once per target (cached for the chunk) and once per
    query, then does brute-force Jaccard comparison.

    Changed from original: accepts a *target chunk* instead of ALL targets,
    so each worker only reconstructs a subset -- dramatically reducing peak memory.
    """
    (query_items, target_items, jaccard_threshold, max_hits_per_query,
     skip_self_same_fasta, skip_self_same_contig) = args

    # Pre-build target MinHash objects for this chunk only
    target_mhs = {}
    target_parsed = {}
    for m_wid, m_mh_hashes, m_mh_ksize, m_mh_scaled in target_items:
        mh_m = MinHash(n=0, ksize=m_mh_ksize, scaled=m_mh_scaled)
        mh_m.add_many(m_mh_hashes)
        if mh_m:
            target_mhs[m_wid] = mh_m
            target_parsed[m_wid] = parse_window_id(m_wid)

    results = []
    for q_wid, q_mh_hashes, q_mh_ksize, q_mh_scaled in query_items:
        q_fa, q_contig, q_s, q_e = parse_window_id(q_wid)

        mh_q = MinHash(n=0, ksize=q_mh_ksize, scaled=q_mh_scaled)
        mh_q.add_many(q_mh_hashes)

        if not mh_q:
            continue

        hits = []
        for m_wid, mh_m in target_mhs.items():
            if m_wid == q_wid:
                continue

            m_fa, m_contig, m_s, m_e = target_parsed[m_wid]

            if skip_self_same_fasta and (m_fa == q_fa):
                continue
            if skip_self_same_contig and (m_fa == q_fa and m_contig == q_contig):
                continue

            j = mh_q.jaccard(mh_m)
            if j < jaccard_threshold:
                continue

            c1 = mh_q.avg_containment(mh_m)
            c2 = mh_m.avg_containment(mh_q)
            hits.append((j, m_fa, m_contig, m_s, m_e, c1, c2))

        if not hits:
            continue
        hits.sort(key=lambda x: x[0], reverse=True)
        hits = hits[:max_hits_per_query]

        for (j, m_fa, m_contig, m_s, m_e, c1, c2) in hits:
            results.append((q_fa, q_contig, q_s, q_e, m_fa, m_contig, m_s, m_e, j, c1, c2))

    return results


# ── Worker-initializer state for fast pairwise search ───────────────────────
# Populated once per worker process by _init_search_worker.
#
# _SEARCH_WORKER_TARGETS          — windows grouped by FASTA label (default path)
#                                    lets same-FASTA groups be skipped in O(1)
# _SEARCH_WORKER_TARGETS_BY_CONTIG — windows grouped by contig name (ANI-filter path)
#                                    lets the worker jump directly to partner contigs
# _SEARCH_WORKER_PARTNER_CONTIGS  — {contig_id: frozenset(partner_contig_ids)} from
#                                    Pass 0 ANI pre-filter; None → no filtering
_SEARCH_WORKER_TARGETS: Dict[str, List[Tuple]] = {}
_SEARCH_WORKER_TARGETS_BY_CONTIG: Dict[str, List[Tuple]] = {}
_SEARCH_WORKER_PARTNER_CONTIGS: Optional[Dict[str, frozenset]] = None


def _init_search_worker(
    serialized: List[Tuple],
    partner_contigs: Optional[Dict[str, frozenset]] = None,
) -> None:
    """Worker initializer: parse & group all target windows once per process.

    Builds two indices from *serialized*:
      - by FASTA label  (for the default full-scan path, O(1) same-FASTA skip)
      - by contig name  (for the ANI-filtered path, O(1) partner lookup)

    partner_contigs, when provided, restricts Phase 2 comparisons to windows
    whose contig is a known high-ANI partner of the query contig.
    """
    global _SEARCH_WORKER_TARGETS, _SEARCH_WORKER_TARGETS_BY_CONTIG, _SEARCH_WORKER_PARTNER_CONTIGS
    tgt_fa: Dict[str, List[Tuple]] = defaultdict(list)
    tgt_contig: Dict[str, List[Tuple]] = defaultdict(list)
    for m_wid, m_hashes, _ksize, _scaled in serialized:
        m_fa, m_contig, m_s, m_e = parse_window_id(m_wid)
        m_set = frozenset(m_hashes)
        if m_set:
            entry = (m_wid, m_set, m_fa, m_contig, m_s, m_e)
            tgt_fa[m_fa].append(entry)
            tgt_contig[m_contig].append(entry)
    _SEARCH_WORKER_TARGETS = dict(tgt_fa)
    _SEARCH_WORKER_TARGETS_BY_CONTIG = dict(tgt_contig)
    _SEARCH_WORKER_PARTNER_CONTIGS = partner_contigs


def _search_query_chunk_np(args: Tuple) -> List[Tuple]:
    """Worker: frozenset-based Jaccard against pre-loaded targets.

    Two execution paths, selected by whether _SEARCH_WORKER_PARTNER_CONTIGS is set:

    Default path (no ANI pre-filter):
        Iterates all targets grouped by FASTA label; same-FASTA groups are
        skipped in O(1).  O(N_targets) per query window.

    ANI-filtered path (partner_contigs provided at init):
        For each query contig, looks up its pre-computed set of high-ANI partner
        contigs and only iterates windows from those contigs.  If a query contig
        has no partners (e.g. it is unique) it is skipped entirely.
        O(n_partner_windows) per query window — typically << O(N_targets) in
        large databases.
    """
    query_items, jaccard_threshold, max_hits_per_query, skip_same_fasta, skip_same_contig = args
    results: List[Tuple] = []
    use_ani_filter = _SEARCH_WORKER_PARTNER_CONTIGS is not None

    for q_wid, q_hashes, _ksize, _scaled in query_items:
        q_fa, q_contig, q_s, q_e = parse_window_id(q_wid)
        q_set = frozenset(q_hashes)
        if not q_set:
            continue
        q_size = len(q_set)
        hits = []

        if use_ani_filter:
            # ── ANI-filtered path ─────────────────────────────────────
            # Look up partner contigs for this query contig.
            # If this contig had no high-ANI partners in Pass 0 it won't
            # appear in the dict at all — skip it entirely.
            partner_contig_ids = _SEARCH_WORKER_PARTNER_CONTIGS.get(q_contig)
            if not partner_contig_ids:
                continue
            for p_contig in partner_contig_ids:
                for m_wid, m_set, m_fa, m_contig, m_s, m_e in \
                        _SEARCH_WORKER_TARGETS_BY_CONTIG.get(p_contig, ()):
                    if skip_same_fasta and m_fa == q_fa:
                        continue
                    # m_contig != q_contig by construction (Pass 0 skips same-FASTA)
                    n_shared = len(q_set & m_set)
                    if n_shared == 0:
                        continue
                    n_union = q_size + len(m_set) - n_shared
                    j = n_shared / n_union
                    if j < jaccard_threshold:
                        continue
                    hits.append((j, m_fa, m_contig, m_s, m_e,
                                  n_shared / q_size, n_shared / len(m_set)))
        else:
            # ── Default full-scan path ────────────────────────────────
            for t_fa, t_items in _SEARCH_WORKER_TARGETS.items():
                if skip_same_fasta and t_fa == q_fa:
                    continue  # skip entire FASTA group — O(1)
                for m_wid, m_set, m_fa, m_contig, m_s, m_e in t_items:
                    if m_wid == q_wid:
                        continue
                    if skip_same_contig and m_fa == q_fa and m_contig == q_contig:
                        continue
                    n_shared = len(q_set & m_set)
                    if n_shared == 0:
                        continue  # early-exit before union math
                    n_union = q_size + len(m_set) - n_shared
                    j = n_shared / n_union
                    if j < jaccard_threshold:
                        continue
                    hits.append((j, m_fa, m_contig, m_s, m_e,
                                  n_shared / q_size, n_shared / len(m_set)))

        if not hits:
            continue
        hits.sort(key=lambda x: x[0], reverse=True)
        for j, m_fa, m_contig, m_s, m_e, c1, c2 in hits[:max_hits_per_query]:
            results.append((q_fa, q_contig, q_s, q_e, m_fa, m_contig, m_s, m_e, j, c1, c2))
    return results


def _find_high_ani_pairs(
    fasta_files: List[str],
    aligned_refs: Optional[set] = None,
    *,
    ksize: int = 21,
    scaled: int = 200,
    ani_threshold: float = 0.90,
) -> Tuple[Dict[str, frozenset], List[Tuple]]:
    """Pass 0: whole-genome MinHash pre-filter to find high-ANI contig pairs.

    One sketch is built per contig (no windowing) at a coarse ksize/scaled that
    is fast to compute and has good Jaccard sensitivity at the target ANI.  Then
    all cross-FASTA contig pairs are tested; only pairs whose Jaccard ≥
    ani_threshold^ksize / (2 - ani_threshold^ksize) are kept.

    Returns
    -------
    partner_map : dict  {contig_id: frozenset(partner_contig_ids)}
        Only contigs that have at least one high-ANI partner appear as keys.
        Partners are always from a *different* FASTA file.
    pair_records : list of (contig_a, fasta_label_a, contig_b, fasta_label_b, jaccard)
        One entry per unique ordered pair (a < b by contig name).  Useful for
        writing a lightweight shared_windows_report.csv without re-running Pass 0.
    """
    p_thresh = ani_threshold ** ksize
    j_thresh = p_thresh / (2.0 - p_thresh)

    # ── sketch each contig once (whole contig, no sliding window) ───────
    # contig_id -> (frozenset_of_hashes, fasta_label)
    contig_sigs: Dict[str, Tuple[frozenset, str]] = {}
    for fp in fasta_files:
        label = os.path.basename(fp)
        try:
            fa = pysam.FastaFile(fp)
        except Exception as exc:
            print(f"  [ANI pass0] WARNING: could not open {fp}: {exc}")
            continue
        for contig in fa.references:
            if aligned_refs is not None and contig not in aligned_refs:
                continue
            seq = fa.fetch(contig)
            mh = MinHash(n=0, ksize=ksize, scaled=scaled)
            mh.add_sequence(seq, force=True)
            if mh.hashes:
                contig_sigs[contig] = (frozenset(mh.hashes), label)
        fa.close()

    n_contigs = len(contig_sigs)
    print(f"  [ANI pass0] Sketched {n_contigs} contig(s) "
          f"(ksize={ksize}, scaled={scaled}); "
          f"ANI threshold={ani_threshold:.2f} → Jaccard threshold={j_thresh:.4f}")

    if n_contigs == 0:
        return {}

    # ── all-vs-all Jaccard across different FASTAs ───────────────────────
    partner_map: Dict[str, set] = defaultdict(set)
    pair_records: List[Tuple] = []  # (contig_a, fa_a, contig_b, fa_b, jaccard)
    contig_list = list(contig_sigs.items())
    n_pairs_checked = 0
    n_pairs_kept = 0
    for i, (cid_a, (hashes_a, fa_a)) in enumerate(contig_list):
        for cid_b, (hashes_b, fa_b) in contig_list[i + 1:]:
            if fa_a == fa_b:
                continue  # skip same-FASTA pairs (handled by skip_self_same_fasta)
            n_pairs_checked += 1
            n_shared = len(hashes_a & hashes_b)
            if n_shared == 0:
                continue
            n_union = len(hashes_a) + len(hashes_b) - n_shared
            j = n_shared / n_union
            if j >= j_thresh:
                partner_map[cid_a].add(cid_b)
                partner_map[cid_b].add(cid_a)
                pair_records.append((cid_a, fa_a, cid_b, fa_b, j))
                n_pairs_kept += 1

    print(f"  [ANI pass0] Checked {n_pairs_checked} cross-FASTA contig pair(s); "
          f"kept {n_pairs_kept} high-ANI pair(s) "
          f"({len(partner_map)} contig(s) have at least one partner)")

    return {k: frozenset(v) for k, v in partner_map.items()}, pair_records


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
    n_jobs: Optional[int] = None,
    aligned_refs: Optional[set] = None,
    ani_prefilter: Optional[float] = 0.90,
) -> None:
    """Sketch windows for each FASTA, then report cross-fasta hits (parallelized).

    Parameters
    ----------
    n_jobs : int or None
        Number of parallel workers.  ``None`` (default) uses all available CPUs.
    aligned_refs : set or None
        When provided, only contigs whose IDs appear in this set are sketched.
        Contigs absent from the set (i.e. no reads aligned to them) are silently
        skipped.  For large reference databases this can reduce sketching time by
        orders of magnitude — a 41k-reference DB where only 200 refs have reads
        means ~99.5% of FASTA sequence is irrelevant and can be skipped.
        Pass ``None`` (default) to sketch all contigs in all FASTAs.
    ani_prefilter : float or None
        When set to a float in (0, 1], a whole-genome sketch pass (Pass 0) is run
        before window sketching.  Each contig is sketched once (no windowing) at
        ksize=21/scaled=200, then all cross-FASTA contig pairs are tested.  Only
        contigs that share a partner with Jaccard ≥ ani_threshold^21/(2-...) are
        window-sketched and compared in Phases 1/2.  Contigs with no high-ANI
        partner are silently skipped — in a 200-organism DB this can reduce Phase
        2 work by 90%+ when only a handful of organisms are genuinely similar.
        Pass ``None`` to disable the pre-filter and sketch+compare everything.
        Default: 0.90 (skip contig pairs with estimated ANI < 90%).
    """
    if n_jobs is None:
        n_jobs = os.cpu_count() or 1

    if aligned_refs is not None:
        print(f"  aligned_refs filter active: will only sketch {len(aligned_refs)} "
              f"contig(s) with BAM alignments (skipping all others)")

    # ── Pass 0: whole-genome ANI pre-filter ─────────────────────────────
    # Build one coarse sketch per contig (no windows) and find all cross-FASTA
    # pairs with estimated ANI ≥ ani_prefilter.  Only those partner contigs are
    # sketched in Phase 1 and compared in Phase 2.  Contigs with no high-ANI
    # partner skip both phases entirely.
    #
    # In a 200-organism DB where only E. coli and Shigella are similar, this
    # reduces Phase 2 work from O(all_windows²) to O(ecoli_windows × shigella_windows).
    _partner_contigs: Optional[Dict[str, frozenset]] = None
    if ani_prefilter is not None and len(fasta_files) > 1:
        _partner_contigs, _ = _find_high_ani_pairs(
            fasta_files,
            aligned_refs=aligned_refs,
            ani_threshold=ani_prefilter,
        )
        if not _partner_contigs:
            print("  [ANI pass0] No high-ANI pairs found; "
                  "falling back to full pairwise search.")
            _partner_contigs = None
        else:
            # Restrict Phase 1 sketching to only contigs that have a partner.
            # Build the union of (query contig IDs ∪ all their partner contig IDs)
            # so that both sides of every high-ANI pair get window-sketched.
            _ani_contig_allowlist: set = set()
            for cid, partners in _partner_contigs.items():
                _ani_contig_allowlist.add(cid)
                _ani_contig_allowlist.update(partners)
            # Intersect with aligned_refs if it was already provided
            if aligned_refs is not None:
                _ani_contig_allowlist &= aligned_refs
            print(f"  [ANI pass0] Restricting Phase 1 sketching to "
                  f"{len(_ani_contig_allowlist)} contig(s) "
                  f"(down from {len(aligned_refs) if aligned_refs is not None else 'all'})")
            # Override aligned_refs with the tighter ANI-based allowlist for Phase 1
            aligned_refs = _ani_contig_allowlist

    _CSV_HEADER = [
        "query_fasta", "query_contig", "q_start", "q_end",
        "match_fasta", "match_contig", "m_start", "m_end",
        "jaccard", "containment_q_in_m", "containment_m_in_q",
    ]

    # ── Phase 1: parallel signature generation ──────────────────────────
    all_sigs: Dict[str, SourmashSignature] = {}

    # Multi-FASTA: parallelize across files (inner sketching is serial to avoid nested pools).
    # Single-FASTA or serial: parallelize across windows inside each file.
    multi_fasta = len(fasta_files) > 1 and n_jobs > 1
    n_jobs_inner = 1 if multi_fasta else n_jobs

    # Pass aligned_refs (frozenset for pickling safety) to each worker
    _aligned_refs_frozen = frozenset(aligned_refs) if aligned_refs is not None else None
    worker_args = [
        (fp, ksize, scaled, window, step, max_windows_per_fasta, n_jobs_inner, _aligned_refs_frozen)
        for fp in fasta_files
    ]

    if not multi_fasta:
        # Serial across files, but each file parallelizes its own windows
        for wa in tqdm(worker_args, desc="Sketching FASTAs", unit="file"):
            all_sigs.update(_generate_sigs_for_fasta(wa))
    else:
        with ProcessPoolExecutor(max_workers=min(n_jobs, len(fasta_files))) as pool:
            futures = {pool.submit(_generate_sigs_for_fasta, wa): wa for wa in worker_args}
            for fut in tqdm(as_completed(futures), total=len(futures),
                           desc="Sketching FASTAs", unit="file"):
                all_sigs.update(fut.result())

    print(f"  Generated {len(all_sigs)} window signatures from {len(fasta_files)} FASTAs")

    if not all_sigs:
        print("WARNING: No window signatures generated; skipping shared-window comparisons "
              "(check window/step/ksize and FASTA content).")
        with open(output_csv, "w", newline="") as out:
            csv.writer(out).writerow(_CSV_HEADER)
        return

    # ── Scale max_hits_per_query with unique contig count ────────────────
    # A fixed max_hits_per_query causes cross-genus hits (e.g. E. coli vs
    # Shigella, Jaccard ~0.22 at ksize=51) to be silently evicted in large
    # databases.  When a DB has 10+ E. coli strains, an E. coli query window's
    # top-12 list fills entirely with intra-genus hits (Jaccard ~0.63), and
    # Shigella never appears — even though it would have appeared in a 2-
    # organism run because there was no intra-genus competition.
    #
    # Fix: derive an effective cap from the number of unique target contigs
    # so that the cap grows proportionally with DB size.  The caller-supplied
    # max_hits_per_query becomes a MINIMUM floor rather than an absolute cap.
    #
    #   n_unique_contigs = number of distinct reference sequences in the DB
    #   effective_cap = max(max_hits_per_query, ceil(n_unique_contigs * 0.10))
    #
    # At 10% of contigs, a 2-organism DB (2 contigs) gives cap=max(N,1)≈N
    # (unchanged), while a 200-organism DB with ~400 contigs gives cap≈40 —
    # enough to surface cross-genus hits that would otherwise be crowded out.
    # The 10% factor keeps memory and output size manageable even for large DBs.
    _n_unique_contigs = len({
        sig.split("::")[1] if "::" in sig else sig
        for sig in (wid for wid, _, _, _ in (
            # re-use the (wid, hashes, ksize, scaled) tuples we're about to build
            (wid, None, None, None) for wid in (s.name if hasattr(s, "name") else k
                                                  for k, s in all_sigs.items())
        ))
    }) if all_sigs else 0
    # Simpler: count distinct contig names from the sig keys
    _n_unique_contigs = len({
        parse_window_id(wid)[1]          # index 1 = contig name
        for wid in all_sigs.keys()
    })
    _effective_max_hits = max(max_hits_per_query, math.ceil(_n_unique_contigs * 0.10))
    if _effective_max_hits != max_hits_per_query:
        print(f"  Scaling max_hits_per_query: {max_hits_per_query} → {_effective_max_hits} "
              f"({_n_unique_contigs} unique contigs × 10% floor)")
    max_hits_per_query = _effective_max_hits

    # ── Phase 2: parallel pairwise search ───────────────────────────────
    # Serialize MinHash data as (wid, hashes, ksize, scaled) tuples for pickling.
    # We convert hashes to a compact tuple instead of list to save memory.
    serialized = []
    for wid, sig in all_sigs.items():
        mh = sig.minhash
        serialized.append((wid, tuple(mh.hashes), mh.ksize, mh.scaled))

    # Free the full SourmashSignature objects now that we have serialized data
    del all_sigs

    n_queries = len(serialized)

    _filter_label = (f"ANI≥{ani_prefilter} pre-filter active"
                     if _partner_contigs is not None else "full scan")
    if n_queries <= 50 or n_jobs <= 1:
        # Serial: initialise worker state in-process then run directly
        _init_search_worker(serialized, partner_contigs=_partner_contigs)
        all_results = _search_query_chunk_np((
            serialized, jaccard_threshold, max_hits_per_query,
            skip_self_same_fasta, skip_self_same_contig,
        ))
    else:
        # ── Initializer-based parallel search ──────────────────────────
        # Each worker receives the full target list ONCE at startup via the
        # initializer (not once per job).  Per-job payload is only the small
        # query chunk.  This eliminates the O(N_tiles × target_data) IPC
        # overhead of the previous tiled approach.
        #
        # partner_contigs (from Pass 0) is also passed to the initializer so
        # each worker knows which target contigs are high-ANI partners.  When
        # set, the worker skips all non-partner contigs in O(1) per query window.
        chunk_size = max(1, -(-n_queries // n_jobs))  # ceiling division
        query_chunks = [
            serialized[i:i + chunk_size]
            for i in range(0, n_queries, chunk_size)
        ]
        actual_workers = min(n_jobs, len(query_chunks))
        print(f"  Searching {n_queries} windows: {len(query_chunks)} query chunks "
              f"across {actual_workers} workers "
              f"(targets pre-loaded via initializer; {_filter_label})")

        all_results = []
        with ProcessPoolExecutor(
            max_workers=actual_workers,
            initializer=_init_search_worker,
            initargs=(serialized, _partner_contigs),
        ) as pool:
            jobs = [
                pool.submit(
                    _search_query_chunk_np,
                    (chunk, jaccard_threshold, max_hits_per_query,
                     skip_self_same_fasta, skip_self_same_contig),
                )
                for chunk in query_chunks
            ]
            for fut in tqdm(as_completed(jobs), total=len(jobs),
                            desc="  Tile search", unit="chunk"):
                all_results.extend(fut.result())

    # Free serialized data before writing results
    del serialized

    # ── Phase 3: write results ──────────────────────────────────────────
    with open(output_csv, "w", newline="") as out:
        w = csv.writer(out)
        w.writerow(_CSV_HEADER)
        for (q_fa, q_contig, q_s, q_e, m_fa, m_contig, m_s, m_e, j, c1, c2) in all_results:
            w.writerow([q_fa, q_contig, q_s, q_e, m_fa, m_contig, m_s, m_e,
                        f"{j:.6f}", f"{c1:.6f}", f"{c2:.6f}"])

    print(f"  Wrote {len(all_results)} shared-window hits to {output_csv}")


@dataclass(frozen=True)
class SharedWindow:
    contig: str
    start: int
    end: int
    alt_contig: str
    alt_start: int
    alt_end: int
    jaccard: float


def load_shared_windows_csv(
    csv_path: str,
    *,
    min_jaccard: float = 1.0,
    skip_same_contig: bool = True,
) -> Dict[str, List[SharedWindow]]:
    """Load window-pair CSV into an index keyed by contig; includes both directions."""
    idx: Dict[str, List[SharedWindow]] = defaultdict(list)
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

            qc = str(row["query_contig"])
            mc = str(row["match_contig"])
            if skip_same_contig and qc == mc:
                continue

            qs, qe = int(row["q_start"]), int(row["q_end"])
            ms, me = int(row["m_start"]), int(row["m_end"])

            idx[qc].append(SharedWindow(qc, qs, qe, mc, ms, me, j))
            idx[mc].append(SharedWindow(mc, ms, me, qc, qs, qe, j))

    for contig in list(idx.keys()):
        idx[contig].sort(key=lambda w: (w.start, w.end, w.alt_contig))
    return dict(idx)


def _windows_overlapping(windows_sorted: List[SharedWindow], pos_start: int, pos_end: int) -> List[SharedWindow]:
    """Return windows that overlap [pos_start, pos_end). windows_sorted must be sorted by start."""
    if not windows_sorted:
        return []
    starts = [w.start for w in windows_sorted]
    i = np.searchsorted(starts, pos_end, side="left")

    hits = []
    for w in reversed(windows_sorted[:i]):
        if w.end <= pos_start:
            break
        if w.start < pos_end and w.end > pos_start:
            hits.append(w)
    return hits


def alignment_alt_contigs(shared_idx: Dict[str, List[SharedWindow]], contig: str, aln_s: int, aln_e: int):
    """Return (alt_contigs, overlapped_bp) for an alignment on contig."""
    windows = shared_idx.get(contig)
    if not windows:
        return set(), 0

    overlaps = _windows_overlapping(windows, aln_s, aln_e)
    if not overlaps:
        return set(), 0

    alt_set = set()
    ovl_bp = 0
    for w in overlaps:
        s = max(aln_s, w.start)
        e = min(aln_e, w.end)
        if s < e:
            alt_set.add(w.alt_contig)
            ovl_bp += (e - s)

    alt_set.discard(contig)
    return alt_set, ovl_bp


def _build_read_sharing_clusters(
    bam_path: str,
    aligned_refs: Optional[set] = None,
    min_shared_reads: int = 1,
) -> Tuple[Dict[str, frozenset], List[Tuple]]:
    """Build cross-reference conflict clusters purely from BAM multi-mapping.

    Any two references that share ≥ min_shared_reads multi-mapped reads are
    considered to be competing for those reads and are placed in the same
    conflict cluster.  No FASTA files are required.

    Parameters
    ----------
    bam_path : str
        Path to the indexed BAM file.
    aligned_refs : set or None
        When provided, only references in this set are considered.  Reads that
        align exclusively to non-aligned refs are silently ignored.
    min_shared_reads : int
        Minimum number of shared reads for a pair to be kept (default 1).

    Returns
    -------
    conflict_map : dict  {ref_id: frozenset(competing_ref_ids)}
        Only refs with at least one competitor appear as keys.
    pair_records : list of (ref_a, ref_b, shared_read_count)
        One entry per unique pair; shared_read_count is the number of reads
        that align to both.  Useful for writing shared_windows_report.csv.
    """
    bam = pysam.AlignmentFile(bam_path, "rb")

    # collect every reference each read aligns to (primary + secondary)
    read_to_refs: Dict[str, set] = defaultdict(set)
    for r in bam:
        if r.is_unmapped or not r.reference_name:
            continue
        ref = r.reference_name
        if aligned_refs is not None and ref not in aligned_refs:
            continue
        read_to_refs[r.query_name].add(ref)
    bam.close()

    # count how many reads each (ref_a, ref_b) pair shares
    pair_counts: Dict[Tuple[str, str], int] = defaultdict(int)
    for refs in read_to_refs.values():
        if len(refs) < 2:
            continue  # singly-mapped read — no conflict
        ref_list = sorted(refs)
        for i, ra in enumerate(ref_list):
            for rb in ref_list[i + 1:]:
                pair_counts[(ra, rb)] += 1

    # filter by min_shared_reads and build outputs
    conflict_map: Dict[str, set] = defaultdict(set)
    pair_records: List[Tuple] = []
    for (ra, rb), count in pair_counts.items():
        if count < min_shared_reads:
            continue
        conflict_map[ra].add(rb)
        conflict_map[rb].add(ra)
        pair_records.append((ra, rb, count))

    return {k: frozenset(v) for k, v in conflict_map.items()}, pair_records


def _count_reads_per_ref(bam_path: str) -> Dict[str, int]:
    """Per-reference read counts; uses BAM index (fast) and falls back to full scan."""
    counts: Dict[str, int] = {}
    bam = pysam.AlignmentFile(bam_path, "rb")
    try:
        for stat in bam.get_index_statistics():
            counts[stat.contig] = int(stat.mapped)
    except Exception:
        bam.reset()
        for r in bam:
            if not r.is_unmapped and r.reference_name:
                counts[r.reference_name] = counts.get(r.reference_name, 0) + 1
    bam.close()
    return counts


def compute_shared_window_stats(
    shared_idx: Dict[str, List[SharedWindow]],
) -> Dict[str, dict]:
    """Per-reference summary of shared-window conflict information.

    Returns a dict keyed by reference/contig name with:
      n_shared_windows   – how many windows on this ref overlap another ref
      n_conflicting_refs – number of distinct other refs sharing windows
      mean_window_jaccard – average jaccard of those windows
      max_window_jaccard  – maximum jaccard seen
      shared_bp          – total bp of windows (windows may overlap; approx)
    """
    stats: Dict[str, dict] = {}
    for contig, windows in shared_idx.items():
        if not windows:
            stats[contig] = {
                "n_shared_windows": 0,
                "n_conflicting_refs": 0,
                "mean_window_jaccard": 0.0,
                "max_window_jaccard": 0.0,
                "shared_bp": 0,
            }
            continue
        alt_contigs = {w.alt_contig for w in windows}
        jaccards = [w.jaccard for w in windows]
        stats[contig] = {
            "n_shared_windows": len(windows),
            "n_conflicting_refs": len(alt_contigs),
            "mean_window_jaccard": statistics.mean(jaccards),
            "max_window_jaccard": max(jaccards),
            "shared_bp": sum(w.end - w.start for w in windows),
        }
    return stats


def build_removed_ids_best_alignment(
    bam_path: str,
    shared_idx: Dict[str, List[SharedWindow]],
    *,
    penalize_weight: float = 50.0,
    as_weight: float = 0.0,
    drop_contigs: Optional[set[str]] = None,
    drop_if_ambiguous: bool = True,
    min_alt_count: int = 1,
    only_primary: bool = False,
    ref_read_counts: Optional[Dict[str, int]] = None,
    ani_boost_weight: float = 10.0,
) -> Dict[str, List[str]]:
    """
    For each read with multiple alignments, keep the alignment(s) with best score.

    Score = MAPQ + as_weight*AS
            - penalize_weight * (# alt contigs overlapping shared windows)
            + ani_boost_weight * ani_dominance

    ani_dominance: for each alt contig sharing a high-ANI window with this contig,
      accumulate max_window_jaccard * log2(our_reads / alt_reads).  Positive when
      the current reference has MORE reads than its competitor (boosted), negative
      when it has fewer (penalised).  This makes removal aggressive in favour of
      the reference with the highest read support when ANI is high.
    """
    if drop_contigs is None:
        drop_contigs = set()

    # Pre-build (contig -> alt_contig -> max_jaccard) index for O(1) lookup per pair.
    pair_max_jaccard: Optional[Dict[str, Dict[str, float]]] = None
    if ref_read_counts is not None and ani_boost_weight > 0.0:
        _pmj: Dict[str, Dict[str, float]] = defaultdict(lambda: defaultdict(float))
        for _c, _wins in shared_idx.items():
            for _w in _wins:
                if _w.jaccard > _pmj[_c][_w.alt_contig]:
                    _pmj[_c][_w.alt_contig] = _w.jaccard
        pair_max_jaccard = {k: dict(v) for k, v in _pmj.items()}

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

        alt_set, _ = alignment_alt_contigs(shared_idx, contig, s, e)
        alt_n = len(alt_set)

        try:
            AS = int(r.get_tag("AS"))
        except Exception:
            AS = 0

        mapq = int(r.mapping_quality)

        # ANI-weighted read-dominance term:
        #   For every competitor that shares a high-ANI window with this contig,
        #   add jaccard * log2(our_count / their_count).  The reference with the
        #   most reads wins decisively when regions are nearly identical.
        ani_dominance = 0.0
        if pair_max_jaccard is not None and alt_set:
            our_count = max(1, ref_read_counts.get(contig, 1))
            pmj_for_contig = pair_max_jaccard.get(contig, {})
            for alt_contig in alt_set:
                alt_count = max(1, ref_read_counts.get(alt_contig, 1))
                log_ratio = math.log2(our_count / alt_count)
                max_j = pmj_for_contig.get(alt_contig, 0.5)
                ani_dominance += max_j * log_ratio

        score = mapq + as_weight * AS - penalize_weight * alt_n + ani_boost_weight * ani_dominance

        per_read[r.query_name].append(
            dict(contig=contig, score=score, mapq=mapq, AS=AS, alt_n=alt_n)
        )

    bam.close()

    removed = defaultdict(list)
    for rid, alns in per_read.items():
        if not alns:
            continue
        best_score = max(a["score"] for a in alns)
        keep_contigs = {a["contig"] for a in alns if a["score"] == best_score}

        for a in alns:
            c = a["contig"]
            ambiguous = (a["alt_n"] >= min_alt_count)

            if c in drop_contigs and drop_if_ambiguous and ambiguous:
                removed[rid].append(c)
                continue

            if c not in keep_contigs:
                removed[rid].append(c)

    return {rid: sorted(set(refs)) for rid, refs in removed.items() if refs}


def report_removed_read_stats(bam_path: str, removed_read_ids: Dict[str, List[str]]) -> Dict[str, float]:
    """Quick summary of removal plan vs BAM."""
    bam = pysam.AlignmentFile(bam_path, "rb")

    total_reads = set()
    total_alns = 0

    removed_reads = set(removed_read_ids.keys())
    removed_alns = 0

    per_contig_total = defaultdict(int)
    per_contig_removed = defaultdict(int)

    for r in bam:
        if r.is_unmapped:
            continue
        total_alns += 1
        total_reads.add(r.query_name)

        contig = r.reference_name
        per_contig_total[contig] += 1

        if r.query_name in removed_read_ids and contig in removed_read_ids[r.query_name]:
            removed_alns += 1
            per_contig_removed[contig] += 1

    bam.close()

    n_reads = len(total_reads)
    n_removed_reads = len(removed_reads)

    print("\n=== Removal Summary ===")
    print(f"Unique reads:             {n_reads}")
    print(f"Reads with ≥1 removal:    {n_removed_reads} ({100.0*n_removed_reads/max(1,n_reads):.2f}%)")
    print(f"Total alignments:         {total_alns}")
    print(f"Alignments removed:       {removed_alns} ({100.0*removed_alns/max(1,total_alns):.2f}%)")

    print("\nPer-contig alignment removal:")
    for contig in sorted(per_contig_total):
        tot = per_contig_total[contig]
        rem = per_contig_removed.get(contig, 0)
        if tot == 0:
            continue
        print(f"  {contig:20s}  {rem:8d}/{tot:8d}  ({100.0*rem/tot:6.2f}%)")

    return {
        "total_reads": float(n_reads),
        "removed_reads": float(n_removed_reads),
        "pct_reads_removed": 100.0 * n_removed_reads / max(1, n_reads),
        "total_alignments": float(total_alns),
        "removed_alignments": float(removed_alns),
        "pct_alignments_removed": 100.0 * removed_alns / max(1, total_alns),
    }


# -----------------------------
# BAM filtering + coverage + metrics
# -----------------------------

def create_filtered_bam(bam_in: pysam.AlignmentFile, output_bam_path: str, alignments_to_skip: Dict[str, List[str]]) -> None:
    """Write BAM excluding alignments whose (read_id, ref) is listed for removal."""
    out = pysam.AlignmentFile(output_bam_path, "wb", template=bam_in)
    bam_in.reset()

    total = 0
    written = 0

    for r in bam_in:
        total += 1
        if r.query_name in alignments_to_skip and r.reference_name in alignments_to_skip[r.query_name]:
            continue
        out.write(r)
        written += 1

    out.close()
    pysam.index(output_bam_path)
    print(f"Alignments processed: {total}")
    print(f"Alignments written:   {written}")


def calculate_breadth_coverage_from_bam(
    bam_fs: pysam.AlignmentFile,
    reflengths: Dict[str, int],
    removed_ids: Dict[str, List[str]],
) -> Dict[str, float]:
    """
    Breadth of coverage (% of reference covered by any kept alignment).
    Interval merge is done on-the-fly per reference.
    """
    cov = {}
    bam_fs.reset()

    for ref, L in reflengths.items():
        total_covered = 0
        cur_s = None
        cur_e = None

        for r in bam_fs.fetch(ref, 0, L):
            if r.is_unmapped or r.reference_start is None or r.reference_end is None:
                continue
            rid = r.query_name
            if rid in removed_ids and ref in removed_ids[rid]:
                continue

            s = max(0, int(r.reference_start))
            e = min(L, int(r.reference_end))
            if s >= e:
                continue

            if cur_s is None:
                cur_s, cur_e = s, e
            elif s <= cur_e + 1:
                cur_e = max(cur_e, e)
            else:
                total_covered += (cur_e - cur_s)
                cur_s, cur_e = s, e

        if cur_s is not None:
            total_covered += (cur_e - cur_s)

        cov[ref] = (100.0 * total_covered / L) if L > 0 else 0.0

    return cov
def compare_metrics(per_ref_stats: Dict[str, dict], reflengths: Dict[str, int],
                    accession_to_taxid: Optional[Dict[str, str]] = None,
                    taxid_to_name: Optional[Dict[str, str]] = None) -> pd.DataFrame:
    # Denominators for RPKM/RPM (total aligned reads across all refs, pre- and post-removal)
    total_reads_all = sum(st.get("total_reads", 0) for st in per_ref_stats.values())
    total_passed_all = sum(st.get("pass_filtered_reads", 0) for st in per_ref_stats.values())

    rows = []
    for ref, st in sorted(per_ref_stats.items()):
        total_reads = st.get("total_reads", 0)
        pass_reads = st.get("pass_filtered_reads", 0)
        removed_reads = total_reads - pass_reads
        delta_pct = (100.0 * removed_reads / total_reads) if total_reads else 0.0

        b0 = st.get("breadth_old", 0.0)
        b1 = st.get("breadth", 0.0)

        ref_len = max(1, reflengths.get(ref, 1))

        # RPKM = reads / (ref_length_kb * millions_mapped)
        rpkm_pre  = (total_reads  / (ref_len / 1000.0)) / max(1e-9, total_reads_all  / 1e6)
        rpkm_post = (pass_reads   / (ref_len / 1000.0)) / max(1e-9, total_passed_all / 1e6)
        rpm_pre   =  total_reads  / max(1e-9, total_reads_all  / 1e6)
        rpm_post  =  pass_reads   / max(1e-9, total_passed_all / 1e6)

        # Resolve taxid for this accession
        taxid = ""
        if accession_to_taxid:
            taxid = accession_to_taxid.get(ref, "")
            if not taxid:
                ref_noversion = re.sub(r"\.\d+$", "", ref)
                taxid = accession_to_taxid.get(ref_noversion, "")

        # Resolve organism name from taxid
        organism = ""
        if taxid and taxid_to_name:
            organism = taxid_to_name.get(str(taxid), "")

        rows.append(
            {
                "Reference": ref,
                "Taxid": taxid,
                "Organism": organism,
                "Reference Length": reflengths.get(ref, 0),
                "TP Original": st.get("TP Original", 0),
                "FP Original": st.get("FP Original", 0),
                "FN Original": st.get("FN Original", 0),
                "TP New": st.get("TP New", 0),
                "FP New": st.get("FP New", 0),
                "FN New": st.get("FN New", 0),
                "Total Reads": total_reads,
                "Pass Filtered Reads": pass_reads,
                "Proportion Aligned": st.get("proportion_aligned", 0.0),
                "Precision": st.get("precision", 0.0),
                "Recall": st.get("recall", 0.0),
                "F1": st.get("f1", 0.0),
                "Δ All": removed_reads,
                "Δ All%": delta_pct,
                "Breadth Original": b0,
                "Breadth New": b1,
                "Δ Breadth": (b1 - b0),
                "Δ^-1 Breadth": (b1 / b0) if b0 else 0.0,
                # RPKM / RPM (pre- and post-removal)
                "RPKM Pre": rpkm_pre,
                "RPKM Post": rpkm_post,
                "RPM Pre": rpm_pre,
                "RPM Post": rpm_post,
                # Shared-window ANI conflict info
                "Shared Windows": st.get("n_shared_windows", 0),
                "Conflicting Refs": st.get("n_conflicting_refs", 0),
                "Mean Shared ANI": st.get("mean_window_jaccard", 0.0),
                "Max Shared ANI": st.get("max_window_jaccard", 0.0),
                "Shared BP": st.get("shared_bp", 0),
            }
        )

    return pd.DataFrame(rows)
# -----------------------------
# Optional: shared-window parameter tuner (kept because determine_conflicts references it)
# -----------------------------

def infer_gt_from_readname(readname: str) -> Optional[str]:
    parts = readname.split("_")
    if len(parts) < 3:
        return None
    return "_".join(parts[:-2])


def evaluate_removed_vs_bam(bam_path: str, removed_read_ids: Dict[str, List[str]]) -> Dict[str, int]:
    """Count TP/FP and how many were removed, where GT is inferred from read name."""
    total_TP = total_FP = TP_removed = FP_removed = 0

    bam = pysam.AlignmentFile(bam_path, "rb")
    for r in bam:
        if r.is_unmapped or r.reference_name is None:
            continue
        gt = infer_gt_from_readname(r.query_name)
        if gt is None:
            continue

        is_tp = (gt == r.reference_name)
        removed_for_ref = (r.query_name in removed_read_ids and r.reference_name in removed_read_ids[r.query_name])

        if is_tp:
            total_TP += 1
            if removed_for_ref:
                TP_removed += 1
        else:
            total_FP += 1
            if removed_for_ref:
                FP_removed += 1
    bam.close()

    return {
        "total_TP": total_TP,
        "total_FP": total_FP,
        "TP_removed": TP_removed,
        "FP_removed": FP_removed,
        "TP_remaining": total_TP - TP_removed,
        "FP_remaining": total_FP - FP_removed,
    }


def tune_shared_window_params(
    bam_path: str,
    fasta_files: List[str],
    *,
    attempts: int = 12,
    param_grid: Optional[dict] = None,
    sampler: str = "random",
    allowed_tp_loss_frac: float = 0.01,
    alpha: float = 1.0,
    min_jaccard: float = 0.8,
    tmp_dir: str = ".shared_tune_tmp",
):
    """
    Try a small sweep over (ksize, scaled, window, step) and pick the lowest:
        score = FP_remaining + alpha * TP_removed
    """
    os.makedirs(tmp_dir, exist_ok=True)

    if param_grid is None:
        param_grid = {
            "ksize": [21, 31, 51],
            "scaled": [2000, 8000, 20000],
            "window": [5000, 10000, 20000],
            "step": [2500, 5000, 10000],
        }

    combos = list(itertools.product(param_grid["ksize"], param_grid["scaled"], param_grid["window"], param_grid["step"]))
    if sampler == "random":
        trials = combos if attempts >= len(combos) else random.sample(combos, attempts)
    else:
        trials = combos[:attempts]

    best = None
    best_score = float("inf")
    all_results = []

    for idx, (ksize, scaled, window, step) in enumerate(trials, start=1):
        out_csv = os.path.join(tmp_dir, f"shared_windows_trial_{idx}.csv")
        if os.path.exists(out_csv):
            os.remove(out_csv)

        try:
            report_shared_windows_across_fastas(
                fasta_files=fasta_files,
                output_csv=out_csv,
                ksize=ksize,
                scaled=scaled,
                window=window,
                step=step,
                jaccard_threshold=min_jaccard,
                max_hits_per_query=40,
                skip_self_same_fasta=False,
                skip_self_same_contig=True,
            )
            shared_idx = load_shared_windows_csv(out_csv, min_jaccard=min_jaccard, skip_same_contig=True)

            removed_read_ids = build_removed_ids_best_alignment(
                bam_path=bam_path,
                shared_idx=shared_idx,
                penalize_weight=23.0,
                as_weight=0.0,
                drop_contigs=set(),
                drop_if_ambiguous=True,
                min_alt_count=1,
                only_primary=False,
            )

            stats = evaluate_removed_vs_bam(bam_path, removed_read_ids)
            score = stats["FP_remaining"] + alpha * stats["TP_removed"]

            result = dict(idx=idx, ksize=ksize, scaled=scaled, window=window, step=step, score=score, stats=stats)
            all_results.append(result)

            if score < best_score:
                best_score = score
                best = result

        except Exception as e:
            all_results.append(dict(idx=idx, ksize=ksize, scaled=scaled, window=window, step=step, error=str(e)))

    if best is None:
        return False, None, None, all_results

    tot_TP = best["stats"]["total_TP"]
    tp_removed = best["stats"]["TP_removed"]
    fp_remaining = best["stats"]["FP_remaining"]
    tp_loss_frac = (tp_removed / tot_TP) if tot_TP else 0.0

    success = (fp_remaining == 0) and (tp_loss_frac <= allowed_tp_loss_frac)
    best_params = {
        "ksize": best["ksize"],
        "scaled": best["scaled"],
        "window": best["window"],
        "step": best["step"],
        "score": best["score"],
        "tp_loss_frac": tp_loss_frac,
        "fp_remaining": fp_remaining,
    }

    return success, best_params, best["stats"], all_results

def build_organism_signatures_from_fastas_ani(
    fasta_paths: List[str],
    accession_to_taxid: Dict[str, str | int],
    taxid_to_desc: Optional[Dict[str | int, str]] = None,
    *,
    # Sketch params (tune for ANI-ish behavior)
    ksize: int = 21,
    scaled: int = 1000,
    max_n_frac: float = 0.05,
    allow_missing_accessions: bool = True,
) -> Tuple[List[Tuple[str, SourmashSignature]], Dict[str, dict]]:
    """
    Build one sourmash signature per organism (taxid) by combining sequences across
    all accessions mapped to that taxid.

    "ANI-ish" defaults:
      - ksize=21, scaled=1000 (denser sketches give a more stable containment_ani)
      - You still must compare these signatures using `.containment_ani(...)`

    Returns:
      org_siglist: list[(taxid_str, signature)]
      stats: dict[taxid_str] -> counters (JSON-serializable)
    """
    taxid_mh: Dict[str, MinHash] = {}
    stats: Dict[str, dict] = defaultdict(lambda: {
        "n_records_added": 0,
        "n_bases_added": 0,
        "n_records_skipped_N": 0,
        "n_records_skipped_short": 0,
        "n_records_unmapped": 0,
        "source_fastas": set(),
        "accessions_seen": set(),
    })

    def _normalize_accession(header: str) -> str:
        # FASTA record name is usually up to first whitespace
        return header.split()[0]

    for fp in fasta_paths:
        fa = pysam.FastaFile(fp)
        fasta_label = os.path.basename(fp)
        for rec_name in fa.references:
            acc = _normalize_accession(rec_name)
            taxid = accession_to_taxid.get(acc)
            if taxid is None:
                if not allow_missing_accessions:
                    stats["_UNMAPPED_"]["n_records_unmapped"] += 1
                    stats["_UNMAPPED_"]["source_fastas"].add(fasta_label)
                    stats["_UNMAPPED_"]["accessions_seen"].add(acc)
                continue

            taxid_str = str(taxid)
            seq = fa.fetch(rec_name).upper()

            stats[taxid_str]["source_fastas"].add(fasta_label)
            stats[taxid_str]["accessions_seen"].add(acc)

            if not seq or len(seq) < ksize:
                stats[taxid_str]["n_records_skipped_short"] += 1
                continue

            if max_n_frac is not None and max_n_frac >= 0:
                n_frac = seq.count("N") / max(1, len(seq))
                if n_frac > max_n_frac:
                    stats[taxid_str]["n_records_skipped_N"] += 1
                    continue

            if taxid_str not in taxid_mh:
                taxid_mh[taxid_str] = MinHash(n=0, ksize=ksize, scaled=scaled)

            taxid_mh[taxid_str].add_sequence(seq, force=True)

            stats[taxid_str]["n_records_added"] += 1
            stats[taxid_str]["n_bases_added"] += len(seq)
        fa.close()

    org_siglist: List[Tuple[str, SourmashSignature]] = []
    for taxid_str, mh in taxid_mh.items():
        desc = taxid_to_desc.get(taxid_str) if taxid_to_desc else None
        name = f"{taxid_str}" if not desc else f"{taxid_str} {desc}"
        org_siglist.append((taxid_str, SourmashSignature(mh, name=name)))

    # JSON-ify sets
    for k in list(stats.keys()):
        stats[k]["source_fastas"] = sorted(stats[k]["source_fastas"])
        stats[k]["accessions_seen"] = sorted(stats[k]["accessions_seen"])

    org_siglist.sort(key=lambda x: x[0])
    return org_siglist, dict(stats)

def _safe_ani(mh1, mh2):
    res = mh1.containment_ani(mh2)
    return float(res.ani) if res and res.ani is not None else 0.0
def organism_ani_matrix_from_sigs(
    org_siglist: List[Tuple[str, SourmashSignature]],
    *,
    symmetrize: str = "mean",   # "mean" | "min" | "max"
    diagonal: float = 1.0,
) -> pd.DataFrame:
    """
    Build an ANI-like matrix using sourmash MinHash.containment_ani.

    For each pair (i,j):
      ani_ij = mh_i.containment_ani(mh_j).ani
      ani_ji = mh_j.containment_ani(mh_i).ani
    then combine by symmetrize policy.
    """
    taxids = [t for t, _ in org_siglist]
    mhs = {t: sig.minhash for t, sig in org_siglist}

    mat = pd.DataFrame(0.0, index=taxids, columns=taxids, dtype=float)

    for i, t1 in enumerate(taxids):
        mh1 = mhs[t1]
        for j, t2 in enumerate(taxids):
            if t1 == t2:
                mat.at[t1, t2] = diagonal
                continue
            mh2 = mhs[t2]
            a12 = _safe_ani(mh1, mh2)
            a21 = _safe_ani(mh2, mh1)

            if symmetrize == "min":
                v = min(a12, a21)
            elif symmetrize == "max":
                v = max(a12, a21)
            else:
                v = 0.5 * (a12 + a21)

            mat.at[t1, t2] = v

    return mat

def finalize_proportional_removal(conflict_groups, bam_fs, fetch_reads_in_region, remove_mode='random', random_seed=None, dominance_protect_ratio: float = 3.0):
    """
    Dominance-aware read removal across conflict groups.

    For each conflict group the reference with the *most* reads (dominant) is
    protected — fewer of its reads are removed.  References with fewer reads
    (minor) lose proportionally more, reflecting the assumption that shared
    regions in a minor genome are cross-mapped from the dominant one.

    Removal budget per reference in a group:
      base     = min_cov (smallest read count in the group — previous behaviour)
      fraction = 1.0 - (ref_reads / max_reads)   # 0 for dominant, ~1 for minor
      to_remove = base + int(fraction * (ref_reads - base))

    KEY CHANGE — dominant protection via read-count ratio:
      When the dominant ref (most reads in the conflict region) has a clear
      lead over the runner-up (dominance_ratio ≥ dominance_protect_ratio),
      it is completely protected from the base removal so it does not bleed
      reads across every conflict group it participates in.  This prevents
      the accumulation of small-but-repeated base removals that artificially
      depresses coverage of the true positive.

      Dominance is determined purely by raw read count in the shared region.

    Parameters
    ----------
    dominance_protect_ratio : float
        Raw-read-count ratio (top / runner-up) above which the dominant
        reference is fully exempt from the ``base = min_cov`` removal in each
        conflict group.  3.0 means the dominant must have ≥3× the reads of
        the runner-up to be fully protected.  Lower values protect more
        aggressively; higher values are more conservative.  Pass via
        ``--dominance_protect_ratio`` on the command line.

    Returns
    -------
    removed_read_ids : dict[str, list[str]]
        read_id → list of reference names where it is removed
    cluster_dominance : dict[str, dict]
        ref → {n_shared_regions, n_cluster_partners, read_rank,
               cluster_removal_frac, mean_mapq_in_conflict}
        Penalty metadata for downstream minhash scoring.
    """
    # Alias so the rest of the function body stays readable
    DOMINANCE_PROTECT_RATIO = dominance_protect_ratio

    if random_seed is not None:
        random.seed(random_seed)
    global global_bam
    removed_read_ids = defaultdict(list)

    # Accumulate per-reference cluster stats across ALL groups
    _ref_shared_regions: Dict[str, int] = defaultdict(int)
    _ref_cluster_partners: Dict[str, set] = defaultdict(set)
    _ref_total_cluster_reads: Dict[str, int] = defaultdict(int)
    _ref_removed_cluster_reads: Dict[str, int] = defaultdict(int)
    # Track mean mapq per ref across all conflict groups for diagnostic output
    _ref_mapq_sum: Dict[str, float] = defaultdict(float)
    _ref_mapq_count: Dict[str, int] = defaultdict(int)

    for group in conflict_groups:
        # 1) Group sub-regions by reference
        ref_subregions = defaultdict(list)
        for (ref, s, e) in group:
            ref_subregions[ref].append((s, e))

        # Track partner counts and shared-region counts
        all_refs_in_group = set(ref_subregions.keys())
        for ref in all_refs_in_group:
            _ref_shared_regions[ref] += len(ref_subregions[ref])
            _ref_cluster_partners[ref] |= (all_refs_in_group - {ref})

        # 2) Gather reads per reference, collecting mapq alongside read IDs
        coverage_by_ref = {}
        reads_by_ref = {}
        mapq_by_ref: Dict[str, float] = {}
        for ref, subregs in ref_subregions.items():
            unified_reads = {}   # read_id → (read_id, ref, mapq)
            for (start, end) in subregs:
                for read in bam_fs.fetch(ref, start, end):
                    rid = read.query_name
                    if rid not in unified_reads:
                        unified_reads[rid] = (rid, ref, int(read.mapping_quality))
            final_reads_with_mapq = list(unified_reads.values())
            # Strip mapq for the removal list (keeps downstream API unchanged)
            final_reads = [(r_id, r_ref) for r_id, r_ref, _ in final_reads_with_mapq]
            coverage_by_ref[ref] = len(final_reads)
            reads_by_ref[ref] = final_reads
            # Mean mapq of reads hitting this ref within conflict regions
            mapqs = [mq for _, _, mq in final_reads_with_mapq]
            mapq_by_ref[ref] = (sum(mapqs) / len(mapqs)) if mapqs else 0.0
            # Accumulate for cluster_dominance metadata
            _ref_mapq_sum[ref] += sum(mapqs)
            _ref_mapq_count[ref] += len(mapqs)
        if not coverage_by_ref:
            continue

        # 3) Read-count-based dominance
        max_cov = max(coverage_by_ref.values())
        min_cov = min(coverage_by_ref.values())

        # Identify dominant ref and compute dominance ratio (top vs runner-up)
        # using raw read counts only.
        sorted_covs = sorted(coverage_by_ref.values(), reverse=True)
        second_cov = sorted_covs[1] if len(sorted_covs) > 1 else 0
        dominance_ratio = max_cov / max(second_cov, 1e-9)
        dominant_ref = max(coverage_by_ref, key=coverage_by_ref.get)

        for ref, cov_count in coverage_by_ref.items():
            rlist = reads_by_ref[ref]
            _ref_total_cluster_reads[ref] += cov_count

            # Fraction of removal: how minor is this ref relative to the
            # dominant?  0 for dominant, ~1 for clear minor.
            if max_cov > 0:
                fraction = 1.0 - (cov_count / max_cov)
            else:
                fraction = 0.0

            # Base removal for the dominant ref.
            # When dominance_ratio >= DOMINANCE_PROTECT_RATIO the dominant ref
            # is fully protected (base=0), preventing accumulation of small
            # per-group losses across many shared-region groups.
            # Linearly interpolate: ratio=1 → full base; ratio≥threshold → 0.
            if ref == dominant_ref and DOMINANCE_PROTECT_RATIO > 1.0:
                protection = min(
                    1.0,
                    max(0.0, (dominance_ratio - 1.0) / (DOMINANCE_PROTECT_RATIO - 1.0))
                )
                base = int(min_cov * (1.0 - protection))
            else:
                base = min_cov

            # Extra removal beyond base, proportional to how minor the ref is
            extra = int(fraction * max(0, cov_count - base))
            to_remove = min(cov_count, base + extra)

            if to_remove >= cov_count:
                for (read_id, r_ref) in rlist:
                    removed_read_ids[read_id].append(r_ref)
            else:
                if remove_mode == 'random':
                    sampled = random.sample(rlist, to_remove)
                else:
                    sampled = rlist[:to_remove]
                for (read_id, r_ref) in sampled:
                    removed_read_ids[read_id].append(r_ref)

            _ref_removed_cluster_reads[ref] += to_remove

    # Build cluster_dominance metadata for downstream minhash penalty
    cluster_dominance: Dict[str, dict] = {}
    for ref in set(_ref_shared_regions.keys()) | set(_ref_cluster_partners.keys()):
        total = _ref_total_cluster_reads.get(ref, 0)
        removed = _ref_removed_cluster_reads.get(ref, 0)
        cnt = _ref_mapq_count.get(ref, 0)
        cluster_dominance[ref] = {
            "n_shared_regions": _ref_shared_regions.get(ref, 0),
            "n_cluster_partners": len(_ref_cluster_partners.get(ref, set())),
            "cluster_removal_frac": (removed / total) if total > 0 else 0.0,
            # Mean mapQ of reads in conflict regions — used as a quality signal
            # in optimize_weights to protect high-confidence true positives.
            "mean_conflict_mapq": (_ref_mapq_sum[ref] / cnt) if cnt > 0 else 0.0,
        }

    return removed_read_ids, cluster_dominance


def compare_organism_signatures_pairwise(
    org_siglist: List[Tuple[str, SourmashSignature]],
    output_csv: str,
) -> Dict[str, Dict[str, float]]:
    """
    Pairwise compare organism (taxid-level) signatures.

    Writes a long-form CSV:
      taxid1,taxid2,jaccard,containment_1_in_2,containment_2_in_1

    Returns dict-of-dict of jaccard.
    """
    # stable ordering
    org_siglist = sorted(org_siglist, key=lambda x: x[0])

    out_j = defaultdict(dict)
    with open(output_csv, "w", newline="") as out:
        w = csv.writer(out)
        w.writerow(["taxid1", "taxid2", "jaccard", "containment_1_in_2", "containment_2_in_1"])

        for i in range(len(org_siglist)):
            t1, s1 = org_siglist[i]
            mh1 = s1.minhash
            out_j[t1][t1] = 1.0

            for j in range(i + 1, len(org_siglist)):
                t2, s2 = org_siglist[j]
                mh2 = s2.minhash

                jac = float(mh1.jaccard(mh2))
                c12 = float(mh1.avg_containment(mh2))
                c21 = float(mh2.avg_containment(mh1))

                out_j[t1][t2] = jac
                out_j[t2][t1] = jac

                w.writerow([t1, t2, jac, c12, c21])

    return dict(out_j)


def _build_organism_sbt(
    org_siglist: List[Tuple[str, SourmashSignature]],
    *,
    ksize: int,
) -> SBT:
    """
    Build an in-memory SBT over organism signatures for fast candidate retrieval.
    Leaf names are taxids (string).
    """
    factory = GraphFactory(ksize=ksize, n_tables=1, starting_size=max(1, len(org_siglist)))
    sbt = SBT(factory)
    for taxid, sig in org_siglist:
        # make sure leaf name is just the key we want back
        leaf = SigLeaf(str(taxid), sig)
        sbt.add_node(leaf)
    return sbt


def compare_region_signatures_to_organisms(
    region_siglist: List[Tuple[str, SourmashSignature]],
    org_siglist: List[Tuple[str, SourmashSignature]],
    output_csv: str,
    *,
    ksize: int,
    min_jaccard_candidate: float = 0.01,
    top_k: int = 10,
) -> Dict[str, List[dict]]:
    """
    Compare merged-region signatures (your earlier bed/merged regions) to taxid-level organism signatures.

    Strategy:
      1) Use SBT with a low Jaccard threshold to get candidates (fast).
      2) For candidates, compute containment in both directions:
           c_region_in_org = region_mh.avg_containment(org_mh)
           c_org_in_region = org_mh.avg_containment(region_mh)
         and also write the candidate Jaccard from the SBT result.

    Writes CSV:
      region,region_ref,region_start,region_end,taxid,jaccard,containment_region_in_organism,containment_organism_in_region

    Returns:
      hits_by_region: dict[region_name] -> list of hit dicts (sorted best-first by containment_region_in_organism)
    """
    # Build lookup for org signatures by taxid for containment calc
    org_by_taxid = {str(t): sig for t, sig in org_siglist}

    # Build SBT index over organisms
    sbt = _build_organism_sbt(org_siglist, ksize=ksize)

    hits_by_region: Dict[str, List[dict]] = defaultdict(list)

    with open(output_csv, "w", newline="") as out:
        w = csv.writer(out)
        w.writerow([
            "region",
            "region_ref", "region_start", "region_end",
            "taxid",
            "jaccard",
            "containment_region_in_organism",
            "containment_organism_in_region",
        ])

        for region_name, region_sig in region_siglist:
            # Parse region if it's in your "ref:start-end" form
            # If not, keep it as raw in CSV.
            region_ref = ""
            region_start = ""
            region_end = ""
            try:
                region_ref, coords = region_name.split(":", 1)
                s, e = coords.split("-", 1)
                region_start = int(s)
                region_end = int(e)
            except Exception:
                region_ref = ""
                region_start = ""
                region_end = ""

            mh_r = region_sig.minhash

            # SBT candidates based on Jaccard
            candidates = []
            for sr in sbt.search(region_sig, threshold=min_jaccard_candidate, best_only=False):
                taxid = str(sr.signature.name)  # leaf name we used
                jac = float(sr.score)
                candidates.append((jac, taxid))

            if not candidates:
                continue

            # Prefer higher Jaccard first as a proxy for “worth computing containments”
            candidates.sort(key=lambda x: x[0], reverse=True)
            candidates = candidates[:top_k]

            # Compute containments for candidates
            for jac, taxid in candidates:
                org_sig = org_by_taxid.get(taxid)
                if org_sig is None:
                    continue
                mh_o = org_sig.minhash

                c_r_in_o = float(mh_r.avg_containment(mh_o))
                c_o_in_r = float(mh_o.avg_containment(mh_r))

                row = {
                    "region": region_name,
                    "region_ref": region_ref,
                    "region_start": region_start,
                    "region_end": region_end,
                    "taxid": taxid,
                    "jaccard": jac,
                    "containment_region_in_organism": c_r_in_o,
                    "containment_organism_in_region": c_o_in_r,
                }
                hits_by_region[region_name].append(row)

            # Sort region hits by containment (region in organism) descending
            hits_by_region[region_name].sort(
                key=lambda d: (d["containment_region_in_organism"], d["jaccard"]),
                reverse=True
            )

            # Write rows (already sorted)
            for row in hits_by_region[region_name]:
                w.writerow([
                    row["region"],
                    row["region_ref"], row["region_start"], row["region_end"],
                    row["taxid"],
                    row["jaccard"],
                    row["containment_region_in_organism"],
                    row["containment_organism_in_region"],
                ])

    return dict(hits_by_region)



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

def generate_ani_matrix(
    fasta_files: List[str],
    output_dir: str,
    matchfile: Optional[str] = None,
    matchfile_accession_col: int = 0,
    matchfile_taxid_col: int = 4,
    matchfile_desc_col: int = 2,
):
    """Generate a taxid-level ANI matrix from FASTA signatures.

    Returns the ANI DataFrame (index and columns are taxid strings) so callers
    can use it directly without re-reading the CSV, or None if generation was
    skipped.  The matrix is also written to
    ``<output_dir>/organism_ani_matrix.csv`` for caching.
    """
    import pandas as pd
    os.makedirs(output_dir, exist_ok=True)
    if matchfile:
        accession_to_taxid, taxid_to_desc, _ = load_matchfile(
            matchfile, matchfile_accession_col, matchfile_taxid_col, matchfile_desc_col
        )
        org_sigs, _ = build_organism_signatures_from_fastas_ani(
            fasta_paths=fasta_files,
            accession_to_taxid=accession_to_taxid,
            taxid_to_desc=taxid_to_desc,
            ksize=31,
            scaled=200,
        )
        ani_df = organism_ani_matrix_from_sigs(org_sigs, symmetrize="mean")
        csv_path = os.path.join(output_dir, "organism_ani_matrix.csv")
        ani_df.to_csv(csv_path)
        print(f"ANI matrix written to {csv_path}")
        return ani_df
    else:
        print("No matchfile provided; skipping ANI matrix generation.")
        return None
# -----------------------------
# Main pipeline
# -----------------------------

def determine_conflicts(
    output_dir: Optional[str] = None,
    input_bam: Optional[str] = None,
    min_threshold: float = 0.4,
    fasta_files: Optional[List[str]] = None,
    use_variance: bool = False,
    sigfile: Optional[str] = None,
    bedfile: Optional[str] = None,
    scaled: int = 100,
    kmer_size: int = 31,
    filtered_bam_create: Optional[str] = None,
    FAST_MODE: bool = True,
    sensitive: bool = False,
    cpu_count: Optional[int] = None,
    jump_threshold: Optional[float] = None,
    gap_allowance: float = 0.1,
    sim_ani_threshold: float = 0.96,
    compare_to_reference_windows: bool = False,
    find_optimal_windows: bool = False,
    accession_to_taxid: Optional[Dict[str, str]] = None,
    taxid_to_name: Optional[Dict[str, str]] = None,
    taxid_removal_stats: bool = False,
    window_size: int = 5000,
    step_size: int = 2500,
    dominance_protect_ratio: float = 3.0,
):
    # only raise error if bedfile is missing AND compare_to_reference_windows is not selected
    if output_dir is None or input_bam is None or (bedfile is None and not compare_to_reference_windows):
        raise ValueError("output_dir, input_bam, and bedfile are required.")
    os.makedirs(output_dir, exist_ok=True)
    print(f"Starting conflict detection pipeline: {time.ctime()}")


    fasta_files = fasta_files or []

    print(f"Starting conflict detection: {time.ctime()}")
    print(f"FASTA inputs: {len(fasta_files)}")

    if not fasta_files:
        print("No FASTA provided; region sketching will skip regions without FASTA sequence.")
    elif sensitive:
        print("Sensitive mode enabled; ignoring FASTA inputs.")
        fasta_files = []
    else:
        print(f"Using {len(fasta_files)} FASTA file(s) for region sketches.")



    bam_fs = pysam.AlignmentFile(input_bam, "rb")
    reflengths = {ref: bam_fs.get_reference_length(ref) for ref in bam_fs.references}

    t0 = time.time()

    # ── Sketch-density sanity check ───────────────────────────────────────
    # With scaled=5000 and window=100000 you get ~100000/5000 = 20 k-mers per
    # sketch.  At 98% ANI (E. coli vs Shigella) the expected Jaccard is ~0.22
    # at ksize=51.  With only 20 k-mers the Jaccard estimator's std deviation
    # is sqrt(0.22*0.78/20) ≈ 0.09, making detection unreliable — and in a
    # large DB, the noisy low-Jaccard estimate gets evicted by the max_hits_per_query
    # cap before it can be reported.  Warn loudly so users know to reduce scaled
    # or increase window_size.
    _est_kmers_per_window = window_size / scaled
    _MIN_SKETCH_KMERS = 100
    if _est_kmers_per_window < _MIN_SKETCH_KMERS:
        print(
            f"WARNING: Estimated k-mers per window sketch = {_est_kmers_per_window:.0f} "
            f"(window_size={window_size} / scaled={scaled}).  "
            f"This is below the recommended minimum of {_MIN_SKETCH_KMERS} and will produce "
            f"unreliable Jaccard estimates, especially for cross-genus pairs "
            f"(e.g. E. coli vs Shigella, expected Jaccard ~0.22 at ksize=51).  "
            f"Consider reducing --scaled (e.g. --scaled 500) or increasing --window_size "
            f"so that window_size / scaled ≥ {_MIN_SKETCH_KMERS}."
        )


    # Resolve the cached CSV path for the shared-window / region-comparison report.
    # Use sigfile as the path if it points to an existing CSV; otherwise fall back
    # to a default inside output_dir.
    if compare_to_reference_windows:
        default_csv = os.path.join(output_dir, "shared_windows_report.csv")
    else:
        default_csv = os.path.join(output_dir, "region_comparisons.csv")

    if sigfile and sigfile.endswith(".csv"):
        # Guard against cross-mode contamination: a shared_windows_report.csv
        # left over from a prior --compare_references run must not be used as
        # the region_comparisons cache when running without --compare_references
        # (the two formats have completely different columns and are not
        # interchangeable).  Similarly, a region_comparisons.csv should not be
        # accepted as the shared-window cache.
        _sigfile_basename = os.path.basename(sigfile)
        _expected_basename = (
            "shared_windows_report.csv" if compare_to_reference_windows
            else "region_comparisons.csv"
        )
        if _sigfile_basename == _expected_basename:
            report_path = sigfile
        else:
            print(
                f"WARNING: --sigfile '{sigfile}' is a '{_sigfile_basename}' but "
                f"the current mode expects '{_expected_basename}'; "
                f"ignoring sigfile as CSV cache and using default path."
            )
            report_path = default_csv
    else:
        report_path = default_csv

    # Optional: compare_to_reference_windows adds shared-window signatures (and later uses alignment-based removal)
    shared_idx = None
    if sigfile and (not os.path.exists(sigfile) or os.path.getsize(sigfile) <= 0):
        sigfile = None  # only use sigfile if it points to an existing file; otherwise ignore it for signatures

    # ── Pre-compute aligned reference set for FASTA sketch filtering ─────
    # Only sketch FASTA contigs that actually have reads in the BAM.  For
    # large reference databases (tens-of-thousands of references) this
    # eliminates the vast majority of sketching work — e.g. a 41k-reference
    # DB where only 200 references have reads means ~99.5% of sequence can
    # be skipped.  _ref_counts_pre is reused later for ANI-weighted dominance.
    _ref_counts_pre = _count_reads_per_ref(input_bam)
    _aligned_refs = {r for r, c in _ref_counts_pre.items() if c > 0}
    print(f"  BAM aligned refs: {len(_aligned_refs)} of {len(reflengths)} total "
          f"references have ≥1 read — sketching only these contigs")

    # ── Standalone read-sharing conflict-cluster pass (BAM only) ─────────
    # When --compare_references is off the full FASTA window-sketch pipeline is
    # skipped entirely, but we still want cross-genus conflict cluster IDs so
    # that the winner/loser minhash adjustment in optimize_weights can group
    # organisms like E. coli and Shigella into the same bucket.
    #
    # This block derives conflict clusters purely from the BAM: any two
    # references that share ≥1 multi-mapped read are placed in the same cluster.
    # No FASTA files are read.  The result is written as a minimal
    # shared_windows_report.csv (window positions = 0; the union-find in
    # match_paths.py only uses query_contig / match_contig).
    #
    # The file is cached — re-runs skip this step if the CSV already exists.
    _sw_path = os.path.join(output_dir, "shared_windows_report.csv")
    if (
        not compare_to_reference_windows
        and not (os.path.exists(_sw_path) and os.path.getsize(_sw_path) > 0)
    ):
        print("Building read-sharing conflict clusters from BAM "
              "(no FASTA comparison; --compare_references is off)...")
        _rs_map, _rs_pairs = _build_read_sharing_clusters(
            input_bam,
            aligned_refs=_aligned_refs,
        )
        if _rs_pairs:
            _sw_header = [
                "query_fasta", "query_contig", "q_start", "q_end",
                "match_fasta", "match_contig", "m_start", "m_end",
                "jaccard", "containment_q_in_m", "containment_m_in_q",
            ]
            with open(_sw_path, "w", newline="") as _sw_out:
                _sw_writer = csv.writer(_sw_out)
                _sw_writer.writerow(_sw_header)
                for (ref_a, ref_b, shared_count) in _rs_pairs:
                    # jaccard placeholder = shared / (total_a + total_b - shared)
                    cnt_a = max(1, _ref_counts_pre.get(ref_a, 1))
                    cnt_b = max(1, _ref_counts_pre.get(ref_b, 1))
                    _j = shared_count / (cnt_a + cnt_b - shared_count)
                    _sw_writer.writerow(["bam", ref_a, 0, 0, "bam", ref_b, 0, 0,
                                         f"{_j:.6f}", f"{_j:.6f}", f"{_j:.6f}"])
            print(f"  {len(_rs_pairs)} competing reference pair(s) from "
                  f"{len(_rs_map)} reference(s) → {_sw_path}")
        else:
            print("  No multi-mapped reads found across references; "
                  "shared_windows_report.csv not written.")

    if compare_to_reference_windows:
        print(f"Shared-window report path: {report_path}")
        if os.path.exists(report_path) and os.path.getsize(report_path) > 0:
            print(f"Reusing cached shared-window report: {report_path}")
        else:
            print("Building shared-window report across FASTAs...")

            if find_optimal_windows:
                success, best_params, _, _ = tune_shared_window_params(
                    bam_path=input_bam,
                    fasta_files=fasta_files,
                    attempts=28,
                    allowed_tp_loss_frac=0.02,
                    alpha=1.0,
                    min_jaccard=sim_ani_threshold,
                    tmp_dir=os.path.join(output_dir, "sbt_tune_tmp"),
                )
                print("Tuner:", success, best_params)
                if success and best_params:
                    report_shared_windows_across_fastas(
                        fasta_files=fasta_files,
                        output_csv=report_path,
                        ksize=best_params["ksize"],
                        scaled=best_params["scaled"],
                        window=best_params["window"],
                        step=best_params["step"],
                        jaccard_threshold=sim_ani_threshold,
                        max_hits_per_query=120,
                        skip_self_same_fasta=False,
                        aligned_refs=_aligned_refs,
                    )
                else:
                    report_shared_windows_across_fastas(
                        fasta_files=fasta_files,
                        output_csv=report_path,
                        ksize=21,
                        scaled=500,
                        window=window_size,
                        step=step_size,
                        jaccard_threshold=sim_ani_threshold,
                        max_hits_per_query=120,
                        skip_self_same_fasta=False,
                        aligned_refs=_aligned_refs,
                    )
            else:
                print("Creating shared FASTA report from scratch")
                # ── Split-threshold strategy ──────────────────────────────────
                # We use two different Jaccard thresholds:
                #
                #   _write_j_thresh  (0.05) — written to the CSV.
                #     Captures cross-genus pairs like E. coli O104:H4 vs Shigella
                #     (ANI ≈ 95%, ksize=31 → J ≈ 0.11).  These pairs are used
                #     ONLY for building the union-find conflict cluster IDs in
                #     match_paths.py, NOT for actual read removal.
                #
                #   _removal_j_thresh (0.20) — used when loading shared_idx for
                #     read-removal decisions.  Only strong window similarity
                #     (intra-Shigella, intra-Klebsiella, etc.) drives removal;
                #     weaker cross-genus links do not.
                #
                # Without this split, organisms separated by ~95% ANI were never
                # placed in the same conflict cluster, so the winner/loser minhash
                # adjustment in optimize_weights.py had no effect between them.
                _write_j_thresh   = 0.05   # write threshold: capture cross-genus links
                _removal_j_thresh = 0.20   # removal threshold: only strong similarity
                sim_ani_threshold = _removal_j_thresh  # used below for shared_idx load
                report_shared_windows_across_fastas(
                    fasta_files=fasta_files,
                    output_csv=report_path,
                    ksize=kmer_size,
                    scaled=scaled,
                    window=window_size,
                    step=step_size,
                    jaccard_threshold=_write_j_thresh,   # low threshold: write more pairs
                    max_hits_per_query=120,
                    skip_self_same_fasta=False,
                    aligned_refs=_aligned_refs,
                )

        # ── BAM supplement: append multi-mapper pairs to shared_windows_report ──
        # Even after the FASTA window comparison, some cross-genus reference
        # pairs may be missed if their window Jaccard falls below _write_j_thresh
        # (e.g., pairs with ANI < 93%).  Reads that multi-map in the BAM provide
        # ground-truth evidence of reference ambiguity.  We scan the BAM and
        # append any new (ref_a, ref_b) pairs not already in the report — with
        # nominal j=0.50 and positions=0 — so the union-find in match_paths.py
        # can group them correctly.  These rows have no effect on read removal
        # (which uses window positions from the FASTA comparison).
        if os.path.exists(report_path) and os.path.getsize(report_path) > 0:
            print("Supplementing shared-window report with BAM multi-mapper pairs...")
            _bam_rs_map, _bam_rs_pairs = _build_read_sharing_clusters(
                input_bam,
                aligned_refs=_aligned_refs,
            )
            if _bam_rs_pairs:
                # Load existing pairs to avoid duplicates
                _existing_pairs: set = set()
                try:
                    with open(report_path, newline="") as _efp:
                        _er = csv.DictReader(_efp)
                        for _erow in _er:
                            _qa = str(_erow.get("query_contig", ""))
                            _ma = str(_erow.get("match_contig", ""))
                            if _qa and _ma:
                                _existing_pairs.add((_qa, _ma))
                                _existing_pairs.add((_ma, _qa))
                except Exception as _e_read:
                    print(f"  Warning: could not pre-scan existing pairs: {_e_read}")

                _new_pairs = [
                    (ra, rb, cnt) for ra, rb, cnt in _bam_rs_pairs
                    if (ra, rb) not in _existing_pairs
                ]
                if _new_pairs:
                    _sw_header = [
                        "query_fasta", "query_contig", "q_start", "q_end",
                        "match_fasta", "match_contig", "m_start", "m_end",
                        "jaccard", "containment_q_in_m", "containment_m_in_q",
                    ]
                    # Check if header row already exists
                    _header_exists = os.path.getsize(report_path) > 0
                    with open(report_path, "a", newline="") as _sw_app:
                        _sw_writer = csv.writer(_sw_app)
                        for (ref_a, ref_b, shared_count) in _new_pairs:
                            cnt_a = max(1, _ref_counts_pre.get(ref_a, 1))
                            cnt_b = max(1, _ref_counts_pre.get(ref_b, 1))
                            _j = shared_count / (cnt_a + cnt_b - shared_count)
                            # Use a nominal jaccard floor of 0.10 so these
                            # rows pass any downstream min_jaccard filter.
                            _j_out = max(0.10, _j)
                            _sw_writer.writerow([
                                "bam", ref_a, 0, 0,
                                "bam", ref_b, 0, 0,
                                f"{_j_out:.6f}", f"{_j_out:.6f}", f"{_j_out:.6f}",
                            ])
                    print(f"  Appended {len(_new_pairs)} BAM-derived pair(s) "
                          f"to {report_path}")
                else:
                    print("  BAM multi-mapper pairs already present in report; "
                          "no supplement needed.")
            else:
                print("  No BAM multi-mapper pairs found; supplement skipped.")

        shared_idx = load_shared_windows_csv(report_path, min_jaccard=sim_ani_threshold, skip_same_contig=True)

    # Removal plan
    if compare_to_reference_windows and shared_idx is not None:
        print("Building removal plan based on shared-window alignments...")
        # _ref_counts_pre already computed above for aligned_refs filtering;
        # reuse it here for ANI-weighted dominance scoring.
        print(f"  Pre-counted reads for {len(_ref_counts_pre)} references (used for ANI-weighted dominance scoring)")
        removed_read_ids = build_removed_ids_best_alignment(
            bam_path=input_bam,
            shared_idx=shared_idx,
            penalize_weight=0,
            as_weight=1.0,
            drop_contigs=set(),
            drop_if_ambiguous=True,
            min_alt_count=2,
            only_primary=False,
            ref_read_counts=_ref_counts_pre,
            ani_boost_weight=10.0,
        )
    else:
        # Parse bedgraph + merge regions
        regions = parse_bed_file(bedfile)
        print(f"Input intervals selected from bam/bed: {len(regions)}")

        if use_variance:
            stat_name = "variance"
            threshold = 0.8 if jump_threshold is None else float(jump_threshold)
        else:
            stat_name = "jump"
            threshold = 1.0 if jump_threshold is None else float(jump_threshold)
        merged_regions = merge_bedgraph_regions(
            regions,
            merging_method=stat_name,
            max_stat_threshold=threshold,
            max_group_size=4_000_000,
            reflengths=reflengths,
            gap_allowance=gap_allowance,
        )
        print(f"Merged regions: {len(merged_regions)} (from {len(regions)}) in {time.time()-t0:.2f}s")

        # Check for cached region_comparisons CSV — skip signature + SBT if it exists
        output_csv = report_path  # uses the resolved path (sigfile CSV or default)
        if os.path.exists(output_csv) and os.path.getsize(output_csv) > 0:
            print(f"Reusing cached region comparisons: {output_csv}")
            sum_comparisons = load_region_comparisons_csv(output_csv)
        else:
            signatures = {}
            # Signatures: load or generate (BAM-only — no FASTA needed)
            if not sigfile or not os.path.exists(sigfile):
                nworkers = cpu_count if cpu_count else max(1, int(os.cpu_count() / 2))
                print(f"Sketching merged regions from BAM only (workers={nworkers})")
                t0 = time.time()
                signatures = create_signatures_from_bam(
                    regions_df=merged_regions,
                    bam_path=input_bam,
                    kmer_size=kmer_size,
                    scaled=scaled,
                    num_workers=nworkers,
                )
                print(f"Signatures: {len(signatures)} in {time.time()-t0:.2f}s")

                sig_dir = os.path.join(output_dir, "signatures")
                single_sigfile = os.path.join(sig_dir, "merged_regions.sig")
                save_signatures_sourmash(signatures, single_sigfile)
            else:
                print(f"Loading signatures from: {sigfile}")
                signatures = rebuild_sig_dict(load_signatures_sourmash(sigfile))

            # Search via parallel SBT (scales to millions of signatures)
            sig_items = list(signatures.items())

            nworkers = cpu_count if cpu_count else max(1, int(os.cpu_count() / 2))
            print(f"Searching regions via parallel SBT (workers={nworkers})...")
            sum_comparisons = search_sbt_parallel(
                sig_items, output_csv,
                min_threshold=min_threshold,
                n_jobs=nworkers,
                kmer_size=kmer_size,
                scaled=scaled,
            )

        # Build conflict groups
        print("Building conflict groups...")
        start_time = time.time()
        conflict_groups = build_conflict_groups(sum_comparisons, min_jaccard=0.0)
        print(f"Conflict groups: {len(conflict_groups)}")
        print(f"Conflict groups length: {len(conflict_groups)} built in {time.time() - start_time:.2f} seconds. Next up is proportion removal")
        start_time = time.time()
        # 2) For each group, dominance-aware read removal
        removed_read_ids, cluster_dominance = finalize_proportional_removal(
            conflict_groups,
            bam_fs,
            fetch_reads_in_region,
            remove_mode='random',
            dominance_protect_ratio=dominance_protect_ratio,
        )
        # Write cluster_dominance JSON for downstream minhash penalty
        import json as _json
        _cd_path = os.path.join(output_dir, "cluster_dominance.json")
        with open(_cd_path, "w") as _cd_fh:
            _json.dump(cluster_dominance, _cd_fh, indent=2)
        print(f"Wrote cluster dominance: {_cd_path} ({len(cluster_dominance)} references)")

    # Write removals table
    failed_path = os.path.join(output_dir, "failed_reads.txt")
    with open(failed_path, "w") as fp:
        for read_id, refs in removed_read_ids.items():
            for ref in refs:
                fp.write(f"{ref}\t{read_id}\n")
    print(f"Wrote removals: {failed_path}")

    # ── Single-pass: breadth (old+new), per-ref accounting, removal stats ──
    # Previously this was 4 separate full BAM scans. Now one pass per reference.
    print(f"Computing breadth + per-ref stats (single pass): {time.ctime()}")
    t_stats = time.time()
    per_ref = defaultdict(lambda: defaultdict(int))
    breadth_old = {}
    breadth_new = {}

    # Pre-convert removed_read_ids values to sets for O(1) lookup
    _removed_sets = {rid: set(refs) for rid, refs in removed_read_ids.items()}

    total_reads_all = set()
    total_alns_all = 0
    removed_alns_all = 0

    bam_fs.reset()
    for ref, L in reflengths.items():
        total = 0
        passed = 0

        # Breadth tracking: two interval merge streams (old=all, new=kept)
        old_s = old_e = None
        new_s = new_e = None
        covered_old = 0
        covered_new = 0

        for r in bam_fs.fetch(ref, 0, L):
            if r.is_unmapped or r.reference_start is None or r.reference_end is None:
                continue

            total += 1
            total_alns_all += 1
            rid = r.query_name
            total_reads_all.add(rid)

            s = max(0, int(r.reference_start))
            e = min(L, int(r.reference_end))

            # -- breadth OLD (all reads) --
            if s < e:
                if old_s is None:
                    old_s, old_e = s, e
                elif s <= old_e + 1:
                    old_e = max(old_e, e)
                else:
                    covered_old += (old_e - old_s)
                    old_s, old_e = s, e

            is_removed = rid in _removed_sets and ref in _removed_sets[rid]
            if is_removed:
                removed_alns_all += 1

            # -- breadth NEW (kept reads) --
            if not is_removed and s < e:
                if new_s is None:
                    new_s, new_e = s, e
                elif s <= new_e + 1:
                    new_e = max(new_e, e)
                else:
                    covered_new += (new_e - new_s)
                    new_s, new_e = s, e

            # -- GT accounting --
            gt = infer_gt_from_readname(rid)
            if gt is not None:
                if gt == r.reference_name:
                    per_ref[ref]["TP Original"] += 1
                else:
                    per_ref[ref]["FP Original"] += 1

                if not is_removed:
                    passed += 1
                    if gt == r.reference_name:
                        per_ref[ref]["TP New"] += 1
                    else:
                        per_ref[ref]["FP New"] += 1
                else:
                    if gt == r.reference_name:
                        per_ref[ref]["FN New"] += 1
            else:
                if not is_removed:
                    passed += 1

        # Flush last intervals
        if old_s is not None:
            covered_old += (old_e - old_s)
        if new_s is not None:
            covered_new += (new_e - new_s)

        breadth_old[ref] = (100.0 * covered_old / L) if L > 0 else 0.0
        breadth_new[ref] = (100.0 * covered_new / L) if L > 0 else 0.0

        per_ref[ref]["total_reads"] = total
        per_ref[ref]["pass_filtered_reads"] = passed

        per_ref[ref]["proportion_aligned"] = (per_ref[ref]["TP Original"] / total) if total else 0.0
        denom = per_ref[ref]["TP Original"] + per_ref[ref]["FP Original"]
        per_ref[ref]["precision"] = (per_ref[ref]["TP Original"] / denom) if denom else 0.0
        per_ref[ref]["recall"] = (per_ref[ref]["TP Original"] / total) if total else 0.0
        pr = per_ref[ref]["precision"] + per_ref[ref]["recall"]
        per_ref[ref]["f1"] = (2 * per_ref[ref]["precision"] * per_ref[ref]["recall"] / pr) if pr else 0.0

        per_ref[ref]["breadth_old"] = breadth_old.get(ref, 0.0)
        per_ref[ref]["breadth"] = breadth_new.get(ref, 0.0)

    # Print removal summary (replaces separate report_removed_read_stats call)
    n_reads = len(total_reads_all)
    n_removed_reads = len(_removed_sets)
    print(f"\n=== Removal Summary ===")
    print(f"Unique reads:             {n_reads}")
    print(f"Reads with ≥1 removal:    {n_removed_reads} ({100.0*n_removed_reads/max(1,n_reads):.2f}%)")
    print(f"Total alignments:         {total_alns_all}")
    print(f"Alignments removed:       {removed_alns_all} ({100.0*removed_alns_all/max(1,total_alns_all):.2f}%)")
    print(f"Single-pass stats completed in {time.time()-t_stats:.2f}s")

    # Attach shared-window conflict stats to per_ref (populated in xlsx as informational columns)
    if shared_idx is not None:
        _sw_stats = compute_shared_window_stats(shared_idx)
        for _ref in list(per_ref.keys()):
            _sw = _sw_stats.get(_ref, {})
            per_ref[_ref]["n_shared_windows"]   = _sw.get("n_shared_windows", 0)
            per_ref[_ref]["n_conflicting_refs"]  = _sw.get("n_conflicting_refs", 0)
            per_ref[_ref]["mean_window_jaccard"] = _sw.get("mean_window_jaccard", 0.0)
            per_ref[_ref]["max_window_jaccard"]  = _sw.get("max_window_jaccard", 0.0)
            per_ref[_ref]["shared_bp"]           = _sw.get("shared_bp", 0)

    comparison_df = compare_metrics(per_ref, reflengths, accession_to_taxid=accession_to_taxid,
                                    taxid_to_name=taxid_to_name)


    try:
        out_xlsx = os.path.join(output_dir, "removal_stats.xlsx")
        # add a "Total" row that is the sum of all the TP/FP/FN columns, empty for others
        total_row = {
            "Reference": "Total",
            "Taxid": "",
            "Organism": "",
            "TP Original": comparison_df["TP Original"].sum(),
            "FP Original": comparison_df["FP Original"].sum(),
            "FN Original": comparison_df["FN Original"].sum(),
            "TP New": comparison_df["TP New"].sum(),
            "FP New": comparison_df["FP New"].sum(),
            "FN New": comparison_df["FN New"].sum(),
            "Total Reads": comparison_df["Total Reads"].sum(),
            "Pass Filtered Reads": comparison_df["Pass Filtered Reads"].sum(),
            "Proportion Aligned": "",
            "Precision": "",
            "Recall": "",
            "F1": "",
            "Δ All": comparison_df["Δ All"].sum(),
            "Δ All%": "",
            "Δ^-1 Breadth": "",
            "Breadth New": "",
            "Breadth Original": "",
            # RPKM/RPM are rates — not meaningful to sum
            "RPKM Pre": "",
            "RPKM Post": "",
            "RPM Pre": "",
            "RPM Post": "",
            # Shared-window counts can be summed (with caveats on double-counting)
            "Shared Windows": comparison_df["Shared Windows"].sum() if "Shared Windows" in comparison_df.columns else "",
            "Conflicting Refs": "",
            "Mean Shared ANI": "",
            "Max Shared ANI": "",
            "Shared BP": comparison_df["Shared BP"].sum() if "Shared BP" in comparison_df.columns else "",
        }
        comparison_df = pd.concat([comparison_df, pd.DataFrame([total_row])], ignore_index=True)
        # put Total at Top of excel sheet
        comparison_df = pd.concat([comparison_df.tail(1), comparison_df.iloc[:-1]], ignore_index=True)

        # ── Build Conflict Groups sheet ───────────────────────────────────
        # Summarise which accession pairs share windows and who won (most reads).
        # This gives users an explicit view of which organisms competed for reads
        # and how the minhash winner was determined, which the per-accession
        # Removal Stats sheet alone cannot convey.
        conflict_groups_rows = []
        if shared_idx is not None and len(shared_idx) > 0:
            # Build lookup: accession → organism name and read counts from comparison_df
            _cdf_data = comparison_df[comparison_df["Reference"] != "Total"].set_index("Reference")
            _seen_pairs: set = set()
            group_id = 1

            for ref_a, windows in shared_idx.items():
                alt_contigs = {w.alt_contig for w in windows}
                for ref_b in sorted(alt_contigs):
                    pair_key = tuple(sorted([ref_a, ref_b]))
                    if pair_key in _seen_pairs:
                        continue
                    _seen_pairs.add(pair_key)

                    # Read counts from comparison_df (fall back to per_ref dict)
                    def _reads(ref):
                        if ref in _cdf_data.index:
                            return int(_cdf_data.loc[ref, "Total Reads"] or 0)
                        return int(per_ref.get(ref, {}).get("numreads", 0))

                    def _org(ref):
                        if ref in _cdf_data.index:
                            return str(_cdf_data.loc[ref, "Organism"] or "")
                        return ""

                    def _delta(ref):
                        if ref in _cdf_data.index:
                            return int(_cdf_data.loc[ref, "Δ All"] or 0)
                        return 0

                    reads_a = _reads(ref_a)
                    reads_b = _reads(ref_b)

                    # Shared-window stats between this pair
                    pair_windows = [w for w in windows if w.alt_contig == ref_b]
                    n_shared = len(pair_windows)
                    mean_jaccard = (sum(w.jaccard for w in pair_windows) / n_shared) if n_shared else 0.0
                    max_jaccard  = max((w.jaccard for w in pair_windows), default=0.0)

                    winner = ref_a if reads_a >= reads_b else ref_b
                    loser  = ref_b if winner == ref_a else ref_a

                    conflict_groups_rows.append({
                        "Group ID":         group_id,
                        "Winner Accession": winner,
                        "Winner Organism":  _org(winner),
                        "Winner Reads":     _reads(winner),
                        "Winner Δ All":     _delta(winner),
                        "Loser Accession":  loser,
                        "Loser Organism":   _org(loser),
                        "Loser Reads":      _reads(loser),
                        "Loser Δ All":      _delta(loser),
                        "Shared Windows":   n_shared,
                        "Mean Jaccard (ANI proxy)": round(mean_jaccard, 4),
                        "Max Jaccard":      round(max_jaccard, 4),
                    })
                    group_id += 1

        conflict_groups_df = pd.DataFrame(conflict_groups_rows) if conflict_groups_rows else pd.DataFrame(
            columns=["Group ID", "Winner Accession", "Winner Organism", "Winner Reads",
                     "Winner Δ All", "Loser Accession", "Loser Organism", "Loser Reads",
                     "Loser Δ All", "Shared Windows", "Mean Jaccard (ANI proxy)", "Max Jaccard"]
        )
        if conflict_groups_df.empty:
            conflict_groups_df = pd.DataFrame([{
                "Note": (
                    "No shared-window pairs detected. This sheet is populated when FASTA files are "
                    "provided (--fasta) so that sourmash window sketching can detect cross-accession "
                    "sequence similarity. For BAM-level multi-mapper competition (e.g. E. coli vs "
                    "Shigella) see the Removal Stats sheet: accessions with high Δ All% competed "
                    "for reads even without shared windows being detected at the sketch level."
                )
            }])

        # Write multi-sheet xlsx
        with pd.ExcelWriter(out_xlsx, engine="openpyxl") as _xw:
            comparison_df.to_excel(_xw, sheet_name="Removal Stats", index=False)
            conflict_groups_df.to_excel(_xw, sheet_name="Conflict Groups", index=False)

        print(f"Wrote: {out_xlsx}")
        print(
            comparison_df[comparison_df["Δ All"] != 0][
                ["Reference", "Δ All%", "Δ^-1 Breadth", "Breadth New", "Breadth Original", "TP New", "TP Original", "FP Original", "FP New"]
            ].to_string(index=False)
        )
    except Exception as e:
        print(f"Stats export failed: {e}")

    # Optional: taxid-aggregated removal stats
    if taxid_removal_stats and "Taxid" in comparison_df.columns:
        try:
            print("trying to match taxids in removal reads xlsx output")
            # Work on copy without the Total row
            _tdf = comparison_df[comparison_df["Reference"] != "Total"].copy()
            _tdf["Taxid"] = _tdf["Taxid"].astype(str)
            # Only aggregate rows that have a non-empty taxid
            _has_taxid = _tdf[_tdf["Taxid"].str.strip().ne("")].copy()
            if not _has_taxid.empty:
                _sum_cols = ["TP Original", "FP Original", "FN Original",
                             "TP New", "FP New", "FN New",
                             "Total Reads", "Pass Filtered Reads", "Δ All",
                             "Reference Length"]
                _grouped = _has_taxid.groupby("Taxid", sort=True)
                taxid_agg = _grouped[_sum_cols].sum()
                # Weighted-average for breadth columns (weighted by Reference Length)
                for bcol in ["Breadth Original", "Breadth New"]:
                    if bcol in _has_taxid.columns:
                        _has_taxid[f"_w_{bcol}"] = _has_taxid[bcol] * _has_taxid["Reference Length"]
                        taxid_agg[bcol] = _grouped[f"_w_{bcol}"].sum() / taxid_agg["Reference Length"].replace(0, 1)
                # Recompute derived columns
                taxid_agg["Δ All%"] = (100.0 * taxid_agg["Δ All"] / taxid_agg["Total Reads"]).fillna(0.0)
                taxid_agg["Δ^-1 Breadth"] = (taxid_agg["Breadth New"] / taxid_agg["Breadth Original"].replace(0, float("nan"))).fillna(0.0)
                taxid_agg["Δ Breadth"] = taxid_agg["Breadth New"] - taxid_agg["Breadth Original"]
                denom = taxid_agg["TP Original"] + taxid_agg["FP Original"]
                taxid_agg["Precision"] = (taxid_agg["TP Original"] / denom.replace(0, float("nan"))).fillna(0.0)
                taxid_agg["Recall"] = (taxid_agg["TP Original"] / taxid_agg["Total Reads"].replace(0, float("nan"))).fillna(0.0)
                pr = taxid_agg["Precision"] + taxid_agg["Recall"]
                taxid_agg["F1"] = (2 * taxid_agg["Precision"] * taxid_agg["Recall"] / pr.replace(0, float("nan"))).fillna(0.0)
                taxid_agg["Proportion Aligned"] = (taxid_agg["TP Original"] / taxid_agg["Total Reads"].replace(0, float("nan"))).fillna(0.0)
                # Collect accessions per taxid for reference
                taxid_agg["Accessions"] = _grouped["Reference"].apply(lambda x: "; ".join(sorted(x)))
                taxid_agg = taxid_agg.reset_index()
                # Map organism name from taxid
                if taxid_to_name:
                    taxid_agg["Organism"] = taxid_agg["Taxid"].map(
                        lambda t: taxid_to_name.get(str(t), ""))
                else:
                    taxid_agg["Organism"] = ""
                # Reorder columns
                col_order = ["Taxid", "Organism", "Accessions", "Reference Length",
                             "TP Original", "FP Original", "FN Original",
                             "TP New", "FP New", "FN New",
                             "Total Reads", "Pass Filtered Reads",
                             "Proportion Aligned", "Precision", "Recall", "F1",
                             "Δ All", "Δ All%",
                             "Breadth Original", "Breadth New", "Δ Breadth", "Δ^-1 Breadth"]
                col_order = [c for c in col_order if c in taxid_agg.columns]
                taxid_agg = taxid_agg[col_order]
                # Add Total row
                taxid_total = {"Taxid": "Total", "Organism": "", "Accessions": ""}
                for c in _sum_cols:
                    if c in taxid_agg.columns:
                        taxid_total[c] = taxid_agg[c].sum()
                taxid_total["Δ All%"] = ""
                taxid_total["Precision"] = ""
                taxid_total["Recall"] = ""
                taxid_total["F1"] = ""
                taxid_total["Proportion Aligned"] = ""
                taxid_total["Δ^-1 Breadth"] = ""
                taxid_total["Δ Breadth"] = ""
                taxid_total["Breadth Original"] = ""
                taxid_total["Breadth New"] = ""
                taxid_agg = pd.concat([pd.DataFrame([taxid_total]), taxid_agg], ignore_index=True)
                out_taxid_xlsx = os.path.join(output_dir, "removal_stats_by_taxid.xlsx")
                taxid_agg.to_excel(out_taxid_xlsx, index=False)
                print(f"Wrote: {out_taxid_xlsx}")
            else:
                print("Skipping taxid removal stats: no taxid mappings found in comparison data.")
        except Exception as e:
            print(f"Taxid removal stats export failed: {e}")

    # Optional: create filtered BAM
    if filtered_bam_create:
        try:
            create_filtered_bam(bam_fs, filtered_bam_create, removed_read_ids)
        except Exception as e:
            print(f"Filtered BAM failed: {e}")
    bam_fs.close()

    return removed_read_ids, comparison_df
