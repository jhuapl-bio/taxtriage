#!/usr/bin/env python3
"""
conflict_regions.py

Minimal cleaned “slate” version of your script.

Goals of this refactor:
- Remove duplicate imports / dead code / globals.
- Separate responsibilities (I/O, sketching, comparisons, removal).
- Provide a small, readable, testable layout you can extend.
- Keep only the essential plumbing + clear TODO hooks.

NOTE:
This is intentionally a *minimal scaffold*—it will run, but many “business logic”
pieces are stubbed with TODOs so you can plug back the exact algorithms you want.
"""

from __future__ import annotations

import argparse
import csv
import logging
import math
import os
import random
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import pysam
from sourmash import MinHash, SourmashSignature, load_file_as_signatures, save_signatures
from sourmash.sbt import SBT, GraphFactory
from sourmash.sbtmh import SigLeaf
from tqdm import tqdm
from utils import load_matchfile

# If you have this in your repo, keep it:
# from utils import load_matchfile


# -----------------------------------------------------------------------------
# Logging
# -----------------------------------------------------------------------------

log = logging.getLogger("conflict_regions")


def setup_logging(verbosity: int) -> None:
    level = logging.WARNING
    if verbosity == 1:
        level = logging.INFO
    elif verbosity >= 2:
        level = logging.DEBUG
    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    )


# -----------------------------------------------------------------------------
# Data structures
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class Region:
    Reference: str
    start: int  # 0-based, inclusive
    end: int    # 0-based, exclusive
    mean_depth: float = 0.0

@dataclass(frozen=True)
class SharedWindow:
    contig: str
    start: int
    end: int
    alt_contig: str
    alt_start: int
    alt_end: int
    jaccard: float


def iter_windows(length: int, window: int, step: int):
    """
    Yields (start,end) half-open windows across [0,length).
    Guarantees end <= length.
    """
    if window <= 0 or step <= 0:
        raise ValueError("window and step must be positive integers")
    if length <= 0:
        return
    s = 0
    while s < length:
        e = min(s + window, length)
        if e > s:
            yield s, e
        if s + step == s:  # paranoia
            break
        s += step


def make_window_regions_from_fastas(
    fasta_paths: Sequence[str],
    *,
    window: int,
    step: int,
    contigs: Optional[Sequence[str]] = None,
    min_len: Optional[int] = None,
) -> List[Region]:
    """
    Returns a list of Region(Reference,start,end) derived ONLY from FASTA contig lengths.
    No comparisons, no signatures, no bedgraph.
    """
    if not fasta_paths:
        raise ValueError("fasta_paths is required to generate window regions")

    regions: List[Region] = []
    # We only need lengths, so any one fasta that contains contig works;
    # but if you have multiple fastas and overlapping contig names, this will generate windows per contig per fasta file.
    # If you instead want a single set of contigs, see note below.
    for fp in fasta_paths:
        fa = pysam.FastaFile(fp)
        try:
            contig_list = list(contigs) if contigs is not None else list(fa.references)
            for ctg in contig_list:
                if ctg not in fa.references:
                    continue
                L = int(fa.get_reference_length(ctg))
                if min_len is not None and L < min_len:
                    continue
                for s, e in iter_windows(L, window, step):
                    regions.append(Region(ctg, s, e))
        finally:
            fa.close()

    return regions


def window_regions_df(regions: Iterable[Region]) -> pd.DataFrame:
    """
    Convenience: convert to the same shape as your merged bed regions DataFrame.
    """
    rows = [{"chrom": r.Reference, "start": r.start, "end": r.end} for r in regions]
    return pd.DataFrame(rows)
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# BED/BEDGRAPH parsing + merging
# -----------------------------------------------------------------------------

def parse_bedgraph(path: str) -> pd.DataFrame:
    """
    Expected columns: Reference, start, end, depth
    """
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["Reference", "start", "end", "depth"],
        dtype={"Reference": "category", "start": "int32", "end": "int32", "depth": "float32"},
    )
    return df


def compute_gini(values: Sequence[float]) -> float:
    """Gini coefficient in [0,1] (0 = perfect equality)."""
    if not values:
        return 0.0
    xs = sorted(float(v) for v in values)
    n = len(xs)
    s = sum(xs)
    if s == 0:
        return 0.0
    cum = 0.0
    cum_sum = 0.0
    for x in xs:
        cum += x
        cum_sum += cum
    return (2 * cum_sum) / (n * s) - (n + 1) / n


def merge_bedgraph_regions(
    intervals: pd.DataFrame,
    *,
    method: str = "jump",                  # "jump" | "variance" | "gini"
    stat_threshold: Optional[float] = None,
    max_group_size: int = 4_000_000,
    max_length: Optional[int] = None,
    gap_allowance_frac: float = 0.1,       # allowed gap = ref_len * frac
    reflengths: Optional[Dict[str, int]] = None,
) -> pd.DataFrame:
    """
    Returns DataFrame with columns: Reference, start, end, mean_depth

    This is a simplified/cleaned version of your merge logic.
    """
    if intervals.empty:
        return pd.DataFrame(columns=["Reference", "start", "end", "mean_depth"])

    intervals = intervals.sort_values(["Reference", "start"]).reset_index(drop=True)
    intervals["depth"] = intervals["depth"].fillna(0)

    # Determine per-Reference allowed gap
    if reflengths:
        allowed_gap = {c: int(reflengths.get(str(c), 0) * gap_allowance_frac) for c in intervals["Reference"].unique()}
    else:
        allowed_gap = {c: math.inf for c in intervals["Reference"].unique()}

    # Default jump threshold if none provided
    if method == "jump" and stat_threshold is None:
        try:
            q = intervals.groupby("Reference", observed=True)["depth"].quantile(0.975).mean()
            stat_threshold = float(math.ceil(q) + 1)
        except Exception:
            stat_threshold = 200.0

    merged: List[Dict[str, object]] = []

    for Reference, grp in intervals.groupby("Reference", observed=True):
        grp = grp.reset_index(drop=True)
        if grp.empty:
            continue

        buf_start = int(grp.at[0, "start"])
        buf_end = int(grp.at[0, "end"])
        buf_depths: List[float] = [float(grp.at[0, "depth"])]
        per_gap = allowed_gap.get(Reference, math.inf)

        for i in range(1, len(grp)):
            new_s = int(grp.at[i, "start"])
            new_e = int(grp.at[i, "end"])
            new_d = float(grp.at[i, "depth"])

            gap = new_s - buf_end
            can_merge = True

            if len(buf_depths) >= max_group_size:
                can_merge = False
            if max_length is not None and (new_e - buf_start) > max_length:
                can_merge = False
            if gap > per_gap:
                can_merge = False

            if can_merge:
                if method == "jump":
                    jump = abs(new_d - buf_depths[-1])
                    if stat_threshold is not None and jump > stat_threshold:
                        can_merge = False
                elif method == "variance":
                    if stat_threshold is not None:
                        v = float(np.var(buf_depths + [new_d]))
                        if v > stat_threshold:
                            can_merge = False
                elif method == "gini":
                    if stat_threshold is not None:
                        g = compute_gini(buf_depths + [new_d])
                        if g > stat_threshold:
                            can_merge = False
                else:
                    raise ValueError(f"Unknown merge method: {method}")

            if can_merge:
                buf_end = new_e
                buf_depths.append(new_d)
            else:
                merged.append(
                    {"Reference": str(Reference), "start": buf_start, "end": buf_end, "mean_depth": float(np.mean(buf_depths))}
                )
                buf_start, buf_end, buf_depths = new_s, new_e, [new_d]

        merged.append({"Reference": str(Reference), "start": buf_start, "end": buf_end, "mean_depth": float(np.mean(buf_depths))})

    return pd.DataFrame(merged)


# -----------------------------------------------------------------------------
# Signature creation
# -----------------------------------------------------------------------------

def region_id(Reference: str, start: int, end: int) -> str:
    # keep your convention (note: many parts of your code assume ref:start-end)
    return f"{Reference}:{start}-{end}"


def build_region_signatures_from_fasta(
    regions: Iterable[Region],
    fasta_paths: Sequence[str],
    *,
    ksize: int,
    scaled: int,
) -> Dict[str, SourmashSignature]:
    """
    Prefer FASTA-based signatures when available (faster & deterministic).
    Supports multiple FASTA paths; first one containing the contig is used.
    """
    fastas = []
    for fp in fasta_paths:
        fastas.append(pysam.FastaFile(fp))

    out: Dict[str, SourmashSignature] = {}
    try:
        for r in regions:
            seq = None
            for fa in fastas:
                if r.Reference in fa.references:
                    seq = fa.fetch(r.Reference, r.start, r.end).upper()
                    break
            if not seq or len(seq) < ksize:
                continue

            mh = MinHash(n=0, ksize=ksize, scaled=scaled)
            mh.add_sequence(seq, force=True)
            rid = region_id(r.Reference, r.start, r.end)
            out[rid] = SourmashSignature(mh, name=rid)
    finally:
        for fa in fastas:
            fa.close()

    return out


def build_region_signatures_from_bam(
    regions: Iterable[Region],
    bam_path: str,
    *,
    ksize: int,
    scaled: int,
    only_primary: bool = True,
) -> Dict[str, SourmashSignature]:
    """
    Build signatures by sketching *read sequences* overlapping each region.
    NOTE: This is slow for lots of regions. Use FASTA if possible.
    """
    out: Dict[str, SourmashSignature] = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for r in tqdm(list(regions), desc="Sketching regions from BAM", miniters=1000):
            mh = MinHash(n=0, ksize=ksize, scaled=scaled)

            try:
                for aln in bam.fetch(r.Reference, r.start, r.end):
                    if aln.is_unmapped:
                        continue
                    if only_primary and (aln.is_secondary or aln.is_supplementary):
                        continue
                    seq = aln.query_sequence
                    if not seq:
                        continue
                    mh.add_sequence(seq, force=True)
            except ValueError:
                # contig not in BAM header, etc.
                continue

            rid = region_id(r.Reference, r.start, r.end)
            out[rid] = SourmashSignature(mh, name=rid)

    return out


def save_signatures_file(signatures: Dict[str, SourmashSignature], out_sig: str) -> None:
    os.makedirs(os.path.dirname(out_sig), exist_ok=True)
    with open(out_sig, "wt") as fp:
        save_signatures(list(signatures.values()), fp)
    log.info("Wrote signatures: %s (%d)", out_sig, len(signatures))


def load_signatures_file(sig_path: str) -> Dict[str, SourmashSignature]:
    sigs = list(load_file_as_signatures(sig_path))
    return {s.name: s for s in sigs}


# -----------------------------------------------------------------------------
# Comparison (SBT)
# -----------------------------------------------------------------------------

def build_sbt_index(
    sigs: Dict[str, SourmashSignature],
    *,
    ksize: int,
) -> SBT:
    factory = GraphFactory(ksize=ksize, n_tables=1, starting_size=max(1, len(sigs)))
    sbt = SBT(factory)
    for name, sig in sigs.items():
        sbt.add_node(SigLeaf(name, sig))
    return sbt


def compare_signatures_sbt(
    sigs: Dict[str, SourmashSignature],
    *,
    ksize: int,
    min_jaccard: float,
    out_csv: str,
) -> Dict[str, List[dict]]:
    """
    Writes a region_comparisons.csv and returns:
      sum_comparisons[ref] -> list[{jaccard,to,s1,e1,s2,e2}]
    """
    sbt = build_sbt_index(sigs, ksize=ksize)

    def parse_region(rid: str) -> Tuple[str, int, int]:
        ref, coords = rid.split(":", 1)
        s, e = coords.split("-", 1)
        return ref, int(s), int(e)

    sum_comparisons: Dict[str, List[dict]] = defaultdict(list)

    with open(out_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([
            "reference1", "start1", "end1",
            "reference2", "start2", "end2",
            "jaccard", "containment_1_in_2", "containment_2_in_1"
        ])

        items = list(sigs.items())
        for rid1, sig1 in tqdm(items, desc="SBT search", miniters=1000):
            r1, s1, e1 = parse_region(rid1)
            mh1 = sig1.minhash

            for sr in sbt.search(sig1, threshold=min_jaccard, best_only=False):
                rid2 = sr.signature.name
                if rid2 == rid1:
                    continue
                r2, s2, e2 = parse_region(rid2)
                if r1 == r2:
                    continue

                j = float(sr.score)
                mh2 = sr.signature.minhash
                c12 = mh1.avg_containment(mh2)
                c21 = mh2.avg_containment(mh1)

                w.writerow([r1, s1, e1, r2, s2, e2, j, c12, c21])

                sum_comparisons[r1].append(dict(jaccard=j, to=r2, s1=s1, e1=e1, s2=s2, e2=e2))
                sum_comparisons[r2].append(dict(jaccard=j, to=r1, s1=s2, e1=e2, s2=s1, e2=e1))

    return sum_comparisons


# -----------------------------------------------------------------------------
# Removal planning (minimal)
# -----------------------------------------------------------------------------

def infer_gt_from_readname(readname: str) -> Optional[str]:
    """
    Your convention: accession_subseq_mutlen etc => drop last two underscore fields.
    """
    parts = readname.split("_")
    if len(parts) < 3:
        return None
    return "_".join(parts[:-2])


def build_removed_ids_gt_biased(
    bam_path: str,
    conflict_groups: List[set[Tuple[str, int, int]]],
    *,
    random_seed: Optional[int] = None,
) -> Dict[str, List[str]]:
    """
    Minimal “GT-biased” removal: remove alignments where ref != inferred GT
    within conflict groups. (You can add ANI/matrix gating later.)
    """
    if random_seed is not None:
        random.seed(random_seed)

    removed: Dict[str, set[str]] = defaultdict(set)

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for group in conflict_groups:
            by_ref: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
            for ref, s, e in group:
                by_ref[ref].append((s, e))

            for ref, segs in by_ref.items():
                for s, e in segs:
                    for aln in bam.fetch(ref, s, e):
                        if aln.is_unmapped:
                            continue
                        if aln.is_secondary or aln.is_supplementary:
                            continue
                        rid = aln.query_name
                        gt = infer_gt_from_readname(rid)
                        if gt is None:
                            continue
                        if gt != ref:
                            removed[rid].add(ref)

    return {rid: sorted(list(refs)) for rid, refs in removed.items()}


def build_conflict_groups(sum_comparisons: Dict[str, List[dict]], *, min_jaccard: float) -> List[set[Tuple[str, int, int]]]:
    """
    Connected components on regions using edges from sum_comparisons.
    """
    graph: Dict[Tuple[str, int, int], set[Tuple[str, int, int]]] = defaultdict(set)

    for ref, conflicts in sum_comparisons.items():
        for c in conflicts:
            if float(c.get("jaccard", 0.0)) < min_jaccard:
                continue
            a = (ref, int(c["s1"]), int(c["e1"]))
            b = (str(c["to"]), int(c["s2"]), int(c["e2"]))
            graph[a].add(b)
            graph[b].add(a)

    visited: set[Tuple[str, int, int]] = set()
    groups: List[set[Tuple[str, int, int]]] = []

    for node in graph:
        if node in visited:
            continue
        stack = [node]
        comp = set()
        visited.add(node)
        while stack:
            x = stack.pop()
            comp.add(x)
            for nb in graph[x]:
                if nb not in visited:
                    visited.add(nb)
                    stack.append(nb)
        groups.append(comp)

    return groups


def write_filtered_bam(
    bam_in_path: str,
    bam_out_path: str,
    removed_read_ids: Dict[str, List[str]],
) -> None:
    """
    Keep alignments unless (read_id in removed_read_ids AND ref in removed_read_ids[read_id]).
    """
    os.makedirs(os.path.dirname(bam_out_path), exist_ok=True)
    with pysam.AlignmentFile(bam_in_path, "rb") as bam_in, pysam.AlignmentFile(bam_out_path, "wb", template=bam_in) as bam_out:
        total = 0
        kept = 0
        for aln in bam_in:
            total += 1
            ref = aln.reference_name
            rid = aln.query_name
            if rid in removed_read_ids and ref in removed_read_ids[rid]:
                continue
            bam_out.write(aln)
            kept += 1

    pysam.index(bam_out_path)
    log.info("Filtered BAM: %s (kept %d/%d alignments)", bam_out_path, kept, total)


def build_support_weighted_matrix(
    region_comparisons_csv: str,
    *,
    hit_threshold: float = 0.10,
    min_hit_windows: int = 1,
    score_col: str = "jaccard",
    symmetrize: bool = True,
    fill_diagonal: float = 1.0,
    output_csv: str | None = None,
) -> pd.DataFrame:
    """
    Build an accession-level similarity matrix from region_comparisons.csv.

    region_comparisons.csv expected columns:
      reference1,start1,end1,reference2,start2,end2,jaccard,containment_1_in_2,containment_2_in_1

    Score definition (directed):
      score(A->B) = (# unique windows of A that have >=1 hit to B with score>=hit_threshold)
                   / (total # unique windows of A observed in the file)

    Then optionally:
      - filter out pairs with < min_hit_windows supporting windows
      - symmetrize by averaging A->B and B->A
      - fill diagonal with fill_diagonal

    Returns:
      pd.DataFrame with accessions as index/columns and values in [0..1].
    """
    df = pd.read_csv(region_comparisons_csv)

    required = {"reference1", "start1", "end1", "reference2", "start2", "end2", score_col}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in {region_comparisons_csv}: {sorted(missing)}")

    # Ensure numeric score
    df[score_col] = pd.to_numeric(df[score_col], errors="coerce")
    df = df.dropna(subset=[score_col])
    if df.empty:
        raise ValueError("No valid rows after parsing score column.")

    # Normalize types
    df["reference1"] = df["reference1"].astype(str)
    df["reference2"] = df["reference2"].astype(str)
    df["start1"] = pd.to_numeric(df["start1"], errors="coerce").astype("Int64")
    df["end1"] = pd.to_numeric(df["end1"], errors="coerce").astype("Int64")

    # Define the query node A and match node B
    df["A"] = df["reference1"]
    df["B"] = df["reference2"]

    # A window is uniquely identified by (A,start1,end1)
    df["A_win"] = (
        df["A"].astype(str) + ":" +
        df["start1"].astype(str) + "-" +
        df["end1"].astype(str)
    )

    # Total windows observed per A (denominator)
    nwin_A = df.groupby("A")["A_win"].nunique().to_dict()

    # Keep only sufficiently strong hits
    df = df[df[score_col] >= hit_threshold].copy()
    if df.empty:
        raise ValueError("No hits left after applying hit_threshold.")

    # Count unique query windows supporting each directed pair (A,B)
    win_hits = (
        df.groupby(["A", "B"])["A_win"]
          .nunique()
          .reset_index(name="hit_windows")
    )

    # Apply minimum support requirement
    if min_hit_windows > 1:
        win_hits = win_hits[win_hits["hit_windows"] >= min_hit_windows].copy()
        if win_hits.empty:
            raise ValueError("No pairs left after applying min_hit_windows.")

    # Directed fraction A->B
    win_hits["frac_A_to_B"] = win_hits.apply(
        lambda r: float(r["hit_windows"]) / float(max(1, nwin_A.get(r["A"], 1))),
        axis=1
    )

    # Nodes (include those that only appear as B as well)
    accs = sorted(set(df["A"]).union(set(df["B"])))
    mat = pd.DataFrame(0.0, index=accs, columns=accs, dtype=float)

    # Fill directed
    for r in win_hits.itertuples(index=False):
        A = r.A
        B = r.B
        v = float(r.frac_A_to_B)
        if v > mat.at[A, B]:
            mat.at[A, B] = v

    # Symmetrize if requested
    if symmetrize:
        mat_sym = mat.copy()
        for i in accs:
            for j in accs:
                if i == j:
                    mat_sym.at[i, j] = fill_diagonal
                else:
                    mat_sym.at[i, j] = 0.5 * (mat.at[i, j] + mat.at[j, i])
        mat = mat_sym
    else:
        # still fill diagonal if present
        for a in accs:
            mat.at[a, a] = fill_diagonal

    if output_csv:
        mat.to_csv(output_csv)

    return mat

# -----------------------------------------------------------------------------
# Orchestration
# -----------------------------------------------------------------------------
def determine_conflicts(
    *,
    output_dir: str,
    bam_path: str,
    bedgraph_path: str,
    fasta_paths: Sequence[str],
    compare_to_reference_windows: Optional[bool] = False,
    ksize: int,
    scaled: int,
    merge_method: str,
    merge_threshold: Optional[float],
    gap_allowance_frac: float,
    min_jaccard: float = 0.8,
    conflict_min_jaccard: float,
    use_bam_for_signatures: bool,
    sig_cache: Optional[str],
    write_filtered: bool,
    random_seed: Optional[int],
    compute_breadth: bool = True,

    window_size: int = 10_000,
    step_size: int = 10_000,
    matchfile: Optional[str] = None,

) -> Tuple[Dict[str, List[str]], pd.DataFrame]:
    os.makedirs(output_dir, exist_ok=True)
    sig_path = sig_cache or os.path.join(output_dir, "signatures", "merged_regions.sig")
    if matchfile:

        accession_to_taxid, taxid_to_desc, taxid_to_accessions = load_matchfile(
            matchfile,
            accession_col="Accession",
            taxid_col="TaxID",
            desc_col="Description",
        )
    acc_mat = build_support_weighted_matrix(
        os.path.join(output_dir, "region_comparisons.csv"),
        hit_threshold=0.8,
        min_hit_windows=5,
        output_csv=os.path.join(output_dir, "accession_support_matrix.csv")
    )
    # Load BAM reference lengths
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        reflengths = {ref: bam.get_reference_length(ref) for ref in bam.references}
    intervals = parse_bedgraph(bedgraph_path)
    if compare_to_reference_windows:
        merged_regions = make_window_regions_from_fastas(
            fasta_paths,
            window=window_size,
            step=step_size,
            min_len=ksize,      # optional: skip contigs shorter than ksize
        )
    else:
        merged_df = merge_bedgraph_regions(
            intervals,
            method=merge_method,
            stat_threshold=merge_threshold,
            gap_allowance_frac=gap_allowance_frac,
            reflengths=reflengths,
        )
        merged_regions = [
            Region(str(r.Reference), int(r.start), int(r.end), float(r.mean_depth))
            for r in merged_df.itertuples(index=False)
        ]
    print("Total list of all Merged Regions", len(merged_regions))
    # Signatures (cache-aware)
    if sig_cache and os.path.exists(sig_path):
        sigs = load_signatures_file(sig_path)
        log.info("Loaded signatures from cache: %s (%d)", sig_path, len(sigs))
    else:
        if fasta_paths and not use_bam_for_signatures:
            sigs = build_region_signatures_from_fasta(
                merged_regions, fasta_paths, ksize=ksize, scaled=scaled
            )
        else:
            sigs = build_region_signatures_from_bam(
                merged_regions, bam_path, ksize=ksize, scaled=scaled
            )
        save_signatures_file(sigs, sig_path)
    # Compare signatures
    comparisons_csv = os.path.join(output_dir, "region_comparisons.csv")
    sum_comparisons = compare_signatures_sbt(
        sigs,
        ksize=ksize,
        min_jaccard=min_jaccard,
        out_csv=comparisons_csv,
    )

    # Conflict groups
    groups = build_conflict_groups(sum_comparisons, min_jaccard=conflict_min_jaccard)
    log.info("Conflict groups: %d", len(groups))

    # Removal plan (GT-biased minimal)
    removed_read_ids = build_removed_ids_gt_biased(
        bam_path,
        groups,
        random_seed=random_seed,
    )
    log.info("Removal plan affects %d reads", len(removed_read_ids))

    # Optionally write filtered BAM
    filtered_bam_path = None
    if write_filtered:
        filtered_bam_path = os.path.join(output_dir, "filtered.bam")
        write_filtered_bam(bam_path, filtered_bam_path, removed_read_ids)

    # Restore your “removal stats” dataframe
    comparison_df = compute_removal_stats_df(
        bam_path=bam_path,
        removed_read_ids=removed_read_ids,
        reflengths=reflengths,
        compute_breadth=compute_breadth,
    )
    comparison_df.to_excel(os.path.join(output_dir, "removal_stats.xlsx"), index=False)

    # Also keep a simple removed-reads CSV (useful for debugging)
    pd.DataFrame(
        [{"read_id": rid, "remove_from_refs": ",".join(refs)} for rid, refs in removed_read_ids.items()]
    ).to_csv(os.path.join(output_dir, "removed_reads.csv"), index=False)
    print_metrics_comparison(comparison_df)
    return removed_read_ids, comparison_df
def print_metrics_comparison(comparison_df):
    """
    Print the same console summary you had before:
      - only rows where Δ All != 0
      - selected columns
    """
    try:
        print("\n=== Metrics Comparison ===\n")

        if comparison_df is None or comparison_df.empty:
            print("(comparison_df is empty)")
            return

        # Ensure Δ All% is numeric if present (optional; harmless if already numeric)
        if "Δ All%" in comparison_df.columns:
            comparison_df["Δ All%"] = pd.to_numeric(comparison_df["Δ All%"], errors="coerce")

        cols = [
            "Reference",
            "Δ All%",
            "Δ^-1 Breadth",
            "Breadth New",
            "Breadth Original",
            "TP New",
            "TP Original",
            "FP Original",
            "FP New",
        ]
        # Keep only columns that exist (prevents KeyError if you tweak schema)
        cols = [c for c in cols if c in comparison_df.columns]

        if "Δ All" in comparison_df.columns:
            view = comparison_df[comparison_df["Δ All"] != 0]
        else:
            # fallback: just print everything if Δ All is missing
            view = comparison_df

        if view.empty:
            print("(No rows where Δ All != 0)")
            return

        print(view[cols].to_string(index=False))

    except Exception as ex:
        print(f"Error while printing metrics comparison to console: {ex}")

def calculate_breadth_coverage_from_bam(
    *,
    bam_path: str,
    reflengths: Dict[str, int],
    removed_read_ids: Dict[str, List[str]],
) -> Dict[str, float]:
    """
    Breadth% per reference, skipping alignments removed for that reference.

    This is alignment-interval breadth (union of aligned intervals).
    """
    breadth = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for ref, L in reflengths.items():
            total_covered = 0
            current_start = None
            current_end = None

            try:
                it = bam.fetch(ref, 0, L)
            except ValueError:
                breadth[ref] = 0.0
                continue

            for aln in it:
                if aln.is_unmapped or aln.reference_start is None or aln.reference_end is None:
                    continue
                rid = aln.query_name
                if rid in removed_read_ids and ref in removed_read_ids[rid]:
                    continue

                s = max(int(aln.reference_start), 0)
                e = min(int(aln.reference_end), L)
                if s >= e:
                    continue

                if current_start is None:
                    current_start, current_end = s, e
                else:
                    if s <= current_end + 1:
                        current_end = max(current_end, e)
                    else:
                        total_covered += (current_end - current_start)
                        current_start, current_end = s, e

            if current_start is not None:
                total_covered += (current_end - current_start)

            breadth[ref] = (100.0 * total_covered / L) if L > 0 else 0.0

    return breadth

def compute_removal_stats_df(
    *,
    bam_path: str,
    removed_read_ids: Dict[str, List[str]],
    reflengths: Dict[str, int],
    compute_breadth: bool = False,
) -> pd.DataFrame:
    """
    Returns a per-reference dataframe similar to your original comparison_df.

    - "Original" counts: all alignments in BAM.
    - "New" counts: alignments that survive the removal plan (read_id not removed for that ref).
    - TP/FP based on GT inferred from read name (infer_gt_from_readname()).
    - FN New: a TP alignment that got removed for that ref.

    Breadth is optional (expensive). If compute_breadth=True, we compute breadth_old and breadth_new.
    """
    stats = defaultdict(lambda: defaultdict(int))

    # Optional breadth
    breadth_old = {}
    breadth_new = {}

    # Iterate alignments once for counts
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for aln in bam:
            if aln.is_unmapped:
                continue
            if aln.reference_name is None:
                continue

            ref = aln.reference_name
            rid = aln.query_name

            gt = infer_gt_from_readname(rid)
            if gt is None:
                # If you prefer, count these separately; for now skip
                continue

            stats[ref]["total_reads"] += 1

            # Original TP/FP for this alignment
            if gt == ref:
                stats[ref]["TP Original"] += 1
            else:
                stats[ref]["FP Original"] += 1

            # "Pass filtered" = not removed for this reference
            removed_for_ref = (rid in removed_read_ids and ref in removed_read_ids[rid])
            if not removed_for_ref:
                stats[ref]["pass_filtered_reads"] += 1
                if gt == ref:
                    stats[ref]["TP New"] += 1
                else:
                    stats[ref]["FP New"] += 1
            else:
                # removed alignment
                if gt == ref:
                    stats[ref]["FN New"] += 1
                else:
                    stats[ref]["TN New"] += 1

    # Optional breadth computation (second pass)
    if compute_breadth:
        # breadth_old: all alignments
        breadth_old = calculate_breadth_coverage_from_bam(
            bam_path=bam_path,
            reflengths=reflengths,
            removed_read_ids={},  # remove nothing
        )
        # breadth_new: skip removed alignments
        breadth_new = calculate_breadth_coverage_from_bam(
            bam_path=bam_path,
            reflengths=reflengths,
            removed_read_ids=removed_read_ids,
        )

    # Build DF
    rows = []
    for ref in sorted(reflengths.keys()):
        s = stats.get(ref, {})
        tp_o = int(s.get("TP Original", 0))
        fp_o = int(s.get("FP Original", 0))
        tp_n = int(s.get("TP New", 0))
        fp_n = int(s.get("FP New", 0))
        fn_n = int(s.get("FN New", 0))

        total = int(s.get("total_reads", 0))
        passed = int(s.get("pass_filtered_reads", 0))

        precision = tp_o / (tp_o + fp_o) if (tp_o + fp_o) > 0 else 0.0
        recall = tp_o / total if total > 0 else 0.0
        f1 = (2 * precision * recall) / (precision + recall) if (precision + recall) > 0 else 0.0

        breadth_old_val = float(breadth_old.get(ref, 0.0)) if compute_breadth else 0.0
        breadth_new_val = float(breadth_new.get(ref, 0.0)) if compute_breadth else 0.0

        delta_reads = passed - total
        delta_all_pct = 100.0 * (delta_reads / total) if total > 0 else 0.0
        delta_breadth = breadth_new_val - breadth_old_val
        delta_breadth_ratio = (breadth_new_val / breadth_old_val) if breadth_old_val > 0 else 0.0

        rows.append({
            "Reference": ref,
            "TP Original": tp_o,
            "FP Original": fp_o,
            "FN Original": 0,  # not tracked in this alignment-wise formulation
            "TP New": tp_n,
            "FP New": fp_n,
            "FN New": fn_n,
            "Total Reads": total,
            "Pass Filtered Reads": passed,
            "Proportion Aligned": (tp_o / total) if total > 0 else 0.0,
            "Precision": precision,
            "Recall": recall,
            "F1": f1,
            "Δ All": delta_reads,
            "Δ All%": delta_all_pct,
            "Breadth Original": breadth_old_val,
            "Breadth New": breadth_new_val,
            "Δ Breadth": delta_breadth,
            "Δ^-1 Breadth": delta_breadth_ratio,
        })

    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Detect ambiguous regions and plan read removals.")
    p.add_argument("--output-dir", required=True)
    p.add_argument("--bam", required=True)
    p.add_argument("--bedgraph", required=True, help="BEDGRAPH-like file: Reference start end depth")

    p.add_argument("--fasta", nargs="*", default=[], help="Optional FASTA(s) for FASTA-based sketches")
    p.add_argument("--use-bam-signatures", action="store_true", help="Sketch from BAM reads instead of FASTA")
    p.add_argument("--sig-cache", default=None, help="Path to .sig cache; if exists, load it")

    p.add_argument("--ksize", type=int, default=31)
    p.add_argument("--scaled", type=int, default=200)

    p.add_argument("--merge-method", choices=["jump", "variance", "gini"], default="jump")
    p.add_argument("--merge-threshold", type=float, default=None, help="Jump/variance/gini threshold (method-dependent)")
    p.add_argument("--gap-allowance-frac", type=float, default=0.1, help="Allowed gap = ref_len * frac")

    p.add_argument("--min-jaccard", type=float, default=0.2, help="SBT search threshold for reporting comparisons")
    p.add_argument("--conflict-min-jaccard", type=float, default=0.0, help="Edge threshold when building conflict groups")

    p.add_argument("--write-filtered-bam", action="store_true")
    p.add_argument("--seed", type=int, default=None)

    p.add_argument("-v", "--verbose", action="count", default=0)
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_argparser().parse_args(argv)
    setup_logging(args.verbose)

    determine_conflicts(
        output_dir=args.output_dir,
        bam_path=args.bam,
        bedgraph_path=args.bedgraph,
        fasta_paths=args.fasta,
        ksize=args.ksize,
        scaled=args.scaled,
        merge_method=args.merge_method,
        merge_threshold=args.merge_threshold,
        gap_allowance_frac=args.gap_allowance_frac,
        min_jaccard=args.min_jaccard,
        use_bam_for_signatures=args.use_bam_signatures,
        sig_cache=args.sig_cache,
        write_filtered=args.write_filtered_bam,
        conflict_min_jaccard=args.conflict_min_jaccard,
        random_seed=args.seed,
    )

    log.info("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
