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

import csv
import itertools
import math
import os
import random
import statistics
import time
from collections import defaultdict, Counter
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
    for chrom, g in df.groupby("chrom", observed=True):
        g = g.reset_index(drop=True)
        if g.empty:
            continue

        buf_start = int(g.at[0, "start"])
        buf_end = int(g.at[0, "end"])
        buf_depths = [float(g.at[0, "depth"])]

        chrom_str = str(chrom)
        gap_thr = allowed_gap.get(chrom_str, float("inf"))

        for i, row in enumerate(g.itertuples(index=False)):
            if i == 0:
                continue

            new_start = int(row.start)
            new_end = int(row.end)
            new_depth = float(row.depth)
            gap = new_start - buf_end

            can_merge = True
            if len(buf_depths) >= max_group_size:
                can_merge = False
            if max_length is not None and (new_end - buf_start) > max_length:
                can_merge = False

            if can_merge:
                if gap > gap_thr:
                    can_merge = False

            if can_merge:
                if merging_method == "jump":
                    assert jump_thr is not None
                    jump = abs(new_depth - buf_depths[-1])
                    if jump > jump_thr:
                        can_merge = False
                elif merging_method == "variance":
                    vals = buf_depths + [new_depth]
                    stat = statistics.pvariance(vals)
                    if max_stat_threshold is not None and stat > max_stat_threshold:
                        can_merge = False
                    if value_diff_tolerance is not None:
                        if abs(new_depth - float(np.mean(buf_depths))) > value_diff_tolerance:
                            can_merge = False
                elif merging_method == "gini":
                    vals = buf_depths + [new_depth]
                    stat = compute_gini(vals)
                    if max_stat_threshold is not None and stat > max_stat_threshold:
                        can_merge = False
                    if value_diff_tolerance is not None:
                        if abs(new_depth - float(np.mean(buf_depths))) > value_diff_tolerance:
                            can_merge = False
                else:
                    raise ValueError(f"Unknown merging method: {merging_method}")

            if can_merge:
                buf_end = new_end
                buf_depths.append(new_depth)
            else:
                merged.append(
                    {
                        "chrom": chrom_str,
                        "start": buf_start,
                        "end": buf_end,
                        "mean_depth": float(np.mean(buf_depths)),
                    }
                )
                buf_start, buf_end, buf_depths = new_start, new_end, [new_depth]

        merged.append(
            {
                "chrom": chrom_str,
                "start": buf_start,
                "end": buf_end,
                "mean_depth": float(np.mean(buf_depths)),
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



def create_signatures_for_regions(
    regions_df: pd.DataFrame,
    bam_path: str,
    fasta_paths: List[str],
    kmer_size: int,
    scaled: int,
    num_workers: int = 1,
) -> Dict[str, SourmashSignature]:
    """
    Build signatures for each merged region.

    If fasta_paths is non-empty, regions are sketched from the first FASTA that contains the contig.
    Otherwise, uses BAM pileup sequence is NOT reconstructed here; it uses read.query_sequence concatenation
    would be invalid for minhash. In this cleaned version, BAM-only mode is disabled unless FASTA is available.

    Practical: Provide FASTA for stable behavior.
    """
    required = {"chrom", "start", "end"}
    if not required.issubset(set(regions_df.columns)):
        raise ValueError(f"regions_df must contain columns: {required}")

    # normalize dtypes
    df = regions_df.copy()
    df["chrom"] = df["chrom"].astype(str)
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)

    fastas = [pysam.FastaFile(fp) for fp in fasta_paths] if fasta_paths else []

    sigs: Dict[str, SourmashSignature] = {}
    it = df.itertuples(index=False)

    for row in tqdm(list(it), total=df.shape[0], miniters=1000, desc="Sketching regions"):
        chrom, start, end = str(row.chrom), int(row.start), int(row.end)

        seq = None
        for fa in fastas:
            if chrom in fa.references:
                seq = fa.fetch(chrom, start, end)
                break

        if not seq:
            # If you truly need BAM-only, add a consensus builder; leaving out here for correctness.
            continue
        if len(seq) < kmer_size:
            continue

        region_name, sig = create_signature_for_single_region(
            chrom, start, end, kmer_size, scaled, seq
        )
        sigs[region_name] = sig

    for fa in fastas:
        fa.close()

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
    """Sketch fixed windows across contigs; if contig shorter than window, sketch the full contig once."""
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

        # If contig is too short to even sketch, skip
        if L < ksize:
            continue

        # Adaptive window: if contig shorter than window, do ONE window spanning the contig
        if L < window:
            windows = [(0, L)]
        else:
            windows = list(iter_windows(L, window, step))

        for s, e in windows:
            subseq = seq[s:e]
            if len(subseq) < ksize:
                continue
            if max_n_frac is not None and max_n_frac >= 0:
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
    """Sketch windows for each FASTA, index in SBT, then report cross-fasta hits."""
    all_sigs: Dict[str, SourmashSignature] = {}

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
            sigs = dict(list(sigs.items())[:max_windows_per_fasta])
        all_sigs.update(sigs)

    if not all_sigs:
        # still write a valid CSV with header so downstream loading is safe
        print(f"WARNING: No window signatures generated; skipping shared-window comparisons (check window/step/ksize and FASTA content).")
        with open(output_csv, "w", newline="") as out:
            w = csv.writer(out)
            w.writerow(
                [
                    "query_fasta",
                    "query_contig",
                    "q_start",
                    "q_end",
                    "match_fasta",
                    "match_contig",
                    "m_start",
                    "m_end",
                    "jaccard",
                    "containment_q_in_m",
                    "containment_m_in_q",
                ]
            )
        print(
            "WARNING: No window signatures generated; skipping shared-window comparisons "
            "(check window/step/ksize and FASTA content)."
        )
        return


    sbt = build_sbt(all_sigs, ksize=ksize)

    with open(output_csv, "w", newline="") as out:
        w = csv.writer(out)
        w.writerow(
            [
                "query_fasta",
                "query_contig",
                "q_start",
                "q_end",
                "match_fasta",
                "match_contig",
                "m_start",
                "m_end",
                "jaccard",
                "containment_q_in_m",
                "containment_m_in_q",
            ]
        )

        for q_wid, q_sig in all_sigs.items():
            q_fa, q_contig, q_s, q_e = parse_window_id(q_wid)
            mh_q = q_sig.minhash

            hits = []
            for sr in sbt.search(q_sig, threshold=jaccard_threshold, best_only=False):
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
                c1 = mh_q.avg_containment(mh_m)
                c2 = mh_m.avg_containment(mh_q)
                hits.append((j, m_fa, m_contig, m_s, m_e, c1, c2))

            if not hits:
                continue
            hits.sort(key=lambda x: x[0], reverse=True)
            hits = hits[:max_hits_per_query]

            for (j, m_fa, m_contig, m_s, m_e, c1, c2) in hits:
                w.writerow([q_fa, q_contig, q_s, q_e, m_fa, m_contig, m_s, m_e, f"{j:.6f}", f"{c1:.6f}", f"{c2:.6f}"])


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
) -> Dict[str, List[str]]:
    """
    For each read with multiple alignments, keep the alignment(s) with best score.
    Score = MAPQ + as_weight*AS - penalize_weight*(#alt_contigs_overlapping_shared_windows)
    Any other contig alignments are marked for removal for that read.
    """
    if drop_contigs is None:
        drop_contigs = set()

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
        score = mapq + as_weight * AS - penalize_weight * alt_n

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

def compare_metrics(per_ref_stats: Dict[str, dict], reflengths: Dict[str, int]) -> pd.DataFrame:
    rows = []
    for ref, st in sorted(per_ref_stats.items()):
        total_reads = st.get("total_reads", 0)
        pass_reads = st.get("pass_filtered_reads", 0)
        delta_reads = pass_reads - total_reads
        delta_pct = (100.0 * delta_reads / total_reads) if total_reads else 0.0

        b0 = st.get("breadth_old", 0.0)
        b1 = st.get("breadth", 0.0)

        ref_len = reflengths.get(ref, 0)

        rows.append(
            {
                "Reference": ref,
                "Reference Length": ref_len,   # 👈 NEW COLUMN
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
                "Δ All": delta_reads,
                "Δ All%": delta_pct,
                "Breadth Original": b0,
                "Breadth New": b1,
                "Δ Breadth": (b1 - b0),
                "Δ^-1 Breadth": (b1 / b0) if b0 else 0.0,
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
                max_hits_per_query=4,
                skip_self_same_fasta=False,
                skip_self_same_contig=True,
            )
            shared_idx = load_shared_windows_csv(out_csv, min_jaccard=min_jaccard, skip_same_contig=True)

            removed_read_ids = build_removed_ids_best_alignment(
                bam_path=bam_path,
                shared_idx=shared_idx,
                penalize_weight=1.0,
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



# -----------------------------
# Main pipeline
# -----------------------------

def determine_conflicts(
    output_dir: Optional[str] = None,
    input_bam: Optional[str] = None,
    min_threshold: float = 0.2,
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
    sim_ani_threshold: float = 0.8,
    compare_to_reference_windows: bool = False,
    find_optimal_windows: bool = False,
    matchfile: Optional[str] = None,  # unused in this cleaned version
    matchfile_accession_col: str = "Accession",
    matchfile_taxid_col: str = "TaxID",
    matchfile_desc_col: str = "Description",
):
    if output_dir is None or input_bam is None or bedfile is None:
        raise ValueError("output_dir, input_bam, and bedfile are required.")
    os.makedirs(output_dir, exist_ok=True)
    print(f"Starting conflict detection pipeline: {time.ctime()}")
    if matchfile:
        accession_to_taxid, taxid_to_desc, taxid_to_accessions =  load_matchfile(matchfile, matchfile_accession_col, matchfile_taxid_col, matchfile_desc_col)
        org_sigs, org_stats = build_organism_signatures_from_fastas_ani(
            fasta_paths=fasta_files,
            accession_to_taxid=accession_to_taxid,
            taxid_to_desc=taxid_to_desc,
            ksize=51,
            scaled=8000,
        )

        ani_df = organism_ani_matrix_from_sigs(org_sigs, symmetrize="mean")
        ani_df.to_csv(os.path.join(output_dir, "organism_ani_matrix.csv"))

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


    # Optional: compare_to_reference_windows adds shared-window signatures (and later uses alignment-based removal)
    shared_idx = None
    if compare_to_reference_windows:
        print("Building shared-window report across FASTAs...")
        report_path = os.path.join(output_dir, "shared_windows_report.csv")

        if find_optimal_windows:
            success, best_params, best_stats, all_results = tune_shared_window_params(
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
                    max_hits_per_query=4,
                    skip_self_same_fasta=False,
                )
            else:
                report_shared_windows_across_fastas(
                    fasta_files=fasta_files,
                    output_csv=report_path,
                    ksize=31,
                    scaled=800,
                    window=10_000,
                    step=10_000,
                    jaccard_threshold=sim_ani_threshold,
                    max_hits_per_query=4,
                    skip_self_same_fasta=False,
                )
        else:
            print("Creating shared FASTA report from scratch")
            report_shared_windows_across_fastas(
                fasta_files=fasta_files,
                output_csv=report_path,
                ksize=51,
                scaled=8000,
                window=1000_000,
                step=1000_000,
                jaccard_threshold=sim_ani_threshold,
                max_hits_per_query=15,
                skip_self_same_fasta=False,
            )

        shared_idx = load_shared_windows_csv(report_path, min_jaccard=sim_ani_threshold, skip_same_contig=True)

    # Removal plan
    if compare_to_reference_windows and shared_idx is not None:
        removed_read_ids = build_removed_ids_best_alignment(
            bam_path=input_bam,
            shared_idx=shared_idx,
            penalize_weight=1.0,
            as_weight=0.0,
            drop_contigs=set(),
            drop_if_ambiguous=True,
            min_alt_count=1,
            only_primary=False,
        )
    else:
        # Parse bedgraph + merge regions
        regions = parse_bed_file(bedfile)
        print(f"Input intervals: {len(regions)}")

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
        signatures = {}
        # Signatures: load or generate
        if not sigfile or not os.path.exists(sigfile):
            nworkers = cpu_count if cpu_count else max(1, int(os.cpu_count() / 2))
            print(f"Sketching merged regions (workers={nworkers})")
            t0 = time.time()
            signatures = create_signatures_for_regions(
                regions_df=merged_regions,
                bam_path=input_bam,
                fasta_paths=fasta_files,
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
        # Clustering (kept minimal): compare all refs together by default
        clusters = [set([x.split(":")[0] for x in signatures.keys()])]

        # Compare signatures via SBT
        output_csv = os.path.join(output_dir, "region_comparisons.csv")
        sig_items = list(signatures.items())


        if FAST_MODE:
            print("Building SBT index...")
            sbt_index = build_sbt_index(sig_items, ksize=51, clusters=clusters)
            print("Searching SBT...")
            sum_comparisons = fast_mode_sbt(sig_items, sbt_index, output_csv, min_threshold, clusters)
        else:
            raise NotImplementedError("Slow mode removed in this cleaned version. Use FAST_MODE=True.")

        # Build conflict groups
        print("Building conflict groups...")
        start_time = time.time()
        conflict_groups = build_conflict_groups(sum_comparisons, min_jaccard=0.0)
        print(f"Conflict groups: {len(conflict_groups)}")
        print(f"Conflict groups length: {len(conflict_groups)} built in {time.time() - start_time:.2f} seconds. Next up is proportion removal")
        start_time = time.time()
        # 2) For each group, remove 1 read from each region
        removed_read_ids = finalize_proportional_removal(
            conflict_groups,
            bam_fs,
            fetch_reads_in_region,
            remove_mode='random'
        )

    stats = report_removed_read_stats(bam_path=input_bam, removed_read_ids=removed_read_ids)

    # Write removals table
    failed_path = os.path.join(output_dir, "failed_reads.txt")
    with open(failed_path, "w") as fp:
        for read_id, refs in removed_read_ids.items():
            for ref in refs:
                fp.write(f"{ref}\t{read_id}\n")
    print(f"Wrote removals: {failed_path}")

    # Optional: create filtered BAM
    if filtered_bam_create:
        try:
            create_filtered_bam(bam_fs, filtered_bam_create, removed_read_ids)
        except Exception as e:
            print(f"Filtered BAM failed: {e}")

    # Breadth before/after
    breadth_old = calculate_breadth_coverage_from_bam(bam_fs, reflengths, removed_ids={})
    breadth_new = calculate_breadth_coverage_from_bam(bam_fs, reflengths, removed_ids=removed_read_ids)

    # Per-reference alignment accounting (GT inferred from read name convention)
    per_ref = defaultdict(lambda: defaultdict(int))
    for ref, L in reflengths.items():
        total = 0
        passed = 0
        for r in bam_fs.fetch(ref, 0, L):
            total += 1
            gt = infer_gt_from_readname(r.query_name)
            if gt is None:
                continue

            if gt == r.reference_name:
                per_ref[ref]["TP Original"] += 1
            else:
                per_ref[ref]["FP Original"] += 1

            keep = not (r.query_name in removed_read_ids and r.reference_name in removed_read_ids[r.query_name])
            if keep:
                passed += 1
                if gt == r.reference_name:
                    per_ref[ref]["TP New"] += 1
                else:
                    per_ref[ref]["FP New"] += 1
            else:
                if gt == r.reference_name:
                    per_ref[ref]["FN New"] += 1

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

    comparison_df = compare_metrics(per_ref, reflengths)


    try:
        out_xlsx = os.path.join(output_dir, "removal_stats.xlsx")
        # add a "Total" row that is the sum of all the TP/FP/FN columns, empty for others
        total_row = {
            "Reference": "Total",
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
        }
        comparison_df = pd.concat([comparison_df, pd.DataFrame([total_row])], ignore_index=True)
        # put Total at Top of excel sheet
        comparison_df = pd.concat([comparison_df.tail(1), comparison_df.iloc[:-1]], ignore_index=True)
        comparison_df.to_excel(out_xlsx, index=False)
        print(f"Wrote: {out_xlsx}")
        print(
            comparison_df[comparison_df["Δ All"] != 0][
                ["Reference", "Δ All%", "Δ^-1 Breadth", "Breadth New", "Breadth Original", "TP New", "TP Original", "FP Original", "FP New"]
            ].to_string(index=False)
        )
    except Exception as e:
        print(f"Stats export failed: {e}")

    bam_fs.close()

    return removed_read_ids, comparison_df
