#!/usr/bin/env python3
"""
ground_truth.py

Compute per-reference and global confusion metrics (OLD vs NEW alignments),
and produce an XLSX report of confusion metrics and remaining false-positive reads.

Usage example:
    python ground_truth.py \
        --input aligned_sorted.bam \
        --match orgs.txt \
        --filtered_bam search_results/filtered.hash.bam \
        --output_xlsx alignment_confusion_report.xlsx \
        --report_confusion_xlsx

The script exposes a `main()` function so it can be imported and called programmatically.
"""
from __future__ import annotations

import argparse
import os
import re
from collections import defaultdict
from typing import Dict, Iterable, Optional, Set, Tuple
from optimize_weights import calculate_scores
import pandas as pd
import pysam
import math


# regex to extract truth accession from read_id like "NC_003310.1_1133_6"
_GT_RE = re.compile(r"^(.*)_\d+_\d+$")


def ground_truth_from_read_id(read_id: str) -> str:
    """Extract ground-truth accession from read_id like NZ_CP076401.1_27_5 -> NZ_CP076401.1"""
    m = _GT_RE.match(read_id or "")
    return m.group(1).strip() if m else (read_id or "").strip()


def _strip_version(ref: str) -> str:
    ref = (ref or "").strip()
    return ref.split(".", 1)[0] if "." in ref else ref


def meta_for_ref(ref: str, ref_meta: Optional[Dict[str, Dict[str, str]]] = None) -> Tuple[Optional[str], Optional[str]]:
    """
    Resolve (taxid, organism) for a reference using:
      1) exact key
      2) stripped-version key
    """
    ref_meta = ref_meta or {}
    if not ref:
        return None, None
    r0 = ref.strip()
    md = ref_meta.get(r0) or ref_meta.get(_strip_version(r0)) or {}
    taxid = md.get("taxid")
    org = md.get("organism")
    taxid = None if taxid in [None, "", "None"] else str(taxid)
    org = None if org in [None, "", "None"] else str(org)
    return taxid, org


def join_taxids(refs: Iterable[str], ref_meta: Optional[Dict[str, Dict[str, str]]] = None) -> str:
    ref_meta = ref_meta or {}
    out = []
    for r in sorted(set(refs)):
        md = ref_meta.get(r) or ref_meta.get(_strip_version(r)) or {}
        t = md.get("taxid")
        if t not in [None, "", "None"]:
            out.append(str(t))
    return ",".join(out)


def join_orgs(refs: Iterable[str], ref_meta: Optional[Dict[str, Dict[str, str]]] = None) -> str:
    ref_meta = ref_meta or {}
    out = []
    for r in sorted(set(refs)):
        md = ref_meta.get(r) or ref_meta.get(_strip_version(r)) or {}
        o = md.get("organism")
        if o not in [None, "", "None"]:
            out.append(str(o))
    return ",".join(out)


def collect_read_alignments(bam_path: str, alignments_to_remove=None):
    """
    Collect per-read aligned references and sequence.

    Returns:
      read_to_refs: dict(read_id -> set(refs aligned))
      read_to_seq:  dict(read_id -> seq)
      all_reads:    set(read_id) includes unmapped too

    alignments_to_remove is optional. If provided, it is assumed to be:
      alignments_to_remove[read_id] = set(refs_to_remove_for_that_read)
    """
    read_to_refs = defaultdict(set)
    read_to_seq = {}
    all_reads = set()

    if not bam_path or not os.path.exists(bam_path):
        raise FileNotFoundError(f"BAM not found: {bam_path}")

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            rid = read.query_name
            all_reads.add(rid)

            # store sequence once per read
            if rid not in read_to_seq and getattr(read, "query_sequence", None):
                read_to_seq[rid] = read.query_sequence

            if read.is_unmapped:
                continue

            ref = read.reference_name

            if alignments_to_remove:
                rm = alignments_to_remove.get(rid)
                if isinstance(rm, set) and ref in rm:
                    continue

            read_to_refs[rid].add(ref)

    return read_to_refs, read_to_seq, all_reads


def build_confusion_dataframe(
    read_to_refs: Dict[str, Set[str]],
    all_reads: Set[str],
    ref_meta: Optional[Dict[str, Dict[str, str]]] = None,
    prefix: str = "NEW_",
) -> pd.DataFrame:
    """
    Per-reference confusion metrics:
      TP/FP/FN/TN, Precision/Recall/F1, TruthSupport, PredSupport
    Adds Reference_TaxID and Reference_Organism.

    prefix: "OLD_" or "NEW_" for side-by-side merge later.
    """
    ref_meta = ref_meta or {}
    read_to_truth = {rid: ground_truth_from_read_id(rid) for rid in all_reads}

    refs = set(read_to_truth.values())
    for rset in read_to_refs.values():
        refs.update(rset)

    total = len(all_reads)
    rows = []

    for ref in sorted(refs):
        tp = fp = fn = 0

        for rid in all_reads:
            truth = read_to_truth[rid]
            aligned = ref in read_to_refs.get(rid, set())

            if truth == ref and aligned:
                tp += 1
            elif truth != ref and aligned:
                fp += 1
            elif truth == ref and not aligned:
                fn += 1

        tn = total - tp - fp - fn
        precision = tp / (tp + fp) if (tp + fp) else 0.0
        recall = tp / (tp + fn) if (tp + fn) else 0.0
        f1 = (2 * precision * recall / (precision + recall)) if (precision + recall) else 0.0

        taxid, org = meta_for_ref(ref, ref_meta)

        rows.append(
            {
                "Reference": ref,
                "Reference_TaxID": taxid,
                "Reference_Organism": org,
                f"{prefix}TP": tp,
                f"{prefix}FP": fp,
                f"{prefix}FN": fn,
                f"{prefix}TN": tn,
                f"{prefix}Precision": precision,
                f"{prefix}Recall": recall,
                f"{prefix}F1": f1,
                f"{prefix}TruthSupport": sum(1 for rid in all_reads if read_to_truth[rid] == ref),
                f"{prefix}PredSupport": sum(1 for rid in all_reads if ref in read_to_refs.get(rid, set())),
            }
        )

    return pd.DataFrame(rows)


def build_false_positive_reads_df(
    read_to_refs_original: Dict[str, Set[str]],
    read_to_refs_filtered: Dict[str, Set[str]],
    read_to_seq: Dict[str, str],
    all_reads: Set[str],
    ref_meta: Optional[Dict[str, Dict[str, str]]] = None,
) -> pd.DataFrame:
    """
    Sheet: remaining false positives after filtering.

    Keeps:
      - truth_ref + truth_taxid + truth_organism
      - OLD aligned refs + OLD aligned taxids + OLD aligned orgs + n_refs_original
      - NEW aligned refs + NEW aligned taxids + n_refs_filtered (NO new org columns)
      - remaining_fp_refs + remaining_fp_taxids + n_remaining_fp_refs
      - sequence
    """
    ref_meta = ref_meta or {}
    rows = []

    for rid in sorted(all_reads):
        truth = ground_truth_from_read_id(rid)
        truth_taxid, truth_org = meta_for_ref(truth, ref_meta)

        refs_f = read_to_refs_filtered.get(rid, set())
        refs_o = read_to_refs_original.get(rid, set())

        # still FP after filtering if any aligned ref != truth remains
        remaining_fp = sorted([r for r in refs_f if r != truth])
        if not remaining_fp:
            continue

        rows.append(
            {
                "read_id": rid,
                "truth_ref": truth,
                "truth_taxid": truth_taxid,
                "truth_organism": truth_org,
                # OLD
                "aligned_refs_original": ",".join(sorted(refs_o)) if refs_o else "",
                "aligned_taxids_original": join_taxids(refs_o, ref_meta) if refs_o else "",
                "aligned_orgs_original": join_orgs(refs_o, ref_meta) if refs_o else "",
                "n_refs_original": len(refs_o),
                # NEW (no org cols)
                "aligned_refs_filtered": ",".join(sorted(refs_f)) if refs_f else "",
                "aligned_taxids_filtered": join_taxids(refs_f, ref_meta) if refs_f else "",
                "n_refs_filtered": len(refs_f),
                # remaining FP (post)
                "remaining_fp_refs": ",".join(remaining_fp),
                "remaining_fp_taxids": join_taxids(remaining_fp, ref_meta) if remaining_fp else "",
                "n_remaining_fp_refs": len(remaining_fp),
                "sequence": read_to_seq.get(rid, ""),
            }
        )

    return pd.DataFrame(rows)


def write_confusion_xlsx(out_xlsx_path: str, confusion_df: pd.DataFrame, fp_reads_df: pd.DataFrame) -> None:
    """
    Writes:
      Sheet 1: confusion_matrix
      Sheet 2: confusion_summary (compact)
      Sheet 3: remaining_false_positives
    """
    # Build compact summary (safe columns)
    summary_columns = [
        "Reference",
        "OLD_TP", "OLD_FP", "OLD_FN", "OLD_TN", "OLD_F1",
        "NEW_TP", "NEW_FP", "NEW_FN", "NEW_TN", "NEW_F1",
        "Reference_TaxID", "Reference_Organism",
    ]
    # Ensure columns exist
    df_copy = confusion_df.copy()
    for col in summary_columns:
        if col not in df_copy.columns:
            df_copy[col] = 0

    numeric_cols = ["OLD_TP", "OLD_FP", "OLD_FN", "OLD_TN", "OLD_F1",
                    "NEW_TP", "NEW_FP", "NEW_FN", "NEW_TN", "NEW_F1"]

    summary_df = df_copy.loc[:, summary_columns].copy()
    for c in numeric_cols:
        summary_df[c] = pd.to_numeric(summary_df[c], errors="coerce").fillna(0)
        if "F1" in c:
            summary_df[c] = summary_df[c].round(4).astype(float)
        else:
            summary_df[c] = summary_df[c].astype(int)

    # Sort summary so TOTAL is last
    summary_df["__is_total"] = summary_df["Reference"].eq("TOTAL")
    summary_df = summary_df.sort_values("__is_total").drop(columns="__is_total").reset_index(drop=True)

    # Write all sheets
    with pd.ExcelWriter(out_xlsx_path, engine="openpyxl") as writer:
        confusion_df.to_excel(writer, sheet_name="confusion_matrix", index=False)
        summary_df.to_excel(writer, sheet_name="confusion_summary", index=False)
        fp_reads_df.to_excel(writer, sheet_name="remaining_false_positives", index=False)


def build_ref_meta_from_match_file(
    matcher_path: str,
    accessioncol: int = 0,
    namecol: int = 2,
    taxcol: int = 4,
    has_header: bool = True,
    sep: str = "\t",
) -> Dict[str, Dict[str, str]]:
    """
    Build ref_meta from your --match TSV file.
    Returns:
      ref_meta[accession] = {"taxid": "...", "organism": "..."}
      and also stores a no-version key (e.g., NC_003310)
    """
    ref_meta: Dict[str, Dict[str, str]] = {}

    if not matcher_path or not os.path.exists(matcher_path):
        return ref_meta

    with open(matcher_path, "r") as f:
        for i, line in enumerate(f):
            line = line.rstrip("\n")
            if not line:
                continue
            if has_header and i == 0:
                continue
            parts = line.split(sep)

            accession = parts[accessioncol].strip() if len(parts) > accessioncol else ""
            name = parts[namecol].strip() if len(parts) > namecol else ""
            taxid = parts[taxcol].strip() if len(parts) > taxcol else ""

            if not accession:
                continue

            md = {
                "taxid": taxid if taxid not in ["", "None", None] else None,
                "organism": name if name not in ["", "None", None] else accession,
            }
            ref_meta[accession] = md
            ref_meta[_strip_version(accession)] = md

    return ref_meta


def compute_global_confusion(read_to_refs: Dict[str, Set[str]], all_reads: Set[str]):
    """
    Compute global TP/FP/FN/TN at the READ level.
    A read is TP if it aligns to its truth ref at least once.
    """
    tp = fp = fn = tn = 0

    for rid in all_reads:
        truth = ground_truth_from_read_id(rid)
        aligned_refs = read_to_refs.get(rid, set())

        if truth in aligned_refs:
            tp += 1
        elif aligned_refs:
            fp += 1
        else:
            fn += 1

    tn = 0  # not defined in this context
    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = (2 * precision * recall / (precision + recall)) if (precision + recall) else 0.0

    return {"TP": tp, "FP": fp, "FN": fn, "TN": tn, "Precision": precision, "Recall": recall, "F1": f1}
def optimize_three_weights(
    bam_path: str,
    match_tsv: str,
    aggregated_stats: dict,
    pathogens: dict,
    sample_name: str,
    sample_type: str,
    total_reads: float,
    aligned_total: int,
    initial=(0.33, 0.33, 0.34),
    lambda_fp_max: float = 50.0,
    lambda_fp_mean: float = 10.0,
    lambda_tp: float = 1.0,
    maxiter: int = 200,
    verbose: bool = True,
):
    """
    Optimize ONLY 3 weights (breadth_weight, gini_weight, minhash_weight) at TAXID level.

    Goal:
      - True-positive taxids should have TASS close to 1
      - False-positive taxids should have TASS as close to 0 as possible (HARSHLY)
        (even if only 1 FP read exists)

    This uses the BAM to derive:
      - truth_taxid per read (from read_id -> truth accession -> accession->taxid)
      - predicted taxids per read (aligned ref accessions -> accession->taxid)
      - TP_taxids: truth taxid that is hit by at least one aligned ref with same taxid
      - FP_taxids: any predicted taxid != truth taxid for any read

    Requires your existing functions in scope:
      - collect_read_alignments(bam_path, alignments_to_remove=None)
      - ground_truth_from_read_id(read_id)
      - calculate_scores(aggregated_stats, pathogens, sample_name, sample_type, total_reads, aligned_total, weights)
    """

    import os
    import math
    import numpy as np
    from collections import defaultdict
    from scipy.optimize import minimize

    # -------- helpers --------
    def _strip_version(acc: str) -> str:
        acc = (acc or "").strip()
        return acc.split(".", 1)[0] if "." in acc else acc

    def _build_acc_to_taxid(match_path: str, accessioncol=0, taxcol=4, has_header=True, sep="\t"):
        """
        Build accession -> taxid mapping.
        Stores BOTH versioned and no-version keys.
        """
        acc_to_taxid = {}
        if not match_path or not os.path.exists(match_path):
            return acc_to_taxid

        with open(match_path, "r") as fh:
            for i, line in enumerate(fh):
                line = line.rstrip("\n")
                if not line:
                    continue
                if has_header and i == 0:
                    continue
                parts = line.split(sep)
                if len(parts) <= max(accessioncol, taxcol):
                    continue
                acc = (parts[accessioncol] or "").strip()
                taxid = (parts[taxcol] or "").strip()
                if not acc or taxid in ["", "None", None]:
                    continue
                taxid = str(taxid)
                acc_to_taxid[acc] = taxid
                acc_to_taxid[_strip_version(acc)] = taxid
        return acc_to_taxid

    def _derive_tp_fp_taxids(bam_path: str, acc_to_taxid: dict):
        """
        Returns:
          tp_taxids: set
          fp_taxids: set
          stats dict with counts for debugging
        """
        read_to_refs, _, all_reads = collect_read_alignments(bam_path, alignments_to_remove=None)

        tp_taxids = set()
        fp_taxids = set()

        n_reads_total = 0
        n_reads_mapped_truth = 0
        n_reads_skipped_no_truth_taxid = 0
        n_alignments_total = 0
        n_alignments_skipped_no_pred_taxid = 0

        for rid in all_reads:
            n_reads_total += 1

            truth_acc = ground_truth_from_read_id(rid)
            truth_taxid = acc_to_taxid.get(truth_acc) or acc_to_taxid.get(_strip_version(truth_acc))
            if not truth_taxid:
                n_reads_skipped_no_truth_taxid += 1
                continue
            truth_taxid = str(truth_taxid)
            n_reads_mapped_truth += 1

            aligned_refs = read_to_refs.get(rid, set())

            # determine if read is TP at taxid-level (any aligned ref has same taxid)
            is_tp_read = False
            for ref_acc in aligned_refs:
                n_alignments_total += 1
                pred_taxid = acc_to_taxid.get(ref_acc) or acc_to_taxid.get(_strip_version(ref_acc))
                if not pred_taxid:
                    n_alignments_skipped_no_pred_taxid += 1
                    continue
                pred_taxid = str(pred_taxid)

                if pred_taxid == truth_taxid:
                    is_tp_read = True
                else:
                    fp_taxids.add(pred_taxid)

            if is_tp_read:
                tp_taxids.add(truth_taxid)

        dbg = dict(
            n_reads_total=n_reads_total,
            n_reads_mapped_truth=n_reads_mapped_truth,
            n_reads_skipped_no_truth_taxid=n_reads_skipped_no_truth_taxid,
            n_alignments_total=n_alignments_total,
            n_alignments_skipped_no_pred_taxid=n_alignments_skipped_no_pred_taxid,
            n_tp_taxids=len(tp_taxids),
            n_fp_taxids=len(fp_taxids),
        )
        return tp_taxids, fp_taxids, dbg

    def _score_by_taxid(weights3):
        """
        weights3: (breadth, gini, minhash) -> dict taxid->tass_score
        """
        bw, gw, mw = weights3
        # enforce sum=1 (robustness)
        s = bw + gw + mw
        if s <= 0:
            bw, gw, mw = 1/3, 1/3, 1/3
        else:
            bw, gw, mw = bw/s, gw/s, mw/s

        weights = {
            # only 3 active
            "breadth_weight": float(bw),
            "gini_coefficient": float(gw),
            "minhash_weight": float(mw),

            # everything else OFF
            "mapq_score": 0.0,
            "disparity_score": 0.0,
            "hmp_weight": 0.0,
            "siblings_score": 0.0,
            "diamond_identity": 0.0,
            "k2_disparity_weight": 0.0,
        }

        final_scores = calculate_scores(
            aggregated_stats=aggregated_stats,
            pathogens=pathogens,
            sample_name=sample_name,
            sample_type=sample_type,
            total_reads=total_reads,
            aligned_total=aligned_total,
            weights=weights,
        )

        out = {}
        for rec in final_scores:
            # IMPORTANT: here 'ref' must be TAXID key (string)
            t = rec.get("ref")
            if t is None:
                continue
            out[str(t)] = float(rec.get("tass_score", 0.0))
        return out

    def _loss(weights3, tp_taxids, fp_taxids):
        """
        Loss that:
          - punishes any FP taxid with high TASS, extremely (max-based)
          - also penalizes average FP TASS
          - also encourages TP taxids to have TASS near 1
        """
        scores = _score_by_taxid(weights3)

        tp_scores = [scores.get(t, 0.0) for t in tp_taxids]
        fp_scores = [scores.get(t, 0.0) for t in fp_taxids]

        # If sets are empty, keep stable behavior but don't make loss constant-zero silently.
        if len(tp_scores) == 0:
            tp_term = 1.0  # strong penalty: "we can't reward TP if none map"
        else:
            tp_term = float(np.mean([(1.0 - s) ** 2 for s in tp_scores]))

        if len(fp_scores) == 0:
            fp_max_term = 0.0
            fp_mean_term = 0.0
        else:
            fp_max = float(np.max(fp_scores))
            fp_max_term = fp_max ** 2
            fp_mean_term = float(np.mean([s ** 2 for s in fp_scores]))

        # Total loss (minimize)
        return (lambda_tp * tp_term) + (lambda_fp_max * fp_max_term) + (lambda_fp_mean * fp_mean_term)

    # -------- build GT taxid sets from BAM --------
    acc_to_taxid = _build_acc_to_taxid(match_tsv, accessioncol=0, taxcol=1, has_header=True, sep="\t")
    if not acc_to_taxid:
        raise ValueError(f"Could not build acc->taxid mapping from match file: {match_tsv}")

    tp_taxids, fp_taxids, dbg = _derive_tp_fp_taxids(bam_path, acc_to_taxid)

    if verbose:
        print("=== TAXID GT sets derived from BAM ===")
        for k, v in dbg.items():
            print(f"{k}: {v}")
        # sanity: show a few
        print("TP taxids sample:", list(sorted(tp_taxids))[:10])
        print("FP taxids sample:", list(sorted(fp_taxids))[:10])

    # -------- optimize with SLSQP --------
    x0 = np.array(initial, dtype=float)
    if np.any(x0 < 0):
        x0 = np.clip(x0, 0, None)
    if x0.sum() <= 0:
        x0 = np.array([1/3, 1/3, 1/3], dtype=float)
    else:
        x0 = x0 / x0.sum()

    bounds = [(0.0, 1.0), (0.0, 1.0), (0.0, 1.0)]
    cons = [{"type": "eq", "fun": lambda w: float(w[0] + w[1] + w[2] - 1.0)}]

    # baseline
    baseline_loss = _loss(x0, tp_taxids, fp_taxids)

    def obj(w):
        # minimize loss
        return _loss(w, tp_taxids, fp_taxids)

    res = minimize(
        obj,
        x0,
        method="SLSQP",
        bounds=bounds,
        constraints=cons,
        options={"maxiter": int(maxiter), "disp": bool(verbose)},
    )

    wopt = res.x
    wopt = np.clip(wopt, 0, 1)
    wopt = wopt / wopt.sum() if wopt.sum() > 0 else np.array([1/3, 1/3, 1/3], dtype=float)

    opt_loss = _loss(wopt, tp_taxids, fp_taxids)

    # compute post scores for reporting
    scores_opt = _score_by_taxid(wopt)
    tp_scores_opt = [scores_opt.get(t, 0.0) for t in tp_taxids]
    fp_scores_opt = [scores_opt.get(t, 0.0) for t in fp_taxids]

    report = dict(
        success=bool(res.success),
        message=str(res.message),
        n_iter=int(getattr(res, "nit", -1)),
        baseline_loss=float(baseline_loss),
        optimized_loss=float(opt_loss),
        initial_weights=dict(breadth_weight=float(x0[0]), gini_weight=float(x0[1]), minhash_weight=float(x0[2])),
        optimized_weights=dict(breadth_weight=float(wopt[0]), gini_weight=float(wopt[1]), minhash_weight=float(wopt[2])),
        tp_taxids=tp_taxids,
        fp_taxids=fp_taxids,
        tp_mean_score=float(np.mean(tp_scores_opt)) if tp_scores_opt else 0.0,
        fp_mean_score=float(np.mean(fp_scores_opt)) if fp_scores_opt else 0.0,
        fp_max_score=float(np.max(fp_scores_opt)) if fp_scores_opt else 0.0,
    )

    if verbose:
        print("\n=== 3-weight TAXID optimization summary ===")
        print("Initial weights:", report["initial_weights"])
        print(f"Baseline loss: {report['baseline_loss']:.6f}")
        print("Optimized weights:", report["optimized_weights"])
        print(f"Optimized loss: {report['optimized_loss']:.6f}")
        print(f"TP mean TASS: {report['tp_mean_score']:.4f}")
        print(f"FP mean TASS: {report['fp_mean_score']:.4f}")
        print(f"FP max  TASS: {report['fp_max_score']:.4f}")

    return report

def parse_args(argv=None):
    p = argparse.ArgumentParser(description="Build confusion metrics and per-read remaining false positives.")
    p.add_argument("--input", "-i", required=True, help="Input BAM (original alignments).")
    p.add_argument("--filtered_bam", required=False, default=None, help="Filtered BAM (post-filter alignments).")
    p.add_argument("--match", required=False, default=None, help="TSV file mapping accession -> name -> taxid (tab-separated).")
    p.add_argument("--accessioncol", type=int, default=0, help="Column index for accession in match file (0-based).")
    p.add_argument("--namecol", type=int, default=2, help="Column index for organism name in match file (0-based).")
    p.add_argument("--taxcol", type=int, default=1, help="Column index for taxid in match file (0-based).")
    p.add_argument("--output_xlsx", required=False, default=None, help="Output XLSX path. Default: <input_dir>/alignment_confusion_report.xlsx")
    p.add_argument("--report_confusion_xlsx", action="store_true", help="Write XLSX report.")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    # ensure BAM exists
    if not os.path.exists(args.input):
        raise FileNotFoundError(f"Input BAM not found: {args.input}")

    out_dir = os.path.dirname(args.input)
    if args.output_xlsx:
        out_xlsx = args.output_xlsx
    else:
        out_xlsx = os.path.join(out_dir, "alignment_confusion_report.xlsx")
    outputdir=os.path.dirname(args.input)
    outputfile = os.path.join(outputdir, "output.txt")
    matchfile = os.path.join(outputdir, "..", "orgs.txt")
    if not args.match:
        args.match = matchfile
    # build ref_meta (optional)
    ref_meta = {}
    if args.match:
        ref_meta = build_ref_meta_from_match_file(
            matcher_path=args.match,
            accessioncol=args.accessioncol,
            namecol=args.namecol,
            taxcol=args.taxcol,
            has_header=True,
            sep="\t",
        )

    # collect original alignments
    read_to_refs_orig, read_to_seq, all_reads = collect_read_alignments(args.input, alignments_to_remove=None)

    if not args.filtered_bam:
        filtered_bam = os.path.join(outputdir, "search_results", "filtered.hash.bam")
        args.filtered_bam = filtered_bam
    # collect filtered alignments if provided; otherwise apply no removal (treated as identical to orig)
    if args.filtered_bam and os.path.exists(args.filtered_bam):
        read_to_refs_filt, _, all_reads2 = collect_read_alignments(args.filtered_bam, alignments_to_remove=None)
        # merge read sets (in case some reads exist only in one of the BAMs)
        all_reads = set(all_reads).union(set(all_reads2))
    else:
        # no separate filtered BAM provided -> treat filtered as same as original (no change)
        read_to_refs_filt = read_to_refs_orig

    # build per-reference confusion tables
    conf_old = build_confusion_dataframe(read_to_refs_orig, all_reads, ref_meta=ref_meta, prefix="OLD_")
    conf_new = build_confusion_dataframe(read_to_refs_filt, all_reads, ref_meta=ref_meta, prefix="NEW_")

    confusion_df = conf_old.merge(
        conf_new.drop(columns=["Reference_TaxID", "Reference_Organism"], errors="ignore"),
        on="Reference",
        how="outer",
    )

    # compute TOTAL global row
    old_global = compute_global_confusion(read_to_refs_orig, all_reads)
    new_global = compute_global_confusion(read_to_refs_filt, all_reads)

    total_row = {
        "Reference": "TOTAL",
        "Reference_TaxID": None,
        "Reference_Organism": "ALL_READS",
        "OLD_TP": old_global["TP"],
        "OLD_FP": old_global["FP"],
        "OLD_FN": old_global["FN"],
        "OLD_TN": old_global["TN"],
        "OLD_Precision": old_global["Precision"],
        "OLD_Recall": old_global["Recall"],
        "OLD_F1": old_global["F1"],
        "NEW_TP": new_global["TP"],
        "NEW_FP": new_global["FP"],
        "NEW_FN": new_global["FN"],
        "NEW_TN": new_global["TN"],
        "NEW_Precision": new_global["Precision"],
        "NEW_Recall": new_global["Recall"],
        "NEW_F1": new_global["F1"],
    }

    confusion_df = pd.concat([confusion_df, pd.DataFrame([total_row])], ignore_index=True)

    # build per-read remaining false positives sheet
    fp_reads_df = build_false_positive_reads_df(
        read_to_refs_orig, read_to_refs_filt, read_to_seq, all_reads, ref_meta=ref_meta
    )

    # write XLSX if requested
    if args.report_confusion_xlsx:
        write_confusion_xlsx(out_xlsx, confusion_df, fp_reads_df)
        print(f"Wrote confusion XLSX report: {out_xlsx}")

    # print confusion_df to stdout
    # pretty print using pandas
    pd.set_option("display.max_rows", None)
    pd.set_option("display.max_colwidth", None)
    print("\n=== Confusion DataFrame (including TOTAL row) ===\n")
    print(confusion_df.to_string(index=False))

    # return results for programmatic use
    return {"confusion_df": confusion_df, "fp_reads_df": fp_reads_df, "ref_meta": ref_meta}


if __name__ == "__main__":
    main()


