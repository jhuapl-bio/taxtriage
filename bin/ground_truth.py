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
import numpy as np
import json

from collections import defaultdict
from typing import Dict, Iterable, Optional, Set, Tuple
import pandas as pd
import pysam
from optimize_weights import compute_tass_score_from_metrics


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
    # sort on the Reference_Organism column alphabetically, keeping TOTAL last
    summary_df = summary_df.sort_values(by=["Reference_Organism"], na_position="last").reset_index(drop=True)
    confusion_df = confusion_df.sort_values(by=["Reference_Organism"]).reset_index(drop=True)
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

def optimize_weights(
    input_bam=None,
    final_json=None,
    accession_to_taxid=None,
    breadth_weight=1/3,
    minhash_weight=1/3,
    gini_weight=1/3,
    alpha=1.0,
    optimize_pos_weight=1.0,
    optimize_neg_weight=1.0,
    optimize_reg=0.0,
    optimize_seed=42,
    optimize_maxiter=100,
    optimize_local_maxiter=50,
    optimize_mode="taxids", # taxids | reads | hybrid
    hybrid_lambda=0.5,
    optimize_report="optimization_report.json",
):
    import json
    import numpy as np

    tp_fp_counts = compute_tp_fp_counts_by_taxid(
        input_bam,
        accession_to_taxid=accession_to_taxid,
        debug_n=10,
    )

    metrics_df = build_metrics_df_from_final_json(final_json, tp_fp_counts)
    accession_col = "taxid"

    start_w = dict(
        breadth_weight=float(breadth_weight),
        minhash_weight=float(minhash_weight),
        gini_weight=float(gini_weight),
    )

    opt = optimize_weights_for_tp_fp(
        metrics_df=metrics_df,
        tp_fp_counts=tp_fp_counts,
        accession_col=accession_col,
        start_weights=start_w,
        alpha=float(alpha),
        pos_weight=float(optimize_pos_weight),
        neg_weight=float(optimize_neg_weight),
        reg_lambda=float(optimize_reg),
        seed=int(optimize_seed),
        maxiter_global=int(optimize_maxiter),
        maxiter_local=int(optimize_local_maxiter),
        optimize_mode=optimize_mode,
        hybrid_lambda=float(hybrid_lambda),
    )

    print("[optimize] status:", opt["status"])
    print("[optimize] mode:", opt.get("optimize_mode"), "hybrid_lambda:", opt.get("hybrid_lambda"))
    print("[optimize] TP reads:", opt.get("tp_total_reads"), "FP reads:", opt.get("fp_total_reads"))
    print("[optimize] TP taxids:", opt.get("tp_total_taxa"), "FP taxids:", opt.get("fp_total_taxa"))

    if opt["status"] != "ok":
        report_obj = {
            "status": opt["status"],
            "config": {
                "alpha": float(alpha),
                "pos_weight": float(optimize_pos_weight),
                "neg_weight": float(optimize_neg_weight),
                "reg_lambda": float(optimize_reg),
                "seed": int(optimize_seed),
                "maxiter_global": int(optimize_maxiter),
                "maxiter_local": int(optimize_local_maxiter),
                "optimize_mode": optimize_mode,
                "hybrid_lambda": float(hybrid_lambda),
                "start_weights": start_w,
            },
            "optimizer": opt,
            "best_weights": start_w,
            "report": [],
        }
        with open(optimize_report, "w") as f:
            json.dump(report_obj, f, indent=2)
        return report_obj

    best_w = opt["weights"]
    breadth_weight = best_w["breadth_weight"]
    minhash_weight = best_w["minhash_weight"]
    gini_weight = best_w["gini_weight"]

    scores = compute_tass_score_from_metrics(
        metrics_df,
        breadth_w=float(breadth_weight),
        minhash_w=float(minhash_weight),
        gini_w=float(gini_weight),
        alpha=float(alpha),
    )
    scores = np.asarray(scores, dtype=float)

    tp = metrics_df["tp_reads"].to_numpy(dtype=int)
    fp = metrics_df["fp_reads"].to_numpy(dtype=int)
    total = tp + fp

    report_rows = []
    for i in range(len(metrics_df)):
        report_rows.append({
            "taxid": str(metrics_df.iloc[i]["taxid"]),
            "name": str(metrics_df.iloc[i].get("name", "")),
            "tp_reads": int(tp[i]),
            "fp_reads": int(fp[i]),
            "total_reads": int(total[i]),
            "fp_fraction": float(fp[i] / total[i]) if total[i] else 0.0,
            "tass_score": float(scores[i]),
            "features": {
                "breadth_log_score": float(metrics_df.iloc[i].get("breadth_log_score", 0.0)),
                "minhash_reduction": float(metrics_df.iloc[i].get("minhash_reduction", 0.0)),
                "gini_coefficient": float(metrics_df.iloc[i].get("gini_coefficient", 0.0)),
            }
        })

    report_rows.sort(key=lambda r: (r["fp_reads"], -r["tp_reads"], r["tass_score"]), reverse=True)

    report_obj = {
        "status": "ok",
        "config": {
            "alpha": float(alpha),
            "pos_weight": float(optimize_pos_weight),
            "neg_weight": float(optimize_neg_weight),
            "reg_lambda": float(optimize_reg),
            "seed": int(optimize_seed),
            "maxiter_global": int(optimize_maxiter),
            "maxiter_local": int(optimize_local_maxiter),
            "optimize_mode": optimize_mode,
            "hybrid_lambda": float(hybrid_lambda),
            "start_weights": start_w,
        },
        "best_weights": {
            "breadth_weight": float(breadth_weight),
            "minhash_weight": float(minhash_weight),
            "gini_weight": float(gini_weight),
        },
        "optimizer": opt,
        "report": report_rows,
    }

    with open(optimize_report, "w") as f:
        json.dump(report_obj, f, indent=2)

    return report_obj



def _clip01(x: np.ndarray, eps: float = 1e-6) -> np.ndarray:
    return np.clip(x, eps, 1.0 - eps)
def optimize_weights_for_tp_fp(
    metrics_df,
    tp_fp_counts: dict[str, tuple[int, int]],
    *,
    accession_col: str,
    start_weights: dict,
    alpha: float,
    pos_weight: float,
    neg_weight: float,
    reg_lambda: float,
    seed: int,
    maxiter_global: int,
    maxiter_local: int,
    optimize_mode: str = "taxids",   # "reads" | "taxids" | "hybrid"
    hybrid_lambda: float = 0.5,
):
    from scipy.optimize import differential_evolution, minimize
    import numpy as np

    accessions = metrics_df[accession_col].astype(str).tolist()
    tp = np.array([tp_fp_counts.get(a, (0, 0))[0] for a in accessions], dtype=float)
    fp = np.array([tp_fp_counts.get(a, (0, 0))[1] for a in accessions], dtype=float)

    y_pos = (tp > 0).astype(float)
    y_neg = (fp > 0).astype(float)

    tp_total_reads = float(tp.sum())
    fp_total_reads = float(fp.sum())
    tp_total_taxa = float(y_pos.sum())
    fp_total_taxa = float(y_neg.sum())

    # mode feasibility checks
    if optimize_mode == "reads":
        if tp_total_reads <= 0 or fp_total_reads <= 0:
            return {"weights": start_weights, "status": "skipped (need both TP and FP reads)"}
    elif optimize_mode == "taxids":
        if tp_total_taxa <= 0 or fp_total_taxa <= 0:
            return {"weights": start_weights, "status": "skipped (need both TP and FP taxids)"}
    elif optimize_mode == "hybrid":
        if (tp_total_reads <= 0 and tp_total_taxa <= 0) or (fp_total_reads <= 0 and fp_total_taxa <= 0):
            return {"weights": start_weights, "status": "skipped (need TP/FP signal)"}
    else:
        raise ValueError(f"Unknown optimize_mode: {optimize_mode}")

    # normalize class weights
    s = pos_weight + neg_weight
    pos_w = pos_weight / s
    neg_w = neg_weight / s

    w0 = np.array([
        start_weights["breadth_weight"],
        start_weights["minhash_weight"],
        start_weights["gini_weight"],
    ], dtype=float)

    w0 = np.clip(w0, 0.0, 1.0)
    w0 = (w0 / w0.sum()) if w0.sum() > 0 else np.array([1/3, 1/3, 1/3], dtype=float)

    def loss_for_w(w: np.ndarray) -> float:
        # simplex constraints (soft)
        if np.any(w < 0) or np.any(w > 1) or abs(w.sum() - 1.0) > 1e-6:
            return 1e6 + 1e6 * (w.sum() - 1.0) ** 2

        scores = compute_tass_score_from_metrics(
            metrics_df,
            breadth_w=float(w[0]),
            minhash_w=float(w[1]),
            gini_w=float(w[2]),
            alpha=alpha,
        )
        p = _clip01(np.asarray(scores, dtype=float))

        # read-weighted terms
        read_pos_term = (-np.sum(tp * np.log(p)) / tp_total_reads) if tp_total_reads > 0 else 0.0
        read_neg_term = (-np.sum(fp * np.log(1.0 - p)) / fp_total_reads) if fp_total_reads > 0 else 0.0

        # taxid-level terms
        tax_pos_term = (-np.sum(y_pos * np.log(p)) / tp_total_taxa) if tp_total_taxa > 0 else 0.0
        tax_neg_term = (-np.sum(y_neg * np.log(1.0 - p)) / fp_total_taxa) if fp_total_taxa > 0 else 0.0

        if optimize_mode == "reads":
            pos_term, neg_term = read_pos_term, read_neg_term
        elif optimize_mode == "taxids":
            pos_term, neg_term = tax_pos_term, tax_neg_term
        else:
            lam = max(0.0, min(1.0, float(hybrid_lambda)))
            pos_term = lam * tax_pos_term + (1.0 - lam) * read_pos_term
            neg_term = lam * tax_neg_term + (1.0 - lam) * read_neg_term

        reg = reg_lambda * float(np.sum((w - w0) ** 2))
        return pos_w * pos_term + neg_w * neg_term + reg

    # --- Global stage (DE) over 2 vars, 3rd implied by sum=1 ---
    # x = [w_breadth, w_minhash], w_gini = 1 - sum(x)
    def unpack_x(x: np.ndarray) -> np.ndarray:
        wb, wm = float(x[0]), float(x[1])
        wg = 1.0 - (wb + wm)
        return np.array([wb, wm, wg], dtype=float)

    def loss_for_x(x: np.ndarray) -> float:
        w = unpack_x(x)
        if np.any(w < 0.0) or np.any(w > 1.0):
            return 1e6 + 1e6 * max(0.0, -w.min())
        return loss_for_w(w)

    bounds = [(0.0, 1.0), (0.0, 1.0)]
    rng = np.random.RandomState(seed)

    de_res = differential_evolution(
        loss_for_x,
        bounds=bounds,
        maxiter=maxiter_global,
        seed=rng,
        polish=False,
        updating="deferred",
        workers=1,
    )

    w_de = unpack_x(de_res.x)
    if np.any(w_de < 0) or abs(w_de.sum() - 1.0) > 1e-6:
        w_de = w0.copy()
    else:
        w_de = np.clip(w_de, 0.0, 1.0)
        w_de = w_de / w_de.sum()

    # --- Local refinement (SLSQP) on 3 vars ---
    cons = [{"type": "eq", "fun": lambda w: np.sum(w) - 1.0}]
    bnds = [(0.0, 1.0)] * 3

    local_res = minimize(
        loss_for_w,
        x0=w_de,
        method="SLSQP",
        bounds=bnds,
        constraints=cons,
        options={"maxiter": maxiter_local, "ftol": 1e-9},
    )

    w_best = np.clip(local_res.x, 0.0, 1.0)
    w_best = w_best / w_best.sum()

    # summaries (both read and taxid views)
    scores_best = compute_tass_score_from_metrics(
        metrics_df,
        breadth_w=float(w_best[0]),
        minhash_w=float(w_best[1]),
        gini_w=float(w_best[2]),
        alpha=alpha,
    )
    scores_best = np.asarray(scores_best, dtype=float)

    tp_mean_reads = float(np.sum(tp * scores_best) / tp_total_reads) if tp_total_reads > 0 else None
    fp_mean_reads = float(np.sum(fp * scores_best) / fp_total_reads) if fp_total_reads > 0 else None
    tp_mean_taxa = float(np.sum(y_pos * scores_best) / tp_total_taxa) if tp_total_taxa > 0 else None
    fp_mean_taxa = float(np.sum(y_neg * scores_best) / fp_total_taxa) if fp_total_taxa > 0 else None

    return {
        "weights": {
            "breadth_weight": float(w_best[0]),
            "minhash_weight": float(w_best[1]),
            "gini_weight": float(w_best[2]),
        },
        "status": "ok",
        "optimize_mode": optimize_mode,
        "hybrid_lambda": float(hybrid_lambda),
        "tp_total_reads": int(tp_total_reads),
        "fp_total_reads": int(fp_total_reads),
        "tp_total_taxa": int(tp_total_taxa),
        "fp_total_taxa": int(fp_total_taxa),
        "tp_mean_score_reads": tp_mean_reads,
        "fp_mean_score_reads": fp_mean_reads,
        "tp_mean_score_taxa": tp_mean_taxa,
        "fp_mean_score_taxa": fp_mean_taxa,
        "loss": float(loss_for_w(w_best)),
        "global_fun": float(de_res.fun),
        "local_success": bool(local_res.success),
        "local_message": str(local_res.message),
    }


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


_ACCESSION_RE = re.compile(r"([A-Z]{1,4}_[A-Z0-9]{3,}\.\d+|[A-Z]{1,4}[0-9]{5,}\.\d+)")

def canonical_accession(s: str) -> Optional[str]:
    """
    Extract canonical accession from read ids / BAM reference names.

    Examples:
      "@NC_003310.1_3323_4" -> "NC_003310.1"
      "NC_003310.1|Monkeypox virus" -> "NC_003310.1"
      "NC_003310.1 some desc" -> "NC_003310.1"
    """
    if not s:
        return None
    s = s.strip()
    if s.startswith("@"):
        s = s[1:]

    m = _ACCESSION_RE.search(s)
    if m:
        return m.group(1)

    # fallback to your original rule
    return s.split("_", 1)[0]


def compute_tp_fp_counts_by_taxid(
    bam_path: str,
    accession_to_taxid: dict,
    *,
    skip_secondary: bool = True,
    skip_supplementary: bool = True,
    debug_n: int = 10,
):
    """
    Count TP/FP reads per *aligned taxid*.

    For each alignment:
      truth_acc = canonical_accession(QNAME)
      aligned_acc = canonical_accession(RNAME)
      truth_taxid = accession_to_taxid.get(truth_acc)
      aligned_taxid = accession_to_taxid.get(aligned_acc)

    If aligned_taxid is missing -> skip (can't attribute)
    If truth_taxid is missing -> count as FP against aligned_taxid (conservative), OR skip if you prefer.

    TP if truth_taxid == aligned_taxid else FP

    Returns:
      dict[taxid_str] = (tp_reads, fp_reads)
    """
    counts = defaultdict(lambda: [0, 0])
    debug = []

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for aln in bam.fetch(until_eof=True):
            if aln.is_unmapped:
                continue
            if skip_secondary and aln.is_secondary:
                continue
            if skip_supplementary and aln.is_supplementary:
                continue

            raw_ref = bam.get_reference_name(aln.reference_id)

            aligned_acc = canonical_accession(raw_ref)
            truth_acc = canonical_accession(aln.query_name)

            if not aligned_acc or not truth_acc:
                continue

            aligned_taxid = accession_to_taxid.get(aligned_acc)
            truth_taxid = accession_to_taxid.get(truth_acc)

            # if we can't map aligned ref -> taxid, we can't assign this alignment
            if aligned_taxid is None:
                continue

            aligned_taxid = str(aligned_taxid)
            truth_taxid = str(truth_taxid) if truth_taxid is not None else None

            if truth_taxid is not None and truth_taxid == aligned_taxid:
                counts[aligned_taxid][0] += 1
            else:
                counts[aligned_taxid][1] += 1
                if len(debug) < debug_n:
                    debug.append((truth_acc, truth_taxid, aligned_acc, aligned_taxid))

    if debug:
        print("[debug] first mismatches (truth_acc, truth_taxid -> aligned_acc, aligned_taxid):")
        for tacc, ttx, aacc, atx in debug:
            print(f"  {tacc} ({ttx})  ->  {aacc} ({atx})")

    return {k: (v[0], v[1]) for k, v in counts.items()}


def build_metrics_df_from_final_json(
    final_json: dict,  # <-- CHANGED: now dict[taxid] -> stats dict
    tp_fp_by_taxid: dict[str, tuple[int, int]],
):
    """
    Build per-taxid metrics_df from a 1D dict final_json and merge TP/FP counts.

    final_json format (new):
      {
        "28871": {"name": "...", "breadth_log_score": ..., "minhash_reduction": ..., ...},
        "2200830": {...},
        ...
      }

    Returns DataFrame with:
      taxid, name, accessions, breadth_log_score, minhash_reduction, disparity_score, gini_coefficient,
      tp_reads, fp_reads, total_reads, fp_fraction
    """
    rows = []

    for taxid_key, stats in (final_json or {}).items():
        if stats is None:
            continue

        taxid = str(taxid_key).strip()
        if not taxid:
            continue

        # Accessions may be present but often empty in your new structure
        accs = stats.get("accessions") or []
        accs = sorted({canonical_accession(a) for a in accs if canonical_accession(a)})

        # Features (match your optimizer targets)
        breadth = float(stats.get("breadth_log_score", 0.0))
        minhash = float(stats.get("minhash_reduction", stats.get("minhash_score", 0.0)))

        # Prefer k2_disparity_score; fall back to disparity if present
        disparity = float(stats.get("k2_disparity_score", stats.get("disparity", 0.0)))

        gini = float(stats.get("gini_coefficient", 0.0))

        tp, fp = tp_fp_by_taxid.get(taxid, (0, 0))
        total = tp + fp

        rows.append({
            "taxid": taxid,
            "name": stats.get("name", ""),
            "accessions": accs,
            "breadth_log_score": breadth,
            "minhash_reduction": minhash,
            "disparity_score": disparity,
            "gini_coefficient": gini,
            "tp_reads": int(tp),
            "fp_reads": int(fp),
            "total_reads": int(total),
            "fp_fraction": (fp / total) if total else 0.0,
        })

    df = pd.DataFrame(rows)
    if df.empty:
        return df

    # Defensive: if taxid somehow duplicated, consolidate
    def union_lists(series):
        s = set()
        for lst in series:
            s.update(lst or [])
        return sorted(s)

    agg = {
        "name": "first",
        "accessions": union_lists,
        "breadth_log_score": "max",
        "minhash_reduction": "max",
        "disparity_score": "max",
        "gini_coefficient": "max",
        "tp_reads": "sum",
        "fp_reads": "sum",
        "total_reads": "sum",
    }

    df = df.groupby("taxid", as_index=False).agg(agg)
    df["fp_fraction"] = df.apply(
        lambda r: (r["fp_reads"] / r["total_reads"]) if r["total_reads"] else 0.0,
        axis=1
    )

    # Helpful ordering for debugging
    df = df.sort_values(["fp_reads", "tp_reads"], ascending=[False, False]).reset_index(drop=True)
    return df

def build_ground_truth_metrics_df(
    bam_path: str,
    *,
    skip_secondary: bool = True,
    skip_supplementary: bool = True,
    debug_n: int = 10,
) -> pd.DataFrame:
    """
    Build per-aligned-reference TP/FP counts using BAM-derived ground truth.

    Truth accession is extracted from QNAME (first accession-looking token).
    Aligned accession is extracted from RNAME (reference name in BAM).
      TP if aligned_accession == truth_accession else FP

    Returns DataFrame:
      accession | tp_reads | fp_reads | total_reads | fp_fraction
    """

    counts = defaultdict(lambda: [0, 0])  # aligned_acc -> [tp, fp]
    debug_examples: list[Tuple[str, str]] = []  # (truth, aligned)

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for aln in bam.fetch(until_eof=True):
            if aln.is_unmapped:
                continue
            if skip_secondary and aln.is_secondary:
                continue
            if skip_supplementary and aln.is_supplementary:
                continue

            # aligned ref name from BAM header
            raw_ref = bam.get_reference_name(aln.reference_id)
            aligned_acc = _canonical_accession(raw_ref)

            # truth accession from read id
            truth_acc = _canonical_accession(aln.query_name)

            # If we can't canonicalize, skip (or you could count separately)
            if not aligned_acc or not truth_acc:
                continue

            if aligned_acc == truth_acc:
                counts[aligned_acc][0] += 1
            else:
                counts[aligned_acc][1] += 1
                if len(debug_examples) < debug_n:
                    debug_examples.append((truth_acc, aligned_acc))

    if debug_examples:
        print("[debug] first mismatches (truth -> aligned):")
        for t, a in debug_examples:
            print(f"  {t}  ->  {a}")

    rows = []
    for accession, (tp, fp) in counts.items():
        total = tp + fp
        rows.append({
            "accession": accession,
            "tp_reads": tp,
            "fp_reads": fp,
            "total_reads": total,
            "fp_fraction": (fp / total) if total else 0.0,
        })

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(["fp_reads", "tp_reads"], ascending=[False, False]).reset_index(drop=True)

    return df

# A reasonably broad accession matcher (NC_, NZ_, CP, etc.), with version if present
_ACCESSION_RE = re.compile(r"([A-Z]{1,4}_[A-Z0-9]{3,}\.\d+|[A-Z]{1,4}[0-9]{5,}\.\d+)")


def _canonical_accession(s: str) -> Optional[str]:
    """
    Extract a canonical accession from a string.

    Examples it should handle:
      - "NC_003310.1_3323_4" -> "NC_003310.1"
      - "@NC_003310.1_3323_4" -> "NC_003310.1"
      - "NC_003310.1|foo bar" -> "NC_003310.1"
      - "NC_003310.1 some desc" -> "NC_003310.1"
    """
    if not s:
        return None

    s = s.strip()
    if s.startswith("@"):
        s = s[1:]

    m = _ACCESSION_RE.search(s)
    if m:
        return m.group(1)

    # fallback: your original rule (left of first underscore)
    return s.split("_", 1)[0]


def build_ground_truth_metrics_df(
    bam_path: str,
    *,
    skip_secondary: bool = True,
    skip_supplementary: bool = True,
    debug_n: int = 10,
) -> pd.DataFrame:
    """
    Build per-aligned-reference TP/FP counts using BAM-derived ground truth.

    Truth accession is extracted from QNAME (first accession-looking token).
    Aligned accession is extracted from RNAME (reference name in BAM).
      TP if aligned_accession == truth_accession else FP

    Returns DataFrame:
      accession | tp_reads | fp_reads | total_reads | fp_fraction
    """

    counts = defaultdict(lambda: [0, 0])  # aligned_acc -> [tp, fp]
    debug_examples: list[Tuple[str, str]] = []  # (truth, aligned)

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for aln in bam.fetch(until_eof=True):
            if aln.is_unmapped:
                continue
            if skip_secondary and aln.is_secondary:
                continue
            if skip_supplementary and aln.is_supplementary:
                continue

            # aligned ref name from BAM header
            raw_ref = bam.get_reference_name(aln.reference_id)
            aligned_acc = _canonical_accession(raw_ref)

            # truth accession from read id
            truth_acc = _canonical_accession(aln.query_name)

            # If we can't canonicalize, skip (or you could count separately)
            if not aligned_acc or not truth_acc:
                continue

            if aligned_acc == truth_acc:
                counts[aligned_acc][0] += 1
            else:
                counts[aligned_acc][1] += 1
                if len(debug_examples) < debug_n:
                    debug_examples.append((truth_acc, aligned_acc))

    if debug_examples:
        print("[debug] first mismatches (truth -> aligned):")
        for t, a in debug_examples:
            print(f"  {t}  ->  {a}")

    rows = []
    for accession, (tp, fp) in counts.items():
        total = tp + fp
        rows.append({
            "accession": accession,
            "tp_reads": tp,
            "fp_reads": fp,
            "total_reads": total,
            "fp_fraction": (fp / total) if total else 0.0,
        })

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(["fp_reads", "tp_reads"], ascending=[False, False]).reset_index(drop=True)

    return df
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


