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
    accession_to_key=None,
    accession_to_subkey=None,
    accession_to_toplevelkey=None,
    breadth_weight=1/3,
    minhash_weight=1/3,
    gini_weight=1/3,
    disparity_weight=0.0,
    hmp_weight=0.0,
    abundance_confidence_weight=0.0,
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
    sampletype="sterile",
    entropy_lambda=0.0,
    min_weight=0.0,
    threshold_pref_lambda=0.05,
    optimize_tp_target=None,
    optimize_tp_target_weight=0.0,
    optimize_tp_target_scope="taxa",
    optimize_youden_weight=0.0,
    optimize_fp_cutoff=None,
    optimize_fp_cutoff_weight=0.0,
    optimize_curve_scope="taxa",
    tp_score_floor=0.7,
    tp_floor_weight=0.0,
    fp_score_ceiling=0.15,
    fp_ceiling_weight=0.0,
    separation_weight=0.0,
    weight_prior=None,
    weight_prior_lambda=0.0,
    plasmid_bonus_weight=0.05,
    youden_min_threshold=None,
    prefer_granularity="subkey",
    platform=None
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

    gran_tp_fp_counts = {"key": tp_fp_counts}
    if accession_to_key and accession_to_key is not accession_to_taxid:
        gran_tp_fp_counts["key"] = compute_tp_fp_counts_by_taxid(
            input_bam,
            accession_to_taxid=accession_to_key,
            debug_n=10,
        )
    if accession_to_subkey:
        gran_tp_fp_counts["subkey"] = compute_tp_fp_counts_by_taxid(
            input_bam,
            accession_to_taxid=accession_to_subkey,
            debug_n=10,
        )
    if accession_to_toplevelkey:
        gran_tp_fp_counts["toplevelkey"] = compute_tp_fp_counts_by_taxid(
            input_bam,
            accession_to_taxid=accession_to_toplevelkey,
            debug_n=10,
        )

    start_w = dict(
        breadth_weight=float(breadth_weight),
        minhash_weight=float(minhash_weight),
        gini_weight=float(gini_weight),
        disparity_weight=float(disparity_weight),
        hmp_weight=float(hmp_weight),
        plasmid_bonus_weight=float(plasmid_bonus_weight),
        abundance_confidence_weight=float(abundance_confidence_weight),
    )

    # ── Helper: run optimizer at a given granularity ─────────────────────────
    def _run_optimize(df, counts, col_name, label):
        """Run optimize_weights_for_tp_fp on a given metrics_df & counts dict."""
        res = optimize_weights_for_tp_fp(
            metrics_df=df,
            tp_fp_counts=counts,
            accession_col=col_name,
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
            entropy_lambda=float(entropy_lambda),
            min_weight=float(min_weight),
            threshold_pref_lambda=float(threshold_pref_lambda),
            tp_target=optimize_tp_target,
            tp_target_weight=float(optimize_tp_target_weight),
            tp_target_scope=str(optimize_tp_target_scope or "taxa"),
            youden_weight=float(optimize_youden_weight),
            fp_cutoff=optimize_fp_cutoff,
            fp_cutoff_weight=float(optimize_fp_cutoff_weight),
            curve_scope=str(optimize_curve_scope or "taxa"),
            tp_score_floor=float(tp_score_floor),
            tp_floor_weight=float(tp_floor_weight),
            fp_score_ceiling=float(fp_score_ceiling),
            fp_ceiling_weight=float(fp_ceiling_weight),
            separation_weight=float(separation_weight),
            weight_prior=weight_prior,
            weight_prior_lambda=float(weight_prior_lambda),
            plasmid_bonus_weight=float(plasmid_bonus_weight),
            abundance_confidence_weight=float(abundance_confidence_weight),
            youden_min_threshold=youden_min_threshold,
        )
        print(f"[optimize:{label}] status: {res['status']}")
        if res["status"] == "ok":
            print(f"[optimize:{label}] TP reads: {res.get('tp_total_reads')}, "
                  f"FP reads: {res.get('fp_total_reads')}")
            print(f"[optimize:{label}] TP taxa: {res.get('tp_total_taxa')}, "
                  f"FP taxa: {res.get('fp_total_taxa')}")
            print(f"[optimize:{label}] weights: {res['weights']}")
        return res

    # ── Primary optimization at taxid level (existing behaviour) ─────────────
    opt = _run_optimize(metrics_df, tp_fp_counts, accession_col, "taxid")
    ## general stats on total taxa expected
    tp_total_taxa = metrics_df[metrics_df['tp_reads']>0]['key'].unique().shape[0]
    fp_total_taxa = metrics_df[(metrics_df['tp_reads']==0) & (metrics_df['fp_reads']>0)]['key'].unique().shape[0]
    total_taxa = metrics_df["key"].unique().shape[0]


    # ── Multi-granularity optimization (key, subkey, toplevelkey) ─────────────
    granularity_results = {}
    granularity_reports = {}
    for gran_field in ("key", "subkey", "toplevelkey"):
        if gran_field not in metrics_df.columns:
            continue
        # Build a grouped df: aggregate features by gran_field, merge TP/FP
        gran_df = metrics_df.copy()
        gran_df[gran_field] = gran_df[gran_field].astype(str)

        # Aggregate: take max of features (consistent with build_metrics_df logic)
        agg_spec = {
            "name": "first",
            "category": "first",
            "key": "first",
            "subkey": "first",
            "toplevelkey": "first",
            "breadth_log_score": "max",
            "minhash_reduction": "max",
            "disparity_score": "max",
            "gini_coefficient": "max",
            "hmp_percentile": "max",
            "plasmid_score": "max",
            "abundance_confidence": "max",
            "has_plasmid": "max",
            "tp_reads": "sum",
            "fp_reads": "sum",
            "total_reads": "sum",
        }
        # Only include columns that exist
        agg_spec = {k: v for k, v in agg_spec.items() if k in gran_df.columns}
        grouped = gran_df.groupby(gran_field, as_index=False).agg(agg_spec)
        grouped["fp_fraction"] = grouped.apply(
            lambda r: (r["fp_reads"] / r["total_reads"]) if r["total_reads"] else 0.0, axis=1
        )

        # Build tp_fp_counts keyed by gran_field value
        gran_counts = {}
        if gran_field in gran_tp_fp_counts:
            src_counts = gran_tp_fp_counts[gran_field]
            for gval in grouped[gran_field].astype(str).tolist():
                gran_counts[gval] = src_counts.get(gval, (0, 0))
            grouped["tp_reads"] = grouped[gran_field].map(lambda v: gran_counts.get(str(v), (0, 0))[0]).astype(int)
            grouped["fp_reads"] = grouped[gran_field].map(lambda v: gran_counts.get(str(v), (0, 0))[1]).astype(int)
            grouped["total_reads"] = grouped["tp_reads"] + grouped["fp_reads"]
            grouped["fp_fraction"] = grouped.apply(
                lambda r: (r["fp_reads"] / r["total_reads"]) if r["total_reads"] else 0.0, axis=1
            )
        else:
            for _, row in grouped.iterrows():
                gval = str(row[gran_field])
                gran_counts[gval] = (int(row["tp_reads"]), int(row["fp_reads"]))

        gran_opt = _run_optimize(grouped, gran_counts, gran_field, gran_field)
        granularity_results[gran_field] = gran_opt

        if gran_opt.get("status") == "ok":
            g_weights = gran_opt.get("weights") or {}
            g_scores = compute_tass_score_from_metrics(
                grouped,
                breadth_w=float(g_weights.get("breadth_weight", breadth_weight)),
                minhash_w=float(g_weights.get("minhash_weight", minhash_weight)),
                gini_w=float(g_weights.get("gini_weight", gini_weight)),
                disparity_w=float(g_weights.get("disparity_weight", disparity_weight)),
                hmp_w=float(g_weights.get("hmp_weight", hmp_weight)),
                alpha=float(alpha),
                plasmid_bonus_w=float(g_weights.get("plasmid_bonus_weight", plasmid_bonus_weight)),
                abundance_confidence_w=float(g_weights.get("abundance_confidence_weight", abundance_confidence_weight)),
            )
            g_scores = np.asarray(g_scores, dtype=float)

            gran_rows = []
            for i in range(len(grouped)):
                gran_rows.append({
                    "group": str(grouped.iloc[i][gran_field]),
                    "name": str(grouped.iloc[i].get("name", "")),
                    "key": str(grouped.iloc[i].get("key", "")),
                    "subkey": str(grouped.iloc[i].get("subkey", "")),
                    "toplevelkey": str(grouped.iloc[i].get("toplevelkey", "")),
                    "microbial_category": str(grouped.iloc[i].get("category", "")),
                    "tp_reads": int(grouped.iloc[i]["tp_reads"]),
                    "fp_reads": int(grouped.iloc[i]["fp_reads"]),
                    "total_reads": int(grouped.iloc[i]["total_reads"]),
                    "fp_fraction": float(grouped.iloc[i]["fp_reads"] / grouped.iloc[i]["total_reads"]) if grouped.iloc[i]["total_reads"] else 0.0,
                    "tass_score": float(g_scores[i]),
                    "features": {
                        "breadth_log_score": float(grouped.iloc[i].get("breadth_log_score", 0.0)),
                        "minhash_reduction": float(grouped.iloc[i].get("minhash_reduction", 0.0)),
                        "disparity_score": float(grouped.iloc[i].get("disparity_score", 0.0)),
                        "hmp_percentile": float(grouped.iloc[i].get("hmp_percentile", 0.0)),
                        "gini_coefficient": float(grouped.iloc[i].get("gini_coefficient", 0.0)),
                        "plasmid_score": float(grouped.iloc[i].get("plasmid_score", 0.0)),
                        "has_plasmid": bool(grouped.iloc[i].get("has_plasmid", False)),
                    },
                })
            gran_rows.sort(key=lambda r: (r["fp_reads"], -r["tp_reads"], r["tass_score"]), reverse=True)
            granularity_reports[gran_field] = gran_rows

    # ── Select best overall result ───────────────────────────────────────────
    # Prefer the user-specified granularity if it succeeded; otherwise fall
    # back to the level with the lowest loss across all granularities.
    all_opts = {"taxid": opt}
    all_opts.update(granularity_results)

    # Try the preferred granularity first
    _pref = prefer_granularity or "subkey"
    if _pref in all_opts and all_opts[_pref].get("status") == "ok":
        best_label = _pref
        best_opt = all_opts[_pref]
        print(f"[optimize] Using preferred granularity: {_pref} "
              f"(loss={best_opt.get('loss', 'N/A')})")
    else:
        # Fallback: pick the level with the lowest loss
        if _pref in all_opts:
            print(f"[optimize] Preferred granularity '{_pref}' did not succeed "
                  f"(status={all_opts[_pref].get('status', '?')}), falling back to best loss")
        best_label = "taxid"
        best_opt = opt
        for label, result in all_opts.items():
            if result.get("status") != "ok":
                continue
            if best_opt.get("status") != "ok":
                best_label = label
                best_opt = result
                continue
            # Lower loss = better
            if result.get("loss", float("inf")) < best_opt.get("loss", float("inf")):
                best_label = label
                best_opt = result

    print(f"\n[optimize] Best granularity: {best_label} (loss={best_opt.get('loss', 'N/A')})")

    if best_opt.get("status") != "ok":
        report_obj = {
            "status": best_opt.get("status", "failed"),
            "config": {
                "alpha": float(alpha),
                "total_fp_taxa": fp_total_taxa,
                "total_tp_taxa": tp_total_taxa,
                "total_taxa": total_taxa,
                "platform": platform,
                "sampletype": sampletype,
                "pos_weight": float(optimize_pos_weight),
                "neg_weight": float(optimize_neg_weight),
                "reg_lambda": float(optimize_reg),
                "entropy_lambda": float(entropy_lambda),
                "min_weight": float(min_weight),
                "threshold_pref_lambda": float(threshold_pref_lambda),
                "tp_target": None if optimize_tp_target is None else float(optimize_tp_target),
                "tp_target_weight": float(optimize_tp_target_weight),
                "tp_target_scope": str(optimize_tp_target_scope or "taxa"),
                "youden_weight": float(optimize_youden_weight),
                "fp_cutoff": None if optimize_fp_cutoff is None else float(optimize_fp_cutoff),
                "fp_cutoff_weight": float(optimize_fp_cutoff_weight),
                "curve_scope": str(optimize_curve_scope or "taxa"),
                "seed": int(optimize_seed),
                "maxiter_global": int(optimize_maxiter),
                "maxiter_local": int(optimize_local_maxiter),
                "optimize_mode": optimize_mode,
                "hybrid_lambda": float(hybrid_lambda),
                "start_weights": start_w,
            },
            "optimizer": best_opt,
            "best_weights": start_w,
            "best_granularity": best_label,
            "granularity_results": {
                k: {"status": v.get("status"), "weights": v.get("weights"), "loss": v.get("loss")}
                for k, v in all_opts.items()
            },
            "granularity_reports": granularity_reports,
            "report": [],
        }
        if optimize_report:
            with open(optimize_report, "w") as f:
                json.dump(report_obj, f, indent=2)
        return report_obj

    best_w = best_opt["weights"]
    breadth_weight = best_w["breadth_weight"]
    minhash_weight = best_w["minhash_weight"]
    gini_weight = best_w["gini_weight"]
    disparity_weight = best_w.get("disparity_weight", 0.0)
    hmp_weight = best_w.get("hmp_weight", 0.0)
    _best_pbw = float(best_w.get("plasmid_bonus_weight", plasmid_bonus_weight))

    scores = compute_tass_score_from_metrics(
        metrics_df,
        breadth_w=float(breadth_weight),
        minhash_w=float(minhash_weight),
        gini_w=float(gini_weight),
        disparity_w=float(disparity_weight),
        hmp_w=float(hmp_weight),
        alpha=float(alpha),
        plasmid_bonus_w=_best_pbw,
        abundance_confidence_w=float(abundance_confidence_weight),
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
            "microbial_category": str(metrics_df.iloc[i].get("category", "")),
            "key": metrics_df.iloc[i].get("key"),
            "subkey": metrics_df.iloc[i].get("subkey"),
            "toplevelkey": metrics_df.iloc[i].get("toplevelkey"),
            "tp_reads": int(tp[i]),
            "fp_reads": int(fp[i]),
            "total_reads": int(total[i]),
            "fp_fraction": float(fp[i] / total[i]) if total[i] else 0.0,
            "tass_score": float(scores[i]),
            "features": {
                "breadth_log_score": float(metrics_df.iloc[i].get("breadth_log_score", 0.0)),
                "minhash_reduction": float(metrics_df.iloc[i].get("minhash_reduction", 0.0)),
                "disparity_score": float(metrics_df.iloc[i].get("disparity_score", 0.0)),
                "hmp_percentile": float(metrics_df.iloc[i].get("hmp_percentile", 0.0)),
                "gini_coefficient": float(metrics_df.iloc[i].get("gini_coefficient", 0.0)),
                "plasmid_score": float(metrics_df.iloc[i].get("plasmid_score", 0.0)),
                "abundance_confidence": float(metrics_df.iloc[i].get("abundance_confidence", 0.0)),
                "has_plasmid": bool(metrics_df.iloc[i].get("has_plasmid", False)),
            }
        })

    report_rows.sort(key=lambda r: (r["fp_reads"], -r["tp_reads"], r["tass_score"]), reverse=True)

    report_obj = {
        "status": "ok",
        "config": {
            "alpha": float(alpha),
            "sampletype": sampletype,
            "pos_weight": float(optimize_pos_weight),
            "neg_weight": float(optimize_neg_weight),
            "reg_lambda": float(optimize_reg),
            "entropy_lambda": float(entropy_lambda),
            "min_weight": float(min_weight),
            "threshold_pref_lambda": float(threshold_pref_lambda),
            "tp_target": None if optimize_tp_target is None else float(optimize_tp_target),
            "tp_target_weight": float(optimize_tp_target_weight),
            "tp_target_scope": str(optimize_tp_target_scope or "taxa"),
            "youden_weight": float(optimize_youden_weight),
            "fp_cutoff": None if optimize_fp_cutoff is None else float(optimize_fp_cutoff),
            "fp_cutoff_weight": float(optimize_fp_cutoff_weight),
            "curve_scope": str(optimize_curve_scope or "taxa"),
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
            "disparity_weight": float(disparity_weight),
            "hmp_weight": float(hmp_weight),
            "plasmid_bonus_weight": _best_pbw,
            "abundance_confidence_weight": float(abundance_confidence_weight),
        },
        "best_granularity": best_label,
        "granularity_results": {
            k: {"status": v.get("status"), "weights": v.get("weights"), "loss": v.get("loss")}
            for k, v in all_opts.items()
        },
        "granularity_reports": granularity_reports,
        "optimizer": best_opt,
        "report": report_rows,
    }

    if optimize_report:
        with open(optimize_report, "w") as f:
            json.dump(report_obj, f, indent=2)

    return report_obj



def _clip01(x: np.ndarray, eps: float = 1e-6) -> np.ndarray:
    return np.clip(x, eps, 1.0 - eps)

def _entropy(w, eps=1e-8):
    """Shannon entropy of weight vector (higher = more balanced)."""
    w_safe = np.clip(w, eps, 1.0)
    return -float(np.sum(w_safe * np.log(w_safe)))

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
    entropy_lambda: float = 0.0,
    min_weight: float = 0.0,
    threshold_pref_lambda: float = 0.0,
    tp_target: float | None = None,
    tp_target_weight: float = 0.0,
    tp_target_scope: str = "taxa",
    youden_weight: float = 0.0,
    fp_cutoff: float | None = None,
    fp_cutoff_weight: float = 0.0,
    curve_scope: str = "taxa",
    tp_score_floor: float = 0.7,
    tp_floor_weight: float = 0.0,
    fp_score_ceiling: float = 0.15,
    fp_ceiling_weight: float = 0.0,
    separation_weight: float = 0.0,
    weight_prior: dict | None = None,
    weight_prior_lambda: float = 0.0,
    plasmid_bonus_weight: float = 0.0,
    abundance_confidence_weight: float = 0.0,
    youden_min_threshold: float | None = None,
):
    from scipy.optimize import differential_evolution, minimize
    import numpy as np

    # Store the abundance_confidence_weight for use in score computations.
    # It's additive (like plasmid_bonus_weight) — not on the simplex.
    _acw = float(abundance_confidence_weight)

    # Minimum allowed Youden J threshold — prevents optimizer from picking
    # unreasonably low cutoffs for sterile/blood sites.
    _youden_min_t = float(youden_min_threshold) if youden_min_threshold is not None else None

    accessions = metrics_df[accession_col].astype(str).tolist()
    tp_raw = np.array([tp_fp_counts.get(a, (0, 0))[0] for a in accessions], dtype=float)
    fp_raw = np.array([tp_fp_counts.get(a, (0, 0))[1] for a in accessions], dtype=float)

    # Filter out accessions with zero TP *and* zero FP reads –
    # they carry no ground-truth signal and skew mean scores / taxa counts.
    has_signal = (tp_raw + fp_raw) >= 1
    if has_signal.sum() == 0:
        return {"weights": start_weights, "status": "skipped (no accessions with TP+FP >= 1)"}

    # Keep only rows that have at least 1 read of ground-truth signal
    keep_idx = np.where(has_signal)[0]
    tp = tp_raw[keep_idx]
    fp = fp_raw[keep_idx]
    metrics_df = metrics_df.iloc[keep_idx].reset_index(drop=True)
    accessions = [accessions[i] for i in keep_idx]

    # Taxid-level labels: TP if any TP reads exist, else negative.
    y_pos = (tp > 0).astype(float)
    y_neg = (tp == 0).astype(float)

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
        start_weights.get("disparity_weight", 0.0),
        start_weights.get("hmp_weight", 0.0),
    ], dtype=float)

    # Starting plasmid bonus weight (optimized separately, not on simplex)
    _pbw0 = float(start_weights.get("plasmid_bonus_weight", plasmid_bonus_weight))
    _pbw_max = 0.20  # upper bound for plasmid bonus

    # Apply min-weight floor to w0 and bounds
    _mw = max(0.0, min(0.15, float(min_weight)))  # cap at 15% per weight
    n_weights = 5

    w0 = np.clip(w0, _mw, 1.0)
    w0 = (w0 / w0.sum()) if w0.sum() > 0 else np.full(n_weights, 1.0 / n_weights)

    # ── Weight prior: quadratic pull toward preferred weight targets ──────
    # weight_prior is a dict like {"breadth_weight": 0.4} meaning "I want
    # breadth to be near 0.4".  weight_prior_lambda controls the strength.
    _weight_names = ["breadth_weight", "minhash_weight", "gini_weight",
                     "disparity_weight", "hmp_weight"]
    _prior_target = None
    _prior_mask = None
    if weight_prior and weight_prior_lambda > 0:
        _prior_target = np.zeros(n_weights, dtype=float)
        _prior_mask = np.zeros(n_weights, dtype=float)
        for wname, wtarget in weight_prior.items():
            if wname in _weight_names:
                idx = _weight_names.index(wname)
                _prior_target[idx] = float(wtarget)
                _prior_mask[idx] = 1.0

    # Max entropy for reference (uniform distribution)
    _max_entropy = _entropy(np.full(n_weights, 1.0 / n_weights))

    # Check if any organisms actually have plasmid data — skip pbw optimization if not
    _has_plasmid_data = ("plasmid_score" in metrics_df.columns
                         and (metrics_df["plasmid_score"] > 0).any())

    def _threshold_pref_penalty(scores: np.ndarray) -> float:
        if threshold_pref_lambda <= 0:
            return 0.0

        # Prefer taxid-level retention when available; fall back to reads.
        if tp_total_taxa > 0 and fp_total_taxa > 0:
            tp_w = y_pos
            fp_w = y_neg
            total_tp = tp_total_taxa
            total_fp = fp_total_taxa
        elif tp_total_reads > 0 and fp_total_reads > 0:
            tp_w = tp
            fp_w = fp
            total_tp = tp_total_reads
            total_fp = fp_total_reads
        else:
            return 0.0

        order = np.argsort(scores)[::-1]
        scores_sorted = scores[order]
        tp_sorted = tp_w[order]
        fp_sorted = fp_w[order]

        tp_ret = np.cumsum(tp_sorted) / (total_tp + 1e-12)
        fp_ret = np.cumsum(fp_sorted) / (total_fp + 1e-12)
        J = tp_ret - fp_ret

        # Apply minimum threshold floor if configured
        if _youden_min_t is not None:
            floor_mask = scores_sorted >= _youden_min_t
            if floor_mask.any():
                best_idx = int(np.where(floor_mask, J, -np.inf).argmax())
            else:
                best_idx = int(np.nanargmax(J))
        else:
            best_idx = int(np.nanargmax(J))
        t_min = float(scores_sorted.min())
        t_max = float(scores_sorted.max())
        norm_t = (float(scores_sorted[best_idx]) - t_min) / (t_max - t_min + 1e-12)

        # Penalize lower best-thresholds so the optimizer favors higher cutoffs.
        return float(threshold_pref_lambda) * (1.0 - norm_t)

    def _retention_curves(scores: np.ndarray, scope: str):
        scope = (scope or "taxa").lower().strip()

        def _curves(tp_w, fp_w, total_tp, total_fp):
            if total_tp <= 0 or total_fp <= 0:
                return None, None, None
            order = np.argsort(scores)[::-1]
            s_sorted = scores[order]
            tp_sorted = tp_w[order]
            fp_sorted = fp_w[order]
            tp_ret = np.cumsum(tp_sorted) / (total_tp + 1e-12)
            fp_ret = np.cumsum(fp_sorted) / (total_fp + 1e-12)
            return tp_ret, fp_ret, s_sorted

        if scope == "reads":
            return _curves(tp, fp, tp_total_reads, fp_total_reads)
        if scope == "hybrid":
            tp_ret_taxa, fp_ret_taxa, s_sorted = _curves(y_pos, y_neg, tp_total_taxa, fp_total_taxa)
            tp_ret_reads, fp_ret_reads, _ = _curves(tp, fp, tp_total_reads, fp_total_reads)
            if tp_ret_taxa is None or tp_ret_reads is None:
                return None, None, None
            lam = max(0.0, min(1.0, float(hybrid_lambda)))
            tp_ret = lam * tp_ret_taxa + (1.0 - lam) * tp_ret_reads
            fp_ret = lam * fp_ret_taxa + (1.0 - lam) * fp_ret_reads
            return tp_ret, fp_ret, s_sorted

        return _curves(y_pos, y_neg, tp_total_taxa, fp_total_taxa)

    def _retention_at_cutoff(scores: np.ndarray, scope: str, cutoff: float):
        scope = (scope or "taxa").lower().strip()

        def _at(tp_w, fp_w, total_tp, total_fp):
            if total_tp <= 0 or total_fp <= 0:
                return None, None
            mask = scores >= cutoff
            tp_ret = float(tp_w[mask].sum() / (total_tp + 1e-12))
            fp_ret = float(fp_w[mask].sum() / (total_fp + 1e-12))
            return tp_ret, fp_ret

        if scope == "reads":
            return _at(tp, fp, tp_total_reads, fp_total_reads)
        if scope == "hybrid":
            tp_tax, fp_tax = _at(y_pos, y_neg, tp_total_taxa, fp_total_taxa)
            tp_rd, fp_rd = _at(tp, fp, tp_total_reads, fp_total_reads)
            if tp_tax is None or tp_rd is None:
                return None, None
            lam = max(0.0, min(1.0, float(hybrid_lambda)))
            return (lam * tp_tax + (1.0 - lam) * tp_rd, lam * fp_tax + (1.0 - lam) * fp_rd)

        return _at(y_pos, y_neg, tp_total_taxa, fp_total_taxa)

    # Mutable container for the current plasmid bonus weight during optimization.
    # Updated by the outer DE/SLSQP loops so loss_for_w always sees the latest value.
    _current_pbw = [_pbw0]

    def loss_for_w(w: np.ndarray) -> float:
        # simplex + min-weight constraints (soft)
        if np.any(w < 0) or np.any(w > 1) or abs(w.sum() - 1.0) > 1e-6:
            return 1e6 + 1e6 * (w.sum() - 1.0) ** 2
        if _mw > 0 and np.any(w < _mw - 1e-8):
            return 1e6 + 1e6 * float(np.sum(np.maximum(0, _mw - w)))

        scores = compute_tass_score_from_metrics(
            metrics_df,
            breadth_w=float(w[0]),
            minhash_w=float(w[1]),
            gini_w=float(w[2]),
            disparity_w=float(w[3]),
            hmp_w=float(w[4]),
            alpha=alpha,
            plasmid_bonus_w=float(_current_pbw[0]),
            abundance_confidence_w=_acw,
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

        # Entropy regularization: penalize weight concentration
        # Low entropy = one weight dominates → add penalty
        # We negate entropy (since we minimize loss) and scale by _max_entropy
        ent_penalty = 0.0
        if entropy_lambda > 0:
            ent = _entropy(w)
            # Penalty = entropy_lambda * (1 - entropy / max_entropy)
            # = 0 when weights are uniform, = entropy_lambda when one weight = 1
            ent_penalty = entropy_lambda * (1.0 - ent / _max_entropy) if _max_entropy > 0 else 0.0

        # Weight prior penalty: pull specified weights toward target values
        prior_penalty = 0.0
        if _prior_target is not None and _prior_mask is not None:
            diffs = (w - _prior_target) * _prior_mask
            prior_penalty = float(weight_prior_lambda) * float(np.sum(diffs ** 2))

        thresh_penalty = _threshold_pref_penalty(p)

        tp_target_penalty = 0.0
        if tp_target is not None and tp_target_weight > 0:
            scope = (tp_target_scope or "taxa").lower().strip()
            tp_mean_taxa = (float(np.sum(y_pos * p) / tp_total_taxa) if tp_total_taxa > 0 else None)
            tp_mean_reads = (float(np.sum(tp * p) / tp_total_reads) if tp_total_reads > 0 else None)

            if scope == "reads":
                tp_mean = tp_mean_reads
            elif scope == "hybrid":
                lam = max(0.0, min(1.0, float(hybrid_lambda)))
                if tp_mean_taxa is None and tp_mean_reads is None:
                    tp_mean = None
                elif tp_mean_taxa is None:
                    tp_mean = tp_mean_reads
                elif tp_mean_reads is None:
                    tp_mean = tp_mean_taxa
                else:
                    tp_mean = lam * tp_mean_taxa + (1.0 - lam) * tp_mean_reads
            else:
                tp_mean = tp_mean_taxa

            if tp_mean is not None:
                gap = max(0.0, float(tp_target) - tp_mean)
                tp_target_penalty = float(tp_target_weight) * (gap ** 2)

        youden_penalty = 0.0
        if youden_weight > 0:
            tp_ret, fp_ret, s_sorted = _retention_curves(p, curve_scope)
            if tp_ret is not None and s_sorted is not None:
                j_vals = tp_ret - fp_ret
                # Apply minimum threshold floor: mask out Youden J values
                # at thresholds below the floor (if configured).  This
                # prevents the optimizer from rewarding very low cutoffs
                # that are impractical for sterile/blood sites.
                if _youden_min_t is not None:
                    floor_mask = s_sorted >= _youden_min_t
                    if floor_mask.any():
                        best_j = float(np.nanmax(j_vals[floor_mask]))
                    else:
                        # All thresholds below floor → penalize heavily
                        best_j = 0.0
                else:
                    best_j = float(np.nanmax(j_vals))
                youden_penalty = float(youden_weight) * (1.0 - best_j)

        fp_cutoff_penalty = 0.0
        if fp_cutoff is not None and fp_cutoff_weight > 0:
            _, fp_ret = _retention_at_cutoff(p, curve_scope, float(fp_cutoff))
            if fp_ret is not None:
                fp_cutoff_penalty = float(fp_cutoff_weight) * float(fp_ret)

        # ── TP Score Floor Hinge Loss ──────────────────────────────────
        # Per-organism quadratic penalty for any TP scoring below tp_score_floor.
        # Pushes ALL individual TP organisms above the floor (e.g. 0.7).
        tp_floor_penalty = 0.0
        if tp_floor_weight > 0 and tp_score_floor is not None:
            tp_mask = y_pos > 0
            if tp_mask.any():
                tp_scores = p[tp_mask]
                shortfalls = np.maximum(0.0, tp_score_floor - tp_scores)
                tp_floor_penalty = float(tp_floor_weight) * float(np.mean(shortfalls ** 2))

        # ── FP Score Ceiling Hinge Loss ────────────────────────────────
        # Per-organism quadratic penalty for any FP scoring above fp_score_ceiling.
        # Pushes ALL individual FP organisms below the ceiling (e.g. 0.15).
        fp_ceiling_penalty = 0.0
        if fp_ceiling_weight > 0 and fp_score_ceiling is not None:
            fp_mask = y_neg > 0
            if fp_mask.any():
                fp_scores = p[fp_mask]
                overshoots = np.maximum(0.0, fp_scores - fp_score_ceiling)
                fp_ceiling_penalty = float(fp_ceiling_weight) * float(np.mean(overshoots ** 2))

        # ── Multi-Threshold Retention Shape Penalty ────────────────────
        # Evaluates the retention curve at multiple thresholds and penalizes
        # deviations from the ideal step-function shape:
        #   - TP retention should stay ≥ 0.95 at all thresholds up to 0.7
        #   - FP retention should drop to ≤ 0.05 at thresholds ≥ 0.15
        separation_penalty = 0.0
        if separation_weight > 0:
            _check_thresholds = [0.1, 0.2, 0.3, 0.5, 0.7]
            shape_loss = 0.0
            for _t in _check_thresholds:
                _tp_ret_t, _fp_ret_t = _retention_at_cutoff(p, curve_scope, float(_t))
                if _tp_ret_t is not None and _fp_ret_t is not None:
                    # TP should be retained at ~1.0 across all thresholds
                    tp_drop = max(0.0, 0.95 - _tp_ret_t)
                    shape_loss += tp_drop ** 2
                    # FP should be ~0.0 at thresholds >= 0.15
                    if _t >= 0.15:
                        fp_excess = max(0.0, _fp_ret_t - 0.05)
                        shape_loss += fp_excess ** 2
            separation_penalty = float(separation_weight) * shape_loss

        return (
            pos_w * pos_term
            + neg_w * neg_term
            + reg
            + ent_penalty
            + thresh_penalty
            + tp_target_penalty
            + youden_penalty
            + fp_cutoff_penalty
            + tp_floor_penalty
            + fp_ceiling_penalty
            + separation_penalty
            + prior_penalty
        )

    # --- Global stage (DE) over 4 simplex vars + optional pbw ---
    # x = [w_breadth, w_minhash, w_gini, w_disparity, (pbw)], w_hmp = 1 - sum(x[:4])
    _upper = 1.0 - (n_weights - 1) * _mw  # max any single weight can be
    _optimize_pbw = _has_plasmid_data  # only optimize pbw if plasmid data exists
    if _optimize_pbw:
        print(f"[optimize] Plasmid data detected — optimizing plasmid_bonus_weight (start={_pbw0:.4f}, max={_pbw_max})")
    else:
        print(f"[optimize] No plasmid data in metrics — plasmid_bonus_weight fixed at {_pbw0:.4f}")

    def unpack_x(x: np.ndarray):
        """Unpack DE vector into (5-weight array, pbw_value)."""
        wb, wm, wg, wd = float(x[0]), float(x[1]), float(x[2]), float(x[3])
        wh = 1.0 - (wb + wm + wg + wd)
        w = np.array([wb, wm, wg, wd, wh], dtype=float)
        if _mw > 0:
            w = np.clip(w, _mw, _upper)
            w = w / w.sum()  # re-normalize after clipping
        pbw_val = float(x[4]) if (_optimize_pbw and len(x) > 4) else _pbw0
        return w, pbw_val

    def loss_for_x(x: np.ndarray) -> float:
        w, pbw_val = unpack_x(x)
        if np.any(w < 0.0) or np.any(w > 1.0):
            return 1e6 + 1e6 * max(0.0, -w.min())
        _current_pbw[0] = pbw_val
        return loss_for_w(w)

    bounds = [(_mw, _upper), (_mw, _upper), (_mw, _upper), (_mw, _upper)]
    if _optimize_pbw:
        bounds.append((0.0, _pbw_max))
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

    w_de, pbw_de = unpack_x(de_res.x)
    if np.any(w_de < 0) or abs(w_de.sum() - 1.0) > 1e-6:
        w_de = w0.copy()
        pbw_de = _pbw0
    else:
        w_de = np.clip(w_de, _mw, _upper)
        w_de = w_de / w_de.sum()

    _current_pbw[0] = pbw_de  # carry forward DE result

    # --- Local refinement (SLSQP) on 5 simplex vars + optional pbw ---
    if _optimize_pbw:
        # 6-variable vector: w[0:5] = simplex weights, w[5] = pbw
        def loss_for_slsqp(x6):
            w5 = x6[:5]
            _current_pbw[0] = float(x6[5])
            return loss_for_w(w5)

        cons_slsqp = [{"type": "eq", "fun": lambda x6: np.sum(x6[:5]) - 1.0}]
        bnds_slsqp = [(_mw, _upper)] * n_weights + [(0.0, _pbw_max)]
        x0_slsqp = np.append(w_de, pbw_de)

        local_res = minimize(
            loss_for_slsqp,
            x0=x0_slsqp,
            method="SLSQP",
            bounds=bnds_slsqp,
            constraints=cons_slsqp,
            options={"maxiter": maxiter_local, "ftol": 1e-9},
        )

        w_best = np.clip(local_res.x[:5], _mw, _upper)
        w_best = w_best / w_best.sum()
        pbw_best = float(np.clip(local_res.x[5], 0.0, _pbw_max))
    else:
        cons = [{"type": "eq", "fun": lambda w: np.sum(w) - 1.0}]
        bnds = [(_mw, _upper)] * n_weights

        local_res = minimize(
            loss_for_w,
            x0=w_de,
            method="SLSQP",
            bounds=bnds,
            constraints=cons,
            options={"maxiter": maxiter_local, "ftol": 1e-9},
        )

        w_best = np.clip(local_res.x, _mw, _upper)
        w_best = w_best / w_best.sum()
        pbw_best = _pbw0  # unchanged

    _current_pbw[0] = pbw_best
    if _optimize_pbw:
        print(f"[optimize] Optimized plasmid_bonus_weight: {_pbw0:.4f} -> {pbw_best:.4f}")

    # summaries (both read and taxid views)
    scores_best = compute_tass_score_from_metrics(
        metrics_df,
        breadth_w=float(w_best[0]),
        minhash_w=float(w_best[1]),
        gini_w=float(w_best[2]),
        disparity_w=float(w_best[3]),
        hmp_w=float(w_best[4]),
        alpha=alpha,
        plasmid_bonus_w=float(pbw_best),
        abundance_confidence_w=_acw,
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
            "disparity_weight": float(w_best[3]),
            "hmp_weight": float(w_best[4]),
            "plasmid_bonus_weight": float(pbw_best),
            "abundance_confidence_weight": _acw,
        },
        "status": "ok",
        "optimize_mode": optimize_mode,
        "hybrid_lambda": float(hybrid_lambda),
        "threshold_pref_lambda": float(threshold_pref_lambda),
        "tp_target": None if tp_target is None else float(tp_target),
        "tp_target_weight": float(tp_target_weight),
        "tp_target_scope": str(tp_target_scope or "taxa"),
        "youden_weight": float(youden_weight),
        "fp_cutoff": None if fp_cutoff is None else float(fp_cutoff),
        "fp_cutoff_weight": float(fp_cutoff_weight),
        "curve_scope": str(curve_scope or "taxa"),
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


_ACCESSION_RE = re.compile(r"([A-Z]{1,4}_[A-Z0-9]{3,}|[A-Z]{1,4}[0-9]{5,})\.\d+")
_ACCESSION_RE2 = re.compile(r"([A-Z]{1,4}-[A-Z0-9]{3,}|[A-Z]{1,4}[0-9]{5,})_\d+")

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

    m2 = _ACCESSION_RE2.search(s)
    if m2:
        return m2.group(1).replace("-", "_")

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
    final_json: dict,
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
        categorical = stats.get("microbial_category", "unknown")
        # Prefer k2_disparity_score; fall back to disparity if present
        disparity = float(stats.get("k2_disparity_score", stats.get("disparity", 0.0)))

        gini = float(stats.get("gini_coefficient", 0.0))
        hmp = float(stats.get("hmp_percentile", 0.0))
        plasmid_sc = float(stats.get("plasmid_score", 0.0))
        abundance_conf = float(stats.get("abundance_confidence", 0.0))

        tp, fp = tp_fp_by_taxid.get(taxid, (0, 0))
        total = tp + fp

        rows.append({
            "taxid": taxid,
            "name": stats.get("name", ""),
            "key": stats.get("key"),
            "subkey": stats.get("subkey"),
            "toplevelkey": stats.get("toplevelkey"),
            "accessions": accs,
            "category": categorical,
            "breadth_log_score": breadth,
            "minhash_reduction": minhash,
            "disparity_score": disparity,
            "gini_coefficient": gini,
            "hmp_percentile": hmp,
            "plasmid_score": plasmid_sc,
            "abundance_confidence": abundance_conf,
            "has_plasmid": bool(stats.get("has_plasmid", False)),
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
        "key": "first",
        "subkey": "first",
        "toplevelkey": "first",
        "accessions": union_lists,
        "category": "first",
        "breadth_log_score": "max",
        "minhash_reduction": "max",
        "disparity_score": "max",
        "gini_coefficient": "max",
        "hmp_percentile": "max",
        "plasmid_score": "max",
        "abundance_confidence": "max",
        "has_plasmid": "max",
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


