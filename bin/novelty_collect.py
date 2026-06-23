#!/usr/bin/env python3
# ##############################################################################################
# # Copyright 2022-25 The Johns Hopkins University Applied Physics Laboratory LLC
# # (header trimmed for the sketch -- copy the full APL header from bin/make_report.py)
# ##############################################################################################
#
# Collect the per-sample NOVELTY_SCORE outputs (*.novelty.summary.tsv + *.novelty.candidates.tsv)
# into downloadable artifacts that the interactive HTML report links to:
#
#   <sample>.novelty.json   per-sample {summary, candidates}
#   <sample>.novelty.xlsx   per-sample workbook (Summary + Candidates sheets)
#   all.novelty.json        combined { "taxtriage_novelty": true, "samples": {<sample>: {...}} }
#   all.novelty.xlsx        combined workbook (all samples on the Summary + Candidates sheets)
#
# The combined all.novelty.json is ALSO the feed make_report.py reads (-n/--novelty) to bake the
# Novelty panel + download links into all.odr.html.
#
# Samples are matched between the summary and candidates files by the 'sample' column (robust to
# how the files are named/staged); the filename stem is used only as a last-resort fallback.

import argparse
import csv
import json
import os
import sys

SUMMARY_COLS = ["sample", "classifier", "gene_mode", "total_reads", "dark_fraction", "highrank_only_fraction",
                "lowident_tail_mass", "z_dark", "z_highrank", "z_lowident",
                "novelty_score", "novelty_flag",
                "ref_aligned", "k2_classified", "mmseqs_assigned",
                "mmseqs_assigned_species", "mmseqs_assigned_highrank",
                "ref_aligned_frac", "k2_frac", "mmseqs_frac"]
CAND_COLS = ["sample", "taxid", "rank", "name", "reads", "frac_of_sample",
             "frac_of_highrank"]

# Columns that should carry through as numbers (the rest stay strings).
_NUM = {"total_reads", "dark_fraction", "highrank_only_fraction", "lowident_tail_mass",
        "z_dark", "z_highrank", "z_lowident", "novelty_score", "reads", "frac_of_sample",
        "frac_of_highrank",
        "ref_aligned", "k2_classified", "mmseqs_assigned",
        "mmseqs_assigned_species", "mmseqs_assigned_highrank",
        "ref_aligned_frac", "k2_frac", "mmseqs_frac"}


def _is_placeholder(path):
    b = os.path.basename(path)
    return (not path) or b.startswith("NO_FILE") or b.startswith("~")


def _coerce(key, val):
    if key in _NUM and val not in (None, ""):
        try:
            f = float(val)
            return int(f) if f.is_integer() else f
        except (TypeError, ValueError):
            return val
    if key in ("novelty_flag", "gene_mode") and val not in (None, ""):
        try:
            return int(val)
        except (TypeError, ValueError):
            return val
    return val


def _read_tsv(path):
    """Return list of dict rows from a TSV (empty list if missing/placeholder/empty)."""
    if _is_placeholder(path) or not os.path.isfile(path):
        return []
    with open(path, newline="") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        rows = []
        for r in rdr:
            rows.append({k: _coerce(k, v) for k, v in r.items()})
        return rows


def _stem(path):
    b = os.path.basename(path)
    for suf in (".novelty.summary.tsv", ".novelty.candidates.tsv"):
        if b.endswith(suf):
            return b[: -len(suf)]
    return b.split(".")[0]


def _sample_of(row, path):
    return (row.get("sample") if row else None) or _stem(path)


def collect(summary_paths, cand_paths):
    """Build {sample: {"summary": {...}|None, "candidates": [..]}} keyed by sample name."""
    samples = {}

    for p in summary_paths:
        rows = _read_tsv(p)
        if not rows:
            continue
        # summary files carry exactly one data row
        row = rows[0]
        s = _sample_of(row, p)
        samples.setdefault(s, {"summary": None, "candidates": []})
        samples[s]["summary"] = row

    for p in cand_paths:
        rows = _read_tsv(p)
        if not rows:
            continue
        s = _sample_of(rows[0], p)
        samples.setdefault(s, {"summary": None, "candidates": []})
        # keep candidate order as written (already sorted by support upstream)
        samples[s]["candidates"].extend(rows)

    return samples


def _write_json(path, obj):
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, ensure_ascii=False, separators=(",", ":"))


def _write_xlsx(path, summary_rows, cand_rows):
    """Write a two-sheet workbook. Uses pandas+openpyxl; degrades to CSV siblings if absent."""
    try:
        import pandas as pd
    except ImportError:
        pd = None

    if pd is not None:
        try:
            sdf = pd.DataFrame(summary_rows, columns=SUMMARY_COLS) if summary_rows \
                else pd.DataFrame(columns=SUMMARY_COLS)
            cdf = pd.DataFrame(cand_rows, columns=CAND_COLS) if cand_rows \
                else pd.DataFrame(columns=CAND_COLS)
            with pd.ExcelWriter(path, engine="openpyxl") as w:
                sdf.to_excel(w, sheet_name="Summary", index=False)
                cdf.to_excel(w, sheet_name="Candidates", index=False)
            return
        except Exception as exc:  # noqa: BLE001 - never fail the whole report over a workbook
            sys.stderr.write(f"[novelty_collect] WARNING: xlsx write failed ({exc}); "
                             f"writing CSV fallbacks for {path}\n")

    base = path[:-5] if path.endswith(".xlsx") else path
    for name, cols, rows in (("summary", SUMMARY_COLS, summary_rows),
                             ("candidates", CAND_COLS, cand_rows)):
        with open(f"{base}.{name}.csv", "w", newline="") as fh:
            wri = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
            wri.writeheader()
            wri.writerows(rows)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="Collect NOVELTY per-sample TSVs into JSON/XLSX.")
    p.add_argument("--summaries", nargs="*", default=[], help="*.novelty.summary.tsv files")
    p.add_argument("--candidates", nargs="*", default=[], help="*.novelty.candidates.tsv files")
    p.add_argument("--outdir", default=".", help="output directory")
    p.add_argument("--combined-prefix", default="all.novelty",
                   help="basename prefix for the combined json/xlsx (default: all.novelty)")
    return p.parse_args(argv)


def main(argv=None):
    a = parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)

    summary_paths = [p for p in a.summaries if not _is_placeholder(p)]
    cand_paths = [p for p in a.candidates if not _is_placeholder(p)]

    samples = collect(summary_paths, cand_paths)

    # ---- per-sample artifacts ----
    for s, blocks in sorted(samples.items()):
        _write_json(os.path.join(a.outdir, f"{s}.novelty.json"), blocks)
        _write_xlsx(os.path.join(a.outdir, f"{s}.novelty.xlsx"),
                    [blocks["summary"]] if blocks.get("summary") else [],
                    blocks.get("candidates", []))

    # ---- combined artifacts ----
    # Surface the novelty backend at the top level so the report can label the tab without
    # digging into a sample. One run = one --novelty method, so the first non-empty wins.
    classifier = ""
    for blocks in samples.values():
        c = (blocks.get("summary") or {}).get("classifier")
        if c:
            classifier = c
            break
    # Gene mode (--novelty_gene): query was Pyrodigal-predicted genes, not whole contigs. Surface
    # at the top level so the report can switch its count-unit labels (contigs -> genes/seqs).
    gene_mode = 0
    for blocks in samples.values():
        if int((blocks.get("summary") or {}).get("gene_mode") or 0):
            gene_mode = 1
            break
    combined = {"taxtriage_novelty": True, "version": "1.0",
                "classifier": classifier, "gene_mode": gene_mode, "samples": samples}
    _write_json(os.path.join(a.outdir, f"{a.combined_prefix}.json"), combined)

    all_summary = [b["summary"] for b in samples.values() if b.get("summary")]
    all_cand = [r for b in samples.values() for r in b.get("candidates", [])]
    _write_xlsx(os.path.join(a.outdir, f"{a.combined_prefix}.xlsx"), all_summary, all_cand)

    sys.stderr.write(f"[novelty_collect] {len(samples)} sample(s); "
                     f"{len(all_cand)} candidate row(s) -> {a.combined_prefix}.json/.xlsx\n")


if __name__ == "__main__":
    main()
