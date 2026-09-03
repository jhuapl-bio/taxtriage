#!/usr/bin/env python3
##############################################################################################
# Copyright 2025 The Johns Hopkins University Applied Physics Laboratory LLC
# All rights reserved.
# Permission is hereby granted, free of charge, to any person obtaining a copy of this
# software and associated documentation files (the "Software"), to deal in the Software
# without restriction, including without limitation the rights to use, copy, modify,
# merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
# PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
# OR OTHER DEALINGS IN THE SOFTWARE.
#

"""
combine_samples_json.py
=======================
Merge all per-sample *.paths.json files produced by ALIGNMENT_PER_SAMPLE into
a single all.samples.json placed in the report folder.

Output format
-------------
{
  "taxtriage_combined": true,
  "version": "1.0",
  "samples": [
    { <full contents of sample_A.paths.json> },
    { <full contents of sample_B.paths.json> },
    ...
  ]
}

This file can be:
  1. Dragged onto any TaxTriage heatmap.html report for instant multi-sample import.
  2. Passed directly to make_report.py as a single -i argument to rebuild the HTML.

Individual per-sample JSONs are still written to the alignment/ folder by
ALIGNMENT_PER_SAMPLE — this script just aggregates them.

Optional extra feeds
--------------------
--novelty <all.novelty.json>       adds  "novelty" + "has_novelty"
--annotate-reports <*.tsv ...>     adds  "prot_data" + "has_prot"

Both are folded into the SAME file so a single drag-and-drop carries the
Novelty and VF/AMR panels too, instead of the user having to locate and drop
report/all.novelty.json and annotate/<sample>.annotate_report.tsv separately.
The annotation rows are emitted in the report's own prot_data column shape
(the same one bin/make_report.py builds), not as a raw copy of the TSV — the
upstream table carries large free-text columns the report never renders.
"""

import argparse
import csv
import glob
import json
import os
import sys

# The `sites` column of an annotate_report can be far larger than csv's default
# 128 KB field cap; without this the reader raises on perfectly valid rows.
try:
    csv.field_size_limit(sys.maxsize)
except (OverflowError, ValueError):  # 32-bit platforms reject sys.maxsize
    csv.field_size_limit(2**31 - 1)


# ── VF/AMR ingestion (mirrors bin/make_report.py) ───────────────────────────
PROT_KEYS = ("genus_summary", "per_gene_hits", "sample_overview", "amr_genes")


def _is_amr_annotation(prop, antibiotics_class, antibiotics, source):
    """Does an annotation row describe Antimicrobial Resistance?

    Kept in step with _is_amr_annotation() in bin/make_report.py: the merged
    report splits VF / Drug Target / Transporter hits into "Per-Gene Hits" and
    AMR hits into "AMR Genes", but annotate_report.tsv mixes them in one table.
    """
    def _s(v):
        return "" if v is None else str(v)

    p = _s(prop).strip().lower()
    if any(w in p for w in ("resist", "amr", "antibiotic", "antimicrobial")):
        return True
    if _s(antibiotics_class).strip() or _s(antibiotics).strip():
        return True
    if any(w in _s(source).strip().lower()
           for w in ("card", "ncbiamr", "resfinder", "argannot", "amrfinder")):
        return True
    return False


def _annot_row(r, sample, pident):
    """Build one report-shaped prot row from a raw annotate_report record.

    Returns (row, is_amr) or None when the row is filtered out by --pident.
    """
    def g(*keys):
        for k in keys:
            v = r.get(k)
            if v is None or v == "":
                continue
            # Empty cells can survive as NaN floats even with dtype=str.
            if isinstance(v, float) and v != v:
                continue
            return v
        return None

    pid_raw = g("pident", "%id")
    try:
        if pid_raw is not None and float(pid_raw) < pident:
            return None
    except (TypeError, ValueError):
        pass

    prop = g("property", "Property")
    abx_class = g("antibiotics_class", "Antibiotics Class")
    abx = g("antibiotics", "Antibiotics")
    source = g("source", "Source")
    row = {
        "Specimen ID":        g("Specimen ID", "sample", "sample_name") or sample,
        "Genus":              g("genus", "Genus") or "Unknown",
        "Species":            g("species", "Species"),
        "Gene":               g("gene_name", "Gene", "Gene Name"),
        "Product":            g("product", "Product"),
        "Property":           prop,
        "Description":        g("function", "product", "Description"),
        "Antibiotics Class":  abx_class,
        "Antibiotics":        abx,
        "Source":             source,
        "Source ID":          g("source_id", "Source ID"),
        "%id":                pid_raw,
        "E-value":            g("evalue", "E-value"),
        "Bitscore":           g("bitscore", "Bitscore"),
        "Reference Organism": g("organism", "Reference Organism"),
        "Level":              g("level", "Level"),
        "taxids": {
            k: v for k, v in {
                "species_taxid": g("species_taxid"),
                "taxon_id":      g("taxon_id"),
                "genus_taxid":   g("genus_taxid"),
            }.items() if v
        },
    }
    # Drop empty keys. The report reads these with `?? ""` fallbacks, so absent
    # and null render identically — but nulls are pure weight in a file this
    # size.
    row = {k: v for k, v in row.items() if v not in (None, "", {})}
    return row, _is_amr_annotation(prop, abx_class, abx, source)


def _read_annot_records(path, name):
    """Yield raw dict records from an annotate report (.tsv/.csv/.txt or .xlsx).

    Excel needs pandas + openpyxl. Both are present in the process container
    (bin/make_report.py already does pd.read_excel there), but the import is
    lazy and failure is non-fatal so this script still runs anywhere.
    """
    low = name.lower()
    if low.endswith((".xlsx", ".xls")):
        try:
            import pandas as pd
        except ImportError:
            print(f"[combine_samples_json] WARNING: {name!r} is Excel but pandas is not "
                  f"available here — pass the .tsv instead; skipping", file=sys.stderr)
            return []
        df = pd.read_excel(path, dtype=str)
        df = df.where(pd.notnull(df), None)
        return df.to_dict(orient="records")

    delim = "," if low.endswith(".csv") else "\t"
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        return list(csv.DictReader(fh, delimiter=delim))


def load_annotate_reports(paths, pident=0):
    """Read standalone annotate_report files into report-shaped prot_data.

    Accepts the .tsv AND the .xlsx form: the pipeline's PROTEINS subworkflow
    emits ANNOTATE_REPORT.out.xlsx (despite the channel being named
    ch_annotate_report_tsv), so an xlsx-only filter silently produced an
    all.odr.json with no VF/AMR block at all.

    Returns {genus_summary, per_gene_hits, sample_overview, amr_genes}.
    """
    out = {k: [] for k in PROT_KEYS}

    expanded = []
    for path in (paths or []):
        path = (path or "").strip()
        if not path:
            continue
        if os.path.isfile(path):
            expanded.append(path)
        else:
            expanded.extend(sorted(glob.glob(path)))

    for path in expanded:
        name = os.path.basename(path)
        if name.startswith("NO_FILE") or name.startswith("~"):
            continue
        if not name.lower().endswith((".tsv", ".csv", ".txt", ".xlsx", ".xls")):
            print(f"[combine_samples_json] WARNING: skipping unrecognised annotate report "
                  f"{name!r}", file=sys.stderr)
            continue
        # Recover the sample name from "<sample>.annotate_report.<ext>".
        sample = name.split(".annotate_report")[0]
        n_before = sum(len(v) for v in out.values())
        try:
            records = _read_annot_records(path, name)
        except Exception as exc:  # noqa: BLE001
            print(f"[combine_samples_json] WARNING: cannot read annotate report {path!r}: {exc}",
                  file=sys.stderr)
            continue
        for r in records:
            built = _annot_row(r, sample, pident)
            if built is None:
                continue
            row, is_amr = built
            out["amr_genes" if is_amr else "per_gene_hits"].append(row)
        n_added = sum(len(v) for v in out.values()) - n_before
        print(f"[combine_samples_json] Annotations {name!r} -> {sample!r} ({n_added} row(s))")

    return out


def load_novelty(path):
    """Read all.novelty.json into the report's novelty payload shape.

    Accepts the combined NOVELTY_COLLECT output
    ({samples: {<s>: {summary, candidates}}, classifier, gene_mode}) and the
    per-sample form ({summary, candidates}). Returns None when unusable.
    """
    path = (path or "").strip()
    if not path or os.path.basename(path).startswith("NO_FILE") or not os.path.isfile(path):
        return None
    try:
        with open(path, encoding="utf-8") as fh:
            data = json.load(fh)
    except Exception as exc:  # noqa: BLE001
        print(f"[combine_samples_json] WARNING: cannot read novelty json {path!r}: {exc}",
              file=sys.stderr)
        return None

    if isinstance(data.get("samples"), dict):
        samples = {s: v for s, v in data["samples"].items() if isinstance(v, dict)}
    elif "summary" in data or "candidates" in data:
        summary = data.get("summary") or {}
        name = summary.get("sample") or os.path.basename(path).split(".novelty")[0]
        samples = {name: {"summary": summary, "candidates": data.get("candidates") or []}}
    else:
        return None

    if not samples:
        return None
    return {
        "samples": samples,
        "classifier": data.get("classifier", "") or "",
        "gene_mode": data.get("gene_mode"),
    }


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Combine per-sample .paths.json files into a single all.samples.json.",
    )
    parser.add_argument(
        "-i", "--input", required=True, nargs="+",
        metavar="FILE",
        help="One or more .paths.json files (globs are expanded).",
    )
    parser.add_argument(
        "-o", "--output", default="all.samples.json",
        metavar="OUTPUT",
        help="Output combined JSON file (default: all.samples.json).",
    )
    parser.add_argument(
        "-n", "--novelty", default=None, metavar="JSON",
        help="Optional: all.novelty.json from NOVELTY_COLLECT. Embedded as a top-level "
             "'novelty' block so the Novelty panel works from this file alone.",
    )
    parser.add_argument(
        "--annotate-reports", nargs="*", default=[], metavar="FILE",
        help="Optional: standalone <sample>.annotate_report.tsv files. Embedded as a "
             "top-level 'prot_data' block so the VF/AMR tab works from this file alone.",
    )
    parser.add_argument(
        "--pident", type=float, default=0.0, metavar="FLOAT",
        help="Drop annotation rows below this %%identity (default: 0, keep all).",
    )
    return parser.parse_args(argv)


def main():
    args = parse_args()

    # Expand globs and collect all paths
    paths = []
    for pattern in args.input:
        pattern = pattern.strip()
        if os.path.isfile(pattern):
            paths.append(pattern)
        else:
            expanded = glob.glob(pattern)
            if expanded:
                paths.extend(sorted(expanded))
            else:
                print(f"[combine_samples_json] WARNING: no files found for {pattern!r}, skipping",
                      file=sys.stderr)

    if not paths:
        print("[combine_samples_json] ERROR: no input files found.", file=sys.stderr)
        sys.exit(1)

    samples = []
    for path in paths:
        try:
            with open(path, "r", encoding="utf-8") as fh:
                data = json.load(fh)
        except Exception as exc:
            print(f"[combine_samples_json] WARNING: failed to parse {path!r}: {exc}",
                  file=sys.stderr)
            continue

        # If somehow a combined file was passed, expand it rather than nesting
        if data.get("taxtriage_combined"):
            samples.extend(data.get("samples", []))
            print(f"[combine_samples_json] Expanded combined file {path!r} "
                  f"({len(data.get('samples', []))} samples)")
        else:
            sample_name = (data.get("metadata") or {}).get("sample_name",
                           os.path.basename(path).split(".")[0])
            samples.append(data)
            print(f"[combine_samples_json] Loaded {path!r} -> {sample_name!r}")

    combined = {
        "taxtriage_combined": True,
        "version": "1.1",
        "samples": samples,
    }

    # ── optional novelty block ────────────────────────────────────────────
    novelty = load_novelty(args.novelty)
    if novelty:
        combined["novelty"] = novelty
        combined["has_novelty"] = True
        print(f"[combine_samples_json] Novelty embedded: {len(novelty['samples'])} sample(s)")

    # ── optional VF/AMR block ─────────────────────────────────────────────
    prot_data = load_annotate_reports(args.annotate_reports, pident=args.pident)
    n_prot = sum(len(v) for v in prot_data.values())
    if n_prot:
        combined["prot_data"] = prot_data
        combined["has_prot"] = True
        print(f"[combine_samples_json] Annotations embedded: {n_prot} row(s) "
              f"({len(prot_data['per_gene_hits'])} VF, {len(prot_data['amr_genes'])} AMR)")

    with open(args.output, "w", encoding="utf-8") as fh:
        json.dump(combined, fh, ensure_ascii=False, allow_nan=False, separators=(",", ":"))

    size_mb = os.path.getsize(args.output) / (1024 * 1024)
    print(f"[combine_samples_json] Written {len(samples)} sample(s) -> {args.output!r} "
          f"({size_mb:.1f} MB)")


if __name__ == "__main__":
    main()
