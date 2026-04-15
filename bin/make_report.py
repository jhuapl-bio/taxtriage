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
make_report.py
==============
Build the all.comparison.report.html from either:
  a) one or more *.paths.json files produced by ALIGNMENT_PER_SAMPLE  (preferred)
  b) a single TSV/XLSX tabular file produced by ORGANISM_MERGE_REPORT (fallback)

Optionally accepts one or more protein-annotation XLSX files produced by the
--annotate_proteins / --annotate_meta steps (NOT use_diamond / get_features).
"""

import argparse
import glob
import json
import math
import os
import re
import sys

import pandas as pd


# ──────────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────────

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Build TaxTriage multi-run comparison HTML report.",
    )
    parser.add_argument(
        "-i", "--input", required=True, nargs="+",
        metavar="FILE",
        help="Input file(s): one or more .paths.json files  OR  a single "
             "TSV/XLSX tabular report (auto-detected by extension).",
    )
    parser.add_argument(
        "-p", "--protein_annotations", nargs="*", default=[],
        metavar="XLSX",
        help="Optional: protein-annotation XLSX file(s) produced by "
             "--annotate_proteins / --annotate_meta. "
             "Do NOT pass files from use_diamond or get_features here.",
    )
    parser.add_argument(
        "-t", "--template",
        metavar="TEMPLATE", default="heatmap.html",
        help="Input HTML template file (default: heatmap.html).",
    )
    parser.add_argument(
        "-m", '--mintass', default=None, type=float,
        help="Minimum TASS score (0-1) for inclusion in the report. "
             "If omitted, the minimum best_threshold across all input JSON files is used "
             "(derived from --cutoff_granularity; defaults to 0.0 if no threshold is found)."
    )
    parser.add_argument(
        "--cutoff_granularity", default="subkey",
        choices=["subkey", "key", "toplevelkey"],
        help="Which best_cutoffs granularity level to read from each JSON when --mintass is not "
             "supplied (default: subkey).",
    )
    parser.add_argument(
        "-pident", '--pident', default=0.0, type=float,
        help="Minimum percent identity (0-100) for inclusion in the report (default: 0.0, i.e. include all). Only used if you have VF/AMR results"
    )
    parser.add_argument(
        "-o", "--output",
        metavar="OUTPUT", default="all.comparison.report.html",
        help="Output HTML file.",
    )
    parser.add_argument(
        "-mc", "--microbial_category", nargs="+",
        default=["all"],
        metavar="CAT",
        help="Microbial category filter (default: Primary). "
             "Accepted values: Primary, Commensal, Opportunistic, Potential, Unknown. "
             "Use 'all' to include every category. "
             "Multiple values are allowed, e.g. -mc Primary Potential",
    )
    return parser.parse_args(argv)


# ──────────────────────────────────────────────────────────────────────────────
# JSON ingestion
# ──────────────────────────────────────────────────────────────────────────────

_TAX_RANKS = ["superkingdom", "phylum", "class", "order", "family", "genus"]


def _extract_best_threshold(paths, granularity="subkey"):
    """Scan metadata of each JSON to find the minimum best_threshold across all files.

    Returns the minimum value (0–1 scale) or None if no threshold is found.
    """
    thresholds = []
    for path in paths:
        path = path.strip()
        if not os.path.isfile(path):
            continue
        try:
            with open(path) as fh:
                data = json.load(fh)
        except Exception:
            continue
        bc = data.get("metadata", {}).get("best_cutoffs", {})
        t = (bc.get(granularity) or {}).get("best_threshold")
        if t is not None:
            thresholds.append(float(t))
    return min(thresholds) if thresholds else None

_VALID_MICROBIAL_CATS = {"Primary", "Commensal", "Opportunistic", "Potential", "Unknown"}


def _resolve_microbial_cats(cat_args):
    """Return a set of accepted microbial categories, or None to mean 'all'."""
    if not cat_args:
        return None
    lowered = [c.strip().lower() for c in cat_args]
    if "all" in lowered:
        return None  # no filtering
    result = set()
    for raw in cat_args:
        canon = raw.strip().capitalize()
        # handle "opportunistic" -> "Opportunistic", title-case normalisation
        # do a case-insensitive lookup against valid set
        matched = next(
            (v for v in _VALID_MICROBIAL_CATS if v.lower() == raw.strip().lower()),
            None,
        )
        if matched:
            result.add(matched)
        else:
            print(
                f"[make_report] WARNING: --microbial_category value {raw!r} is not recognised "
                f"(valid: {sorted(_VALID_MICROBIAL_CATS)} or 'all'); ignoring.",
                file=sys.stderr,
            )
    return result or None  # if everything was invalid, fall back to no filter


def _flatten_organism(org, sample_name, sample_type, total_reads):
    """
    Flatten one organism entry (any hierarchy level) into a flat dict
    suitable for the tabular view and all plots.
    """
    strain_reads = float(org.get("numreads", 0) or 0)
    pct = strain_reads / max(1, total_reads) * 100.0
    covered = int(org.get("covered_bases", 0) or 0)
    genome_len = int(org.get("length", 0) or 0)
    breadth_pct = round(covered / genome_len * 100, 2) if genome_len > 0 else 0.0
    tass = float(org.get("tass_score", 0) or 0)

    tax = org.get("taxonomy", {})

    return {
        "Specimen ID":         sample_name,
        "Sample Type":         sample_type,
        "Detected Organism":   org.get("name", "Unknown"),
        "Taxonomic ID #":      str(org.get("key", "")),
        "Subkey":              str(org.get("subkey", org.get("key", ""))),
        "Microbial Category":  org.get("microbial_category", "Unknown"),
        "Ann Class":           org.get("annClass", ""),
        "IsAnnotated":         "Yes" if org.get("is_annotated", "No") == "Yes" else "No",
        "High Consequence":    bool(org.get("high_cons", False)),
        "Status":              org.get("status", ""),
        "TASS Score":          round(tass * 100, 1),
        "# Reads Aligned":     int(strain_reads),
        "% Reads":             round(pct, 4),
        "Coverage":            round((org.get("coverage", 0) or 0) * 100, 1),
        "Covered Bases":       covered,
        "Genome Length (bp)":  genome_len,
        "Breadth %":           breadth_pct,
        "Mean Depth":          round(float(org.get("meandepth", 0) or 0), 2),
        "Gini Coefficient":    round(float(org.get("gini_coefficient", 0) or 0), 3),
        "Mean MapQ":           round(float(org.get("meanmapq", 0) or 0), 1),
        "Mean BaseQ":          round(float(org.get("meanbaseq", 0) or 0), 1),
        "Minhash Score":       round(float(org.get("minhash_reduction", 0) or 0), 3),
        "Breadth Score":       round(float(org.get("breadth_log_score", 0) or 0), 3),
        "MapQ Score":          round(float(org.get("mapq_score", 0) or 0), 3),
        "Disparity Score":     round(float(org.get("disparity", 0) or 0), 3),
        "Diamond Identity":    round(float(org.get("diamond_identity", 0) or 0), 1),
        "K2 Reads":            int(org.get("k2_reads", 0) or 0),
        "RPM":                 round(float(org.get("rpm", 0) or 0), 2),
        "RPKM":                round(float(org.get("rpkm", 0) or 0), 4),
        "Passes Threshold":    bool(org.get("passes_threshold", False)),
        # taxonomy
        "Superkingdom":        tax.get("superkingdom", ""),
        "Phylum":              tax.get("phylum", ""),
        "Class":               tax.get("class", ""),
        "Order":               tax.get("order", ""),
        "Family":              tax.get("family", ""),
        "Genus":               tax.get("genus", ""),
    }


def _iter_organisms(json_data, sample_name, mintass=0, microbial_cats=None):
    """Yield flat organism dicts from a parsed paths JSON.

    microbial_cats: set of accepted microbial_category strings, or None to allow all.
    """
    meta = json_data.get("metadata", {})
    sample_type = meta.get("sample_type", "unknown")
    total_reads = int(meta.get("total_reads", 1) or 1)

    def _cat_ok(org):
        if microbial_cats is None:
            return True
        return org.get("microbial_category", "Unknown") in microbial_cats

    for grp in json_data.get("organisms", []):
        for sk_m in grp.get("members", []):
            for strain in sk_m.get("members", []):
                if float(strain.get("tass_score", 0) or 0) >= mintass and _cat_ok(strain):
                    yield _flatten_organism(strain, sample_name, sample_type, total_reads)
            # if no nested members, use subkey level directly
            if not sk_m.get("members"):
                if float(sk_m.get("tass_score", 0) or 0) >= mintass and _cat_ok(sk_m):
                    yield _flatten_organism(sk_m, sample_name, sample_type, total_reads)
        # if no members at all, use group level
        if not grp.get("members"):
            if float(grp.get("tass_score", 0) or 0) >= mintass and _cat_ok(grp):
                yield _flatten_organism(grp, sample_name, sample_type, total_reads)


def load_json_inputs(paths, mintass=0, microbial_cats=None):
    """Return flat organism rows, per-sample metadata, and per-organism contig data."""
    rows = []
    sample_meta = {}
    # contig_data: dict keyed by "<sample>||<organism_name>||<taxon_id>"
    # value: {contigs: [...], depth_histogram: {...}}
    contig_data = {}

    for path in paths:
        path = path.strip()
        if not os.path.isfile(path):
            expanded = glob.glob(path)
            if not expanded:
                print(f"[make_report] WARNING: cannot find {path!r}, skipping", file=sys.stderr)
                continue
            for p in expanded:
                rows_, sm_, cd_ = load_json_inputs([p], mintass=mintass, microbial_cats=microbial_cats)
                rows.extend(rows_)
                sample_meta.update(sm_)
                contig_data.update(cd_)
            continue

        try:
            with open(path) as fh:
                data = json.load(fh)
        except Exception as exc:
            print(f"[make_report] WARNING: failed to parse {path}: {exc}", file=sys.stderr)
            continue

        meta = data.get("metadata", {})
        sample_name = meta.get("sample_name", os.path.basename(path).split(".")[0])
        sample_meta[sample_name] = meta

        for row in _iter_organisms(data, sample_name, mintass=mintass, microbial_cats=microbial_cats):
            rows.append(row)

        # Extract per-contig and depth-histogram data from each strain
        _STRIP = {'members', 'subkey', 'key', 'toplevelkey'}
        for grp in data.get("organisms", []):
            for sk_m in grp.get("members", []):
                if float(sk_m.get('tass_score', 0) or 0) < mintass:
                    continue
                if microbial_cats is not None and sk_m.get('microbial_category', 'Unknown') not in microbial_cats:
                    continue
                for strain in sk_m.get("members", []):
                    if float(strain.get('tass_score', 0) or 0) < mintass:
                        continue
                    if microbial_cats is not None and strain.get('microbial_category', 'Unknown') not in microbial_cats:
                        continue
                    _contigs = strain.get("contigs")
                    _dhist   = strain.get("depth_histogram")
                    if _contigs or _dhist:
                        _key = f"{sample_name}||{strain.get('name','')}||{strain.get('key','')}"
                        contig_data[_key] = {
                            "sample":          sample_name,
                            "organism":        strain.get("name", "Unknown"),
                            "taxon_id":        str(strain.get("key", "")),
                            "contigs":         [{k: v for k, v in c.items() if k not in _STRIP} for c in (_contigs or [])],
                            "depth_histogram": _dhist or {},
                        }

    return rows, sample_meta, contig_data


# ──────────────────────────────────────────────────────────────────────────────
# TSV / XLSX fallback ingestion
# ──────────────────────────────────────────────────────────────────────────────

def load_tabular_input(path, mintass=0, microbial_cats=None):
    ext = os.path.splitext(path)[1].lower()
    if ext in (".xlsx", ".xls"):
        df = pd.read_excel(path, dtype=str)
    else:
        df = pd.read_csv(path, sep="\t", dtype=str)
    df = df[df['TASS Score'].astype(float) >= mintass].copy()
    if microbial_cats is not None and 'Microbial Category' in df.columns:
        df = df[df['Microbial Category'].isin(microbial_cats)].copy()
    df.columns = df.columns.str.strip()
    if "Detected Organism" in df.columns:
        df["Detected Organism"] = df["Detected Organism"].str.replace("°", "", regex=False).str.strip()
    if "Index" in df.columns:
        df = df.drop(columns=["Index"], errors="ignore")
    if "index" in df.columns:
        df = df.drop(columns=["index"], errors="ignore")
    df = df.where(pd.notnull(df), None)
    return df.to_dict(orient="records"), {}


# ──────────────────────────────────────────────────────────────────────────────
# Protein annotation ingestion  (only --annotate_proteins / --annotate_meta)
# ──────────────────────────────────────────────────────────────────────────────

def load_protein_annotations(paths, pident=0):
    """
    Read one or more protein-annotation XLSX files (sheets: Genus Summary,
    Per-Gene Hits, Sample Overview, AMR Genes) and return a dict:
      {
        "genus_summary":   [row, ...],
        "per_gene_hits":   [row, ...],
        "sample_overview": [row, ...],
        "amr_genes":       [row, ...],
      }
    Rows from multiple files are concatenated.
    """
    out = {
        "genus_summary": [],
        "per_gene_hits": [],
        "sample_overview": [],
        "amr_genes": [],
    }
    sheet_map = {
        "Genus Summary":   "genus_summary",
        "Per-Gene Hits":   "per_gene_hits",
        "Sample Overview": "sample_overview",
        "AMR Genes":       "amr_genes",
    }

    for path in (paths or []):
        path = path.strip()
        if not path or not os.path.isfile(path):
            continue
        try:
            wb = pd.read_excel(path, sheet_name=None, dtype=str)
        except Exception as exc:
            print(f"[make_report] WARNING: cannot read {path}: {exc}", file=sys.stderr)
            continue
        for sheet_name, key in sheet_map.items():
            if sheet_name in wb:
                df = wb[sheet_name].where(pd.notnull(wb[sheet_name]), None)
                # filter where %id is >= pident (if that column exists)
                if '%id' in df.columns:
                    df = df[df['%id'].astype(float) >= pident].copy()
                out[key].extend(df.to_dict(orient="records"))
    return out


# ──────────────────────────────────────────────────────────────────────────────
# JSON sanitisation helpers
# ──────────────────────────────────────────────────────────────────────────────

def _sanitize(obj):
    if isinstance(obj, dict):
        return {k: _sanitize(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_sanitize(v) for v in obj]
    try:
        if pd.isna(obj):
            return None
    except Exception:
        pass
    if isinstance(obj, float):
        return obj if math.isfinite(obj) else None
    if hasattr(obj, "item"):       # numpy scalar
        return obj.item()
    return obj


# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    # ── detect input mode ─────────────────────────────────────────────────────
    is_json_mode = all(
        f.strip().endswith(".json")
        for f in args.input
        if f.strip()
    )

    microbial_cats = _resolve_microbial_cats(args.microbial_category)
    if microbial_cats is None:
        print("[make_report] Microbial category filter: all")
    else:
        print(f"[make_report] Microbial category filter: {sorted(microbial_cats)}")

    # ── resolve mintass ────────────────────────────────────────────────────────
    if args.mintass is not None:
        mintass = args.mintass
        print(f"[make_report] TASS threshold: {mintass:.4f} (0–1) = {mintass*100:.1f} (from --mintass)")
    elif is_json_mode:
        mintass = _extract_best_threshold(args.input, granularity=args.cutoff_granularity)
        if mintass is not None:
            print(f"[make_report] TASS threshold: {mintass:.4f} (0–1) = {mintass*100:.1f} "
                  f"(min best_threshold across JSONs, granularity={args.cutoff_granularity!r})")
        else:
            mintass = 0.0
            print("[make_report] TASS threshold: 0.0 (no best_cutoffs found in JSON files)")
    else:
        mintass = 0.0
        print("[make_report] TASS threshold: 0.0 (tabular input, no --mintass supplied)")

    if is_json_mode:
        rows, sample_meta, contig_data = load_json_inputs(
            args.input, mintass=mintass, microbial_cats=microbial_cats
        )
        print(f"[make_report] Loaded {len(rows)} organism rows from "
              f"{len(args.input)} JSON file(s); {len(contig_data)} organisms have contig data")
    else:
        contig_data = {}
        if len(args.input) > 1:
            print("[make_report] WARNING: multiple non-JSON inputs given; "
                  "using only the first.", file=sys.stderr)
        rows, sample_meta = load_tabular_input(args.input[0], mintass, microbial_cats=microbial_cats)
        print(f"[make_report] Loaded {len(rows)} rows from tabular file "
              f"{args.input[0]!r}")

    # ── derive column lists ────────────────────────────────────────────────────
    all_cols = list(rows[0].keys()) if rows else []
    numeric_cols = []
    if rows:
        for col in all_cols:
            vals = [r[col] for r in rows if r.get(col) is not None]
            if vals and all(
                isinstance(v, (int, float)) or
                (isinstance(v, str) and _is_numeric_str(v))
                for v in vals[:50]
            ):
                numeric_cols.append(col)

    # ── protein annotations ───────────────────────────────────────────────────
    prot_data = load_protein_annotations(args.protein_annotations, pident=args.pident)
    has_prot = any(len(v) > 0 for v in prot_data.values())
    print(f"[make_report] Protein annotations loaded: {has_prot} "
          f"({sum(len(v) for v in prot_data.values())} total rows)")

    # ── build bootstrap payload ───────────────────────────────────────────────
    payload = _sanitize({
        "records":          rows,
        "all_cols":         all_cols,
        "numeric_cols":     numeric_cols,
        "sample_meta":      sample_meta,
        "prot_data":        prot_data,
        "has_prot":         has_prot,
        "contig_data":      list(contig_data.values()),   # list of organism contig objects
        "mintass":          round(mintass * 100, 1), # 0–100 scale for UI filter-min
    })

    bootstrap_json = json.dumps(payload, ensure_ascii=False, allow_nan=False)
    bootstrap_json = bootstrap_json.replace("</", "<\\/")

    # ── render template ───────────────────────────────────────────────────────
    with open(args.template, "r", encoding="utf-8") as fh:
        tpl = fh.read()

    # Replace whatever is currently inside <script id="BOOTSTRAP">…</script>
    # (works whether the template has the __BOOTSTRAP_JSON__ placeholder or
    #  pre-embedded dummy data from a previous development build).
    html, n_replaced = re.subn(
        r'(<script id="BOOTSTRAP" type="application/json">).*?(</script>)',
        r'\g<1>' + bootstrap_json + r'\2',
        tpl,
        count=1,
        flags=re.DOTALL,
    )
    if not n_replaced:
        # Fallback: simple string replace for the original placeholder
        html = tpl.replace("__BOOTSTRAP_JSON__", bootstrap_json)

    with open(args.output, "w", encoding="utf-8") as fh:
        fh.write(html)

    print(f"[make_report] Written: {args.output}")


def _is_numeric_str(s):
    try:
        float(s)
        return True
    except (ValueError, TypeError):
        return False


if __name__ == "__main__":
    main()
