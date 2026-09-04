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
import csv
import datetime
import glob
import json
import math
import os
import re
import sys
from collections import defaultdict

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
        "--annotate_reports", nargs="*", default=[],
        metavar="TSV",
        help="Optional: standalone per-sample annotation report TSV(s) produced by "
             "annotate_report.py (published under <outdir>/annotate/<sample>.annotate_report.tsv). "
             "These carry de-novo / unaligned VF-AMR hits for samples that have NO reference "
             "alignment, so their annotation is otherwise lost from the merged report. Only "
             "samples NOT already present in --protein_annotations are supplemented (no double-count).",
    )
    parser.add_argument(
        "--mintass", default=1.0, type=float,
        help="Minimum TASS score for inclusion in the report (default: 1.0, i.e. include all). "
             "This is a hard filter; organisms below this threshold will be excluded entirely. "
             "Note that the UI filter slider is pre-populated from best_cutoffs in the input data, "
             "so you can use that to set a more conservative default while still allowing users to see all organisms if they wish."
    )
    parser.add_argument(
        "-c", "--min_conf", default=None, type=float,
        help="Explicit TASS confidence cutoff (0-100). When set, this overrides the "
             "auto-computed best_cutoffs for every sample, so the global filter slider "
             "AND every per-sample-type slider in the UI default to this value instead "
             "of the thresholds-JSON-derived recommendation. Does not hard-filter data "
             "(use --mintass for that) -- it only changes the sliders' starting position."
    )
    parser.add_argument(
        "-t", "--template",
        metavar="TEMPLATE", default="heatmap.html",
        help="Input HTML template file (default: heatmap.html).",
    )
    parser.add_argument(
        "--offline_report", action="store_true",
        help="Download the report's CDN libraries (d3, xlsx, jspdf, Leaflet, "
             "Font Awesome) at build time and embed them inline so the report "
             "opens with no internet access. Requires network at build time.",
    )
    parser.add_argument(
        "--offline_report_files", default=None, metavar="DIR",
        help="Directory containing local copies of the CDN libraries (and any "
             "fonts/marker images they reference). When given, those files are "
             "embedded inline instead of downloading — a fully offline build. "
             "Takes precedence over --offline_report.",
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
        "-x", "--output_xlsx",
        metavar="XLSX", default=None,
        help="Optional: also write the merged data to this XLSX path. NOTE: this is "
             "NOT what -p/--protein_annotations does -- that flag takes XLSX INPUT "
             "file(s) to read, it does not produce output. Writes a 'Detections' "
             "sheet (same rows/columns as the HTML report's table), a 'Metadata' "
             "sheet when --metadata was loaded, and 'Genus Summary' / 'Per-Gene "
             "Hits' / 'Sample Overview' / 'AMR Genes' sheets when "
             "--protein_annotations / --annotate_reports were loaded.",
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
    parser.add_argument(
        "-n", "--novelty", default=None, metavar="JSON",
        help="Optional: combined all.novelty.json produced by NOVELTY_COLLECT. When given, the "
             "report shows a Novelty Detection panel (per-sample score/flag + candidate taxa).",
    )
    parser.add_argument(
        "--novelty-downloads", nargs="*", default=[], metavar="FILE",
        help="Optional: per-sample + combined novelty JSON/XLSX files to expose as download links. "
             "Only their basenames are used (the files are published alongside the report).",
    )
    parser.add_argument(
        "--pathogens", default=None, metavar="CSV",
        help="Optional: pathogen reference sheet (assets/pathogen_sheet.csv). When given, novelty "
             "candidates and VF/AMR hits are cross-referenced against it so the report can flag "
             "listed pathogens that have NO reference alignment (taxid, then name, then genus rollup).",
    )
    parser.add_argument(
        "--metadata", nargs="*", default=[], metavar="FILE",
        help="Optional: sample metadata CSV/TSV/XLSX, the same file you would drag onto the "
             "report's Metadata & Mapping tab. Must have a 'sample' (or 'sample_name') column; "
             "every other column is carried through as-is. A 'specimen' column groups a "
             "specimen's DNA/RNA libraries so the report can merge them. Rows whose sample id "
             "is not in the run are reported and skipped (see --metadata-add-unmatched).",
    )
    parser.add_argument(
        "--metadata-add-unmatched", action="store_true",
        help="With --metadata, also create entries for metadata rows whose sample id is not in "
             "the run. Off by default: an unmatched id is nearly always an id-format mismatch, "
             "and inventing samples hides that.",
    )
    parser.add_argument(
        "--metadata-sheet", default=None, metavar="NAME",
        help="With an XLSX --metadata file, the worksheet to read. Defaults to a sheet named "
             "'metadata' (case-insensitive) if present, else the first sheet.",
    )
    parser.add_argument(
        "--rename-samples", default=None, metavar="FILE",
        help="Optional: a 2-column CSV/TSV/XLSX mapping original sample names (1st "
             "column) to display names (2nd column). Applied after all data is loaded "
             "-- every occurrence of a matched sample name (records, metadata, contig "
             "data, protein/VF-AMR annotations, novelty) is renamed before the HTML is "
             "written. A header row is auto-detected and skipped if present.",
    )
    # ── Sample-QC flag defaults ──────────────────────────────────────────────
    # These seed the report's whole-SAMPLE rule set (the "Sample QC / Flags"
    # panel). They never drop data: every sample is still written to the report,
    # it is just marked -- or, with --flag-action hide, additionally removed from
    # the views. Users can edit every rule live in the report.
    flg = parser.add_argument_group("sample QC flags (report defaults)")
    flg.add_argument(
        "--flag-min-reads", default=None, type=float, metavar="N",
        help="Flag a sample whose total_reads is below N.",
    )
    flg.add_argument(
        "--flag-min-aligned-reads", default=None, type=float, metavar="N",
        help="Flag a sample whose aligned_reads is below N.",
    )
    flg.add_argument(
        "--flag-min-organisms", default=None, type=float, metavar="N",
        help="Flag a sample with fewer than N DISTINCT organisms scoring at or above "
             "--flag-organism-tass.",
    )
    flg.add_argument(
        "--flag-organism-tass", default=None, type=float, metavar="TASS",
        help="TASS cutoff used by --flag-min-organisms (default: --min_conf if given, else 75).",
    )
    flg.add_argument(
        "--flag-min-detections", default=None, type=float, metavar="N",
        help="Flag a sample with fewer than N detections passing their own threshold.",
    )
    flg.add_argument(
        "--flag-metadata", default=None, metavar="SPEC",
        help="Metadata criteria, ';'-separated. Each is 'field:op:value' (or the shorthand "
             "'field=value', which means equals). Operators: == != contains !contains regex "
             "empty !empty < <= > >=. Example: "
             "\"sample_type:==:nasal;host_disease:contains:influenza;site:!empty:\". "
             "The field is looked up in the run metadata first, then in the sample metadata.",
    )
    flg.add_argument(
        "--flag-logic", default="any", choices=["any", "all"],
        help="Flag a sample when ANY rule matches (default) or only when ALL of them do.",
    )
    flg.add_argument(
        "--flag-action", default="flag", choices=["flag", "hide"],
        help="What happens to a matching sample: 'flag' marks it in the Heatmap / Table / "
             "Metadata & Mapping / Summary tabs but leaves it fully visible (default); "
             "'hide' also removes it from every chart and table (reversible in the report).",
    )
    flg.add_argument(
        "--flag-view", default="all", choices=["all", "hide", "only"],
        help="Which samples the report OPENS on: 'all' (default; a rule with action 'hide' "
             "still hides its own matches), 'hide' (every flagged sample hidden) or 'only' "
             "(show the flagged samples and hide everything that passed). Switchable in the "
             "report's sidebar.",
    )
    flg.add_argument(
        "--flag-exclude-taxids", default=None, metavar="IDS",
        help="Taxids (or organism names) that never count toward a detection figure -- the "
             "distinct-organism counts, the detection counts and any aggregated column. "
             "Defaults to '9606' (human), so --flag-min-organisms means the same thing on a "
             "run with dehosting and one without. Comma- or space-separated; pass '' (empty) "
             "to count everything, including host.",
    )
    flg.add_argument(
        "--flag-missing", action="store_true",
        help="Treat a missing / blank value as a match. Off by default, so a rule can never "
             "flag a sample purely because the field was never populated.",
    )
    flg.add_argument(
        "--flag-rules", default=None, metavar="JSON",
        help="A JSON file holding the full rule list -- either a bare list of rule objects or "
             "{\"logic\":…, \"missing_fails\":…, \"rules\":[…]}. Each rule is "
             "{source, field, op, value[, agg, tass, action]} with source one of "
             "meta | derived | runmeta | data. Replaces every other --flag-* criterion.",
    )
    parser.add_argument(
        "--vfamr-taxids", default=None, metavar="TSV",
        help="Optional: bvbrc specialty-gene reference TSV "
             "(assets/bvbrc_specialty_genes_with_sequences_taxids_and_sites.tsv). Used to recover "
             "each VF/AMR hit's taxids from its Source ID so pathogen matching can key on taxid "
             "instead of the merged sheet's (often mis-parsed) Genus/Species text.",
    )
    parser.add_argument(
        "--insilico_params", default=None, metavar="JSON",
        help="Optional: JSON file of the in-silico subsampling run parameters "
             "(mode, series counts, replicates, seed, sim_nreads, iss_model, ...). "
             "Populates the provenance panel on the In-Silico suite tab. When absent, "
             "the params are inferred from the subsample sample names/metadata.",
    )
    parser.add_argument(
        "--insilico_json", nargs="*", default=[], metavar="JSON",
        help="Optional: per-dataset .paths.json file(s) for the in-silico subsample datasets "
             "(from ALIGNMENT_PER_SAMPLE_INSILICO). These are used ONLY to build the In-Silico "
             "suite tab (expected-vs-reality + dilution-series LoD); they are NOT added to the "
             "main multi-run heatmap/table so they don't skew cross-sample views.",
    )
    parser.add_argument(
        "--insilico_manifests", nargs="*", default=[], metavar="TSV",
        help="Optional: *_subsample_manifest.tsv file(s) produced by SUBSAMPLE_INSILICO. "
             "Provide authoritative target/actual read counts, total master reads and per-dataset "
             "seeds for the In-Silico suite tab. When absent, target counts are parsed from the "
             "subsample sample names.",
    )
    return parser.parse_args(argv)


# ──────────────────────────────────────────────────────────────────────────────
# Sample QC flags  (whole-sample rules baked into the report as defaults)
# ──────────────────────────────────────────────────────────────────────────────
# The report evaluates these client-side (assets/src/js/41_sample_flags.js), so
# all we do here is translate the --flag-* CLI surface into the rule objects that
# module understands and hand them over in the bootstrap payload. Nothing is
# filtered out of the data: a rule only decides how a sample is PRESENTED.
#
#   rule = {source, field, op, value, agg, tass, action}
#     source  meta     SAMPLE_META metric (total_reads, aligned_reads, platform…)
#             derived  computed from the detections (organisms above TASS, …)
#             runmeta  a metadata column (run metadata first, sample metadata as
#                      a fallback)
#             data     a numeric detection column, aggregated per sample by `agg`
#     action  flag  -> marked in the report   |   hide -> also removed from views

_FLAG_OPS = {
    "==", "!=", "contains", "!contains", "regex",
    "empty", "!empty", "<", "<=", ">", ">=",
}
_FLAG_SOURCES = {"meta", "derived", "runmeta", "data"}


def _parse_metadata_flag_specs(spec, action):
    """Parse --flag-metadata into runmeta rules.

    Accepts 'field:op:value' and the shorthand 'field=value' (equals), separated
    by ';'. A malformed clause is reported and skipped rather than raising --
    a typo in one criterion must not cost the user the whole report.
    """
    rules = []
    for raw in str(spec).split(";"):
        clause = raw.strip()
        if not clause:
            continue
        field = op = value = None
        if ":" in clause:
            parts = clause.split(":", 2)
            if len(parts) == 3 and parts[1].strip() in _FLAG_OPS:
                field, op, value = parts[0].strip(), parts[1].strip(), parts[2]
            elif len(parts) >= 2 and parts[1].strip() in ("empty", "!empty"):
                field, op, value = parts[0].strip(), parts[1].strip(), ""
        if field is None and "=" in clause:
            field, value = clause.split("=", 1)
            field, op = field.strip(), "=="
        if not field or op not in _FLAG_OPS:
            print(f"[make_report] WARNING: cannot parse --flag-metadata clause {clause!r}; "
                  f"expected 'field:op:value' with op in {sorted(_FLAG_OPS)}, or 'field=value'",
                  file=sys.stderr)
            continue
        rules.append({
            "source": "runmeta",
            "field": field,
            "op": op,
            "value": "" if value is None else str(value).strip(),
            "action": action,
        })
    return rules


def _normalize_flag_rule(rule, default_action):
    """Coerce one rule dict from a --flag-rules JSON file. Returns None if unusable."""
    if not isinstance(rule, dict):
        return None
    source = str(rule.get("source", "")).strip()
    field = str(rule.get("field", "")).strip()
    op = str(rule.get("op", "")).strip()
    if source not in _FLAG_SOURCES or not field or op not in _FLAG_OPS:
        print(f"[make_report] WARNING: skipping malformed flag rule {rule!r}", file=sys.stderr)
        return None
    out = {
        "on": rule.get("on", True) is not False,
        "source": source,
        "field": field,
        "op": op,
        "value": "" if rule.get("value") is None else str(rule.get("value")),
        "action": "hide" if str(rule.get("action", default_action)) == "hide" else "flag",
    }
    if source == "data":
        out["agg"] = str(rule.get("agg") or "max")
    if rule.get("tass") is not None:
        out["tass"] = float(rule["tass"])
    # Minimum classifier reads for the unsupported_k2_organisms field, the way
    # `tass` qualifies the distinct-organism counts.
    if rule.get("k2min") is not None:
        out["k2min"] = float(rule["k2min"])
    return out


def _flag_value(v):
    """Render a threshold for the rule payload without a spurious ".0".

    The --flag-* thresholds parse as float (so 0.5 works), which turns a plain
    count like 2000000000 into "2000000000.0" — noise in the report's rule
    editor and in every "actual:" line it prints.
    """
    if isinstance(v, float) and v.is_integer():
        return str(int(v))
    return str(v)


#: Host, by default. See --flag-exclude-taxids.
FLAG_EXCLUDE_DEFAULT = ["9606"]


def _parse_flag_exclude(spec):
    """Split --flag-exclude-taxids into a list; None means "use the default".

    An explicitly empty string is honoured as "exclude nothing" -- that is how a
    user asks for host to be counted like any other organism.
    """
    if spec is None:
        return list(FLAG_EXCLUDE_DEFAULT)
    out, seen = [], set()
    for tok in re.split(r"[,;\s]+", str(spec)):
        tok = tok.strip()
        if not tok or tok.lower() in seen:
            continue
        seen.add(tok.lower())
        out.append(tok)
    return out


def build_sample_flag_config(args):
    """Turn the --flag-* arguments into the report's default rule set.

    Returns None when nothing was requested, so a run with no flag parameters
    ships a report identical to what it produced before this feature existed.
    """
    action = getattr(args, "flag_action", "flag") or "flag"
    logic = getattr(args, "flag_logic", "any") or "any"
    missing = bool(getattr(args, "flag_missing", False))
    rules = []

    # A rules file replaces every individual criterion (the module docs and the
    # nextflow param both say so) -- it is the escape hatch for rule sets that
    # the flat CLI surface cannot express.
    rules_path = getattr(args, "flag_rules", None)
    if rules_path:
        try:
            with open(rules_path) as fh:
                blob = json.load(fh)
        except Exception as exc:
            print(f"[make_report] WARNING: could not read --flag-rules {rules_path}: {exc}",
                  file=sys.stderr)
            blob = None
        if isinstance(blob, dict):
            logic = blob.get("logic", logic)
            missing = bool(blob.get("missing_fails", missing))
            # A rules file may also carry the two whole-report settings.
            if blob.get("view") is not None:
                args.flag_view = blob.get("view")
            elif blob.get("hide_all"):
                args.flag_view = "hide"
            if blob.get("exclude_taxids") is not None:
                ex = blob.get("exclude_taxids")
                args.flag_exclude_taxids = ",".join(str(x) for x in ex) if isinstance(ex, list) else str(ex)
            raw_rules = blob.get("rules") or []
        elif isinstance(blob, list):
            raw_rules = blob
        else:
            raw_rules = []
        for r in raw_rules:
            norm = _normalize_flag_rule(r, action)
            if norm:
                rules.append(norm)
    else:
        if args.flag_min_reads is not None:
            rules.append({"source": "meta", "field": "total_reads", "op": "<",
                          "value": _flag_value(args.flag_min_reads), "action": action})
        if args.flag_min_aligned_reads is not None:
            rules.append({"source": "meta", "field": "aligned_reads", "op": "<",
                          "value": _flag_value(args.flag_min_aligned_reads), "action": action})
        if args.flag_min_organisms is not None:
            # The organism count is only meaningful next to a score cutoff; fall
            # back to --min_conf so "organisms above TASS" means the same thing
            # here as the slider the report opens with.
            tass = args.flag_organism_tass
            if tass is None:
                tass = args.min_conf if args.min_conf is not None else 75.0
            rules.append({"source": "derived", "field": "unique_taxids_above_tass", "op": "<",
                          "value": _flag_value(args.flag_min_organisms), "tass": float(tass),
                          "action": action})
        if args.flag_min_detections is not None:
            rules.append({"source": "derived", "field": "passing_detections", "op": "<",
                          "value": _flag_value(args.flag_min_detections), "action": action})
        if args.flag_metadata:
            rules.extend(_parse_metadata_flag_specs(args.flag_metadata, action))

    if not rules:
        return None
    view = getattr(args, "flag_view", "all") or "all"
    if view not in ("all", "hide", "only"):
        view = "all"
    return {
        "enabled": True,
        "logic": "all" if str(logic) == "all" else "any",
        "missing_fails": missing,
        "view": view,
        # Kept in step with `view` so a report built by this version still opens
        # correctly in anything that only knows the old boolean.
        "hide_all": view == "hide",
        "exclude_taxids": _parse_flag_exclude(getattr(args, "flag_exclude_taxids", None)),
        "rules": rules,
    }

# ──────────────────────────────────────────────────────────────────────────────
# Sample metadata (CSV / TSV / XLSX)
# ──────────────────────────────────────────────────────────────────────────────
# Build-time equivalent of dragging a metadata file onto the report's
# Metadata & Mapping tab. The parsing rules below deliberately mirror
# _rowToMetaRecord / _applyMetaRecords in assets/src/js/28_meta_csv.js so the
# same file produces the same report either way — if you change one, change the
# other, and scripts/test_metadata_cli.py will tell you if they drift.
#
# Rules (from the browser implementation):
#   • headers are lower-cased and trimmed; spaces are NOT converted
#   • the key column is "sample" or "sample_name"; whitespace → underscores
#   • every other column is carried through untouched
#   • "", "null", "na" (any case) become None
#   • latitude / longitude / depth / salinity are parsed as floats
_META_NUM_FIELDS = {"latitude", "longitude", "depth", "salinity"}
_META_NULL_STRINGS = {"", "null", "na"}


def _meta_read_table(path, sheet=None):
    """Read a metadata file into a list of dicts with lower-cased headers."""
    ext = os.path.splitext(path)[1].lower()
    if ext in (".xlsx", ".xls", ".xlsm"):
        # Match the browser: prefer a sheet literally named "metadata".
        if sheet is None:
            try:
                names = pd.ExcelFile(path).sheet_names
            except Exception:
                names = []
            sheet = next((n for n in names if str(n).strip().lower() == "metadata"), 0)
        df = pd.read_excel(path, sheet_name=sheet, dtype=str)
    elif ext in (".tsv", ".tab"):
        df = pd.read_csv(path, sep="\t", dtype=str)
    else:
        df = pd.read_csv(path, dtype=str)
    df.columns = [str(c).strip().lower() for c in df.columns]
    return df.where(pd.notna(df), None).to_dict("records")


def _meta_row_to_record(row):
    """One parsed row → a metadata record, or None when it has no sample id."""
    raw = row.get("sample") or row.get("sample_name") or ""
    sample = re.sub(r"\s+", "_", str(raw).strip())
    if not sample or sample.lower() in _META_NULL_STRINGS:
        return None
    rec = {"sample_name": sample}
    for k, v in row.items():
        if k in ("sample", "sample_name"):
            continue
        sv = "" if v is None else str(v).strip()
        if sv.lower() in _META_NULL_STRINGS:
            rec[k] = None
        elif k in _META_NUM_FIELDS:
            try:
                rec[k] = float(sv)
            except ValueError:
                rec[k] = None
        else:
            rec[k] = sv
    return rec


def load_sample_metadata(paths, sample_meta, known_samples=None, add_unmatched=False, sheet=None):
    """Merge metadata file(s) into sample_meta, in place.

    `known_samples` supplements sample_meta's keys when deciding whether a row
    matches the run. This matters for the tabular input path, where
    load_tabular_input() returns an EMPTY sample_meta — the samples are real,
    they're just only named in the data rows, and without this every metadata
    row would be reported as unmatched.

    Returns a summary dict for reporting. Unmatched rows are the failure mode
    worth surfacing: the file parses perfectly but its sample ids don't line up
    with the run, which otherwise shows up only as a report with no metadata.
    """
    summary = {"files": 0, "rows": 0, "matched": 0, "unmatched": [], "added": 0,
               "columns": set(), "specimens": 0}
    if not paths:
        return summary

    known = set(sample_meta.keys()) | set(known_samples or ())
    for path in paths:
        if not path or not str(path).strip():
            continue
        try:
            rows = _meta_read_table(path, sheet=sheet)
        except Exception as exc:
            print(f"[make_report] ERROR reading metadata {path!r}: {exc}", file=sys.stderr)
            continue
        summary["files"] += 1
        for row in rows:
            rec = _meta_row_to_record(row)
            if rec is None:
                continue
            summary["rows"] += 1
            name = rec["sample_name"]
            if name not in known:
                if not add_unmatched:
                    summary["unmatched"].append(name)
                    continue
                known.add(name)
                summary["added"] += 1
            else:
                summary["matched"] += 1
            target = sample_meta.setdefault(name, {"sample_name": name})
            for k, v in rec.items():
                if k == "sample_name":
                    continue
                summary["columns"].add(k)
                # A None from the metadata file means "blank cell". Don't let it
                # erase a value the pipeline already established for this sample.
                if v is None and target.get(k) not in (None, ""):
                    continue
                target[k] = v

    # How many multi-sample specimens the grouping column actually produced —
    # the number that tells you whether merge will do anything in the report.
    groups = {}
    for name, meta in sample_meta.items():
        spec = None
        for field in ("specimen", "specimen_id", "specimen id", "specimenid",
                      "specimen_group", "specimen group"):
            val = meta.get(field)
            if val is not None and str(val).strip():
                spec = str(val).strip()
                break
        groups.setdefault(spec or name, []).append(name)
    summary["specimens"] = sum(1 for members in groups.values() if len(members) > 1)
    return summary


# ──────────────────────────────────────────────────────────────────────────────
# Sample rename map (CSV / TSV / XLSX) — original name -> display name
# ──────────────────────────────────────────────────────────────────────────────
_RENAME_HEADER_STRINGS = {
    "old", "old_name", "original", "original_name", "sample", "sample_name",
    "new", "new_name", "renamed", "rename", "display_name", "display",
}


def load_sample_rename_map(path):
    """Read a 2-column CSV/TSV/XLSX into an {original_name: new_name} dict.

    The first column is the sample name as it appears in the input JSON/TSV
    (metadata.sample_name), the second is the desired display name. Extra
    columns are ignored. A header row is auto-detected (either cell matching
    a common header word) and skipped; otherwise every row is treated as data.
    """
    mapping = {}
    if not path:
        return mapping
    p = path.strip()
    if not p:
        return mapping
    if not os.path.isfile(p):
        print(f"[make_report] WARNING: --rename-samples file not found: {p}", file=sys.stderr)
        return mapping

    ext = os.path.splitext(p)[1].lower()
    try:
        if ext in (".xlsx", ".xls", ".xlsm"):
            df = pd.read_excel(p, dtype=str, header=None)
        elif ext in (".tsv", ".tab"):
            df = pd.read_csv(p, sep="\t", dtype=str, header=None)
        else:
            df = pd.read_csv(p, dtype=str, header=None)
    except Exception as exc:  # noqa: BLE001
        print(f"[make_report] WARNING: cannot read --rename-samples file {p}: {exc}", file=sys.stderr)
        return mapping

    if df.shape[1] < 2:
        print(f"[make_report] WARNING: --rename-samples file {p} needs at least 2 columns; ignoring.",
              file=sys.stderr)
        return mapping

    df = df.iloc[:, :2]
    df.columns = ["old", "new"]

    # Skip an obvious header row.
    if len(df) and (
        str(df.iloc[0]["old"]).strip().lower() in _RENAME_HEADER_STRINGS
        or str(df.iloc[0]["new"]).strip().lower() in _RENAME_HEADER_STRINGS
    ):
        df = df.iloc[1:]

    for _, row in df.iterrows():
        old = "" if row["old"] is None else str(row["old"]).strip()
        new = "" if row["new"] is None else str(row["new"]).strip()
        if not old or old.lower() in _META_NULL_STRINGS:
            continue
        if not new or new.lower() in _META_NULL_STRINGS:
            continue
        mapping[old] = new
    return mapping


def apply_sample_renames(rename_map, rows, sample_meta, contig_data, prot_data, novelty_data):
    """Rename sample identifiers everywhere they appear, in place. Returns a summary dict."""
    summary = {"applied": set(), "unmatched": []}
    if not rename_map:
        return {"applied": 0, "unmatched": []}

    def _new(name):
        if name in rename_map:
            summary["applied"].add(name)
            return rename_map[name]
        return name

    # Detection rows: "Specimen ID"
    for r in rows:
        old = r.get("Specimen ID")
        if old in rename_map:
            r["Specimen ID"] = _new(old)

    # Per-sample metadata, keyed by sample name.
    new_sample_meta = {}
    for name, meta in sample_meta.items():
        new_name = _new(name)
        meta = dict(meta)
        meta["sample_name"] = new_name
        if new_name in new_sample_meta:
            new_sample_meta[new_name].update(meta)
        else:
            new_sample_meta[new_name] = meta
    sample_meta.clear()
    sample_meta.update(new_sample_meta)

    # Contig data: keys are "<sample>||<organism>||<taxon_id>", entries carry "sample".
    new_contig_data = {}
    for key, entry in contig_data.items():
        old_sample = entry.get("sample", "")
        new_sample = _new(old_sample)
        entry = dict(entry)
        entry["sample"] = new_sample
        _, _, rest = key.partition("||")
        new_key = f"{new_sample}||{rest}" if rest else new_sample
        new_contig_data[new_key] = entry
    contig_data.clear()
    contig_data.update(new_contig_data)

    # Protein/VF-AMR annotation rows carry "Specimen ID" (or occasionally "Sample").
    for lst in prot_data.values():
        for r in lst:
            for field in ("Specimen ID", "Sample", "sample"):
                if field in r and r[field] in rename_map:
                    r[field] = _new(r[field])

    # Novelty payload: {"samples": {<sample>: {...}}}.
    if novelty_data and isinstance(novelty_data.get("samples"), dict):
        new_samples = {}
        for name, sdata in novelty_data["samples"].items():
            new_samples[_new(name)] = sdata
        novelty_data["samples"] = new_samples

    summary["unmatched"] = [old for old in rename_map if old not in summary["applied"]]
    summary["applied"] = len(summary["applied"])
    return summary


# ──────────────────────────────────────────────────────────────────────────────
# JSON ingestion
# ──────────────────────────────────────────────────────────────────────────────

_TAX_RANKS = ["domain", "kingdom", "phylum", "class", "order", "family", "genus"]


def _collect_best_cutoffs(sample_meta):
    """Aggregate best_cutoffs across all loaded samples.

    For each granularity level (key/subkey/toplevelkey), takes the minimum
    best_threshold so the UI starts at the most conservative recommended cutoff.
    Returns a dict in the same shape as a single sample's best_cutoffs, or None.
    """
    levels = ("key", "subkey", "toplevelkey")
    aggregated = {}
    for level in levels:
        thresholds = []
        for meta in sample_meta.values():
            bc = meta.get("best_cutoffs") or {}
            t = (bc.get(level) or {}).get("best_threshold")
            if t is not None:
                thresholds.append(float(t))
        if thresholds:
            # Use the sample whose best_threshold is the minimum as the representative
            best_meta = min(
                (m for m in sample_meta.values() if (m.get("best_cutoffs") or {}).get(level)),
                key=lambda m: float((m.get("best_cutoffs", {}).get(level) or {}).get("best_threshold") or 9999),
            )
            aggregated[level] = dict((best_meta.get("best_cutoffs") or {}).get(level) or {})
            aggregated[level]["best_threshold"] = min(thresholds)
    return aggregated or None

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


def _flatten_organism(org, sample_name, sample_type, total_reads,
                      species_parent=None, genus_parent=None, level="Strain"):
    """
    Flatten one organism entry (any hierarchy level) into a flat dict
    suitable for the tabular view and all plots.

    species_parent / genus_parent: the subkey (species) and toplevelkey (genus)
    group objects this organism rolls up into. Their tass_score is attached as
    Species TASS / Genus TASS so the UI can show that a strain failing its own
    threshold may still be detected at the species or genus level (LCA-aware
    rollup). For a row that IS already the species/genus level, the parent of
    that same level points to itself.
    """
    strain_reads = float(org.get("numreads", 0) or 0)
    pct = strain_reads / max(1, total_reads) * 100.0
    covered = int(org.get("covered_bases", 0) or 0)
    genome_len = int(org.get("length", 0) or 0)
    breadth_pct = round(min(100.0, covered / genome_len * 100), 2) if genome_len > 0 else 0.0
    tass = float(org.get("tass_score", 0) or 0)

    tax = org.get("taxonomy", {})

    # ── ANI annotation (set by match_paths.py when the ANI matrix is enabled) ──
    # Presence of the 'high_ani_matches' key — even as an empty list — signals
    # that ANI was computed for this run. Its absence means the data predates
    # ANI support, so ANI-dependent views fall back to an "out of date /
    # unsupported" state for this sample.
    _ani_annotated = 'high_ani_matches' in org
    _ani_list = [
        {"key": str(m.get("key", "")), "ani_pct": round(float(m.get("ani_pct", 0) or 0), 2)}
        for m in (org.get('high_ani_matches') or []) if isinstance(m, dict)
    ]

    # Parent-level TASS (species = subkey, genus = toplevelkey). Fall back to the
    # organism's own TASS when a parent level is absent so the rollup never
    # under-reports the row itself.
    _species_src = species_parent if species_parent is not None else org
    _genus_src = genus_parent if genus_parent is not None else (species_parent if species_parent is not None else org)
    species_tass = float(_species_src.get("tass_score", tass) or 0)
    genus_tass = float(_genus_src.get("tass_score", tass) or 0)
    species_name = str(_species_src.get("name", org.get("name", "")) or "")
    genus_name = str(_genus_src.get("name", "") or tax.get("genus", "") or "")

    return {
        "Specimen ID":         sample_name,
        "Sample Type":         sample_type,
        "Detected Organism":   org.get("name", "Unknown"),
        "TASS Score":          round(tass, 100),  # 4th col — 0–100 scale for display
        "Taxonomic ID #":      str(org.get("key", "")),
        "Subkey":              str(org.get("subkey", org.get("key", ""))),
        "Microbial Category":  org.get("microbial_category", "Unknown"),
        "Ann Class":           org.get("annClass", ""),
        "IsAnnotated":         "Yes" if org.get("is_annotated", "No") == "Yes" else "No",
        "High Consequence":    bool(org.get("high_cons", False)),
        "Mol Type":            org.get("mol_type", ""),
        "Status":              org.get("status", ""),
        "# Reads Aligned":     int(strain_reads),
        "% Reads":             round(pct, 4),
        "Coverage":            round(min(100.0, (org.get("coverage", 0) or 0) * 100), 1),
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
        # MicrobeRT (mmbert) classifier probability + model name. mmbert is a
        # 0–1 probability; surfaced here as a 0–100 % to match the PDF report.
        # Kept None when absent so the column renders blank rather than 0.
        "MicrobeRT Probability": (round(float(org.get("mmbert")) * 100, 2)
                                  if org.get("mmbert") not in (None, "")
                                  else None),
        "MicrobeRT Model":     org.get("mmbert_model") or "",
        "K2 Reads":            int(org.get("k2_reads", 0) or 0),
        "RPM":                 round(float(org.get("rpm", 0) or 0), 2),
        "RPKM":                round(float(org.get("rpkm", 0) or 0), 4),
        "Passes Threshold":    bool(org.get("passes_threshold", False)),
        # ANI annotation: capability flag + list of high-ANI partner taxa
        # ({key, ani_pct}). Consumed by the cross-sample Feature Compare view and
        # by client-side capability detection (absence ⇒ ANI unsupported).
        "ANI Annotated":       _ani_annotated,
        "High ANI Matches":    _ani_list,
        # Taxonomic rollup level of this row: "Strain" (key), "Species" (subkey),
        # or "Genus" (toplevelkey). Lets the UI switch the view granularity and
        # surface a species/genus summary row when its children fail their own
        # threshold. Strain rows are the default view.
        "Level":               level,
        # Parent-level rollup TASS — lets the UI show a strain that fails its own
        # threshold but is still detected at the species (subkey) or genus
        # (toplevelkey) level. Pass/fail itself is computed client-side against
        # the active threshold, so only the scores are emitted here.
        "Species TASS":        round(species_tass, 4),
        "Genus TASS":          round(genus_tass, 4),
        "Species Name":        species_name,
        "Genus Name":          genus_name,
        # taxonomy
        "Kingdom":             tax.get("kingdom", ""),
        "Domain":             tax.get("domain", ""),
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

    # grp = toplevelkey (genus) group; sk_m = subkey (species) group; strain = key.
    #
    # Leaf nodes (the organisms actually shown in the default view) are emitted as
    # Level="Strain". For every group that has children we ALSO emit a summary row
    # one level up — Level="Species" for a subkey group with strain members, and
    # Level="Genus" for a genus group with members — carrying that group's own
    # aggregate TASS/coverage. These summary rows let the UI roll up to the
    # species/genus level and surface a species hit even when every child strain
    # falls below the active cutoff. They are hidden in the default Strain view
    # unless promoted by the rollup.
    for grp in json_data.get("organisms", []):
        grp_has_members = bool(grp.get("members"))
        for sk_m in grp.get("members", []):
            sk_has_members = bool(sk_m.get("members"))
            for strain in sk_m.get("members", []):
                if float(strain.get("tass_score", 0) or 0) >= mintass and _cat_ok(strain):
                    yield _flatten_organism(strain, sample_name, sample_type, total_reads,
                                            species_parent=sk_m, genus_parent=grp,
                                            level="Strain")
            if sk_has_members:
                # species summary row over its strain members
                if float(sk_m.get("tass_score", 0) or 0) >= mintass and _cat_ok(sk_m):
                    yield _flatten_organism(sk_m, sample_name, sample_type, total_reads,
                                            species_parent=sk_m, genus_parent=grp,
                                            level="Species")
            else:
                # no nested strains: the subkey node is itself the leaf
                if float(sk_m.get("tass_score", 0) or 0) >= mintass and _cat_ok(sk_m):
                    yield _flatten_organism(sk_m, sample_name, sample_type, total_reads,
                                            species_parent=sk_m, genus_parent=grp,
                                            level="Strain")
        if grp_has_members:
            # genus summary row over its species/strain members
            if float(grp.get("tass_score", 0) or 0) >= mintass and _cat_ok(grp):
                yield _flatten_organism(grp, sample_name, sample_type, total_reads,
                                        species_parent=grp, genus_parent=grp,
                                        level="Genus")
        else:
            # no members at all: the genus node is itself the leaf
            if float(grp.get("tass_score", 0) or 0) >= mintass and _cat_ok(grp):
                yield _flatten_organism(grp, sample_name, sample_type, total_reads,
                                        species_parent=grp, genus_parent=grp,
                                        level="Strain")


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

        # ── Combined all.samples.json: expand into per-sample loads ──────────
        if data.get("taxtriage_combined") and isinstance(data.get("samples"), list):
            print(f"[make_report] Detected combined JSON {path!r} with "
                  f"{len(data['samples'])} sample(s); expanding...")
            for sample_data in data["samples"]:
                s_meta = sample_data.get("metadata", {})
                s_name = s_meta.get("sample_name",
                                     os.path.basename(path).split(".")[0])
                sample_meta[s_name] = s_meta
                for row in _iter_organisms(sample_data, s_name,
                                           mintass=mintass,
                                           microbial_cats=microbial_cats):
                    rows.append(row)
                _STRIP = {'members', 'subkey', 'key', 'toplevelkey'}
                for grp in sample_data.get("organisms", []):
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
                            _bhist   = strain.get("breadth_histogram")
                            if _contigs or _dhist or _bhist:
                                _key = f"{s_name}||{strain.get('name','')}||{strain.get('key','')}"
                                _cd_entry = {
                                    "sample":          s_name,
                                    "organism":        strain.get("name", "Unknown"),
                                    "taxon_id":        str(strain.get("key", "")),
                                    "contigs":         [{k: v for k, v in c.items() if k not in _STRIP} for c in (_contigs or [])],
                                    "depth_histogram": _dhist or {},
                                }
                                if _bhist:
                                    _cd_entry["breadth_histogram"] = _bhist
                                contig_data[_key] = _cd_entry
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
                    _bhist   = strain.get("breadth_histogram")
                    if _contigs or _dhist or _bhist:
                        _key = f"{sample_name}||{strain.get('name','')}||{strain.get('key','')}"
                        _cd_entry = {
                            "sample":          sample_name,
                            "organism":        strain.get("name", "Unknown"),
                            "taxon_id":        str(strain.get("key", "")),
                            "contigs":         [{k: v for k, v in c.items() if k not in _STRIP} for c in (_contigs or [])],
                            "depth_histogram": _dhist or {},
                        }
                        if _bhist:
                            _cd_entry["breadth_histogram"] = _bhist
                        contig_data[_key] = _cd_entry

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

def load_novelty(novelty_json, download_paths):
    """
    Read the combined all.novelty.json (NOVELTY_COLLECT) and return
      (novelty_payload, downloads)
    where novelty_payload = {"samples": {<sample>: {"summary": {...}, "candidates": [...]}}}
    and downloads = [{"label": <sample or "All samples">, "kind": "json"|"xlsx",
                      "filename": <basename>}, ...] for the report's download links.
    Returns ({"samples": {}}, []) when nothing usable is provided.
    """
    payload = {"samples": {}, "classifier": ""}
    if novelty_json:
        p = novelty_json.strip()
        if p and os.path.isfile(p):
            try:
                with open(p, encoding="utf-8") as fh:
                    data = json.load(fh)
                payload["samples"] = data.get("samples", {}) or {}
                payload["classifier"] = data.get("classifier", "") or ""
            except Exception as exc:  # noqa: BLE001
                print(f"[make_report] WARNING: cannot read novelty json {p}: {exc}",
                      file=sys.stderr)

    # Build download link descriptors from the staged files (basenames only — the files are
    # published next to all.odr.html, so relative links resolve).
    downloads = []
    seen = set()
    for path in (download_paths or []):
        name = os.path.basename((path or "").strip())
        if not name or name in seen or name.startswith("NO_FILE") or name.startswith("~"):
            continue
        if not (name.endswith(".json") or name.endswith(".xlsx")):
            continue
        seen.add(name)
        kind = "xlsx" if name.endswith(".xlsx") else "json"
        if name.startswith("all.novelty"):
            label = "All samples (combined)"
        else:
            # strip the .novelty.json / .novelty.xlsx suffix to recover the sample name
            label = name.split(".novelty.")[0]
        downloads.append({"label": label, "kind": kind, "filename": name})

    # Stable order: combined first, then per-sample alphabetical, xlsx before json within a group.
    downloads.sort(key=lambda d: (not d["filename"].startswith("all.novelty"),
                                  d["label"].lower(), d["kind"]))
    return payload, downloads


# ──────────────────────────────────────────────────────────────────────────────
# Pathogen reference sheet (assets/pathogen_sheet.csv)
# ──────────────────────────────────────────────────────────────────────────────

# Severity ranking used when a genus rolls up many species with mixed status.
_PATH_SEVERITY = {"primary": 4, "opportunistic": 3, "potential": 2, "commensal": 1}

# Ranks for which a first-token / name genus rollup is meaningful. Deliberately
# excludes "no rank" so phage / environmental / clade entries (e.g. "Bacillus
# phage SPBc2") are NOT mis-attributed to a host genus. Exact taxid/name matches
# still apply to every candidate regardless of rank.
_PATH_ROLLUP_RANKS = {"genus", "species", "subspecies", "strain", "serotype", "serovar"}


def load_pathogens(path):
    """
    Read the pathogen reference sheet into compact lookups the report can use to
    flag novelty candidates and VF/AMR hits that correspond to listed pathogens
    even when nothing aligned. Returns:

      {
        "by_taxid": {"<taxid>": {n,c,s,hc,g}, ...},   # exact NCBI taxid
        "by_name":  {"<lower name>": {n,c,s,hc,g}, ...},
        "by_genus": {"<lower genus>": {c,hc,n}, ...},  # genus rollup (most-severe class)
      }

    where each record carries: n=name, c=general_classification (lower),
    s=status, hc=high_consequence (bool), g=genus (lower).
    Returns empty lookups when no sheet is supplied or it cannot be read.
    """
    out = {"by_taxid": {}, "by_name": {}, "by_genus": {}}
    if not path:
        return out
    p = path.strip()
    if not p or not os.path.isfile(p):
        if p:
            print(f"[make_report] WARNING: pathogen sheet not found: {p}", file=sys.stderr)
        return out

    by_taxid, by_name, by_genus = out["by_taxid"], out["by_name"], out["by_genus"]
    try:
        with open(p, newline="", encoding="utf-8-sig") as fh:
            for r in csv.DictReader(fh):
                taxid = (r.get("taxid") or "").strip()
                name = (r.get("name") or "").strip()
                cls = (r.get("general_classification") or "").strip().lower()
                status = (r.get("status") or "").strip()
                hc = (r.get("high_consequence") or "").strip().upper() in ("TRUE", "1", "YES", "Y")
                genus = (r.get("genus") or "").strip()
                gl = genus.lower()
                rec = {"n": name, "c": cls, "s": status, "hc": hc, "g": gl}
                if taxid:
                    by_taxid.setdefault(taxid, rec)
                if name:
                    by_name.setdefault(name.lower(), rec)
                if gl:
                    g = by_genus.get(gl)
                    if g is None:
                        by_genus[gl] = {"c": cls, "hc": hc, "n": 1}
                    else:
                        g["n"] += 1
                        if _PATH_SEVERITY.get(cls, 0) > _PATH_SEVERITY.get(g["c"], 0):
                            g["c"] = cls
                        g["hc"] = g["hc"] or hc
    except Exception as exc:  # noqa: BLE001
        print(f"[make_report] WARNING: cannot read pathogen sheet {p}: {exc}", file=sys.stderr)
    return out


def _match_pathogen(taxid, name, rank, paths):
    """
    Resolve a (taxid, name, rank) triple against the pathogen lookups.
    Order: exact taxid -> exact name -> genus rollup (genus-rank uses the name;
    species/strain uses the first name token). Returns a compact descriptor or None.
    """
    by_taxid, by_name, by_genus = paths["by_taxid"], paths["by_name"], paths["by_genus"]
    taxid = str(taxid or "").strip()
    name = str(name or "").strip()
    rank = str(rank or "").strip().lower()

    if taxid and taxid in by_taxid:
        h = by_taxid[taxid]
        return {"match": "taxid", "name": h["n"], "class": h["c"], "status": h["s"], "hc": bool(h["hc"])}
    if name and name.lower() in by_name:
        h = by_name[name.lower()]
        return {"match": "name", "name": h["n"], "class": h["c"], "status": h["s"], "hc": bool(h["hc"])}
    if name and by_genus and rank in _PATH_ROLLUP_RANKS:
        if rank == "genus":
            gkey = name.lower()
        else:
            # First-token rollup for species/strain. Skip phage/host-virus/satellite
            # species (e.g. "Escherichia phage D6") so they are NOT mis-attributed to
            # the host bacterial genus; an exact taxid/name match would already win.
            low = name.lower()
            gkey = "" if any(w in low for w in (" phage", "prophage", " virus", " satellite")) \
                else name.split()[0].lower()
        if gkey and gkey in by_genus:
            g = by_genus[gkey]
            return {"match": "genus", "name": gkey.capitalize(), "class": g["c"],
                    "status": "", "hc": bool(g["hc"]), "genus_n": g.get("n", 0)}
    return None


def annotate_novelty_pathogens(novelty_data, paths):
    """
    Stamp each novelty candidate that corresponds to a listed pathogen with a
    `pathogen` descriptor (in place). No-op when no pathogen sheet was loaded.
    Returns the number of candidates flagged.
    """
    if not (paths["by_taxid"] or paths["by_name"] or paths["by_genus"]):
        return 0
    n = 0
    for _sname, sdata in (novelty_data.get("samples") or {}).items():
        for c in (sdata.get("candidates") or []):
            hit = _match_pathogen(c.get("taxid"), c.get("name"), c.get("rank"), paths)
            if hit:
                c["pathogen"] = hit
                n += 1
    return n


# ──────────────────────────────────────────────────────────────────────────────
# VF/AMR taxid resolution (bvbrc specialty-gene reference) + pathogen stamping
# ──────────────────────────────────────────────────────────────────────────────

def load_source_taxids(path):
    """
    Read the bvbrc specialty-gene reference TSV into a lookup
        source_id -> {"species_taxid": .., "taxon_id": .., "genus_taxid": ..}
    so VF/AMR hits (which carry a Source ID but lose their taxids in the merged
    Per-Gene Hits / AMR Genes sheets) can be re-associated with a canonical taxid.
    Returns {} when no file is supplied / readable.
    """
    out = {}
    if not path:
        return out
    p = path.strip()
    if not p or not os.path.isfile(p):
        if p:
            print(f"[make_report] WARNING: VF/AMR source-taxid TSV not found: {p}", file=sys.stderr)
        return out
    try:
        with open(p, newline="", encoding="utf-8-sig") as fh:
            rdr = csv.DictReader(fh, delimiter="\t")
            for row in rdr:
                sid = (row.get("source_id") or "").strip()
                if not sid or sid in out:
                    continue
                out[sid] = {
                    "species_taxid": (row.get("species_taxid") or "").strip(),
                    "taxon_id": (row.get("taxon_id") or "").strip(),
                    "genus_taxid": (row.get("genus_taxid") or "").strip(),
                }
    except Exception as exc:  # noqa: BLE001
        print(f"[make_report] WARNING: cannot read VF/AMR source-taxid TSV {p}: {exc}", file=sys.stderr)
    return out


def _match_pathogen_taxids(taxids, name, genus, paths):
    """
    Pathogen-sheet match for a VF/AMR hit. Tries, in order: each candidate taxid
    (species_taxid → taxon_id → genus_taxid) against the sheet's taxid index, then
    the species name (exact, then 'Genus species' two-token), then the genus name.
    Returns a compact descriptor (carrying the canonical sheet genus, so genus
    rollup works even when the merged sheet's Genus column is mis-parsed) or None.
    """
    by_taxid, by_name, by_genus = paths["by_taxid"], paths["by_name"], paths["by_genus"]
    for t in taxids:
        t = str(t or "").strip()
        if t and t in by_taxid:
            h = by_taxid[t]
            return {"match": "taxid", "name": h["n"], "class": h["c"], "status": h["s"],
                    "hc": bool(h["hc"]), "genus": h["g"], "taxid": t}
    nl = str(name or "").strip().lower()
    if nl and nl in by_name:
        h = by_name[nl]
        return {"match": "name", "name": h["n"], "class": h["c"], "status": h["s"],
                "hc": bool(h["hc"]), "genus": h["g"]}
    if nl:
        two = " ".join(nl.split()[:2])
        if two and two in by_name:
            h = by_name[two]
            return {"match": "name", "name": h["n"], "class": h["c"], "status": h["s"],
                    "hc": bool(h["hc"]), "genus": h["g"]}
    gl = str(genus or "").strip().lower()
    if gl and gl in by_genus:
        g = by_genus[gl]
        return {"match": "genus", "name": gl.capitalize(), "class": g["c"], "status": "",
                "hc": bool(g["hc"]), "genus": gl, "genus_n": g.get("n", 0)}
    return None


def annotate_protein_pathogens(prot_data, source_taxids, paths):
    """
    Stamp each VF/AMR row (Per-Gene Hits + AMR Genes) with a `pathogen` descriptor
    and the taxids resolved from its Source ID (via the bvbrc reference). Lets the
    report flag listed pathogens by canonical taxid instead of relying on the merged
    sheet's (often mis-parsed) Genus/Species text. Returns the number of rows flagged.
    """
    if not (paths["by_taxid"] or paths["by_name"] or paths["by_genus"]):
        return 0
    n = 0
    for key in ("per_gene_hits", "amr_genes"):
        for r in (prot_data.get(key) or []):
            sid = str(r.get("Source ID") or r.get("source_id") or "").strip()
            st = source_taxids.get(sid, {})
            taxids = [st.get("species_taxid"), st.get("taxon_id"), st.get("genus_taxid")]
            if any(taxids):
                r["taxids"] = {k: v for k, v in st.items() if v}
            hit = _match_pathogen_taxids(
                taxids, r.get("Species") or r.get("species"),
                r.get("Genus") or r.get("genus"), paths,
            )
            if hit:
                r["pathogen"] = hit
                n += 1
    return n


# Keys of prot_data that actually carry VF/AMR annotation.
#   PROT_HIT_KEYS        — per-hit rows; a sample present here genuinely has hits.
#   PROT_ANNOTATION_KEYS — the above plus the genus rollup; any of these being
#                          non-empty means the VF/AMR tab has something to draw.
# "sample_overview" is deliberately excluded from both: create_report.py writes
# that sheet for every run (it is the TASS organism table), so treating it as
# annotation makes has_prot always True and marks every sample as covered.
PROT_HIT_KEYS = ("per_gene_hits", "amr_genes")
PROT_ANNOTATION_KEYS = PROT_HIT_KEYS + ("genus_summary",)


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


def _is_amr_annotation(prop, antibiotics_class, antibiotics, source):
    """Heuristic: does an annotation row describe Antimicrobial Resistance?

    The merged report keeps VF / Drug Target / Transporter hits in "Per-Gene Hits"
    and AMR hits in "AMR Genes". The standalone annotate_report.tsv mixes them in
    one table, so route each row by its property / antibiotics fields / source.
    """
    def _s(v):
        if v is None:
            return ""
        if isinstance(v, float) and math.isnan(v):
            return ""
        return str(v)
    p = _s(prop).strip().lower()
    if any(w in p for w in ("resist", "amr", "antibiotic", "antimicrobial")):
        return True
    if _s(antibiotics_class).strip() or _s(antibiotics).strip():
        return True
    if any(w in _s(source).strip().lower() for w in ("card", "ncbiamr", "resfinder", "argannot", "amrfinder")):
        return True
    return False


def load_standalone_annotations(paths, covered_samples, pident=0):
    """
    Build supplemental VF/AMR annotation rows from standalone annotate_report.tsv
    files (annotate_report.py output) for samples NOT already represented in the
    merged protein-annotation XLSX. Samples with no reference alignment never get
    an organism hierarchy, so their de-novo annotation is otherwise dropped.

    Returns a dict shaped like load_protein_annotations()
    ({genus_summary, per_gene_hits, sample_overview, amr_genes}). Rows are emitted
    in the same column shape the report's prot_data expects, so they merge cleanly
    and still get pathogen-stamped by annotate_protein_pathogens().
    """
    out = {"genus_summary": [], "per_gene_hits": [], "sample_overview": [], "amr_genes": []}
    covered = {str(s) for s in (covered_samples or set())}
    # genus rollup accumulator: (sample, genus, property) -> {genes:set, ids:[%id], evals:[]}
    roll = {}

    expanded = []
    for path in (paths or []):
        path = (path or "").strip()
        if not path:
            continue
        if os.path.isfile(path):
            expanded.append(path)
        else:
            expanded.extend(glob.glob(path))

    for path in expanded:
        name = os.path.basename(path)
        if name.startswith("NO_FILE") or name.startswith("~"):
            continue
        # Recover the sample name from "<sample>.annotate_report.tsv/.xlsx".
        sample = name.split(".annotate_report")[0]
        if not sample or sample in covered:
            continue  # already supplied by the merged XLSX — avoid double counting
        try:
            if path.endswith(".xlsx") or path.endswith(".xls"):
                df = pd.read_excel(path, dtype=str)
            else:
                df = pd.read_csv(path, sep="\t", dtype=str)
        except Exception as exc:  # noqa: BLE001
            print(f"[make_report] WARNING: cannot read annotate report {path}: {exc}", file=sys.stderr)
            continue
        df = df.where(pd.notnull(df), None)
        for r in df.to_dict(orient="records"):
            def g(*keys):
                for k in keys:
                    v = r.get(k)
                    if v is None or v == "":
                        continue
                    # Empty cells can survive as NaN floats despite dtype=str.
                    if isinstance(v, float) and math.isnan(v):
                        continue
                    return v
                return None
            pid_raw = g("pident", "%id")
            try:
                if pid_raw is not None and float(pid_raw) < pident:
                    continue
            except (TypeError, ValueError):
                pass
            prop = g("property", "Property")
            abx_class = g("antibiotics_class", "Antibiotics Class")
            abx = g("antibiotics", "Antibiotics")
            source = g("source", "Source")
            gene = g("gene_name", "Gene", "Gene Name")
            genus = g("genus", "Genus") or "Unknown"
            species = g("species", "Species")
            row = {
                "Specimen ID":       sample,
                "Genus":             genus,
                "Species":           species,
                "Gene":              gene,
                "Product":           g("product", "Product"),
                "Property":          prop,
                "Description":       g("function", "product", "Description"),
                "Antibiotics Class": abx_class,
                "Antibiotics":       abx,
                "Source":            source,
                "Source ID":         g("source_id", "Source ID"),
                "%id":               pid_raw,
                "E-value":           g("evalue", "E-value"),
                "Bitscore":          g("bitscore", "Bitscore"),
                "Reference Organism": g("organism", "Reference Organism"),
                "Level":             g("level", "Level"),
                "taxids": {
                    k: v for k, v in {
                        "species_taxid": g("species_taxid"),
                        "taxon_id":      g("taxon_id"),
                        "genus_taxid":   g("genus_taxid"),
                    }.items() if v
                },
            }
            if _is_amr_annotation(prop, abx_class, abx, source):
                out["amr_genes"].append(row)
            else:
                out["per_gene_hits"].append(row)
            # genus rollup
            rk = (sample, genus, prop or "")
            acc = roll.setdefault(rk, {"genes": set(), "ids": [], "evals": []})
            if gene:
                acc["genes"].add(str(gene))
            try:
                acc["ids"].append(float(pid_raw))
            except (TypeError, ValueError):
                pass
            ev = g("evalue", "E-value")
            if ev is not None:
                acc["evals"].append(str(ev))

    for (sample, genus, prop), acc in roll.items():
        ids = acc["ids"]
        out["genus_summary"].append({
            "Sample":         sample,
            "Genus":          genus,
            "Property":       prop,
            "Genes":          ", ".join(sorted(acc["genes"])) if acc["genes"] else None,
            "# Hits":         str(len(acc["evals"]) or len(acc["genes"]) or 1),
            "Best %id":       (str(round(max(ids), 1)) if ids else None),
            "Avg %id":        (str(round(sum(ids) / len(ids), 1)) if ids else None),
            "Avg E-value":    (acc["evals"][0] if acc["evals"] else None),
            "Median E-value": (acc["evals"][len(acc["evals"]) // 2] if acc["evals"] else None),
        })

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
# In-Silico subsampling suite  (spike-in / dilution-series: expected vs reality)
# ──────────────────────────────────────────────────────────────────────────────

# Subsample dataset sample ids look like:
#   <parent>_insilico_<iss|nanosim>_ss_<mode>_c<count>_r<rep>   (synthetic)
#   <parent>_background_ss_<mode>_c<count>_r<rep>               (natural background)
_SS_ID_RE = re.compile(
    r'^(?P<parent>.+?)_'
    r'(?:insilico_(?P<isstok>iss|nanosim)|background)'
    r'_ss_(?P<mode>consistent|randomized)_c(?P<count>\d+)_r(?P<rep>\d+)$'
)


def _load_insilico_params(path):
    """Load the optional in-silico params JSON. Returns {} on any problem."""
    if not path:
        return {}
    try:
        with open(path) as fh:
            data = json.load(fh)
        return data if isinstance(data, dict) else {}
    except Exception as exc:
        print(f"[make_report] WARNING: could not read --insilico_params {path!r}: {exc}",
              file=sys.stderr)
        return {}


def _load_insilico_manifests(paths):
    """Parse *_subsample_manifest.tsv files → {dataset_id: {row fields}}."""
    out = {}
    for p in paths or []:
        try:
            with open(p) as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    did = (row.get("dataset_id") or "").strip()
                    if did:
                        out[did] = row
        except Exception as exc:
            print(f"[make_report] WARNING: could not read manifest {p!r}: {exc}",
                  file=sys.stderr)
    return out


def _f1(precision, recall):
    return (2 * precision * recall / (precision + recall)) if (precision + recall) > 0 else 0.0


def build_insilico_suite(rows, sample_meta, params_file=None, manifest_files=None,
                         detect_threshold=None):
    """
    Assemble the In-Silico subsampling suite payload: run parameters plus, for
    each (parent sample, platform) group, a per-dataset expected-vs-reality table
    and a per-organism dilution-series / limit-of-detection view.

    "Reality" comes from the subsample datasets already present in `rows` (each
    subsample flows through the pipeline as its own sample). An organism counts as
    *detected* when its TASS score clears `detect_threshold` (the same recommended
    cutoff the report uses — best_cutoffs.subkey.best_threshold — NOT the raw
    per-row passes_threshold flag, which is computed client-side and is unset in
    the JSON). "Expected" is the truth set of organisms detected at full depth; its
    composition is scaled to each dataset's target read count. Authoritative
    target/actual/total counts are taken from the subsample manifest when supplied.

    Rows are collapsed to a SINGLE taxonomic level (Species when available, else
    Strain) so the same reads are not triple-counted across strain/species/genus
    rollup rows.

    Returns None when no subsample datasets are present (feature off).
    """
    thr = float(detect_threshold) if detect_threshold not in (None, "") else 0.0

    # ── Identify subsample datasets from the sample-id pattern ────────────────
    datasets = {}   # sname -> {parent, plat, mode, count, rep}
    sample_names = set(r.get("Specimen ID") for r in rows)
    sample_names.update(sample_meta.keys())
    for sname in sample_names:
        if not sname:
            continue
        m = _SS_ID_RE.match(str(sname))
        if m:
            datasets[sname] = {
                "parent": m.group("parent"),
                "plat": m.group("isstok") or "background",
                "mode": m.group("mode"),
                "count": int(m.group("count")),
                "rep": int(m.group("rep")),
            }
    if not datasets:
        return None

    # Pairing (read pairs vs reads) is best taken from the sample's platform
    # metadata; fall back to the platform token (iss = paired Illumina).
    def _is_paired(sname, plat):
        p = str((sample_meta.get(sname) or {}).get("platform", "")).upper()
        if p:
            return p == "ILLUMINA"
        return plat == "iss"

    manifests = _load_insilico_manifests(manifest_files)

    # ── Choose ONE taxonomic level to avoid triple-counting rollup rows ───────
    subsample_rows = [r for r in rows if r.get("Specimen ID") in datasets]
    _levels = {r.get("Level") for r in subsample_rows}
    level_pick = "Species" if "Species" in _levels else ("Strain" if "Strain" in _levels else None)

    # ── Observed organisms per dataset (taxid -> observed metrics) ────────────
    # Subsampling selects records: for paired-end an index is a read PAIR, so the
    # target/actual counts are in pairs, while the aligner counts each mate
    # separately (≈ 2× reads). Normalise observed read counts back to records
    # (÷2 for paired) so they are directly comparable to the target read counts.
    obs = defaultdict(dict)
    for r in subsample_rows:
        if level_pick is not None and r.get("Level") != level_pick:
            continue
        sn = r.get("Specimen ID")
        tid = str(r.get("Taxonomic ID #", "") or "")
        if not tid:
            continue
        rpr = 2 if _is_paired(sn, datasets[sn]["plat"]) else 1
        tass = float(r.get("TASS Score", 0) or 0)
        obs[sn][tid] = {
            "name": r.get("Detected Organism", "Unknown"),
            "reads": int(round(int(r.get("# Reads Aligned", 0) or 0) / rpr)),
            "tass": tass,
            "passes": tass >= thr,   # detection vs the report's recommended TASS cutoff
            "category": r.get("Microbial Category", "Unknown"),
            # Lineage labels so the report can roll the suite up to Species or
            # Genus, and so a Strain-level detection can still be matched to a
            # series that was collapsed to Species.
            "species": str(r.get("Species Name", "") or ""),
            "genus": str(r.get("Genus Name", "") or ""),
        }

    # ── Group datasets by (parent, platform) ──────────────────────────────────
    groups = defaultdict(list)
    for sname, d in datasets.items():
        groups[(d["parent"], d["plat"])].append((sname, d))

    suite_groups = []
    all_counts = set()
    all_modes = set()
    max_rep = 1

    for (parent, plat), items in sorted(groups.items()):
        items.sort(key=lambda x: (x[1]["count"], x[1]["rep"]))
        mode = items[0][1]["mode"]
        all_modes.add(mode)

        # Truth set = organisms DETECTED at full depth (deepest dataset), i.e. TASS
        # clears the cutoff. Composition (expected read fraction) is taken from their
        # full-depth reads. Organisms present only as low-depth noise are excluded.
        deepest_count = max(d["count"] for _, d in items)
        deepest_snames = [sn for sn, d in items if d["count"] == deepest_count]
        exp_reads_by_tid = defaultdict(float)
        name_by_tid, cat_by_tid, sp_by_tid, gen_by_tid = {}, {}, {}, {}
        for sn in deepest_snames:
            for tid, ov in obs.get(sn, {}).items():
                if not ov["passes"]:
                    continue   # only organisms confidently detected at full depth
                exp_reads_by_tid[tid] += ov["reads"]
                name_by_tid[tid] = ov["name"]
                cat_by_tid[tid] = ov["category"]
                sp_by_tid[tid] = ov["species"]
                gen_by_tid[tid] = ov["genus"]
        # Fallback: if the cutoff excluded everything at full depth, use all organisms
        # present at full depth so the tab still shows the pool rather than nothing.
        if not exp_reads_by_tid:
            for sn in deepest_snames:
                for tid, ov in obs.get(sn, {}).items():
                    exp_reads_by_tid[tid] += ov["reads"]
                    name_by_tid[tid] = ov["name"]
                    cat_by_tid[tid] = ov["category"]
                    sp_by_tid[tid] = ov["species"]
                    gen_by_tid[tid] = ov["genus"]
        total_deep = sum(exp_reads_by_tid.values()) or 1.0
        expected_fraction = {tid: v / total_deep for tid, v in exp_reads_by_tid.items()}
        expected_set = set(expected_fraction)

        # Per-dataset (count × rep) expected-vs-reality rows.
        dataset_rows = []
        for sname, d in items:
            o = obs.get(sname, {})
            observed_total = sum(v["reads"] for v in o.values())
            detected_set = {tid for tid, v in o.items() if v["passes"]}
            tp = len(expected_set & detected_set)
            fp = len(detected_set - expected_set)
            fn = len(expected_set - detected_set)
            precision = tp / (tp + fp) if (tp + fp) else 0.0
            recall = tp / (tp + fn) if (tp + fn) else 0.0
            man = manifests.get(sname, {})
            dataset_rows.append({
                "id": sname,
                "replicate": d["rep"],
                "target_count": int(man.get("target_count") or d["count"]),
                "actual_count": int(man["actual_count"]) if man.get("actual_count") else d["count"],
                "total_master_reads": int(man["total_master_reads"]) if man.get("total_master_reads") else None,
                "seed": man.get("seed"),
                "observed_total_reads": observed_total,
                "n_detected": len(detected_set),
                "tp": tp, "fp": fp, "fn": fn,
                "precision": round(precision, 4),
                "recall": round(recall, 4),
                "f1": round(_f1(precision, recall), 4),
            })
            all_counts.add(d["count"])
            max_rep = max(max_rep, d["rep"])

        # Per-count aggregation across replicates (for the LoD chart).
        counts_sorted = sorted({d["count"] for _, d in items})
        by_count = defaultdict(list)
        for sname, d in items:
            by_count[d["count"]].append(sname)

        organisms = []
        for tid, frac in sorted(expected_fraction.items(), key=lambda kv: -kv[1]):
            series = []
            lod = None
            for c in counts_sorted:
                reps = by_count[c]
                obs_reads_vals, tass_vals, det_flags = [], [], []
                for sn in reps:
                    ov = obs.get(sn, {}).get(tid)
                    obs_reads_vals.append(ov["reads"] if ov else 0)
                    tass_vals.append(ov["tass"] if ov else 0.0)
                    det_flags.append(bool(ov and ov["passes"]))
                mean_obs = sum(obs_reads_vals) / len(obs_reads_vals)
                mean_tass = sum(tass_vals) / len(tass_vals)
                det_rate = sum(1 for f in det_flags if f) / len(det_flags)
                detected = det_rate >= 0.5
                series.append({
                    "count": c,
                    "expected_reads": round(frac * c, 1),
                    "observed_reads": round(mean_obs, 1),
                    "tass": round(mean_tass, 2),
                    "detection_rate": round(det_rate, 3),
                    "detected": detected,
                    "n_reps": len(reps),
                })
                if detected and lod is None:
                    lod = c
            organisms.append({
                "taxid": tid,
                "name": name_by_tid.get(tid, "Unknown"),
                "category": cat_by_tid.get(tid, "Unknown"),
                # Lineage for the report's Species/Genus rollups and for matching
                # a Strain-level detection against a Species-level series.
                "species": sp_by_tid.get(tid, ""),
                "genus": gen_by_tid.get(tid, ""),
                "expected_fraction": round(frac, 6),
                "lod_count": lod,
                "series": series,
            })

        # Group pairing: use the first dataset's platform metadata.
        _grp_paired = _is_paired(items[0][0], plat)
        suite_groups.append({
            "parent": parent,
            "platform": plat,
            # The single level the series rows were collapsed to (Species when
            # available, else Strain). The report labels its plots with this and
            # uses it to explain Strain->Species matches.
            "level": level_pick or "Strain",
            "source": "background" if plat == "background" else "insilico",
            "mode": mode,
            "n_datasets": len(items),
            "counts": counts_sorted,
            # Unit that target/observed counts are expressed in. Paired-end
            # subsamples read pairs; observed reads are normalised to pairs above.
            "read_unit": "read pairs" if _grp_paired else "reads",
            "datasets": dataset_rows,
            "organisms": organisms,
        })

    # ── Parameters (explicit file overrides inferred) ─────────────────────────
    inferred = {
        "mode": "/".join(sorted(all_modes)) if all_modes else None,
        "series_counts": sorted(all_counts),
        "replicates": max_rep,
        "detection_threshold": round(thr, 2),
    }
    params = dict(inferred)
    params.update({k: v for k, v in _load_insilico_params(params_file).items() if v is not None})

    return {
        "enabled": True,
        "params": params,
        "groups": suite_groups,
    }


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

    # All organisms with TASS > 0 are always included; the UI pre-populates its
    # filter slider from best_cutoffs baked into the BOOT payload at load time.
    mintass = args.mintass
    print("[make_report] TASS threshold: 0.0 (all organisms included; UI filter set from best_cutoffs)")

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

    # ── optional sample metadata file ─────────────────────────────────────────
    # Merged into sample_meta BEFORE run_metadata_records is derived below, so
    # the extra columns flow into the report's metadata table, map and specimen
    # grouping exactly as they would had the file been dropped on the tab.
    if args.metadata:
        # Sample ids as they appear in the loaded detections — the authoritative
        # list of what this run actually contains.
        _run_samples = {
            str(r.get("Specimen ID")) for r in rows if r.get("Specimen ID") not in (None, "")
        }
        meta_summary = load_sample_metadata(
            args.metadata, sample_meta,
            known_samples=_run_samples,
            add_unmatched=args.metadata_add_unmatched,
            sheet=args.metadata_sheet,
        )
        cols = sorted(meta_summary["columns"])
        print(f"[make_report] Metadata: {meta_summary['rows']} row(s) from "
              f"{meta_summary['files']} file(s); {meta_summary['matched']} matched this run"
              + (f", {meta_summary['added']} added" if meta_summary["added"] else ""))
        if cols:
            print(f"[make_report] Metadata columns: {', '.join(cols)}")
        if meta_summary["specimens"]:
            print(f"[make_report] Specimen grouping: {meta_summary['specimens']} "
                  f"multi-sample specimen(s) — merge will be available in the report")
        elif any(c.startswith("specimen") for c in cols):
            print("[make_report] WARNING: a specimen column was read but no two samples "
                  "share a specimen value, so nothing will merge.", file=sys.stderr)
        if meta_summary["unmatched"]:
            _u = meta_summary["unmatched"]
            print(f"[make_report] WARNING: {len(_u)} metadata row(s) did not match any sample "
                  f"in this run and were skipped: {', '.join(_u[:5])}"
                  + (" …" if len(_u) > 5 else ""), file=sys.stderr)
            if not meta_summary["matched"]:
                print("[make_report] WARNING: NO metadata row matched this run — check that the "
                      "'sample' column uses the same ids as the pipeline.", file=sys.stderr)

    # ── derive column lists ────────────────────────────────────────────────────
    # These fields are carried on each record for client-side analysis (the
    # Feature Compare view + capability detection) but are NOT human-displayable
    # table columns — 'High ANI Matches' is a nested list — so keep them out of
    # the column picker / detections table.
    _NON_DISPLAY_COLS = {"High ANI Matches", "ANI Annotated"}
    all_cols = [c for c in (rows[0].keys() if rows else []) if c not in _NON_DISPLAY_COLS]
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
    # NOTE: "sample_overview" is the TASS organism table, NOT VF/AMR annotation —
    # create_report.py writes that sheet unconditionally, even when annotation was
    # disabled or matched nothing. Only the annotation-bearing sheets may be used
    # to decide "does this run have VF/AMR data?" or "is this sample covered?";
    # counting sample_overview made has_prot True for every run (rendering an empty
    # VF/AMR tab) and marked every sample as already covered (silently discarding
    # the --annotate_reports fallback below).
    prot_data = load_protein_annotations(args.protein_annotations, pident=args.pident)
    has_prot = any(len(prot_data.get(k, [])) > 0 for k in PROT_ANNOTATION_KEYS)
    print(f"[make_report] Protein annotations loaded: {has_prot} "
          f"({sum(len(v) for v in prot_data.values())} total rows; "
          f"{sum(len(prot_data.get(k, [])) for k in PROT_ANNOTATION_KEYS)} annotation rows)")

    # ── standalone annotation reports (de-novo / unaligned samples) ────────────
    # Samples with no reference alignment never get an organism hierarchy, so their
    # VF/AMR annotation is absent from the merged XLSX. Supplement prot_data from the
    # per-sample annotate_report.tsv files for any sample not already covered.
    if args.annotate_reports:
        # A sample counts as "covered" only if the merged XLSX carried actual VF/AMR
        # rows for it. Appearing in sample_overview means nothing — every sample in
        # the run is listed there regardless of whether annotation matched.
        covered = set()
        for _k in PROT_HIT_KEYS:
            for _r in prot_data.get(_k, []):
                _s = _r.get("Specimen ID") or _r.get("Sample") or _r.get("sample")
                if _s not in (None, ""):
                    covered.add(str(_s))
        supp = load_standalone_annotations(args.annotate_reports, covered, pident=args.pident)
        _n_supp = sum(len(v) for v in supp.values())
        if _n_supp:
            for _k in prot_data:
                prot_data[_k].extend(supp.get(_k, []))
            has_prot = any(len(prot_data.get(k, [])) > 0 for k in PROT_ANNOTATION_KEYS)
            _supp_samples = sorted({
                (r.get("Specimen ID") or r.get("Sample") or r.get("sample"))
                for k in ("per_gene_hits", "amr_genes") for r in supp.get(k, [])
            } - {None, ""})
            print(f"[make_report] Supplemented annotation for {len(_supp_samples)} unaligned "
                  f"sample(s) from standalone reports: {_supp_samples} ({_n_supp} rows)")
        else:
            print("[make_report] No standalone annotation rows to supplement "
                  "(all samples already covered or files empty)")

    # ── novelty detection (reference-free LCA) ────────────────────────────────
    novelty_data, novelty_downloads = load_novelty(args.novelty, args.novelty_downloads)
    has_novelty = bool(novelty_data.get("samples"))
    print(f"[make_report] Novelty loaded: {has_novelty} "
          f"({len(novelty_data.get('samples', {}))} sample(s), "
          f"{len(novelty_downloads)} download link(s))")

    # ── pathogen reference cross-reference ────────────────────────────────────
    # Load the pathogen sheet, pre-flag novelty candidates that are listed
    # pathogens, and ship the lookups so the report can also flag VF/AMR hits
    # (genus/species) that have no reference alignment.
    pathogens = load_pathogens(args.pathogens)
    has_pathogens = bool(pathogens["by_taxid"] or pathogens["by_name"] or pathogens["by_genus"])
    n_flagged = annotate_novelty_pathogens(novelty_data, pathogens) if has_pathogens else 0
    # Resolve VF/AMR taxids from the bvbrc reference (Source ID -> taxids) and stamp
    # each VF/AMR row with its pathogen match by canonical taxid.
    source_taxids = load_source_taxids(args.vfamr_taxids) if has_pathogens else {}
    n_prot_flagged = annotate_protein_pathogens(prot_data, source_taxids, pathogens) if has_pathogens else 0
    print(f"[make_report] Pathogen sheet loaded: {has_pathogens} "
          f"({len(pathogens['by_taxid'])} taxid / {len(pathogens['by_genus'])} genus entries; "
          f"{n_flagged} novelty candidate(s), {n_prot_flagged} VF/AMR row(s) flagged; "
          f"{len(source_taxids)} bvbrc source-id taxids)")

    # ── optional sample rename map ────────────────────────────────────────────
    # Applied last, after every data source (records, metadata, contig data,
    # protein/VF-AMR annotations, novelty) has been loaded and merged, so a
    # single rename touches all of them consistently before the HTML is built.
    if args.rename_samples:
        rename_map = load_sample_rename_map(args.rename_samples)
        if rename_map:
            rn_summary = apply_sample_renames(
                rename_map, rows, sample_meta, contig_data, prot_data, novelty_data
            )
            print(f"[make_report] Sample rename: {len(rename_map)} mapping(s) loaded from "
                  f"{args.rename_samples!r}; {rn_summary['applied']} applied")
            if rn_summary["unmatched"]:
                _u = rn_summary["unmatched"]
                print(f"[make_report] WARNING: {len(_u)} rename entry(ies) matched no sample "
                      f"in this run: {', '.join(_u[:5])}" + (" …" if len(_u) > 5 else ""),
                      file=sys.stderr)
        else:
            print(f"[make_report] WARNING: --rename-samples {args.rename_samples!r} produced no "
                  f"usable mappings", file=sys.stderr)

    # ── extract run-level metadata for map / metadata panel ──────────────────
    # These fields are part of the fixed pipeline/sample identity — not run metadata.
    _META_PIPELINE_KEYS = {
        "sample_name", "sample_type", "platform", "workflow_revision", "commit_id",
        "total_reads", "aligned_reads", "total_organism_reads", "num_species_groups",
        "num_keys", "num_subkeys", "num_toplevelkeys", "minmapq", "mapq_breadth_power",
        "weights", "min_conf_applied", "best_cutoffs", "best_cutoffs_source",
        "preferred_granularity", "control_type", "negative_controls_used",
        "positive_controls_used", "control_fold_threshold", "missing_positive_controls",
        "insilico_controls_used", "insilico_simulator_types", "missing_insilico_controls",
        "missing_insilico_by_type",
    }
    run_metadata_records = []
    for sname, smeta in sample_meta.items():
        rec = {"sample_name": sname}
        # Include ALL metadata fields that are not fixed pipeline/identity fields
        for k, v in smeta.items():
            if k not in _META_PIPELINE_KEYS and k != "sample_name":
                rec[k] = v
        # Only add to records if at least one field has a non-null, non-empty value
        if any(v is not None and v != "" for k, v in rec.items() if k != "sample_name"):
            run_metadata_records.append(rec)

    has_geo = any(
        r.get("latitude") is not None and r.get("longitude") is not None
        for r in run_metadata_records
    )

    # ── derive global pipeline revision / commit ─────────────────────────────
    # Mirrors the ODR PDF logic: pick first non-empty value; "local" is valid
    # for commit_id; "NA" is treated as absent for both fields.
    def _first_val(vals, exclude_upper=("NA", "NULL", "NONE", "")):
        for v in vals:
            sv = str(v).strip() if v is not None else ""
            if sv.upper() not in exclude_upper:
                return sv
        return None

    _all_revisions = [m.get("workflow_revision") for m in sample_meta.values()]
    _all_commits   = [m.get("commit_id")         for m in sample_meta.values()]

    # revision: None when all NA → display as "Not Specified or Local Build"
    pipeline_revision = _first_val(_all_revisions)
    # commit: "local" is a valid value (not filtered); None only if all NA/empty
    pipeline_commit   = _first_val(_all_commits)

    # When commit is None or "local", try git to get the actual hash
    if not pipeline_commit or pipeline_commit.lower() == "local":
        try:
            import subprocess as _sp
            _git_hash = _sp.check_output(
                ["git", "rev-parse", "HEAD"], stderr=_sp.DEVNULL, cwd=os.path.dirname(__file__)
            ).decode().strip()
            if _git_hash:
                pipeline_commit = _git_hash
        except Exception:
            pass
        if not pipeline_commit:
            pipeline_commit = "local"

    print(f"[make_report] Pipeline revision: {pipeline_revision or 'Not Specified or Local Build'}  commit: {pipeline_commit}")

    # ── explicit --min_conf override ──────────────────────────────────────────
    # If the user passed --min_conf, stamp it into every sample's best_cutoffs
    # (all granularities) so BOTH the global filter slider (which reads the
    # aggregated best_cutoffs below) AND the per-sample-type sliders (which read
    # each sample's own best_cutoffs directly, see _defaultTassForType in
    # 03_sample_sidebar.js) default to this value instead of the auto-computed
    # thresholds-JSON recommendation.
    if args.min_conf is not None:
        for smeta in sample_meta.values():
            smeta["best_cutoffs"] = {
                level: {"best_threshold": args.min_conf}
                for level in ("key", "subkey", "toplevelkey")
            }
            smeta["best_cutoffs_source"] = "user-specified (--min_conf)"
            smeta["min_conf_applied"] = args.min_conf
        print(f"[make_report] --min_conf={args.min_conf} overriding all sample/type slider defaults")

    # ── in-silico subsampling suite (spike-in / dilution-series) ──────────────
    # The subsample datasets are kept OUT of the main heatmap/table (report.nf feeds
    # only non-control samples to -i). Load their JSONs here on a SEPARATE track so
    # the suite tab can be built without adding them to the cross-sample records.
    insil_rows, insil_meta = [], {}
    _insil_json = [f for f in (args.insilico_json or [])
                   if f and f.strip().endswith(".json") and not os.path.basename(f).startswith("NO_FILE")]
    if _insil_json:
        try:
            insil_rows, insil_meta, _ = load_json_inputs(
                _insil_json, mintass=mintass, microbial_cats=microbial_cats
            )
            print(f"[make_report] Loaded {len(insil_rows)} in-silico subsample organism rows "
                  f"from {len(_insil_json)} dataset JSON(s)")
        except Exception as exc:
            print(f"[make_report] WARNING: could not load --insilico_json inputs: {exc}",
                  file=sys.stderr)
    # Union with main rows is harmless — the suite only picks rows whose sample id
    # matches the subsample pattern, so non-subsample rows are ignored. This also
    # lets the suite work if a user passes subsample JSONs straight to -i.
    _suite_rows = rows + insil_rows
    _suite_meta = dict(sample_meta); _suite_meta.update(insil_meta)
    # Detection cutoff for the suite: use the SAME recommended TASS threshold the
    # report defaults to (best_cutoffs.subkey → key), so "detected" here matches
    # what the user sees elsewhere. Fall back to the --mintass hard filter.
    _suite_bc = _collect_best_cutoffs(_suite_meta) or {}
    _suite_thr = ((_suite_bc.get("subkey") or {}).get("best_threshold")
                  or (_suite_bc.get("key") or {}).get("best_threshold"))
    if _suite_thr is None:
        _suite_thr = mintass if (mintass and mintass > 1) else 0.0
    insilico_suite = build_insilico_suite(
        _suite_rows, _suite_meta,
        params_file=args.insilico_params,
        manifest_files=args.insilico_manifests,
        detect_threshold=_suite_thr,
    )
    if insilico_suite:
        _ng = len(insilico_suite["groups"])
        _nd = sum(len(g["datasets"]) for g in insilico_suite["groups"])
        print(f"[make_report] In-Silico suite: {_ng} group(s), {_nd} subsample dataset(s)")
    else:
        print("[make_report] In-Silico suite: no subsample datasets detected (tab hidden)")

    # ── collect best_cutoffs for UI pre-population ────────────────────────────
    best_cutoffs_payload = _collect_best_cutoffs(sample_meta)
    if best_cutoffs_payload:
        thresh = (best_cutoffs_payload.get("subkey") or {}).get("best_threshold")
        print(f"[make_report] UI will pre-set TASS filter to: {thresh} (from best_cutoffs.subkey)")
    else:
        print("[make_report] No best_cutoffs found; UI TASS filter will default to 0")

    # ── sample QC flag defaults ───────────────────────────────────────────────
    sample_flags = build_sample_flag_config(args)
    if sample_flags:
        _acts = {r["action"] for r in sample_flags["rules"]}
        print(f"[make_report] Sample QC: {len(sample_flags['rules'])} default rule(s), "
              f"logic={sample_flags['logic']}, action(s)={'/'.join(sorted(_acts))}")
        for _r in sample_flags["rules"]:
            _t = f" (TASS {_r['tass']:g})" if _r.get("tass") is not None else ""
            _a = f" agg={_r['agg']}" if _r.get("agg") else ""
            print(f"[make_report]   - {_r['source']}.{_r['field']}{_t}{_a} "
                  f"{_r['op']} {_r['value']} -> {_r['action']}")
    else:
        print("[make_report] Sample QC: no default rules (no --flag-* criteria given); "
              "the report's Sample QC panel will start empty.")

    # ── build bootstrap payload ───────────────────────────────────────────────
    payload = _sanitize({
        "records":               rows,
        "all_cols":              all_cols,
        "numeric_cols":          numeric_cols,
        "sample_meta":           sample_meta,
        "prot_data":             prot_data,
        "has_prot":              has_prot,
        "contig_data":           list(contig_data.values()),   # list of organism contig objects
        "best_cutoffs":          best_cutoffs_payload,         # for UI filter pre-population
        "run_metadata_records":  run_metadata_records,         # per-sample run metadata
        "has_geo":               has_geo,                      # true if any sample has lat/lon
        "novelty":               novelty_data,                 # {samples: {<s>: {summary, candidates}}}
        "has_novelty":           has_novelty,                  # true if any novelty sample present
        "novelty_downloads":     novelty_downloads,            # [{label, kind, filename}] for links
        "pathogens":             pathogens,                    # {by_taxid, by_name, by_genus} pathogen lookups
        "has_pathogens":         has_pathogens,                # true if a pathogen sheet was loaded
        "sample_flags":          sample_flags,                 # default whole-sample QC rules (None when unconfigured)
        "report_generated_at":   datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "pipeline_revision":     pipeline_revision,            # global branch/tag or "local"
        "pipeline_commit":       pipeline_commit,              # global commit hash or None
        "insilico_suite":        insilico_suite,               # spike-in/dilution suite or None
        "has_insilico_suite":    bool(insilico_suite),         # true when subsample datasets present
    })

    bootstrap_json = json.dumps(payload, ensure_ascii=False, allow_nan=False, separators=(',', ':'))

    # Build inline JS instead of writing a separate heatmap_boot.js file
    bootstrap_script = f"<script>\nwindow.HEATMAP_BOOT = {bootstrap_json};\n</script>"

    # ── render template ───────────────────────────────────────────────────────
    # heatmap.html is a thin shell that references its CSS/JS as external parts
    # under assets/src/. inline_template() folds those back in so the report is
    # one self-contained file. (See bin/report_template.py.)
    # offline_report / offline_report_files fold the CDN libraries inline so the
    # report opens with no internet access. Default leaves them as CDN links.
    from report_template import inline_template
    tpl = inline_template(
        args.template,
        offline=args.offline_report,
        offline_dir=args.offline_report_files,
    )

    # Use a function replacement (not a plain string) so backslashes / `\g`-like
    # sequences inside the JSON payload are inserted verbatim and never treated
    # as regex backreferences.
    _repl = lambda m: bootstrap_script

    # ── Strip the pages.js demo loader (active or commented) ──────────────────
    # The real dataset is injected inline below. pages.js is not copied next to
    # the report, so a leftover tag 404s in the browser console.
    html = re.sub(
        r'[ \t]*<!--\s*<script[^>]+src=["\']pages\.js["\'][^>]*>\s*</script>\s*-->[ \t]*\n?',
        '', tpl, flags=re.IGNORECASE,
    )
    html = re.sub(
        r'[ \t]*<script[^>]+src=["\']pages\.js["\'][^>]*>\s*</script>[ \t]*\n?',
        '', html, flags=re.IGNORECASE,
    )

    # ── Inject the inline bootstrap at the heatmap_boot.js anchor ─────────────
    # Handle the commented placeholder, an active tag, or an existing BOOTSTRAP
    # block. A *commented* anchor MUST be matched together with its <!-- … -->
    # so the boot script is not left commented out — that would blank the report.
    anchor_patterns = [
        r'<!--\s*<script[^>]+src=["\']heatmap_boot\.js["\'][^>]*>\s*</script>\s*-->',
        r'<script[^>]+src=["\']heatmap_boot\.js["\'][^>]*>\s*</script>',
        r'<!--\s*<script id=["\']BOOTSTRAP["\'][^>]*>.*?</script>\s*-->',
        r'<script id=["\']BOOTSTRAP["\'][^>]*>.*?</script>',
    ]
    n_replaced = 0
    for pat in anchor_patterns:
        html, n_replaced = re.subn(pat, _repl, html, count=1, flags=re.DOTALL | re.IGNORECASE)
        if n_replaced:
            break

    if not n_replaced and "__BOOTSTRAP_SCRIPT__" in html:
        # Fallback: explicit placeholder token.
        html = html.replace("__BOOTSTRAP_SCRIPT__", bootstrap_script)
        n_replaced = 1

    if not n_replaced:
        raise SystemExit(
            "[make_report] ERROR: no heatmap_boot.js / BOOTSTRAP anchor found in template"
        )

    # ── Guard: the boot script must not have landed inside an HTML comment ────
    _idx = html.find("window.HEATMAP_BOOT =")
    if _idx != -1:
        _open = html.rfind("<!--", 0, _idx)
        _close = html.rfind("-->", 0, _idx)
        if _open != -1 and _open > _close:
            raise SystemExit(
                "[make_report] ERROR: bootstrap script was inserted inside an HTML "
                "comment — the report would be blank. Check the template anchor."
            )

    with open(args.output, "w", encoding="utf-8") as fh:
        fh.write(html)

    print(f"[make_report] Written: {args.output}")

    # ── optional combined XLSX export ─────────────────────────────────────────
    # Separate from --protein_annotations (-p), which is XLSX *input*. This is
    # the only flag that makes make_report.py *write* an XLSX file.
    if args.output_xlsx:
        def _sheet_df(records, cols=None):
            if not records:
                return pd.DataFrame(columns=cols or [])
            df = pd.DataFrame(records)
            if cols:
                # Keep report column order; include any extra keys at the end.
                extra = [c for c in df.columns if c not in cols]
                df = df.reindex(columns=list(cols) + extra)
            return df

        _xlsx_sheets = [("Detections", _sheet_df(rows, all_cols))]
        if run_metadata_records:
            _xlsx_sheets.append(("Metadata", _sheet_df(run_metadata_records)))
        if has_prot:
            for _key, _sheet_name in (
                ("genus_summary",   "Genus Summary"),
                ("per_gene_hits",   "Per-Gene Hits"),
                ("sample_overview", "Sample Overview"),
                ("amr_genes",       "AMR Genes"),
            ):
                _recs = prot_data.get(_key, [])
                if _recs:
                    _xlsx_sheets.append((_sheet_name, _sheet_df(_recs)))

        with pd.ExcelWriter(args.output_xlsx, engine="openpyxl") as writer:
            for _sheet_name, _df in _xlsx_sheets:
                # Excel sheet names are capped at 31 characters.
                _df.to_excel(writer, sheet_name=_sheet_name[:31], index=False)

        print(f"[make_report] Written: {args.output_xlsx} "
              f"({len(_xlsx_sheets)} sheet(s): {', '.join(n for n, _ in _xlsx_sheets)})")


def _is_numeric_str(s):
    try:
        float(s)
        return True
    except (ValueError, TypeError):
        return False


if __name__ == "__main__":
    main()
