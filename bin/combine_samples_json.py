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

Merging finished runs
---------------------
Inputs may also be already-combined files (report/all.odr.json from earlier
runs). Their samples are expanded (never nested) and their embedded
"novelty" / "prot_data" blocks are carried over, so several runs can be merged
into one file:

    combine_samples_json.py -i runA/report/all.odr.json runB/report/all.odr.json \
        -o merged.odr.json

Sample names must be unique (the report keys samples by name). When two inputs
share a name, --on-duplicate decides: rename (default; "<name>__<label>",
label = --labels entry or the run folder name), skip, keep, or error. Renames
are applied everywhere the name appears (metadata, organism rows, control
links, VF/AMR "Specimen ID", novelty). --tag-all applies the label to every
sample. A top-level "merged_from" list records each input.

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

from microbial_category import canonicalize_categories

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
    parser.add_argument(
        "--on-duplicate", choices=("rename", "skip", "keep", "error"), default="rename",
        help="What to do when a sample name was already loaded from an earlier input: "
             "rename it to <name>__<label> (default), skip it, keep both unchanged "
             "(the report will collapse them), or stop with an error.",
    )
    parser.add_argument(
        "--labels", nargs="+", default=None, metavar="LABEL",
        help="One label per -i argument (in order) used for renaming. Default: the run "
             "folder name (the parent of report/ for .../report/all.odr.json, otherwise "
             "the file's folder name).",
    )
    parser.add_argument(
        "--tag-all", action="store_true",
        help="Rename every sample to <name>__<label>, not only duplicates.",
    )
    return parser.parse_args(argv)


# ── multi-run merging helpers ───────────────────────────────────────────────
# Keys whose string value (or list of strings) names a sample. Renames are only
# applied there, so an organism that happens to share a sample's name is safe.
_SAMPLE_REF_KEYS = {"sample_name", "sample_id", "parent_id", "Specimen ID", "sample"}


def _is_sample_ref_key(key):
    return key in _SAMPLE_REF_KEYS or key.endswith("_controls_used") or key.endswith("_controls")


def _rename_refs(obj, mapping):
    """Rename sample references in place (see _SAMPLE_REF_KEYS)."""
    if not mapping:
        return obj
    if isinstance(obj, dict):
        for k, v in obj.items():
            if _is_sample_ref_key(k):
                if isinstance(v, str) and v in mapping:
                    obj[k] = mapping[v]
                    continue
                if isinstance(v, list) and any(isinstance(x, str) for x in v):
                    obj[k] = [mapping.get(x, x) if isinstance(x, str) else x for x in v]
                    continue
            _rename_refs(v, mapping)
    elif isinstance(obj, list):
        for v in obj:
            _rename_refs(v, mapping)
    return obj


def _default_label(path):
    """Run folder name: <run>/report/all.odr.json -> <run>."""
    d = os.path.dirname(os.path.abspath(path))
    if os.path.basename(d) == "report":
        d = os.path.dirname(d)
    return os.path.basename(d) or os.path.basename(path).split(".")[0]


def _sample_name(data, path):
    return (data.get("metadata") or {}).get("sample_name") or os.path.basename(path).split(".")[0]


def main():
    args = parse_args()

    if args.labels is not None and len(args.labels) != len(args.input):
        print(f"[combine_samples_json] ERROR: --labels has {len(args.labels)} value(s) "
              f"but -i has {len(args.input)}.", file=sys.stderr)
        sys.exit(1)

    # Expand globs; each file remembers the label of the -i argument it came from
    paths = []
    for idx, pattern in enumerate(args.input):
        pattern = pattern.strip()
        hits = [pattern] if os.path.isfile(pattern) else sorted(glob.glob(pattern))
        if not hits:
            print(f"[combine_samples_json] WARNING: no files found for {pattern!r}, skipping",
                  file=sys.stderr)
        for h in hits:
            label = args.labels[idx] if args.labels else _default_label(h)
            paths.append((h, label))

    if not paths:
        print("[combine_samples_json] ERROR: no input files found.", file=sys.stderr)
        sys.exit(1)

    samples = []
    seen = set()
    merged_from = []
    carried_prot = {k: [] for k in PROT_KEYS}
    carried_nov = None
    n_combined_inputs = 0

    for path, label in paths:
        try:
            with open(path, "r", encoding="utf-8") as fh:
                data = json.load(fh)
        except Exception as exc:
            print(f"[combine_samples_json] WARNING: failed to parse {path!r}: {exc}",
                  file=sys.stderr)
            continue

        is_combined = bool(data.get("taxtriage_combined"))
        n_combined_inputs += is_combined
        in_samples = data.get("samples", []) if is_combined else [data]

        # Decide each sample's final name for this input
        mapping, keep_idx, skipped = {}, [], []
        local = set()
        for j, sd in enumerate(in_samples):
            name = _sample_name(sd, path)
            new = name
            if args.tag_all or name in seen or name in local:
                if name in seen and not args.tag_all:
                    if args.on_duplicate == "error":
                        print(f"[combine_samples_json] ERROR: sample {name!r} in {path!r} "
                              f"was already loaded (use --on-duplicate rename|skip|keep).",
                              file=sys.stderr)
                        sys.exit(1)
                    if args.on_duplicate == "skip":
                        skipped.append(name)
                        continue
                    if args.on_duplicate == "keep":
                        print(f"[combine_samples_json] WARNING: duplicate sample {name!r} "
                              f"kept as-is; the report will treat them as one sample.",
                              file=sys.stderr)
                        keep_idx.append(j)
                        local.add(name)
                        continue
                new = f"{name}__{label}"
                k = 2
                while new in seen or new in local:
                    new = f"{name}__{label}_{k}"
                    k += 1
            if new != name:
                mapping[name] = new
            local.add(new)
            keep_idx.append(j)

        for j in keep_idx:
            sd = in_samples[j]
            old = _sample_name(sd, path)
            _rename_refs(sd, mapping)
            if old in mapping:
                sd.setdefault("metadata", {})["sample_name"] = mapping[old]
                sd["metadata"]["original_sample_name"] = old
                sd["metadata"]["merged_from_run"] = label
            samples.append(sd)
        seen |= local

        # Carry over blocks embedded in an already-combined input
        if is_combined:
            emb = data.get("prot_data")
            if isinstance(emb, dict):
                for k in PROT_KEYS:
                    for r in emb.get(k) or []:
                        if not isinstance(r, dict):
                            continue
                        if r.get("Specimen ID") in skipped:
                            continue
                        carried_prot[k].append(_rename_refs(r, mapping))
            nov = data.get("novelty")
            if isinstance(nov, dict) and isinstance(nov.get("samples"), dict):
                if carried_nov is None:
                    carried_nov = {"samples": {}, "classifier": nov.get("classifier", "") or "",
                                   "gene_mode": nov.get("gene_mode")}
                for sname, sv in nov["samples"].items():
                    if sname in skipped:
                        continue
                    carried_nov["samples"][mapping.get(sname, sname)] = _rename_refs(sv, mapping)

        merged_from.append({
            "file": os.path.abspath(path),
            "label": label,
            "samples": len(keep_idx),
            "renamed": mapping,
            "skipped": skipped,
        })
        what = f"combined file ({len(in_samples)} sample(s))" if is_combined else "sample"
        print(f"[combine_samples_json] Loaded {path!r} [{label}] {what}: kept {len(keep_idx)}"
              + (f", renamed {mapping}" if mapping else "")
              + (f", skipped {skipped}" if skipped else ""))

    combined = {
        "taxtriage_combined": True,
        "version": "1.1",
        "samples": samples,
    }
    # Provenance only when merging finished runs, so the pipeline's own
    # all.odr.json (built from per-sample .paths.json) is unchanged.
    if n_combined_inputs:
        combined["merged_from"] = merged_from

    # ── optional novelty block (carried over + --novelty) ─────────────────
    novelty = load_novelty(args.novelty)
    if carried_nov and carried_nov["samples"]:
        if novelty:
            carried_nov["samples"].update(novelty["samples"])
            carried_nov["classifier"] = novelty.get("classifier") or carried_nov["classifier"]
        novelty = carried_nov
    if novelty:
        combined["novelty"] = novelty
        combined["has_novelty"] = True
        print(f"[combine_samples_json] Novelty embedded: {len(novelty['samples'])} sample(s)")

    # ── optional VF/AMR block (carried over + --annotate-reports) ─────────
    prot_data = load_annotate_reports(args.annotate_reports, pident=args.pident)
    for k in PROT_KEYS:
        prot_data[k] = carried_prot[k] + prot_data[k]
    n_prot = sum(len(v) for v in prot_data.values())
    if n_prot:
        combined["prot_data"] = prot_data
        combined["has_prot"] = True
        print(f"[combine_samples_json] Annotations embedded: {n_prot} row(s) "
              f"({len(prot_data['per_gene_hits'])} VF, {len(prot_data['amr_genes'])} AMR)")

    # Blank, unrecognised or legacy ("Unknown") categories -> "Uncategorized".
    canonicalize_categories(combined)
    with open(args.output, "w", encoding="utf-8") as fh:
        json.dump(combined, fh, ensure_ascii=False, allow_nan=False, separators=(",", ":"))

    size_mb = os.path.getsize(args.output) / (1024 * 1024)
    print(f"[combine_samples_json] Written {len(samples)} sample(s) -> {args.output!r} "
          f"({size_mb:.1f} MB)")


if __name__ == "__main__":
    main()
