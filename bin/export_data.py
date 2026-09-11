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
export_data.py
==============
Combined, multi-tab data export for TaxTriage — the pipeline-side twin of the
report's "Export Data" panel (assets/src/js/47_export_data.js).

The interactive report can export ONE rendered table at a time. This module
builds the SAME catalog of datasets the report's Export Data panel offers —
detections, per-sample summary, cross-sample organism rollup, coverage, VF/AMR,
novelty, run metadata, geography, in-silico — and writes them out as a
multi-sheet workbook, per-dataset CSV/TSVs, and/or a single wide table joined on
Specimen ID x Organism. No browser required.

KEEP IN SYNC with assets/src/js/47_export_data.js: the dataset ids and the column
headers are the contract between the two, so a spreadsheet produced by the
pipeline lines up with one a user exported by hand from the report.

Two entry points:

  library   from export_data import write_exports
            write_exports(payload, outdir, formats=..., datasets=...)
            Called by bin/make_report.py when --export_data is given, using the
            same bootstrap payload that is embedded in the HTML.

  CLI       export_data.py -i all.odr.html -o export/
            export_data.py -i boot.json    -o export/ --formats xlsx,csv --datasets detections,organism_summary
            Re-exports from any already-built report (the payload is read back
            out of window.HEATMAP_BOOT), so an old run can be exported without
            re-running the pipeline.
"""

import argparse
import csv
import datetime
import json
import math
import os
import re
import sys
from collections import OrderedDict, defaultdict

# ──────────────────────────────────────────────────────────────────────────────
# Small helpers (mirrors of the JS side)
# ──────────────────────────────────────────────────────────────────────────────

# Join shapes for the wide table. Mirrors TT_JOIN in 47_export_data.js.
JOIN_NONE = None
JOIN_SAMPLE = "sample"
JOIN_ORG = "organism"
JOIN_SAMPLE_ORG = "sample+organism"


def _num(v):
    """float(v) or None — never raises, and never lets NaN/inf through."""
    if v is None or isinstance(v, bool):
        return None
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return f if math.isfinite(f) else None


def _round(v, dp=3):
    f = _num(v)
    if f is None:
        return ""
    return round(f, dp)


def _yn(v):
    """Booleans become Yes/No so a spreadsheet cell reads the same as the report."""
    if isinstance(v, bool):
        return "Yes" if v else "No"
    return v


def _truthy(v):
    if isinstance(v, bool):
        return v
    if v is None:
        return False
    s = str(v).strip().lower()
    return s in ("1", "true", "yes", "y", "t")


def _flatten_scalars(obj, prefix=""):
    """Scalar-only view of a metadata dict.

    Nested objects are dropped (no cell representation); lists of primitives are
    joined with "; ". Mirrors _ttScalars() in the JS module so the Run Metadata
    and Pipeline Sample Metadata sheets carry the same columns either way.
    """
    out = OrderedDict()
    for k, v in (obj or {}).items():
        if v is None:
            continue
        key = prefix + str(k)
        if isinstance(v, (list, tuple)):
            if v and all(not isinstance(x, (dict, list, tuple)) for x in v):
                out[key] = "; ".join("" if x is None else str(x) for x in v)
            continue
        if isinstance(v, dict):
            continue
        out[key] = _yn(v)
    return out


def _key_union(rows):
    """Union of keys across rows, in first-seen order."""
    cols, seen = [], set()
    for r in rows or []:
        for k in r:
            if k not in seen:
                seen.add(k)
                cols.append(k)
    return cols


class Table(object):
    """One exportable dataset: ordered columns plus dict rows."""

    __slots__ = ("columns", "rows")

    def __init__(self, columns, rows):
        self.columns = list(columns) if columns else _key_union(rows)
        self.rows = rows or []

    def __len__(self):
        return len(self.rows)

    def aoa(self):
        out = [list(self.columns)]
        for r in self.rows:
            out.append([("" if r.get(c) is None else r.get(c)) for c in self.columns])
        return out


# ──────────────────────────────────────────────────────────────────────────────
# Row filtering (the pipeline-side equivalent of the report's sidebar filters)
# ──────────────────────────────────────────────────────────────────────────────

def _cutoff_from(block, keys=("subkey", "key", "toplevelkey")):
    """Pull best_threshold out of a best_cutoffs block, preferring the same
    granularity the report's slider pre-populates from (subkey)."""
    if not isinstance(block, dict):
        return None
    for k in keys:
        sub = block.get(k)
        if isinstance(sub, dict):
            v = _num(sub.get("best_threshold"))
            if v is not None:
                return v
    return _num(block.get("best_threshold"))


def sample_thresholds(payload):
    """Recommended TASS cutoff per sample, plus the run-wide fallback.

    `Passes Threshold` is computed CLIENT-side in the report (it is emitted as
    False for every record in the JSON), so a pipeline-side export has to
    recompute it. This uses exactly what the report's slider defaults to: the
    sample's own best_cutoffs.subkey.best_threshold, falling back to the run-wide
    best_cutoffs payload, falling back to 0 (everything passes).
    """
    default = _cutoff_from(payload.get("best_cutoffs")) or 0.0
    per_sample = {}
    for sname, meta in (payload.get("sample_meta") or {}).items():
        v = _cutoff_from((meta or {}).get("best_cutoffs"))
        if v is not None:
            per_sample[sname] = v
    return per_sample, default


def _threshold_for(sample, ctx):
    return ctx["thresholds"].get(sample, ctx["default_threshold"])


def _passes(row, ctx):
    """Does this detection clear its sample's recommended cutoff?

    Honours an explicit truthy `Passes Threshold` when a future build starts
    emitting one; otherwise recomputes it (see sample_thresholds).
    """
    if _truthy(row.get("Passes Threshold")):
        return True
    thr = _threshold_for(str(row.get("Specimen ID") or ""), ctx)
    return (_num(row.get("TASS Score")) or 0.0) >= thr


def filter_records(records, min_tass=None, level=None, passing_only=False,
                   high_consequence_only=False, samples=None, thresholds=None,
                   default_threshold=0.0):
    """Apply the coarse filters the report's sidebar offers.

    Deliberately a subset: text search, kingdom toggles and specimen merge are
    interactive state with no pipeline equivalent. `level` collapses the
    Strain/Species/Genus rollup rows to one taxonomic level, which is what keeps
    a pipeline export from triple-counting the same reads.
    """
    out = []
    keep_samples = set(samples) if samples else None
    thr = _num(min_tass)
    for r in records or []:
        if level and (r.get("Level") or "Strain") != level:
            continue
        if keep_samples is not None and r.get("Specimen ID") not in keep_samples:
            continue
        if thr is not None and (_num(r.get("TASS Score")) or 0.0) < thr:
            continue
        if passing_only:
            thr = (thresholds or {}).get(str(r.get("Specimen ID") or ""), default_threshold)
            if not _truthy(r.get("Passes Threshold")) and (_num(r.get("TASS Score")) or 0.0) < thr:
                continue
        if high_consequence_only and not _truthy(r.get("High Consequence")):
            continue
        out.append(r)
    return out


# ──────────────────────────────────────────────────────────────────────────────
# Dataset builders
# ──────────────────────────────────────────────────────────────────────────────

def _bd_detections(ctx):
    cols = list(ctx["all_cols"] or _key_union(ctx["records"]))
    # `Passes Threshold` ships as False for every record (the report computes it
    # in the browser), so append the recomputed verdict and the cutoff it used
    # rather than silently handing the user a column of Nos.
    extra = ["TASS Cutoff", "Passes Cutoff"]
    rows = []
    for r in ctx["records"]:
        row = OrderedDict((c, _yn(r.get(c, ""))) for c in cols)
        row["TASS Cutoff"] = _round(_threshold_for(str(r.get("Specimen ID") or ""), ctx), 2)
        row["Passes Cutoff"] = "Yes" if _passes(r, ctx) else "No"
        rows.append(row)
    return Table(cols + extra, rows)


def _bd_sample_summary(ctx):
    cols = ["Specimen ID", "Specimen Group", "Sample Type", "Platform", "Total Reads",
            "Aligned Reads", "TASS Cutoff", "# Detections", "# Passing Cutoff",
            "# Distinct Organisms", "# High Consequence", "Max TASS Score", "Top Organism",
            "QC Flag"]
    agg = OrderedDict()
    for r in ctx["records"]:
        s = str(r.get("Specimen ID") or "")
        if not s:
            continue
        e = agg.setdefault(s, {"n": 0, "pass": 0, "hc": 0, "orgs": set(),
                               "max": 0.0, "top": ""})
        e["n"] += 1
        e["orgs"].add(r.get("Detected Organism") or "")
        t = _num(r.get("TASS Score")) or 0.0
        if t > e["max"]:
            e["max"] = t
            e["top"] = r.get("Detected Organism") or ""
        if _passes(r, ctx):
            e["pass"] += 1
        if _truthy(r.get("High Consequence")):
            e["hc"] += 1
    rows = []
    for s in sorted(agg):
        e = agg[s]
        meta = (ctx["sample_meta"].get(s) or {})
        rows.append(OrderedDict([
            ("Specimen ID", s),
            ("Specimen Group", _specimen_of(s, ctx)),
            ("Sample Type", meta.get("sample_type") or ""),
            ("Platform", meta.get("platform") or ""),
            ("Total Reads", meta.get("total_reads", "") if meta.get("total_reads") is not None else ""),
            ("Aligned Reads", meta.get("aligned_reads", "") if meta.get("aligned_reads") is not None else ""),
            ("TASS Cutoff", _round(_threshold_for(s, ctx), 2)),
            ("# Detections", e["n"]),
            ("# Passing Cutoff", e["pass"]),
            ("# Distinct Organisms", len(e["orgs"])),
            ("# High Consequence", e["hc"]),
            ("Max TASS Score", _round(e["max"])),
            ("Top Organism", e["top"]),
            # Whole-sample QC verdicts are evaluated live in the report (the
            # pipeline only seeds the default rules), so this column exists for
            # header parity with a browser export and stays blank here.
            ("QC Flag", ""),
        ]))
    return Table(cols, rows)


_SPECIMEN_FIELDS = ("specimen", "specimen_id", "specimen id", "specimenid",
                    "specimen_group", "specimen group")


def _specimen_of(sample, ctx):
    """Specimen a sample belongs to — samplesheet metadata, else the sample
    itself. Mirrors specimenOf() in 01_global_state.js (minus live UI overrides,
    which do not exist outside the browser)."""
    meta = ctx["sample_meta"].get(sample) or {}
    for f in _SPECIMEN_FIELDS:
        v = meta.get(f)
        if v is not None and str(v).strip():
            return str(v).strip()
    return sample


def _median(vals):
    if not vals:
        return 0.0
    s = sorted(vals)
    n = len(s)
    mid = n // 2
    return s[mid] if n % 2 else (s[mid - 1] + s[mid]) / 2.0


_ANI_STR_RE = re.compile(r"^(.+?)\(([\d.]+)%?\)$")


def _ani_matches_for(row):
    """High-ANI partners as [(taxid, pct)].

    Accepts both shapes the pipeline emits: a JSON list of {key, ani_pct}, and
    the serialized "taxid(pct%);taxid(pct%)" string. Mirrors _aniMatchesFor() in
    19_cross_sample_organism.js.
    """
    v = (row or {}).get("High ANI Matches")
    if isinstance(v, (list, tuple)):
        out = []
        for m in v:
            key = str((m or {}).get("key") or "")
            if key:
                out.append((key, _num((m or {}).get("ani_pct")) or 0.0))
        return out
    if isinstance(v, str) and v.strip():
        out = []
        for part in v.split(";"):
            part = part.strip()
            if not part:
                continue
            m = _ANI_STR_RE.match(part)
            out.append((m.group(1).strip(), _num(m.group(2)) or 0.0) if m else (part, 0.0))
        return out
    return []


def _ani_groups(rows_by_taxid, records):
    """Union-find over organisms that share a high-ANI edge, restricted to
    organisms actually present in this view. Returns taxid -> (group, size).
    Mirrors the ANI grouping at the tail of _xsAggregate()."""
    present = set(rows_by_taxid)
    parent = {t: t for t in present}

    def find(x):
        while parent.get(x, x) != x:
            parent[x] = parent.get(parent[x], parent[x])
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for r in records:
        taxid = str(r.get("Taxonomic ID #") or "")
        if taxid not in present:
            continue
        for partner, _pct in _ani_matches_for(r):
            if partner in present:
                union(taxid, partner)

    sizes = defaultdict(int)
    group_of = {}
    for t in present:
        g = find(t)
        group_of[t] = g
        sizes[g] += 1
    return {t: (g, sizes[g]) for t, g in group_of.items()}


def _bd_organism_summary(ctx):
    """Run-wide per-organism rollup. Mirrors _xsAggregate() in
    19_cross_sample_organism.js at its default aggregation (max per specimen)."""
    cols = ["Taxonomic ID #", "Detected Organism", "Microbial Category", "High Consequence",
            "# Specimens Passing", "# Specimens Below Cutoff", "# Specimens Detected",
            "# Specimens Total", "Prevalence %", "Detected Prevalence %",
            "Mean TASS", "Median TASS", "Min TASS", "Max TASS",
            "Mean Coverage", "Median Coverage", "Min Coverage", "Max Coverage",
            "Total Reads Aligned", "ANI Group", "ANI Group Size", "Specimens"]

    # Denominators: every specimen in the run (not just the filtered view), so
    # prevalence means the same thing here as in the report.
    all_specimens = set()
    for r in ctx["records_all"]:
        s = r.get("Specimen ID")
        if s:
            all_specimens.add(_specimen_of(str(s), ctx))
    total = len(all_specimens) or 1

    # Detected-anywhere counts come from the UNFILTERED records, matching the
    # report's "Pass / Below / Total" breakdown.
    det_by_org = defaultdict(set)
    for r in ctx["records_all"]:
        key = str(r.get("Taxonomic ID #") or "") or (r.get("Detected Organism") or "")
        if key:
            det_by_org[key].add(_specimen_of(str(r.get("Specimen ID") or ""), ctx))

    by_org = OrderedDict()
    for r in ctx["records"]:
        spec = _specimen_of(str(r.get("Specimen ID") or ""), ctx)
        if not spec:
            continue
        taxid = str(r.get("Taxonomic ID #") or "")
        key = taxid or (r.get("Detected Organism") or "")
        if not key:
            continue
        e = by_org.get(key)
        if e is None:
            e = {"taxid": taxid or "-", "name": r.get("Detected Organism") or key,
                 "cat": r.get("Microbial Category") or "Unknown", "hc": False,
                 "tass": {}, "cov": {}, "seen": set(), "reads": 0}
            by_org[key] = e
        t = _num(r.get("TASS Score")) or 0.0
        c = _num(r.get("Coverage")) or 0.0
        # "Detected" = seen at all; "Passing" = clears the sample's cutoff. The
        # report draws the same distinction in its Pass / Below / Total column.
        e["seen"].add(spec)
        if _passes(r, ctx):
            e["tass"][spec] = max(e["tass"].get(spec, 0.0), t)
            e["cov"][spec] = max(e["cov"].get(spec, 0.0), c)
        e["reads"] += int(_num(r.get("# Reads Aligned")) or 0)
        if _truthy(r.get("High Consequence")):
            e["hc"] = True

    ani = _ani_groups({e["taxid"]: e for e in by_org.values() if e["taxid"] != "-"},
                      ctx["records"])
    rows = []
    for key, e in by_org.items():
        tass = list(e["tass"].values())
        cov = list(e["cov"].values())
        n_pass = len(tass)
        n_det = max(len(det_by_org.get(key, ())), len(e["seen"]), n_pass)
        rows.append(OrderedDict([
            ("Taxonomic ID #", e["taxid"]),
            ("Detected Organism", e["name"]),
            ("Microbial Category", e["cat"]),
            ("High Consequence", "Yes" if e["hc"] else "No"),
            ("# Specimens Passing", n_pass),
            ("# Specimens Below Cutoff", max(0, n_det - n_pass)),
            ("# Specimens Detected", n_det),
            ("# Specimens Total", total),
            ("Prevalence %", _round(min(100.0, n_pass / total * 100.0), 2)),
            ("Detected Prevalence %", _round(min(100.0, n_det / total * 100.0), 2)),
            ("Mean TASS", _round(sum(tass) / len(tass) if tass else 0)),
            ("Median TASS", _round(_median(tass))),
            ("Min TASS", _round(min(tass) if tass else 0)),
            ("Max TASS", _round(max(tass) if tass else 0)),
            ("Mean Coverage", _round(sum(cov) / len(cov) if cov else 0)),
            ("Median Coverage", _round(_median(cov))),
            ("Min Coverage", _round(min(cov) if cov else 0)),
            ("Max Coverage", _round(max(cov) if cov else 0)),
            ("Total Reads Aligned", e["reads"]),
            ("ANI Group", ani.get(e["taxid"], ("", ""))[0]),
            ("ANI Group Size", ani.get(e["taxid"], ("", ""))[1]),
            ("Specimens", "; ".join(sorted(e["tass"]))),
        ]))
    rows.sort(key=lambda r: (-r["# Specimens Passing"], str(r["Detected Organism"])))
    return Table(cols, rows)


def _bd_coverage(ctx):
    rows = []
    keep = ctx["view_keys"]
    keep_samples = ctx["view_samples"]
    for cd in ctx["contig_data"]:
        sample = cd.get("sample") or ""
        taxid = str(cd.get("taxon_id") or "")
        if keep is not None and (sample, taxid) not in keep and sample not in keep_samples:
            continue
        contigs = cd.get("contigs") or []
        length = covered = reads = 0
        depth_sum = 0.0
        for c in contigs:
            l = _num(c.get("length")) or 0
            length += l
            covered += _num(c.get("covered_bases")) or 0
            reads += _num(c.get("reads")) or 0
            depth_sum += (_num(c.get("mean_depth")) or 0.0) * l
        dh = cd.get("depth_histogram") or {}
        rows.append(OrderedDict([
            ("Specimen ID", sample),
            ("Detected Organism", cd.get("organism") or ""),
            ("Taxonomic ID #", taxid),
            ("Contigs", len(contigs)),
            ("Genome Length (bp)", int(length) or ""),
            ("Covered Bases", int(covered) or ""),
            ("Breadth %", _round(covered / length * 100.0) if length else ""),
            ("Mean Depth", _round(depth_sum / length) if length else ""),
            ("# Reads Aligned", int(reads) or ""),
            ("Bases 0x", dh.get("0x", "")),
            ("Bases 1-5x", dh.get("1-5x", "")),
            ("Bases 5-10x", dh.get("5-10x", "")),
            ("Bases 10-50x", dh.get("10-50x", "")),
            ("Bases >50x", dh.get(">50x", "")),
        ]))
    return Table(_key_union(rows), rows)


def _bd_contigs(ctx):
    rows = []
    keep_samples = ctx["view_samples"]
    for cd in ctx["contig_data"]:
        sample = cd.get("sample") or ""
        if keep_samples and sample not in keep_samples:
            continue
        for c in (cd.get("contigs") or []):
            dh = c.get("depth_histogram") or {}
            rows.append(OrderedDict([
                ("Specimen ID", sample),
                ("Detected Organism", cd.get("organism") or ""),
                ("Taxonomic ID #", str(cd.get("taxon_id") or "")),
                ("Contig", c.get("name") or ""),
                ("Length (bp)", c.get("length", "")),
                ("# Reads Aligned", c.get("reads", "")),
                ("Mean Depth", c.get("mean_depth", "")),
                ("Covered Bases", c.get("covered_bases", "")),
                ("Coverage", c.get("coverage", "")),
                ("Bases 0x", dh.get("0x", "")),
                ("Bases 1-5x", dh.get("1-5x", "")),
                ("Bases 5-10x", dh.get("5-10x", "")),
                ("Bases 10-50x", dh.get("10-50x", "")),
                ("Bases >50x", dh.get(">50x", "")),
            ]))
    return Table(_key_union(rows), rows)


def _prot_table(ctx, key):
    rows = [OrderedDict((k, _yn(v)) for k, v in (r or {}).items())
            for r in (ctx["prot_data"].get(key) or [])]
    return Table(_key_union(rows), rows)


def _bd_novelty_summary(ctx):
    rows = []
    for sname, sdata in (ctx["novelty"].get("samples") or {}).items():
        summary = (sdata or {}).get("summary") or {}
        if not summary:
            continue
        row = OrderedDict([("Specimen ID", sname)])
        row.update(_flatten_scalars(summary))
        rows.append(row)
    return Table(_key_union(rows), rows)


def _bd_novelty_candidates(ctx):
    rows = []
    for sname, sdata in (ctx["novelty"].get("samples") or {}).items():
        for c in ((sdata or {}).get("candidates") or []):
            row = OrderedDict([("Specimen ID", sname)])
            row.update(_flatten_scalars(c))
            pathogen = c.get("pathogen")
            if isinstance(pathogen, dict):
                row["Pathogen Match"] = pathogen.get("name") or pathogen.get("organism") or "Yes"
                if pathogen.get("taxid") is not None:
                    row["Pathogen Taxid"] = pathogen.get("taxid")
            rows.append(row)
    return Table(_key_union(rows), rows)


def _bd_run_metadata(ctx):
    """Per-sample run metadata.

    The report seeds an empty row for every sample that has no metadata record
    (so a user can type one in on the Run Metadata tab). Do the same here, or a
    pipeline export would silently be missing the samples whose metadata was
    never supplied -- and would not line up row-for-row with a browser export.
    """
    rows = [_flatten_scalars(m) for m in ctx["run_metadata"]]
    have = {r.get("sample_name") for r in rows}
    seen = set()
    extra = []
    for r in ctx["records_all"]:
        sid = str(r.get("Specimen ID") or "")
        if sid and sid not in have and sid not in seen:
            seen.add(sid)
            extra.append(OrderedDict([("sample_name", sid)]))
    extra.sort(key=lambda r: r["sample_name"])
    rows.extend(extra)
    return Table(_key_union(rows), rows)


def _bd_sample_metadata(ctx):
    rows = []
    for s in sorted(ctx["sample_meta"]):
        row = OrderedDict([("Specimen ID", s)])
        row.update(_flatten_scalars(ctx["sample_meta"][s]))
        rows.append(row)
    return Table(_key_union(rows), rows)


def _bd_geo(ctx):
    rows = []
    for m in ctx["run_metadata"]:
        lat, lon = _num(m.get("latitude")), _num(m.get("longitude"))
        if lat is None or lon is None:
            continue
        rows.append(OrderedDict([
            ("Specimen ID", m.get("sample_name") or m.get("sample_id") or ""),
            ("Latitude", lat),
            ("Longitude", lon),
            ("Location", m.get("location") or ""),
            ("Country", m.get("sample_origin_country") or ""),
            ("State/Province", m.get("sample_origin_state_province_territory") or ""),
            ("Environmental Site", m.get("environmental_site") or ""),
            ("Collection Time", m.get("collection_time") or ""),
            ("Run ID", m.get("run_id") or ""),
        ]))
    return Table(_key_union(rows), rows)


def _insilico_head(g):
    return OrderedDict([
        ("Parent Sample", g.get("parent") or ""),
        ("Platform", g.get("platform") or ""),
        ("Series Kind", g.get("series_kind") or "depth"),
        ("Level", g.get("level") or ""),
        ("Read Unit", g.get("read_unit") or "reads"),
    ])


def _bd_insilico_datasets(ctx):
    rows = []
    for g in ((ctx["insilico"] or {}).get("groups") or []):
        for d in (g.get("datasets") or []):
            row = _insilico_head(g)
            row.update(OrderedDict([
                ("Dataset ID", d.get("id", "")),
                ("Replicate", d.get("replicate", "")),
                ("Target Count", d.get("target_count", "")),
                ("Actual Count", d.get("actual_count", "")),
                ("Total Master Reads", d.get("total_master_reads", "")),
                ("Seed", d.get("seed", "")),
                ("Observed Total Reads", d.get("observed_total_reads", "")),
                ("# Detected", d.get("n_detected", "")),
                ("TP", d.get("tp", "")),
                ("FP", d.get("fp", "")),
                ("FN", d.get("fn", "")),
                ("Precision", d.get("precision", "")),
                ("Recall", d.get("recall", "")),
                ("F1", d.get("f1", "")),
            ]))
            rows.append(row)
    return Table(_key_union(rows), rows)


def _bd_insilico_lod(ctx):
    rows = []
    for g in ((ctx["insilico"] or {}).get("groups") or []):
        for o in (g.get("organisms") or []):
            for s in (o.get("series") or []):
                row = _insilico_head(g)
                row.update(OrderedDict([
                    ("Taxonomic ID #", str(o.get("taxid") or "")),
                    ("Detected Organism", o.get("name") or ""),
                    ("Microbial Category", o.get("category") or ""),
                    ("Species Name", o.get("species") or ""),
                    ("Genus Name", o.get("genus") or ""),
                    ("Expected Fraction", o.get("expected_fraction", "")),
                    ("LoD Count", o.get("lod_count") if o.get("lod_count") is not None else ""),
                    ("Series Count", s.get("count", "")),
                    ("Expected Reads", s.get("expected_reads", "")),
                    ("Observed Reads", s.get("observed_reads", "")),
                    ("TASS Score", s.get("tass", "")),
                    ("Detection Rate", s.get("detection_rate", "")),
                    ("Detected", "Yes" if s.get("detected") else "No"),
                    ("# Replicates", s.get("n_reps", "")),
                ]))
                rows.append(row)
    return Table(_key_union(rows), rows)


# ──────────────────────────────────────────────────────────────────────────────
# The catalog  (ids + labels MUST match TT_EXPORT_DATASETS in the JS module)
# ──────────────────────────────────────────────────────────────────────────────

DATASETS = [
    dict(id="detections", label="Detections", tab="Summary / Table",
         join=JOIN_SAMPLE_ORG, default=True, build=_bd_detections),
    dict(id="sample_summary", label="Sample Summary", tab="Summary",
         join=JOIN_SAMPLE, default=True, build=_bd_sample_summary),
    dict(id="organism_summary", label="Cross-Sample Organisms", tab="Explore",
         join=JOIN_ORG, default=True, build=_bd_organism_summary),
    dict(id="coverage", label="Coverage Summary", tab="Coverage / Histogram",
         join=JOIN_SAMPLE_ORG, default=True, build=_bd_coverage),
    dict(id="contigs", label="Coverage by Contig", tab="Coverage / Histogram",
         join=JOIN_NONE, default=False, build=_bd_contigs),
    dict(id="vfamr_hits", label="VF/AMR Per-Gene Hits", tab="VF/AMR",
         join=JOIN_NONE, default=True, build=lambda c: _prot_table(c, "per_gene_hits")),
    dict(id="vfamr_genus", label="VF/AMR Genus Summary", tab="VF/AMR",
         join=JOIN_NONE, default=True, build=lambda c: _prot_table(c, "genus_summary")),
    dict(id="vfamr_amr", label="AMR Genes", tab="VF/AMR",
         join=JOIN_NONE, default=True, build=lambda c: _prot_table(c, "amr_genes")),
    dict(id="novelty_summary", label="Novelty Summary", tab="Novelty",
         join=JOIN_SAMPLE, default=True, build=_bd_novelty_summary),
    dict(id="novelty_candidates", label="Novelty Candidates", tab="Novelty",
         join=JOIN_NONE, default=True, build=_bd_novelty_candidates),
    dict(id="run_metadata", label="Run Metadata", tab="Run Metadata",
         join=JOIN_SAMPLE, default=True, build=_bd_run_metadata),
    dict(id="sample_metadata", label="Pipeline Sample Metadata", tab="Run Metadata",
         join=JOIN_SAMPLE, default=False, build=_bd_sample_metadata),
    dict(id="geo", label="Sample Geography", tab="Map",
         join=JOIN_SAMPLE, default=False, build=_bd_geo),
    dict(id="insilico_datasets", label="In-Silico Datasets", tab="In-Silico",
         join=JOIN_NONE, default=True, build=_bd_insilico_datasets),
    dict(id="insilico_lod", label="In-Silico LoD Series", tab="In-Silico",
         join=JOIN_NONE, default=True, build=_bd_insilico_lod),
]

DATASET_IDS = [d["id"] for d in DATASETS]
_BY_ID = {d["id"]: d for d in DATASETS}


# ──────────────────────────────────────────────────────────────────────────────
# Context + build
# ──────────────────────────────────────────────────────────────────────────────

def make_context(payload, min_tass=None, level=None, passing_only=False,
                 high_consequence_only=False, samples=None):
    """Normalise a bootstrap payload into the dict every builder reads."""
    records_all = payload.get("records") or []
    thresholds, default_threshold = sample_thresholds(payload)
    records = filter_records(records_all, min_tass=min_tass, level=level,
                             passing_only=passing_only,
                             high_consequence_only=high_consequence_only,
                             samples=samples, thresholds=thresholds,
                             default_threshold=default_threshold)
    filtered = len(records) != len(records_all)
    view_keys = None
    view_samples = set()
    if filtered:
        view_keys = {(r.get("Specimen ID"), str(r.get("Taxonomic ID #") or "")) for r in records}
        view_samples = {r.get("Specimen ID") for r in records}
    return {
        "payload": payload,
        "records": records,
        "records_all": records_all,
        "all_cols": [c for c in (payload.get("all_cols") or [])
                     if c not in ("High ANI Matches", "ANI Annotated")],
        "sample_meta": payload.get("sample_meta") or {},
        "prot_data": payload.get("prot_data") or {},
        "contig_data": payload.get("contig_data") or [],
        "run_metadata": payload.get("run_metadata_records") or [],
        "novelty": payload.get("novelty") or {"samples": {}},
        "insilico": payload.get("insilico_suite") or None,
        "thresholds": thresholds,
        "default_threshold": default_threshold,
        "filtered": filtered,
        "view_keys": view_keys,
        "view_samples": view_samples,
    }


def build_tables(ctx, dataset_ids=None, drop_empty=True):
    """Build each requested dataset. Returns an OrderedDict id -> Table.

    A dataset the run carries no data for is dropped rather than written as an
    empty sheet, so an export of a run without VF/AMR does not hand the user
    four blank tabs to wonder about.
    """
    ids = list(dataset_ids) if dataset_ids else [d["id"] for d in DATASETS if d["default"]]
    out = OrderedDict()
    for did in ids:
        spec = _BY_ID.get(did)
        if not spec:
            print(f"[export_data] WARNING: unknown dataset '{did}' (known: "
                  f"{', '.join(DATASET_IDS)})", file=sys.stderr)
            continue
        try:
            table = spec["build"](ctx)
        except Exception as exc:  # noqa: BLE001 - one bad dataset must not sink the export
            print(f"[export_data] WARNING: dataset '{did}' failed: {exc}", file=sys.stderr)
            continue
        if drop_empty and not len(table):
            continue
        out[did] = table
    return out


def wide_join(tables, ctx):
    """Join the built tables onto the detections backbone, one row per
    Specimen ID x Organism. Mirrors _ttWideJoin() in the JS module.

    Datasets with no join key (several rows per organism — per-gene hits,
    novelty candidates, in-silico series) cannot be widened; their labels come
    back in `skipped` so the caller can say so.
    """
    backbone = tables.get("detections") or _bd_detections(ctx)
    columns = list(backbone.columns)
    rows = [OrderedDict(r) for r in backbone.rows]
    skipped = []
    key_cols = {"Specimen ID", "Taxonomic ID #", "Detected Organism"}

    def _key(row, join):
        s = str(row.get("Specimen ID") or "")
        o = str(row.get("Taxonomic ID #") or row.get("Detected Organism") or "")
        if join == JOIN_SAMPLE:
            return s
        if join == JOIN_ORG:
            return o
        return s + "" + o

    for did, table in tables.items():
        if did == "detections":
            continue
        spec = _BY_ID[did]
        if not spec["join"]:
            skipped.append(spec["label"])
            continue
        if not len(table):
            continue
        index = {}
        for r in table.rows:
            index.setdefault(_key(r, spec["join"]), r)
        out_cols = [c for c in table.columns if c not in key_cols]
        prefixed = [spec["label"] + " · " + c for c in out_cols]
        columns.extend(prefixed)
        for row in rows:
            hit = index.get(_key(row, spec["join"]))
            for c, pc in zip(out_cols, prefixed):
                row[pc] = "" if not hit else ("" if hit.get(c) is None else hit.get(c))
    return Table(columns, rows), skipped


def stacked(tables):
    """Every table in one, with a leading Dataset column and the union of
    columns — the single-CSV shape for datasets that share no join key."""
    columns = ["Dataset"]
    seen = {"Dataset"}
    rows = []
    for did, table in tables.items():
        label = _BY_ID[did]["label"]
        for c in table.columns:
            if c not in seen:
                seen.add(c)
                columns.append(c)
        for r in table.rows:
            row = OrderedDict([("Dataset", label)])
            row.update(r)
            rows.append(row)
    return Table(columns, rows)


def manifest_table(tables, ctx, opts):
    """Provenance sheet: what was exported, from which run, under what filters."""
    p = ctx["payload"]
    rows = [
        ["TaxTriage combined data export"],
        ["Generated", datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")],
        ["Report built", p.get("report_generated_at") or ""],
        ["Pipeline revision", p.get("pipeline_revision") or ""],
        ["Pipeline commit", p.get("pipeline_commit") or ""],
        ["Source", opts.get("source") or ""],
        [],
        ["Filters applied"],
        ["Min TASS score", "" if opts.get("min_tass") is None else opts.get("min_tass")],
        ["Taxonomic level", opts.get("level") or "all levels"],
        ["Passing only", "Yes" if opts.get("passing_only") else "No"],
        ["High consequence only", "Yes" if opts.get("high_consequence_only") else "No"],
        ["Samples", "; ".join(opts["samples"]) if opts.get("samples") else "all"],
        ["Detections kept", f"{len(ctx['records'])} of {len(ctx['records_all'])}"],
        ["Recommended TASS cutoff", ctx["default_threshold"]],
        ["Per-sample cutoffs", "; ".join(f"{k}={v}" for k, v in sorted(ctx["thresholds"].items())) or "none"],
        [],
        ["Datasets included", "Tab", "Rows"],
    ]
    for did, table in tables.items():
        rows.append([_BY_ID[did]["label"], _BY_ID[did]["tab"], len(table)])
    # A Table needs dict rows; this sheet is free-form, so hand back the AoA
    # directly and let the writers special-case it.
    return rows


# ──────────────────────────────────────────────────────────────────────────────
# Writers
# ──────────────────────────────────────────────────────────────────────────────

_SHEET_BAD = re.compile(r"[\\/\?\*\[\]:]")


def _sheet_name(label, used):
    base = _SHEET_BAD.sub("-", str(label or "Sheet"))[:31] or "Sheet"
    name, n = base, 2
    while name in used:
        suffix = f"_{n}"
        n += 1
        name = base[:31 - len(suffix)] + suffix
    used.add(name)
    return name


def _write_csv(table, path, delimiter=","):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter=delimiter, quoting=csv.QUOTE_MINIMAL, lineterminator="\r\n")
        for row in table.aoa():
            w.writerow(row)
    return path


def _write_xlsx(sheets, path):
    """sheets: list of (name, aoa). openpyxl is already a make_report.py
    dependency (pandas' Excel engine), so no new requirement."""
    try:
        from openpyxl import Workbook
    except ImportError:
        raise SystemExit("[export_data] ERROR: openpyxl is required for XLSX output "
                         "(pip install openpyxl), or use --formats csv")
    wb = Workbook()
    wb.remove(wb.active)
    used = set()
    for name, aoa in sheets:
        ws = wb.create_sheet(_sheet_name(name, used))
        for row in aoa:
            ws.append(list(row))
        # Freeze the header row of real data sheets so a wide export stays
        # readable when scrolled.
        if len(aoa) > 1:
            ws.freeze_panes = "A2"
    wb.save(path)
    return path


def write_exports(payload, outdir, formats=("xlsx",), datasets=None, prefix="taxtriage",
                  min_tass=None, level=None, passing_only=False,
                  high_consequence_only=False, samples=None, delimiter=",",
                  source=""):
    """Build and write the combined export. Returns the list of paths written.

    formats (any combination):
      xlsx     one workbook, an Export Info sheet plus one sheet per dataset
      wide     the Sample x Organism join — .xlsx sheet and .csv
      csv      one CSV per dataset
      stacked  every dataset in one CSV with a leading Dataset column
    """
    os.makedirs(outdir, exist_ok=True)
    ctx = make_context(payload, min_tass=min_tass, level=level,
                       passing_only=passing_only,
                       high_consequence_only=high_consequence_only,
                       samples=samples)
    tables = build_tables(ctx, datasets)
    if not tables:
        print("[export_data] nothing to export (no dataset produced rows)", file=sys.stderr)
        return []

    opts = dict(min_tass=min_tass, level=level, passing_only=passing_only,
                high_consequence_only=high_consequence_only, samples=samples,
                source=source)
    info = manifest_table(tables, ctx, opts)
    formats = [f.strip().lower() for f in formats if f and f.strip()]
    written = []
    ext = "tsv" if delimiter == "\t" else "csv"

    if "xlsx" in formats:
        sheets = [("Export Info", info)]
        sheets += [(_BY_ID[did]["label"], t.aoa()) for did, t in tables.items()]
        written.append(_write_xlsx(sheets, os.path.join(outdir, f"{prefix}.combined.xlsx")))

    if "wide" in formats:
        table, skipped = wide_join(tables, ctx)
        if skipped:
            print(f"[export_data] wide join: {', '.join(skipped)} have several rows per "
                  f"organism and were left out of the joined table", file=sys.stderr)
        written.append(_write_csv(table, os.path.join(outdir, f"{prefix}.wide.{ext}"), delimiter))
        written.append(_write_xlsx([("Export Info", info), ("Combined", table.aoa())],
                                   os.path.join(outdir, f"{prefix}.wide.xlsx")))

    if "csv" in formats:
        for did, table in tables.items():
            written.append(_write_csv(table, os.path.join(outdir, f"{prefix}.{did}.{ext}"), delimiter))

    if "stacked" in formats:
        written.append(_write_csv(stacked(tables),
                                  os.path.join(outdir, f"{prefix}.stacked.{ext}"), delimiter))

    print(f"[export_data] {len(tables)} dataset(s) "
          f"({', '.join(f'{did}={len(t)}' for did, t in tables.items())})")
    for p in written:
        print(f"[export_data] Written: {p}")
    return written


# ──────────────────────────────────────────────────────────────────────────────
# Payload loading (report HTML or raw JSON)
# ──────────────────────────────────────────────────────────────────────────────

_BOOT_ANCHOR = "window.HEATMAP_BOOT ="


def load_payload(path):
    """Read a bootstrap payload from a built report HTML or a JSON file."""
    with open(path, encoding="utf-8") as fh:
        text = fh.read()
    stripped = text.lstrip()
    if stripped.startswith("{"):
        return json.loads(text)
    idx = text.find(_BOOT_ANCHOR)
    if idx == -1:
        raise SystemExit(f"[export_data] ERROR: no {_BOOT_ANCHOR} payload found in {path} "
                         "(is this a TaxTriage report?)")
    start = idx + len(_BOOT_ANCHOR)
    end = text.find("\n", start)
    raw = (text[start:end] if end != -1 else text[start:]).strip().rstrip(";").strip()
    try:
        return json.loads(raw)
    except json.JSONDecodeError:
        # Older builds pretty-printed the payload across several lines: fall back
        # to a brace scan from the first '{'.
        brace = text.find("{", start)
        block = None
        depth, in_str, esc = 0, False, False
        for i in range(brace, len(text)):
            ch = text[i]
            if in_str:
                if esc:
                    esc = False
                elif ch == "\\":
                    esc = True
                elif ch == '"':
                    in_str = False
                continue
            if ch == '"':
                in_str = True
            elif ch == "{":
                depth += 1
            elif ch == "}":
                depth -= 1
                if depth == 0:
                    block = text[brace:i + 1]
                    break
        if block is not None:
            try:
                return json.loads(block)
            except json.JSONDecodeError:
                pass
        raise SystemExit(
            f"[export_data] ERROR: the payload in {path} is not strict JSON.\n"
            "  Pipeline-built reports embed the payload as one JSON line and export fine.\n"
            "  A hand-maintained demo report (e.g. assets/pages.js, examples/) stores it as a\n"
            "  pretty-printed JavaScript object literal with unquoted keys, which cannot be\n"
            "  read back. Point -i at a report produced by the pipeline, or at the\n"
            "  bootstrap JSON itself."
        )


# ──────────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────────

def parse_args(argv=None):
    ap = argparse.ArgumentParser(
        description="Export TaxTriage report data across multiple tabs into one "
                    "spreadsheet / CSV set. Reads a built report HTML or a raw "
                    "bootstrap JSON.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Datasets: " + ", ".join(DATASET_IDS),
    )
    ap.add_argument("-i", "--input", default=None,
                    help="Built report (all.odr.html / all.comparison.report.html) "
                         "or a bootstrap JSON payload.")
    ap.add_argument("-o", "--outdir", default="taxtriage_export",
                    help="Directory to write into (created if needed).")
    ap.add_argument("--prefix", default="taxtriage",
                    help="File-name prefix for the written files (default: taxtriage).")
    ap.add_argument("--formats", default="xlsx",
                    help="Comma-separated: xlsx (one sheet per dataset), wide "
                         "(Sample x Organism join, written as both .xlsx and .csv), "
                         "csv (one file per dataset), stacked (all datasets in one "
                         "CSV with a Dataset column). Default: xlsx.")
    ap.add_argument("--datasets", default=None,
                    help="Comma-separated dataset ids to include, or 'all'. "
                         "Default: every dataset the run carries data for.")
    ap.add_argument("--list-datasets", action="store_true",
                    help="Print the dataset catalog and exit.")
    ap.add_argument("--min-tass", type=float, default=None,
                    help="Drop detections scoring below this TASS value.")
    ap.add_argument("--level", default=None,
                    choices=["Strain", "Species", "Genus"],
                    help="Keep only this taxonomic level (avoids counting the same "
                         "reads at strain, species and genus).")
    ap.add_argument("--passing-only", action="store_true",
                    help="Keep only detections flagged as passing the threshold.")
    ap.add_argument("--high-consequence-only", action="store_true",
                    help="Keep only high-consequence organisms.")
    ap.add_argument("--samples", default=None,
                    help="Comma-separated Specimen IDs to restrict the export to.")
    ap.add_argument("--delimiter", default=",",
                    help="Delimiter for CSV output: ',' (default), '\\t', ';' or '|'.")
    return ap.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    if args.list_datasets:
        width = max(len(d["id"]) for d in DATASETS)
        for d in DATASETS:
            star = "*" if d["default"] else " "
            print(f"{star} {d['id']:<{width}}  {d['tab']:<22}  {d['label']}")
        print("\n* = included by default")
        return 0

    if not args.input:
        ap_err = "[export_data] ERROR: -i/--input is required (or use --list-datasets)"
        raise SystemExit(ap_err)
    payload = load_payload(args.input)
    datasets = None
    if args.datasets and args.datasets.strip().lower() != "all":
        datasets = [d.strip() for d in args.datasets.split(",") if d.strip()]
    elif args.datasets:
        datasets = DATASET_IDS
    delimiter = "\t" if args.delimiter in ("\\t", "\t", "tab") else args.delimiter

    write_exports(
        payload, args.outdir,
        formats=[f for f in args.formats.split(",")],
        datasets=datasets,
        prefix=args.prefix,
        min_tass=args.min_tass,
        level=args.level,
        passing_only=args.passing_only,
        high_consequence_only=args.high_consequence_only,
        samples=[s.strip() for s in args.samples.split(",")] if args.samples else None,
        delimiter=delimiter,
        source=os.path.basename(args.input),
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
