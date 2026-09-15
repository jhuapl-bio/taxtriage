import re, io, sys
p = sys.argv[1]
s = open(p, encoding="utf-8").read()

# ── 1. threshold helpers, inserted before filter_records ───────────────────────
anchor = "def filter_records(records"
helpers = '''def _cutoff_from(block, keys=("subkey", "key", "toplevelkey")):
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


'''
assert s.count(anchor) == 1
s = s.replace(anchor, helpers + anchor, 1)

# ── 2. filter_records takes a ctx-less threshold map for passing_only ─────────
s = s.replace('''def filter_records(records, min_tass=None, level=None, passing_only=False,
                   high_consequence_only=False, samples=None):''',
'''def filter_records(records, min_tass=None, level=None, passing_only=False,
                   high_consequence_only=False, samples=None, thresholds=None,
                   default_threshold=0.0):''')
s = s.replace('''        if passing_only and not _truthy(r.get("Passes Threshold")):
            continue''',
'''        if passing_only:
            thr = (thresholds or {}).get(str(r.get("Specimen ID") or ""), default_threshold)
            if not _truthy(r.get("Passes Threshold")) and (_num(r.get("TASS Score")) or 0.0) < thr:
                continue''')

# ── 3. detections gains the cutoff columns ───────────────────────────────────
s = s.replace('''def _bd_detections(ctx):
    cols = list(ctx["all_cols"] or _key_union(ctx["records"]))
    rows = []
    for r in ctx["records"]:
        rows.append(OrderedDict((c, _yn(r.get(c, ""))) for c in cols))
    return Table(cols, rows)''',
'''def _bd_detections(ctx):
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
    return Table(cols + extra, rows)''')

# ── 4. sample_summary: passing via _passes + cutoff column ───────────────────
s = s.replace('''    cols = ["Specimen ID", "Specimen Group", "Sample Type", "Platform", "Total Reads",
            "Aligned Reads", "# Detections", "# Passing Threshold", "# Distinct Organisms",
            "# High Consequence", "Max TASS Score", "Top Organism"]''',
'''    cols = ["Specimen ID", "Specimen Group", "Sample Type", "Platform", "Total Reads",
            "Aligned Reads", "TASS Cutoff", "# Detections", "# Passing Cutoff",
            "# Distinct Organisms", "# High Consequence", "Max TASS Score", "Top Organism"]''')
s = s.replace('''        if _truthy(r.get("Passes Threshold")):
            e["pass"] += 1''',
'''        if _passes(r, ctx):
            e["pass"] += 1''')
s = s.replace('''            ("Aligned Reads", meta.get("aligned_reads", "") if meta.get("aligned_reads") is not None else ""),
            ("# Detections", e["n"]),
            ("# Passing Threshold", e["pass"]),''',
'''            ("Aligned Reads", meta.get("aligned_reads", "") if meta.get("aligned_reads") is not None else ""),
            ("TASS Cutoff", _round(_threshold_for(s, ctx), 2)),
            ("# Detections", e["n"]),
            ("# Passing Cutoff", e["pass"]),''')

# ── 5. organism_summary: split passing vs detected on the cutoff ─────────────
s = s.replace('''        e["tass"][spec] = max(e["tass"].get(spec, 0.0), t)
        e["cov"][spec] = max(e["cov"].get(spec, 0.0), c)''',
'''        # "Detected" = seen at all; "Passing" = clears the sample's cutoff. The
        # report draws the same distinction in its Pass / Below / Total column.
        e["seen"].add(spec)
        if _passes(r, ctx):
            e["tass"][spec] = max(e["tass"].get(spec, 0.0), t)
            e["cov"][spec] = max(e["cov"].get(spec, 0.0), c)''')
s = s.replace('''                 "tass": {}, "cov": {}, "reads": 0}''',
'''                 "tass": {}, "cov": {}, "seen": set(), "reads": 0}''')
s = s.replace('''        n_det = max(len(det_by_org.get(key, ())), n_pass)''',
'''        n_det = max(len(det_by_org.get(key, ())), len(e["seen"]), n_pass)''')

# ── 6. make_context wires the thresholds in ──────────────────────────────────
s = s.replace('''    records_all = payload.get("records") or []
    records = filter_records(records_all, min_tass=min_tass, level=level,
                             passing_only=passing_only,
                             high_consequence_only=high_consequence_only,
                             samples=samples)''',
'''    records_all = payload.get("records") or []
    thresholds, default_threshold = sample_thresholds(payload)
    records = filter_records(records_all, min_tass=min_tass, level=level,
                             passing_only=passing_only,
                             high_consequence_only=high_consequence_only,
                             samples=samples, thresholds=thresholds,
                             default_threshold=default_threshold)''')
s = s.replace('''        "filtered": filtered,''',
'''        "thresholds": thresholds,
        "default_threshold": default_threshold,
        "filtered": filtered,''')

# ── 7. manifest records the cutoff used ──────────────────────────────────────
s = s.replace('''        ["Detections kept", f"{len(ctx['records'])} of {len(ctx['records_all'])}"],''',
'''        ["Detections kept", f"{len(ctx['records'])} of {len(ctx['records_all'])}"],
        ["Recommended TASS cutoff", ctx["default_threshold"]],
        ["Per-sample cutoffs", "; ".join(f"{k}={v}" for k, v in sorted(ctx["thresholds"].items())) or "none"],''')

open(p, "w", encoding="utf-8").write(s)
print("patched")
