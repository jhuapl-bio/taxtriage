#!/usr/bin/env python3
"""
mask_compare.py - before / after TASS for region-masked scoring.

match_paths.py --mask_compare re-scores the same BAM WITHOUT the region mask
(a child match_paths.py run: by default run to completion first so peak memory
stays that of one run; --mask_compare_parallel runs it alongside) and folds the unmasked numbers into the masked JSON:

  * every organism entry (genus / species / strain) gets
        "mask_compare": {tass_unmasked, tass_masked, tass_delta,
                         reads_unmasked, reads_masked,
                         breadth_unmasked, breadth_masked,
                         gini_unmasked, gini_masked, present_unmasked}
  * metadata.masking.compare summarises the sample: how many strains moved,
    mean / median change, biggest gain and loss, organisms that only appear
    without the mask (all their reads were masked) and organisms only scored
    with it.

The unmasked JSON itself is kept (default <outdir>/unmasked/<name>.paths.json)
so the stand-alone masking report (mask_report.py --baseline) can use it too.

Can also be run on its own to merge an existing unmasked JSON:
    mask_compare.py --masked S.paths.json --unmasked S.unmasked.paths.json -o S.paths.json
"""
import argparse
import json
import os
import shutil
import statistics
import subprocess
import sys
import time

# match_paths.py options that must NOT reach the unmasked child run: the mask
# itself, the compare switches, and every option naming an output path (the
# child writes into its own directory).
_DROP_OPTS = {
    "--mask_bed", "--mask_exclude_bed", "--mask_mode", "--mask_read_overlap",
    "--mask_untargeted", "--mask_pad", "--mask_min_region", "--mask_penalty_power",
    "--mask_write_bam", "--mask_regions_json_max",
    "--mask_compare", "--mask_compare_parallel", "--mask_compare_out", "--mask_baseline_json",
    "-o", "--output", "--output_dir", "--confusion_xlsx", "--report_confusion_xlsx",
    "--report_metrics", "--config",
}

_DELTA_EPS = 0.5  # |TASS change| below this counts as "unchanged"


def _is_opt(tok):
    if not tok.startswith("-") or tok == "-":
        return False
    try:
        float(tok)
        return False  # negative number = a value
    except ValueError:
        return True


def child_argv(parser, argv, out_json, out_dir):
    """argv for the unmasked re-run: drop mask / output options (and their values)."""
    acts = parser._option_string_actions
    keep, i = [], 0
    while i < len(argv):
        tok = argv[i]
        name = tok.split("=", 1)[0]
        act = acts.get(name) if _is_opt(tok) else None
        n_vals = 0
        if act is not None and "=" not in tok:
            na = act.nargs
            if na is None:
                n_vals = 1
            elif na == 0:
                n_vals = 0
            elif isinstance(na, int):
                n_vals = na
            else:  # '?', '+', '*'
                j = i + 1
                while j < len(argv) and not _is_opt(argv[j]):
                    j += 1
                n_vals = (j - i - 1) if na != "?" else min(1, j - i - 1)
        if act is not None and (name in _DROP_OPTS or any(o in _DROP_OPTS for o in act.option_strings)):
            i += 1 + n_vals
            continue
        keep.extend(argv[i:i + 1 + n_vals])
        i += 1 + n_vals
    return keep + ["-o", out_json, "--output_dir", out_dir]


def default_out(output):
    d = os.path.dirname(output) or "."
    return os.path.join(d, "unmasked", os.path.basename(output).replace(".tsv", ".json"))


def launch(script, parser, argv, out_json):
    """Start the unmasked child run in the background. Returns a handle or None."""
    out_dir = os.path.dirname(out_json) or "."
    work = os.path.join(out_dir, "search_results")
    os.makedirs(work, exist_ok=True)
    cmd = [sys.executable, script] + child_argv(parser, argv, out_json, work)
    log_path = os.path.splitext(out_json)[0] + ".log"
    log = open(log_path, "w")
    env = dict(os.environ, TT_MASK_COMPARE_CHILD="1")
    print(f"[mask_compare] scoring the unmasked baseline -> {out_json}")
    proc = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, env=env)
    return {"proc": proc, "log": log, "log_path": log_path, "out": out_json, "work": work,
            "t0": time.time()}


def wait(handle):
    """Wait for the child; return the unmasked JSON path or None on failure."""
    if not handle:
        return None
    rc = handle["proc"].wait()
    handle["log"].close()
    shutil.rmtree(handle["work"], ignore_errors=True)
    dt = time.time() - handle["t0"]
    if rc != 0 or not os.path.exists(handle["out"]):
        print(f"[mask_compare] WARNING: unmasked baseline run failed (exit {rc}); "
              f"before/after TASS not added. Log: {handle['log_path']}", file=sys.stderr)
        try:
            with open(handle["log_path"]) as fh:
                tail = fh.readlines()[-15:]
            sys.stderr.write("".join(tail))
        except OSError:
            pass
        return None
    print(f"[mask_compare] unmasked baseline finished in {dt:.0f}s")
    return handle["out"]


# ── merge ────────────────────────────────────────────────────────────────────
def _f(v, d=0.0):
    try:
        return float(v)
    except (TypeError, ValueError):
        return d


def _breadth(e):
    ln = _f(e.get("length"))
    return round(100.0 * _f(e.get("covered_bases")) / ln, 3) if ln > 0 else 0.0


def _walk(entries, parent=(), depth=0):
    for e in entries or []:
        if not isinstance(e, dict):
            continue
        k = parent + (f"{e.get('key', e.get('toplevelkey', ''))}|{e.get('name', '')}",)
        yield k, depth, e
        yield from _walk(e.get("members"), k, depth + 1)


_LEVELS = {0: "genus", 1: "species", 2: "strain"}


def merge(masked, unmasked, source="inline"):
    """Fold unmasked TASS into the masked JSON (in place). Returns the summary dict."""
    base = {k: (d, e) for k, d, e in _walk(unmasked.get("organisms"))}
    seen = set()
    deltas, rows = [], []
    n_new = 0
    for k, d, e in _walk(masked.get("organisms")):
        seen.add(k)
        b = base.get(k, (None, None))[1]
        tm = round(_f(e.get("tass_score")), 4)
        tu = round(_f(b.get("tass_score")), 4) if b is not None else 0.0
        e["mask_compare"] = {
            "tass_unmasked": tu,
            "tass_masked": tm,
            "tass_delta": round(tm - tu, 4),
            "reads_unmasked": _f(b.get("numreads")) if b is not None else 0.0,
            "reads_masked": _f(e.get("numreads")),
            "breadth_unmasked": _breadth(b) if b is not None else 0.0,
            "breadth_masked": _breadth(e),
            "gini_unmasked": round(_f(b.get("gini_coefficient")), 4) if b is not None else None,
            "gini_masked": round(_f(e.get("gini_coefficient")), 4),
            "present_unmasked": b is not None,
        }
        if d == 2:
            if b is None:
                n_new += 1
            deltas.append(tm - tu)
            rows.append((tm - tu, e.get("name", ""), tu, tm))
    dropped = []
    for k, (d, e) in base.items():
        if k in seen:
            continue
        # only report the highest level that disappeared
        if any(k[:i] in base and k[:i] not in seen for i in range(1, len(k))):
            continue
        dropped.append({
            "level": _LEVELS.get(d, str(d)),
            "name": e.get("name", ""),
            "key": str(e.get("key", e.get("toplevelkey", ""))),
            "tass_unmasked": round(_f(e.get("tass_score")), 4),
            "reads_unmasked": _f(e.get("numreads")),
        })
    dropped.sort(key=lambda x: -x["tass_unmasked"])
    rows.sort()
    up = sum(1 for x in deltas if x >= _DELTA_EPS)
    down = sum(1 for x in deltas if x <= -_DELTA_EPS)
    summ = {
        "source": source,
        "strains_compared": len(deltas),
        "strains_up": up,
        "strains_down": down,
        "strains_unchanged": len(deltas) - up - down,
        "strains_only_masked": n_new,
        "mean_delta": round(statistics.fmean(deltas), 4) if deltas else 0.0,
        "median_delta": round(statistics.median(deltas), 4) if deltas else 0.0,
        "max_gain": ({"name": rows[-1][1], "tass_unmasked": rows[-1][2], "tass_masked": rows[-1][3],
                      "delta": round(rows[-1][0], 4)} if rows and rows[-1][0] > 0 else None),
        "max_loss": ({"name": rows[0][1], "tass_unmasked": rows[0][2], "tass_masked": rows[0][3],
                      "delta": round(rows[0][0], 4)} if rows and rows[0][0] < 0 else None),
        "dropped": dropped[:500],
        "n_dropped": len(dropped),
        "unchanged_eps": _DELTA_EPS,
    }
    md = masked.setdefault("metadata", {})
    if isinstance(md.get("masking"), dict):
        md["masking"]["compare"] = summ
    else:
        md["mask_compare"] = summ
    return summ


def describe(summ):
    return (f"[mask_compare] {summ['strains_compared']} strains: {summ['strains_up']} up, "
            f"{summ['strains_down']} down, {summ['strains_unchanged']} unchanged "
            f"(|change| < {summ['unchanged_eps']}); mean change {summ['mean_delta']:+.2f}; "
            f"{summ['n_dropped']} organism(s) only scored without the mask")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--masked", required=True, help="masked match_paths JSON")
    ap.add_argument("--unmasked", required=True, help="unmasked match_paths JSON (same sample)")
    ap.add_argument("-o", "--output", required=True, help="merged JSON (may equal --masked)")
    a = ap.parse_args()
    with open(a.masked) as fh:
        m = json.load(fh)
    with open(a.unmasked) as fh:
        u = json.load(fh)
    print(describe(merge(m, u, source=os.path.basename(a.unmasked))))
    with open(a.output, "w") as fh:
        json.dump(m, fh, separators=(",", ":"), indent=4)


if __name__ == "__main__":
    main()
