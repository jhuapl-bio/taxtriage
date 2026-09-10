#!/usr/bin/env python3
"""
Normalise a user-supplied spike-in sheet into the TSV the pipeline consumes.

The sheet says WHICH organisms to spike into the background, HOW MANY reads of
each, and HOW MANY replicate datasets to build at each level. It is deliberately
forgiving about column naming and file format (csv / tsv / txt / xlsx), because
it is a sheet people hand-edit.

Accepted columns (case-insensitive, first match wins):

    accession    accession | assembly | assembly_accession | nuccore | reference | genome | id
    count        count | reads | n_reads | nreads | n | spike | spike_count | amount
    replicates   replicates | replicate | reps | n_reps | repeats        (optional, default 1)
    level        level | series | group | tier                           (optional)
    name         name | organism | label | description                   (optional, cosmetic)

Series semantics
----------------
Each row is one organism at one spike level. Rows are grouped into LEVELS, and
one dataset is built per (level, replicate):

  * With no `level` column, the level IS the count — the common case, where every
    organism is spiked at the same set of amounts:

        accession,count,replicates
        GCF_014621545.1,100,3        -> level 100:  MPXV 100 reads
        GCF_014621545.1,1000,3       -> level 1000: MPXV 1000 reads

  * With a `level` column, organisms in the same level may be spiked at DIFFERENT
    amounts, which is how you build a realistic mixed panel:

        accession,count,replicates,level
        GCF_014621545.1,500,3,low     -> level "low": MPXV 500 + VACV 100
        GCF_000859985.2,100,3,low
        GCF_014621545.1,5000,3,high   -> level "high": MPXV 5000 + VACV 1000
        GCF_000859985.2,1000,3,high

    A non-numeric level label is mapped to the level's TOTAL spiked reads for the
    dataset id (c<total>), because the report's dataset-id grammar only carries an
    integer there and the total is the meaningful x-axis value.

Replicates for a level = the maximum `replicates` seen on its rows.

Output (TSV, one row per accession x level):
    level_key  level_count  accession  count  replicates  name
and, with --accessions, the distinct accession list one per line.
"""

import argparse
import csv
import os
import re
import sys
from collections import OrderedDict, defaultdict

# ── column synonyms ──────────────────────────────────────────────────────────
COLS = {
    "accession": ["accession", "assembly", "assembly_accession", "nuccore",
                  "reference", "genome", "id", "acc"],
    "count": ["count", "reads", "n_reads", "nreads", "n", "spike",
              "spike_count", "amount", "spiked_reads"],
    "replicates": ["replicates", "replicate", "reps", "n_reps", "repeats", "rep"],
    "level": ["level", "series", "group", "tier", "mix"],
    "name": ["name", "organism", "label", "description", "taxon"],
}

# GCF_/GCA_ assemblies and nuccore accessions (NC_045512.2, CP012345, U00096.3…)
ACC_RE = re.compile(r"^(?:GC[AF]_\d+\.\d+|[A-Z]{1,4}_?\d{5,}(?:\.\d+)?)$", re.I)


def norm_header(h):
    return re.sub(r"[^a-z0-9]+", "_", str(h or "").strip().lower()).strip("_")


def read_rows(path):
    """Return a list of dicts from csv/tsv/txt/xlsx, keys normalised."""
    ext = os.path.splitext(path)[1].lower()
    if ext in (".xlsx", ".xlsm", ".xltx"):
        try:
            from openpyxl import load_workbook
        except ImportError:
            raise SystemExit(
                "ERROR: reading an .xlsx spike-in sheet needs openpyxl. "
                "Install it, or save the sheet as .csv / .tsv."
            )
        wb = load_workbook(path, read_only=True, data_only=True)
        ws = wb[wb.sheetnames[0]]
        rows = []
        header = None
        for r in ws.iter_rows(values_only=True):
            if r is None:
                continue
            vals = ["" if v is None else str(v).strip() for v in r]
            if not any(vals):
                continue
            if header is None:
                header = [norm_header(v) for v in vals]
                continue
            rows.append({header[i] if i < len(header) else f"col{i}": vals[i]
                         for i in range(len(vals))})
        wb.close()
        return rows

    with open(path, newline="", encoding="utf-8-sig") as fh:
        sample = fh.read(8192)
        fh.seek(0)
        # Sniffing beats guessing from the extension: people save .csv as TSV.
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",\t;|")
            delim = dialect.delimiter
        except Exception:
            delim = "\t" if ext in (".tsv", ".tab") else ","
        reader = csv.reader(fh, delimiter=delim)
        rows, header = [], None
        for raw in reader:
            vals = [str(v).strip() for v in raw]
            if not any(vals):
                continue
            if vals[0].startswith("#"):
                continue
            if header is None:
                header = [norm_header(v) for v in vals]
                continue
            rows.append({header[i] if i < len(header) else f"col{i}": vals[i]
                         for i in range(len(vals))})
        return rows


def pick(row, key):
    for cand in COLS[key]:
        if cand in row and str(row[cand]).strip() != "":
            return str(row[cand]).strip()
    return ""


def as_int(v, what, lineno):
    s = str(v).strip().replace(",", "").replace("_", "")
    # tolerate "1e4" and "1000.0" from spreadsheet cells
    try:
        f = float(s)
    except ValueError:
        raise SystemExit(f"ERROR: row {lineno}: {what} {v!r} is not a number")
    if f != int(f):
        raise SystemExit(f"ERROR: row {lineno}: {what} {v!r} must be a whole number")
    return int(f)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="Spike-in sheet (csv/tsv/xlsx)")
    ap.add_argument("--output", required=True, help="Normalised TSV to write")
    ap.add_argument("--accessions", default=None,
                    help="Also write the distinct accession list here, one per line")
    ap.add_argument("--default-replicates", type=int, default=1,
                    help="Replicates for rows with no replicates column (default 1)")
    args = ap.parse_args()

    rows = read_rows(args.input)
    if not rows:
        raise SystemExit(f"ERROR: no data rows found in {args.input!r}")

    header_seen = set(rows[0].keys())
    if not any(c in header_seen for c in COLS["accession"]):
        raise SystemExit(
            "ERROR: the spike-in sheet needs an accession column "
            f"(one of: {', '.join(COLS['accession'])}). Found: {', '.join(sorted(header_seen))}"
        )
    if not any(c in header_seen for c in COLS["count"]):
        raise SystemExit(
            "ERROR: the spike-in sheet needs a count column "
            f"(one of: {', '.join(COLS['count'])}). Found: {', '.join(sorted(header_seen))}"
        )

    # ── parse + validate ─────────────────────────────────────────────────────
    parsed = []
    for i, row in enumerate(rows, start=2):   # +1 header, +1 for 1-based
        acc = pick(row, "accession")
        if not acc:
            continue
        cnt = pick(row, "count")
        if not cnt:
            raise SystemExit(f"ERROR: row {i}: accession {acc!r} has no count")
        count = as_int(cnt, "count", i)
        if count <= 0:
            raise SystemExit(f"ERROR: row {i}: count for {acc!r} must be > 0")
        reps_s = pick(row, "replicates")
        reps = as_int(reps_s, "replicates", i) if reps_s else args.default_replicates
        if reps <= 0:
            raise SystemExit(f"ERROR: row {i}: replicates for {acc!r} must be > 0")
        lvl = pick(row, "level")
        if not ACC_RE.match(acc):
            # A warning, not an error: NCBI accession shapes change, and the
            # download step is the real authority on whether it resolves.
            print(f"[parse_spikein] WARNING: row {i}: {acc!r} does not look like an "
                  f"assembly or nuccore accession; passing it through anyway",
                  file=sys.stderr)
        parsed.append({
            "accession": acc,
            "count": count,
            "replicates": reps,
            "level": lvl,
            "name": pick(row, "name"),
            "lineno": i,
        })

    if not parsed:
        raise SystemExit(f"ERROR: {args.input!r} has a header but no usable rows")

    # ── group into levels ────────────────────────────────────────────────────
    # No level column -> the count IS the level, so each amount forms its own
    # dataset. With a level column, rows sharing a label form one mixed dataset.
    groups = OrderedDict()
    for p in parsed:
        key = p["level"] if p["level"] else str(p["count"])
        groups.setdefault(key, []).append(p)

    # One accession must not appear twice in the same level — that would silently
    # double its spike amount.
    for key, members in groups.items():
        seen = defaultdict(list)
        for m in members:
            seen[m["accession"].upper()].append(m["lineno"])
        for acc, lines in seen.items():
            if len(lines) > 1:
                raise SystemExit(
                    f"ERROR: accession {acc} appears {len(lines)} times in level "
                    f"{key!r} (rows {', '.join(map(str, lines))}). Give each organism "
                    "one row per level, or split them into different levels."
                )

    out_rows = []
    used_level_counts = {}
    for key, members in groups.items():
        total = sum(m["count"] for m in members)
        reps = max(m["replicates"] for m in members)
        # The dataset id carries an integer level; use the label when it is numeric
        # (the level IS that amount) and the level's total spiked reads otherwise.
        level_count = int(key) if re.fullmatch(r"\d+", str(key)) else total
        if level_count in used_level_counts and used_level_counts[level_count] != key:
            raise SystemExit(
                f"ERROR: levels {used_level_counts[level_count]!r} and {key!r} both map to "
                f"c{level_count} in the dataset id. Rename one, or change its counts so the "
                "level totals differ."
            )
        used_level_counts[level_count] = key
        for m in members:
            out_rows.append({
                "level_key": key,
                "level_count": level_count,
                "accession": m["accession"],
                "count": m["count"],
                "replicates": reps,
                "name": m["name"],
            })

    out_rows.sort(key=lambda r: (r["level_count"], r["accession"]))

    cols = ["level_key", "level_count", "accession", "count", "replicates", "name"]
    with open(args.output, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in out_rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    if args.accessions:
        seen = OrderedDict()
        for r in out_rows:
            seen.setdefault(r["accession"], 1)
        with open(args.accessions, "w") as fh:
            for a in seen:
                fh.write(a + "\n")

    n_levels = len(groups)
    n_acc = len({r["accession"] for r in out_rows})
    n_datasets = sum(max(m["replicates"] for m in members) for members in groups.values())
    print(f"[parse_spikein] {n_acc} accession(s) x {n_levels} level(s) "
          f"-> {n_datasets} dataset(s); levels: "
          + ", ".join(f"{k}(c{used_level_counts_inv})" for k, used_level_counts_inv in
                      ((k, next(r['level_count'] for r in out_rows if r['level_key'] == k))
                       for k in groups)),
          file=sys.stderr)


if __name__ == "__main__":
    main()
