#!/usr/bin/env python3
"""
Convert assets/pathogen_sheet.csv from the TaxTriage pipeline repo into the
compact JSON payload consumed by docs/pathogen-sheet.md.

Usage:
    python scripts/build_pathogen_data.py [--csv PATH] [--out PATH]

Output format is columnar-ish: a `cols` header list plus `rows` as arrays,
which avoids repeating 20 key names 1700 times (~40% smaller than a list of
objects, before gzip).

Known upstream data issues are normalised here rather than in the browser so
the fix is visible and reviewable. The source CSV is never modified.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path

COLUMNS = [
    "name",
    "taxid",
    "general_classification",
    "alternative_names",
    "pathogenic_sites",
    "commensal_sites",
    "status",
    "high_consequence",
    "pathology",
    "host_organism",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "mol_type",
    "reference",
    "Additional references",
    "assembly_accession",
]

# Multi-valued columns: comma-separated lists that must be split for faceting.
LIST_COLUMNS = {"pathogenic_sites", "commensal_sites", "alternative_names"}

# Upstream typos -> canonical value. Reported at build time so they can be
# fixed at the source eventually.
VALUE_FIXES = {
    "status": {
        "estbalished": "established",
        "etsablished": "established",
    },
}


def clean(value: str) -> str:
    return (value or "").strip()


def split_list(value: str) -> list[str]:
    return [p.strip() for p in clean(value).split(",") if p.strip()]


def main() -> int:
    here = Path(__file__).resolve().parent.parent
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", type=Path, default=here / "assets" / "pathogen_sheet.csv")
    ap.add_argument("--out", type=Path, default=here / "docs" / "data" / "pathogen_sheet.json")
    args = ap.parse_args()

    if not args.csv.is_file():
        print(f"error: CSV not found at {args.csv}", file=sys.stderr)
        print("hint: pass --csv /path/to/taxtriage/assets/pathogen_sheet.csv", file=sys.stderr)
        return 1

    with args.csv.open(newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        missing = [c for c in COLUMNS if c not in (reader.fieldnames or [])]
        if missing:
            print(f"error: CSV is missing expected columns: {missing}", file=sys.stderr)
            return 1
        raw_rows = list(reader)

    applied: Counter[str] = Counter()
    rows: list[list] = []

    for r in raw_rows:
        out_row = []
        for col in COLUMNS:
            val = clean(r.get(col, ""))
            fixes = VALUE_FIXES.get(col)
            if fixes and val in fixes:
                applied[f"{col}: {val!r} -> {fixes[val]!r}"] += 1
                val = fixes[val]
            if col in LIST_COLUMNS:
                out_row.append(split_list(val))
            elif col == "high_consequence":
                out_row.append(val.upper() == "TRUE")
            else:
                out_row.append(val)
        rows.append(out_row)

    # Drop fully-empty rows (the sheet has a stray blank entry).
    rows = [r for r in rows if r[COLUMNS.index("name")]]

    payload = {
        "cols": COLUMNS,
        "rows": rows,
        "source": "https://github.com/jhuapl-bio/taxtriage/blob/main/assets/pathogen_sheet.csv",
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, separators=(",", ":")), encoding="utf-8")

    size_kb = args.out.stat().st_size / 1024
    print(f"wrote {len(rows)} rows ({size_kb:.0f} KB) to {args.out}")
    if applied:
        print("normalised upstream values:")
        for k, n in sorted(applied.items()):
            print(f"  {k}  ({n} row{'s' if n != 1 else ''})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
