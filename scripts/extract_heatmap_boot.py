#!/usr/bin/env python3
"""
Reverse the bootstrap payload out of a generated TaxTriage report.

bin/make_report.py embeds its dataset in the HTML as:

    <script>
    window.HEATMAP_BOOT = {...};
    </script>

and scripts/inline_boot_json.py embeds the same payload for GitHub Pages as:

    <script id="BOOTSTRAP" type="application/json">{...}</script>

This script pulls either form back out (also works on assets/heatmap_boot.js
and assets/pages.js) and writes the original JSON.

Usage
-----
    python scripts/extract_heatmap_boot.py all.odr.html                 # -> all.odr.boot.json
    python scripts/extract_heatmap_boot.py all.odr.html -o boot.json
    python scripts/extract_heatmap_boot.py all.odr.html --stdout --minify
    python scripts/extract_heatmap_boot.py all.odr.html --info          # keys + sizes, no write
    python scripts/extract_heatmap_boot.py all.odr.html --key records -o records.json
    python scripts/extract_heatmap_boot.py all.odr.html --key records --csv records.csv
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

# ── extraction ─────────────────────────────────────────────────────────────

_BOOTSTRAP_RE = re.compile(
    r'<script[^>]*\bid=["\']BOOTSTRAP["\'][^>]*>(.*?)</script>',
    re.DOTALL | re.IGNORECASE,
)
_ASSIGN_RE = re.compile(r'(?:window|globalThis|self)\s*\.\s*HEATMAP_BOOT\s*=\s*')


def _match_braces(text: str, start: int) -> int:
    """
    Return the index just past the object literal that begins at ``text[start]``
    (which must be '{'). String literals, escapes and comments are skipped so
    braces inside data values do not confuse the scan.
    """
    if text[start] != "{":
        raise ValueError(f"expected '{{' at offset {start}, found {text[start]!r}")

    depth = 0
    i = start
    n = len(text)
    quote = None          # active string delimiter, or None
    while i < n:
        c = text[i]

        if quote:
            if c == "\\":
                i += 2
                continue
            if c == quote:
                quote = None
            i += 1
            continue

        if c in "\"'`":
            quote = c
        elif c == "/" and i + 1 < n and text[i + 1] == "/":
            i = text.find("\n", i)
            if i == -1:
                break
            continue
        elif c == "/" and i + 1 < n and text[i + 1] == "*":
            j = text.find("*/", i + 2)
            if j == -1:
                break
            i = j + 2
            continue
        elif c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
            if depth == 0:
                return i + 1
        i += 1

    raise ValueError("unbalanced braces — could not find the end of the payload")


def extract_payload_text(text: str) -> tuple[str, str]:
    """
    Return ``(raw_text, source_kind)`` for the boot payload found in ``text``.

    ``source_kind`` is 'bootstrap-json' (inlined <script id="BOOTSTRAP">) or
    'heatmap-boot-assign' (window.HEATMAP_BOOT = {...}).
    """
    m = _BOOTSTRAP_RE.search(text)
    if m and m.group(1).strip().startswith("{"):
        return m.group(1).strip(), "bootstrap-json"

    # Every `window.HEATMAP_BOOT` hit is tried: the report also contains reader
    # code such as `const BOOT = window.HEATMAP_BOOT || {};`, which is not an
    # assignment and is skipped by the `=` requirement in the regex, but a
    # defensive loop costs nothing.
    for am in _ASSIGN_RE.finditer(text):
        brace = text.find("{", am.end())
        if brace == -1 or text[am.end():brace].strip():
            continue  # not an object literal assignment
        end = _match_braces(text, brace)
        return text[brace:end], "heatmap-boot-assign"

    raise SystemExit(
        "ERROR: no bootstrap payload found. Expected a "
        '<script id="BOOTSTRAP"> block or a `window.HEATMAP_BOOT = {...}` assignment.'
    )


def parse_payload(raw: str) -> dict:
    """Parse the payload as JSON; fall back to Node for JS-literal quirks."""
    try:
        return json.loads(raw)
    except json.JSONDecodeError as exc:
        node = shutil.which("node")
        if not node:
            raise SystemExit(
                f"ERROR: payload is not strict JSON ({exc}) and Node.js is not on "
                "PATH to evaluate it as a JS literal."
            )
        script = "process.stdout.write(JSON.stringify((0,eval)(require('fs').readFileSync(0,'utf8'))))"
        res = subprocess.run(
            [node, "-e", script],
            input="(" + raw + ")",
            capture_output=True,
            text=True,
        )
        if res.returncode != 0:
            raise SystemExit(f"ERROR: node failed to evaluate the payload:\n{res.stderr}")
        return json.loads(res.stdout)


def extract(path: Path) -> tuple[dict, str]:
    text = path.read_text(encoding="utf-8", errors="replace")
    raw, kind = extract_payload_text(text)
    return parse_payload(raw), kind


# ── reporting helpers ──────────────────────────────────────────────────────

def describe(payload: dict) -> str:
    lines = []
    for key, val in payload.items():
        if isinstance(val, list):
            desc = f"list[{len(val)}]"
            if val and isinstance(val[0], dict):
                desc += f"  first-record keys: {len(val[0])}"
        elif isinstance(val, dict):
            desc = f"dict[{len(val)}]  keys: " + ", ".join(list(val)[:6])
            if len(val) > 6:
                desc += ", …"
        else:
            desc = repr(val)
            if len(desc) > 70:
                desc = desc[:67] + "…"
        lines.append(f"  {key:<24} {desc}")
    return "\n".join(lines)


def write_csv(records: list, out: Path) -> None:
    import csv

    if not records or not isinstance(records[0], dict):
        raise SystemExit("ERROR: --csv requires a list of records (e.g. --key records)")
    cols: list[str] = []
    seen = set()
    for rec in records:
        for k in rec:
            if k not in seen:
                seen.add(k)
                cols.append(k)
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for rec in records:
            w.writerow({k: ("" if v is None else v) for k, v in rec.items()})


# ── cli ────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("html", type=Path, help="report .html (or heatmap_boot.js / pages.js)")
    ap.add_argument("-o", "--output", type=Path,
                    help="output JSON path (default: <input>.boot.json)")
    ap.add_argument("--key", help="dump only this top-level key (e.g. records, sample_meta)")
    ap.add_argument("--csv", type=Path, help="also write the selected records as CSV")
    ap.add_argument("--minify", action="store_true",
                    help="write compact JSON exactly as make_report.py serialized it")
    ap.add_argument("--indent", type=int, default=2, help="indent for pretty output (default 2)")
    ap.add_argument("--stdout", action="store_true", help="print JSON instead of writing a file")
    ap.add_argument("--info", action="store_true", help="summarize the payload; write nothing")
    args = ap.parse_args()

    if not args.html.exists():
        sys.exit(f"ERROR: not found: {args.html}")

    payload, kind = extract(args.html)

    if args.info:
        print(f"{args.html}  ({kind})")
        print(describe(payload))
        return

    data = payload
    if args.key:
        if args.key not in payload:
            sys.exit(
                f"ERROR: key {args.key!r} not in payload. Available: "
                + ", ".join(payload)
            )
        data = payload[args.key]

    if args.minify:
        text = json.dumps(data, ensure_ascii=False, allow_nan=False, separators=(",", ":"))
    else:
        text = json.dumps(data, ensure_ascii=False, allow_nan=False, indent=args.indent)

    if args.stdout:
        sys.stdout.write(text + "\n")
    else:
        out = args.output or args.html.with_suffix("").with_suffix(".boot.json")
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(text + "\n", encoding="utf-8")
        print(f"✓ {out}  ({len(text):,} chars, source: {kind})")

    if args.csv:
        write_csv(data, args.csv)
        print(f"✓ {args.csv}  ({len(data):,} rows)")


if __name__ == "__main__":
    main()
