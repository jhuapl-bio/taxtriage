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
"""

import argparse
import glob
import json
import os
import sys


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
        "version": "1.0",
        "samples": samples,
    }

    with open(args.output, "w", encoding="utf-8") as fh:
        json.dump(combined, fh, ensure_ascii=False, allow_nan=False, separators=(",", ":"))

    print(f"[combine_samples_json] Written {len(samples)} sample(s) -> {args.output!r}")


if __name__ == "__main__":
    main()
