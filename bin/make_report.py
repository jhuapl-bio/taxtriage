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

import pandas as pd
import json
import argparse
# Sanitize payload so json.dumps(..., allow_nan=False) won't fail.
import math


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and transform a tabular samplesheet.",
        epilog="Example: python check_samplesheet.py samplesheet.csv samplesheet.valid.csv",
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        default="all.organisms.report.txt",
        metavar="INPUT",
        help="TXT PATHS TaxTriage File to process",
    )
    parser.add_argument(
        "-t",
        "--template",
        metavar="TEMPLATE",
        default="heatmap.html",
        help="Input HTML template file",
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="OUTPUT",
        default="heatmap.personalized.html",
        help="Output HTML file",
    )


    return parser.parse_args(argv)

def main():
  args = parse_args()

  template = args.template

  # 1) Read TSV into DataFrame
  df = pd.read_csv(args.input, sep='\t')
  # df = pd.read_csv('all.organisms.txt', sep='\t')
  df['Detected Organism'] = df['Detected Organism'].str.replace('°', '', regex=False)
  # trim whitespace from Detected Organism column
  df['Detected Organism'] = df['Detected Organism'].str.strip()
  # drop Index nad indx cols
  df = df.drop(columns=['Index', 'index'])

  # 2) Prepare JSON and column lists
  # convert any NaN to None for JSON serialization
  df = df.where(pd.notnull(df), None)
  records = df.to_dict(orient='records')
  all_cols = df.columns.tolist()
  numeric_cols = df.select_dtypes(include='number').columns.tolist()

  # Build payload
  payload = {
      "records": records,
      "all_cols": all_cols,
      "numeric_cols": list(numeric_cols),
  }



  def sanitize_for_json(obj):
      """Recursively replace NaN/Inf/pd.NA with None so JSON serialization with
      allow_nan=False succeeds."""
      if isinstance(obj, dict):
          return {k: sanitize_for_json(v) for k, v in obj.items()}
      if isinstance(obj, list):
          return [sanitize_for_json(v) for v in obj]
      # pandas/numpy missing values
      try:
          if pd.isna(obj):
              return None
      except Exception:
          pass
      # floats that are not finite
      if isinstance(obj, float):
          return obj if math.isfinite(obj) else None
      return obj

  payload_clean = sanitize_for_json(payload)

  # Strict JSON (no NaN/Infinity) + safe for embedding in <script>
  bootstrap_json = json.dumps(payload_clean, ensure_ascii=False, allow_nan=False)
  bootstrap_json = bootstrap_json.replace("</", "<\\/")  # avoid </script> closing the tag

  # 2) Read template and substitute
  with open(template, "r", encoding="utf-8") as f:
      tpl = f.read()

  html = tpl.replace("__BOOTSTRAP_JSON__", bootstrap_json)

  # 3) Write the customized one
  output_path = args.output # "heatmap.personalized.html"
  with open(output_path, "w", encoding="utf-8") as f:
      f.write(html)


  print(f"Generated interactive heatmap with controls: {output_path}")

if __name__ == "__main__":
    main()
