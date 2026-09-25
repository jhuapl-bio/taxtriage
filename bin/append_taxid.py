#!/usr/bin/env python3

##############################################################################################
# Copyright 2024 The Johns Hopkins University Applied Physics Laboratory LLC
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
import argparse
import pandas as pd
from pathlib import Path

try:
    # Shared stdlib-only helper (bin/ncbi_taxid.py). Available alongside this
    # script in the staged bin/ directory at runtime.
    from ncbi_taxid import fetch_taxids
except Exception:  # noqa: BLE001 - keep working even if helper is unavailable
    fetch_taxids = None

def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Map GCF accession to taxids using input and reference files.",
        epilog="Example: python3 map_gcf_to_taxid.py -i input.txt -r reference.txt -o output.tsv -c taxid",
    )
    parser.add_argument(
        "-i",
        "--file_in",
        metavar="FILE_IN",
        type=Path,
        required=True,
        help="Tabular input txt file",
    )
    parser.add_argument(
        "-r",
        "--ref_file",
        metavar="REF_FILE",
        type=Path,
        nargs='+',
        required=True,
        help="Reference file with GCF accession and taxid",
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="OUTPUT",
        type=Path,
        required=False,
        help="Output file to write the mapping",
    )
    parser.add_argument(
        "-c",
        "--column",
        metavar="COLUMN",
        type=str,
        required=False,
        default="taxid",
        help="Column name in the reference file to map to (e.g., taxid, species_taxid)",
    )
    parser.add_argument(
        "-b",
        "--ncbi-backup",
        dest="ncbi_backup",
        action="store_true",
        default=False,
        help=(
            "For accessions with no GCF/GCA -> taxid match in the reference "
            "assembly summaries, query NCBI E-utilities by accession to recover "
            "the taxid and write it into the map file."
        ),
    )
    parser.add_argument(
        "-e",
        "--email",
        metavar="EMAIL",
        type=str,
        required=False,
        default=None,
        help="Email for NCBI E-utilities querying (used with --ncbi-backup).",
    )
    parser.add_argument(
        "--api-key",
        dest="api_key",
        metavar="API_KEY",
        type=str,
        required=False,
        default=None,
        help="Optional NCBI API key for E-utilities querying.",
    )
    parser.add_argument(
        "--ncbi-batch-size",
        dest="ncbi_batch_size",
        metavar="N",
        type=int,
        required=False,
        default=100,
        help=(
            "Number of accessions resolved per E-utilities request pair during "
            "the NCBI backup lookup. Lower this if NCBI returns HTTP 429."
        ),
    )
    parser.add_argument(
        "--custom-map",
        dest="custom_map",
        metavar="CUSTOM_MAP",
        type=Path,
        required=False,
        default=None,
        help=(
            "TSV file with columns 'accession' and 'taxid' for manually mapping "
            "non-standard accessions absent from assembly summaries and NCBI. "
            "The 'accession' column must match the Acc (nuccore accession) in the input."
        ),
    )
    return parser.parse_args(argv)

def read_input_file(input_file):
    return pd.read_csv(input_file, sep='\t', header=None, names=["Acc", "Assembly", "Organism_Name", "Description"])

def read_reference_file(ref_file, wanted=None, keep_columns=None):
    """Read an NCBI assembly_summary_*.txt (refseq or genbank) robustly.

    The file is strictly tab-delimited, but occasional rows contain a stray
    extra (or missing) tab, so pandas' C parser aborts with e.g.
    "Error tokenizing data. C error: Expected 38 fields ..., saw 39".
    Fields such as asm_submitter also contain multiple consecutive spaces
    which must never be treated as delimiters, so parse manually: skip the
    leading comment line, take the column names from the header line, then
    split every data line on tab ONLY and pad/collapse each row to exactly the
    header width so no accession is lost.

    MEMORY: assembly_summary_genbank.txt is ~1.5 GB / ~2.7M rows. Materialising
    every row as a list of 38 Python strings (then a DataFrame) needs tens of GB
    and gets the task OOM-killed. So the file is STREAMED and only rows whose
    accession is in `wanted` are kept (all rows when `wanted` is None), and only
    `keep_columns` are retained (all columns when None).
    """
    with open(ref_file, 'r', newline='') as fh:
        # First line is a free-text comment ("# See ftp ..."); skip it.
        fh.readline()
        header_line = fh.readline().rstrip('\r\n')
        columns = header_line.split('\t')
        ncol = len(columns)
        acc_col = '#assembly_accession' if '#assembly_accession' in columns else columns[0]
        if keep_columns:
            out_cols = [acc_col] + [c for c in keep_columns if c in columns and c != acc_col]
        else:
            out_cols = columns
        idx = [columns.index(c) for c in out_cols]

        rows = []
        for line in fh:
            # Cheap pre-filter on the accession before splitting the whole line.
            if wanted is not None:
                acc = line.split('\t', 1)[0]
                if acc not in wanted:
                    continue
            if not line.strip():
                continue
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) < ncol:
                fields = fields + [''] * (ncol - len(fields))
            elif len(fields) > ncol:
                fields = fields[:ncol - 1] + ['\t'.join(fields[ncol - 1:])]
            rows.append([fields[i] for i in idx])

    df = pd.DataFrame(rows, columns=out_cols)
    if acc_col != '#assembly_accession':
        df = df.rename(columns={acc_col: '#assembly_accession'})
    return df

def _clean_taxid(value):
    """Render a taxid as a clean string ('2697049' not '2697049.0', '' for NaN)."""
    if value is None:
        return ""
    text = str(value).strip()
    if text == "" or text.lower() == "nan":
        return ""
    # pandas may coerce integer taxids to float (e.g. 2697049.0) when the
    # column contains NaNs; strip the spurious trailing decimal.
    if text.endswith(".0"):
        text = text[:-2]
    return text


def map_gcf_to_taxid(input_df, ref_df, column):
    # Create a dictionary from the reference file for mapping
    gcf_to_taxid = ref_df.set_index('#assembly_accession')[column].to_dict()

    # Map the GCF to the desired column (taxid or specified column), keeping the
    # result as clean object-typed strings so NaNs don't coerce taxids to float
    # and so the NCBI backup can assign string taxids without dtype conflicts.
    input_df['Mapped_Value'] = input_df['Assembly'].map(gcf_to_taxid).apply(_clean_taxid).astype(object)

    return input_df

def apply_custom_map(mapped_df, custom_map_file):
    """Apply a user-supplied accession->taxid mapping for non-standard accessions.

    Reads a TSV with columns 'accession' and 'taxid'. For any row where
    Mapped_Value is still empty AND the Acc appears in the custom map, fills in
    the taxid. Rows already resolved by the assembly summary lookup are untouched.

    Expected file format (tab-separated, with header):
        accession\\ttaxid
        OR833055.1\\t2697049
        CUSTOM_SEQ_001\\t12345
    """
    try:
        custom_df = pd.read_csv(custom_map_file, sep='\t', dtype=str)
    except Exception as e:
        print(f"WARNING: Could not read custom map file '{custom_map_file}': {e}. Skipping.")
        return mapped_df

    required_cols = {"accession", "taxid"}
    missing_cols = required_cols - set(custom_df.columns.str.lower())
    if missing_cols:
        print(
            f"WARNING: Custom map file is missing required column(s): {missing_cols}. "
            "Expected a TSV with headers 'accession' and 'taxid'. Skipping."
        )
        return mapped_df

    custom_df.columns = custom_df.columns.str.lower()
    custom_lookup = custom_df.set_index("accession")["taxid"].apply(_clean_taxid).to_dict()

    def _is_missing(value):
        if value is None:
            return True
        text = str(value).strip()
        return text == "" or text.lower() == "nan"

    missing_mask = mapped_df["Mapped_Value"].apply(_is_missing)
    resolved = 0
    for idx in mapped_df.index[missing_mask]:
        acc = str(mapped_df.at[idx, "Acc"]).strip()
        if acc in custom_lookup and custom_lookup[acc]:
            mapped_df.at[idx, "Mapped_Value"] = custom_lookup[acc]
            resolved += 1

    print(f"Custom map resolved {resolved} of {int(missing_mask.sum())} unmapped accession(s).")
    return mapped_df


def backfill_missing_taxids(mapped_df, email=None, api_key=None, batch_size=None):
    """Fill empty Mapped_Value taxids by querying NCBI per accession.

    Accessions that were not present in any assembly summary (e.g. a local
    reference FASTA such as OR833055.1) have no GCF -> taxid mapping, leaving
    Mapped_Value blank. Query NCBI directly by the nuccore accession (the 'Acc'
    column) to recover the taxid so it can be passed downstream in the map file.
    """
    if fetch_taxids is None:
        print("NCBI backup requested but ncbi_taxid helper is unavailable; skipping.")
        return mapped_df

    def _is_missing(value):
        if value is None:
            return True
        text = str(value).strip()
        return text == "" or text.lower() == "nan"

    missing_mask = mapped_df["Mapped_Value"].apply(_is_missing)
    n_missing = int(missing_mask.sum())
    if n_missing == 0:
        return mapped_df

    # Collect the unique accessions first and resolve them in batched, rate
    # limited requests rather than two E-utilities calls per row.
    wanted = []
    for idx in mapped_df.index[missing_mask]:
        acc = str(mapped_df.at[idx, "Acc"]).strip()
        if acc and acc.lower() != "nan":
            wanted.append(acc)
    unique_acc = list(dict.fromkeys(wanted))

    print(
        f"Attempting NCBI taxid backup lookup for {n_missing} unmapped accession(s) "
        f"({len(unique_acc)} unique)..."
    )
    cache = fetch_taxids(
        unique_acc, email=email, api_key=api_key, batch_size=batch_size
    )

    resolved = 0
    for idx in mapped_df.index[missing_mask]:
        acc = str(mapped_df.at[idx, "Acc"]).strip()
        taxid = cache.get(acc)
        if taxid:
            mapped_df.at[idx, "Mapped_Value"] = taxid
            resolved += 1
    print(f"NCBI taxid backup resolved {resolved} of {n_missing} accession(s).")
    return mapped_df


def main(argv=None):
    args = parse_args(argv)

    input_df = read_input_file(args.file_in)
    # Only the assemblies actually present in the input are needed, so stream
    # the (potentially multi-GB) summaries and keep just those rows + column.
    wanted = set(input_df['Assembly'].dropna().astype(str).str.strip())
    ref_df = pd.concat(
        [read_reference_file(f, wanted=wanted, keep_columns=[args.column]) for f in args.ref_file],
        ignore_index=True,
    )
    # RefSeq is passed first; keep its row if an accession appears twice so the
    # set_index().to_dict() lookup is deterministic.
    ref_df = ref_df.drop_duplicates(subset='#assembly_accession', keep='first')

    mapped_df = map_gcf_to_taxid(input_df, ref_df, args.column)
    if args.custom_map:
        mapped_df = apply_custom_map(mapped_df, args.custom_map)
    if args.ncbi_backup:
        mapped_df = backfill_missing_taxids(
            mapped_df,
            email=args.email,
            api_key=args.api_key,
            batch_size=args.ncbi_batch_size,
        )
    # Write the output to the specified file
    if args.output:
        mapped_df[['Acc', 'Assembly',  'Organism_Name', 'Description', 'Mapped_Value']].to_csv(args.output, sep='\t', index=False)
        print(f"Output written to {args.output}")
    else:
        print("Output sample, no output file specified:")
        print(mapped_df[['Acc', 'Assembly', 'Organism_Name', 'Description', 'Mapped_Value']].head())
if __name__ == "__main__":
    main()
