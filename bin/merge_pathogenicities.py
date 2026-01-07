#!/usr/bin/env python3
"""
merge_pathogenicities.py

Merge DIAMOND BLAST outfmt 6 TSV with metadata TSV (AMR or VFs),
optionally keep best hit per contig+group, and output:
  - summary TSV (default)
  - optional per-hit TSV (off by default)

Key features:
- modes: amr | vfs | auto
- default: best-hit filter per (qseqid, group_col) ON
- can disable with --disable-besthit-filter
- drops bulky sequence columns from metadata and outputs
"""

import argparse
import sys
import pandas as pd


BLAST6_DEFAULT_COLS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore"
]

# Add any other large columns you want removed here
DEFAULT_DROP_SEQ_COLS = [
    "sequence", "sequence_protein", "seq", "protein_seq", "dna_seq",
    "nucleotide_sequence", "amino_acid_sequence"
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Merge AMR/VF BLAST hits with metadata and produce summaries.")

    p.add_argument("--blast", required=True, help="DIAMOND BLAST outfmt 6 TSV (no header).")
    p.add_argument("--meta", required=True, help="Metadata TSV.")
    p.add_argument("--out", required=True, help="Output summary TSV (default output).")

    p.add_argument("--mode", choices=["auto", "amr", "vfs"], default="auto",
                   help="Mode controls join key extraction defaults (default: auto).")

    # BLAST input
    p.add_argument("--blast-cols", default=",".join(BLAST6_DEFAULT_COLS),
                   help="Comma-separated BLAST column names (default: standard outfmt 6, 12 cols).")
    p.add_argument("--blast-key-col", default="sseqid",
                   help="BLAST column to extract join key from (default: sseqid).")

    # Overrides for join behavior
    p.add_argument("--split-delim", default=None,
                   help="Delimiter used to extract join key from BLAST key field (override).")
    p.add_argument("--meta-key", default=None,
                   help="Metadata column name used as join key (override).")

    # Best-hit filter behavior
    p.add_argument("--disable-besthit-filter", action="store_true",
                   help="Disable default best-hit filtering per contig+group.")
    p.add_argument("--group-hit-col", default=None,
                   help=("Column used to group hits for best-hit filtering (override). "
                         "Defaults: vfs='VFC ID' if present else 'vf_label'; amr='aroid'."))

    # Summary behavior
    p.add_argument("--group-col", default=None,
                   help=("Column to group summary by (override). "
                         "Default: organism if present, else species, else taxonid."))
    p.add_argument("--top-n", type=int, default=0,
                   help="If >0, keep only top N groups by hit count.")

    # Optional per-hit output
    p.add_argument("--per-hit-out", default=None,
                   help="If set, write per-hit merged TSV (after optional besthit filtering).")

    # Drop bulky columns
    p.add_argument("--drop-cols", default=",".join(DEFAULT_DROP_SEQ_COLS),
                   help="Comma-separated columns to drop from metadata/outputs if present.")
    p.add_argument(
        "--besthit-scope",
        choices=["none", "contig", "contig_taxon", "contig_vfc", "contig_label"],
        default="contig_vfc",
        help="Best-hit filtering scope (default depends on mode)."
    )
    p.add_argument("--taxon-col", default="taxonid", help="Taxon column for contig_taxon scope (default: taxonid).")
    p.add_argument("--vfc-col", default="VFC ID", help="VFC column for contig_vfc scope (default: VFC ID).")
    p.add_argument("--label-col", default="vf_label", help="Label column for contig_label scope (default: vf_label).")

    return p.parse_args()


def read_blast(path: str, colnames: list[str]) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, names=colnames, dtype=str)
    # numeric coercion for sorting/filtering
    for c in ("evalue", "bitscore"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def read_meta(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    df.columns = df.columns.str.strip()
    return df


def drop_cols_if_present(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    present = [c for c in cols if c in df.columns]
    return df.drop(columns=present, errors="ignore")


def choose_defaults(mode: str, meta_cols: list[str]) -> dict:
    """
    Returns dict with:
      meta_key, split_delim, group_hit_col, summary_group_col
    """
    def pick_summary_group():
        for c in ["organism", "Organism", "species", "Species", "taxonid", "TaxonID", "taxid"]:
            if c in meta_cols:
                return c
        return None

    if mode == "amr":
        return {
            "meta_key": "Protein Accession" if "Protein Accession" in meta_cols else None,
            "split_delim": "_",
            "group_hit_col": "aroid" if "aroid" in meta_cols else ("AMR Gene Family" if "AMR Gene Family" in meta_cols else None),
            "summary_group_col": pick_summary_group(),
        }

    if mode == "vfs":
        # VFDB-ish defaults
        return {
            "meta_key": "VFG" if "VFG" in meta_cols else None,
            "split_delim": "|",
            "group_hit_col": "VFC ID" if "VFC ID" in meta_cols else ("vf_label" if "vf_label" in meta_cols else None),
            "summary_group_col": pick_summary_group(),
        }

    # auto: infer from presence of VF columns vs AMR columns
    if "VFG" in meta_cols:
        return choose_defaults("vfs", meta_cols)
    if "Protein Accession" in meta_cols:
        return choose_defaults("amr", meta_cols)

    # fallback
    return {"meta_key": None, "split_delim": "_", "group_hit_col": None, "summary_group_col": pick_summary_group()}


def unique_join(series: pd.Series, sep: str = ";") -> str:
    vals = series.dropna().astype(str)
    if vals.empty:
        return ""
    return sep.join(pd.unique(vals))


def besthit_filter(df, scope, taxon_col="taxonid", vfc_col="VFC ID", label_col="vf_label"):
    if scope == "none":
        return df

    if scope == "contig":
        keys = ["qseqid"]
    elif scope == "contig_taxon":
        if taxon_col not in df.columns:
            print(f"WARNING: taxon-col '{taxon_col}' not found; skipping besthit.", file=sys.stderr)
            return df
        keys = ["qseqid", taxon_col]
    elif scope == "contig_vfc":
        if vfc_col not in df.columns:
            print(f"WARNING: vfc-col '{vfc_col}' not found; skipping besthit.", file=sys.stderr)
            return df
        keys = ["qseqid", vfc_col]
    elif scope == "contig_label":
        if label_col not in df.columns:
            print(f"WARNING: label-col '{label_col}' not found; skipping besthit.", file=sys.stderr)
            return df
        keys = ["qseqid", label_col]
    else:
        return df

    # Ensure numeric
    df["evalue"] = pd.to_numeric(df["evalue"], errors="coerce")
    df["bitscore"] = pd.to_numeric(df["bitscore"], errors="coerce")

    df_sorted = df.sort_values(
        by=keys + ["evalue", "bitscore"],
        ascending=[True]*len(keys) + [True, False],
        na_position="last",
        kind="mergesort",
    )
    return df_sorted.drop_duplicates(subset=keys, keep="first")


def main() -> int:
    args = parse_args()

    blast_cols = [c.strip() for c in args.blast_cols.split(",") if c.strip()]
    blast_df = read_blast(args.blast, blast_cols)

    if args.blast_key_col not in blast_df.columns:
        print(f"ERROR: blast-key-col '{args.blast_key_col}' not in BLAST columns.", file=sys.stderr)
        print("BLAST columns:", blast_df.columns.tolist(), file=sys.stderr)
        return 2

    meta_df = read_meta(args.meta)

    # Drop bulky sequence columns from metadata early
    drop_cols = [c.strip() for c in args.drop_cols.split(",") if c.strip()]
    meta_df = drop_cols_if_present(meta_df, drop_cols)

    defaults = choose_defaults(args.mode, meta_df.columns.tolist())

    meta_key = args.meta_key if args.meta_key else defaults["meta_key"]
    split_delim = args.split_delim if args.split_delim else defaults["split_delim"]
    group_hit_col = args.group_hit_col if args.group_hit_col else defaults["group_hit_col"]
    summary_group_col = args.group_col if args.group_col else defaults["summary_group_col"]

    if not meta_key or meta_key not in meta_df.columns:
        print("ERROR: Could not determine metadata join key column. Provide --meta-key.", file=sys.stderr)
        print("Metadata columns:", meta_df.columns.tolist(), file=sys.stderr)
        return 2

    if not group_hit_col:
        print("ERROR: Could not determine group-hit-col for besthit filtering. Provide --group-hit-col.", file=sys.stderr)
        print("Metadata columns:", meta_df.columns.tolist(), file=sys.stderr)
        return 2

    if not summary_group_col:
        print("ERROR: Could not determine summary-group-col. Provide --summary-group-col.", file=sys.stderr)
        return 2

    # Build join key
    blast_df["_join_key"] = blast_df[args.blast_key_col].astype(str).str.split(split_delim, n=1).str[0]
    meta_df["_join_key"] = meta_df[meta_key].astype(str)

    merged = blast_df.merge(meta_df, on="_join_key", how="left")

    # Best-hit filtering default ON
    if not args.disable_besthit_filter:
        merged = besthit_filter(merged,
            scope=args.besthit_scope,
            taxon_col=args.taxon_col,
            vfc_col=args.vfc_col,
            label_col=args.label_col,
        )

    # Clean up helper + ensure no bulky cols remain
    merged = merged.drop(columns=["_join_key"], errors="ignore")
    merged = drop_cols_if_present(merged, drop_cols)

    # Optional per-hit output
    if args.per_hit_out:
        merged.to_csv(args.per_hit_out, sep="\t", index=False)
        print(f"Wrote per-hit TSV: {args.per_hit_out}", file=sys.stderr)

    # Summary output
    # Build useful summary fields depending on mode / available columns
    candidates = [
        ("contigs", ("qseqid", lambda s: len(pd.unique(s).tolist()))),
        ("hit_count", ("qseqid", "count")),
        ("besthit_groups", (group_hit_col, lambda s: len(pd.unique(s).tolist()) if group_hit_col in s else 0)),
    ]

    # add common VF/AMR descriptors if present
    for col, outname in [
        ("vf_label", "vf_labels"),
        ("VFG", "vfg_ids"),
        ("VFC ID", "vfc_ids"),
        ("Protein Accession", "protein_accessions"),
        ("AMR Gene Family", "amr_gene_families"),
        ("Drug Class", "drug_classes"),
        ("Resistance Mechanism", "resistance_mechanisms"),
        ("aroid", "aroids"),
    ]:
        if col in merged.columns:
            candidates.append((outname, (col, lambda s: unique_join(s))))

    agg_kwargs = {name: spec for name, spec in candidates}

    summary = (
        merged.groupby(summary_group_col, dropna=False)
        .agg(**agg_kwargs)
        .reset_index()
        .rename(columns={summary_group_col: "organism"})
        .sort_values("hit_count", ascending=False)
    )

    if args.top_n and args.top_n > 0:
        summary = summary.head(args.top_n)

    summary.to_csv(args.out, sep="\t", index=False)

    print(f"Wrote summary TSV: {args.out}", file=sys.stderr)
    print(f"Mode: {args.mode} | join: meta_key={meta_key}, split='{split_delim}' | besthit group={group_hit_col}", file=sys.stderr)
    print(f"Rows: BLAST={len(blast_df)} META={len(meta_df)} MERGED(post)={len(merged)} SUMMARY={len(summary)}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
