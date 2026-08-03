#!/usr/bin/env python3

from __future__ import annotations

import argparse
import time
import urllib.error
import urllib.request
from pathlib import Path


ASSEMBLY_SUMMARY_REFSEQ_URL = (
    "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/"
    "assembly_summary_refseq.txt"
)

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------
# General utility functions
# ---------------------------------------------------------------------

def numeric(series: pd.Series) -> pd.Series:
    """Convert an NCBI text column to numeric values."""
    return pd.to_numeric(
        series.replace({"na": np.nan, "": np.nan}),
        errors="coerce",
    )


def normalize_accession(value) -> str | None:
    """Normalize an assembly accession read from a CSV file."""
    if pd.isna(value):
        return None

    text = str(value).strip()
    return text or None


# ---------------------------------------------------------------------
# Assembly-summary processing
# ---------------------------------------------------------------------

def read_assembly_summary(path: Path) -> pd.DataFrame:
    """Read an NCBI assembly_summary_refseq.txt file."""
    header = None

    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            possible_header = line.rstrip("\n").removeprefix("#").lstrip()

            if possible_header.startswith("assembly_accession\t"):
                header = possible_header.split("\t")
                break

    if header is None:
        raise ValueError(
            f"Could not locate the assembly-summary header in {path}"
        )

    df = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        names=header,
        dtype=str,
        low_memory=False,
    )

    required = {
        "assembly_accession",
        "assembly_level",
        "refseq_category",
        "genome_size",
        "contig_count",
        "replicon_count",
        "seq_rel_date",
        "asm_submitter",
        "ftp_path",
    }

    missing = required.difference(df.columns)
    if missing:
        raise ValueError(
            "The assembly summary is missing required columns: "
            + ", ".join(sorted(missing))
        )

    for column in [
        "genome_size",
        "contig_count",
        "replicon_count",
    ]:
        df[column] = numeric(df[column])

    df["assembly_accession"] = (
        df["assembly_accession"].astype("string").str.strip()
    )

    return df


# ---------------------------------------------------------------------
# Assembly-report processing
# ---------------------------------------------------------------------

def assembly_report_url(ftp_path: str | None) -> str | None:
    """Construct an NCBI assembly-report URL from an assembly FTP path."""
    if ftp_path is None or pd.isna(ftp_path):
        return None

    ftp_path = str(ftp_path).strip()
    if not ftp_path or ftp_path == "na":
        return None

    https_path = ftp_path.replace("ftp://", "https://")
    assembly_directory = https_path.rstrip("/").split("/")[-1]

    return (
        f"{https_path.rstrip('/')}/"
        f"{assembly_directory}_assembly_report.txt"
    )


def download_file(
    url: str,
    output_path: Path,
    timeout: int = 60,
    retries: int = 3,
    force: bool = False,
) -> None:
    """Download an NCBI metadata file with caching and retry handling."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if (
        output_path.exists()
        and output_path.stat().st_size > 0
        and not force
    ):
        return

    request = urllib.request.Request(
        url,
        headers={"User-Agent": "assembly-metadata-enrichment/1.0"},
    )

    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                output_path.write_bytes(response.read())
            return
        except (
            urllib.error.HTTPError,
            urllib.error.URLError,
            TimeoutError,
        ):
            if attempt == retries:
                raise
            time.sleep(2 ** attempt)


def read_assembly_report(path: Path) -> pd.DataFrame:
    """Read the sequence table from an NCBI assembly report."""
    column_names = None
    records: list[list[str]] = []

    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")

            if line.startswith("# Sequence-Name"):
                column_names = line.removeprefix("#").lstrip().split("\t")
                continue

            if line.startswith("#") or not line:
                continue

            if column_names is not None:
                fields = line.split("\t")
                if len(fields) == len(column_names):
                    records.append(fields)

    if column_names is None:
        raise ValueError(
            f"Could not find sequence columns in assembly report {path}"
        )

    report = pd.DataFrame(records, columns=column_names)

    required = {"Sequence-Role", "Sequence-Length"}
    missing = required.difference(report.columns)
    if missing:
        raise ValueError(
            f"Assembly report {path} is missing columns: "
            + ", ".join(sorted(missing))
        )

    report["Sequence-Length"] = pd.to_numeric(
        report["Sequence-Length"],
        errors="coerce",
    )

    return report


def summarize_sequence_roles(
    report: pd.DataFrame,
    assembly_size: float,
) -> dict[str, float | int]:
    """Calculate unplaced and unlocalized scaffold counts and percentages."""
    roles = report["Sequence-Role"].fillna("").str.strip().str.lower()

    unplaced = report.loc[roles.eq("unplaced-scaffold")]
    unlocalized = report.loc[roles.eq("unlocalized-scaffold")]

    unplaced_bp = float(unplaced["Sequence-Length"].sum())
    unlocalized_bp = float(unlocalized["Sequence-Length"].sum())

    if np.isfinite(assembly_size) and assembly_size > 0:
        unplaced_percentage = 100.0 * unplaced_bp / assembly_size
        unlocalized_percentage = 100.0 * unlocalized_bp / assembly_size
    else:
        unplaced_percentage = np.nan
        unlocalized_percentage = np.nan

    return {
        "unplaced_scaffold_count": int(len(unplaced)),
        "unplaced_scaffold_percentage": unplaced_percentage,
        "unlocalized_scaffold_count": int(len(unlocalized)),
        "unlocalized_scaffold_percentage": unlocalized_percentage,
    }


def get_assembly_report_metrics(
    assembly_row: pd.Series,
    report_directory: Path,
) -> dict[str, object]:
    """Download one assembly report and return placement metrics."""
    accession = str(assembly_row["assembly_accession"])
    url = assembly_report_url(assembly_row.get("ftp_path"))

    default: dict[str, object] = {
        "unplaced_scaffold_count": np.nan,
        "unplaced_scaffold_percentage": np.nan,
        "unlocalized_scaffold_count": np.nan,
        "unlocalized_scaffold_percentage": np.nan,
    }

    if url is None:
        return default

    report_path = report_directory / f"{accession}_assembly_report.txt"

    try:
        download_file(url, report_path)
        report = read_assembly_report(report_path)
        return summarize_sequence_roles(
            report=report,
            assembly_size=float(assembly_row.get("genome_size", np.nan)),
        )
    except Exception:
        return default


# ---------------------------------------------------------------------
# Pathogen-sheet enrichment
# ---------------------------------------------------------------------

def enrich_pathogen_sheet(
    pathogen_sheet: pd.DataFrame,
    assemblies: pd.DataFrame,
    report_directory: Path,
) -> pd.DataFrame:
    """Add metadata for each assembly accession already in the pathogen sheet."""
    output = pathogen_sheet.copy()

    metadata_columns = [
        "assembly_size_bp",
        "assembly_level",
        "refseq_category",
        "contig_count",
        "replicon_count",
        "unplaced_scaffold_count",
        "unplaced_scaffold_percentage",
        "unlocalized_scaffold_count",
        "unlocalized_scaffold_percentage",
        "assembly_rel_date",
        "assembly_submitter",
    ]

    # Replace existing versions of these columns rather than creating suffixes.
    output = output.drop(
        columns=[column for column in metadata_columns if column in output.columns]
    )

    assembly_lookup = (
        assemblies.drop_duplicates(subset=["assembly_accession"], keep="first")
        .set_index("assembly_accession", drop=False)
    )

    metadata_rows: list[dict[str, object]] = []

    # Cache report-derived values when an accession occurs in multiple rows.
    report_cache: dict[str, dict[str, object]] = {}

    empty_metadata = {
        "assembly_size_bp": np.nan,
        "assembly_level": pd.NA,
        "refseq_category": pd.NA,
        "contig_count": np.nan,
        "replicon_count": np.nan,
        "unplaced_scaffold_count": np.nan,
        "unplaced_scaffold_percentage": np.nan,
        "unlocalized_scaffold_count": np.nan,
        "unlocalized_scaffold_percentage": np.nan,
        "assembly_rel_date": pd.NA,
        "assembly_submitter": pd.NA,
    }

    for _, pathogen in output.iterrows():
        accession = normalize_accession(pathogen.get("assembly_accession"))

        if accession is None or accession not in assembly_lookup.index:
            metadata_rows.append(empty_metadata.copy())
            continue

        assembly = assembly_lookup.loc[accession]

        if accession not in report_cache:
            report_cache[accession] = get_assembly_report_metrics(
                assembly_row=assembly,
                report_directory=report_directory,
            )

        metadata_rows.append({
            "assembly_size_bp": assembly.get("genome_size", np.nan),
            "assembly_level": assembly.get("assembly_level", pd.NA),
            "refseq_category": assembly.get("refseq_category", pd.NA),
            "contig_count": assembly.get("contig_count", np.nan),
            "replicon_count": assembly.get("replicon_count", np.nan),
            **report_cache[accession],
            "assembly_rel_date": assembly.get("seq_rel_date", pd.NA),
            "assembly_submitter": assembly.get("asm_submitter", pd.NA),
        })

    metadata_df = pd.DataFrame(metadata_rows, index=output.index)
    return pd.concat([output, metadata_df], axis=1)


# ---------------------------------------------------------------------
# Command-line entry point
# ---------------------------------------------------------------------

def parse_arguments() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent

    parser = argparse.ArgumentParser(
        description=(
            "Add RefSeq assembly metadata and sequence-placement metrics to "
            "the assembly accessions already listed in a pathogen sheet."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--pathogen-sheet",
        type=Path,
        default="/Users/tas1/Documents/APHL/taxtriage/assets/pathogen_sheet.csv",
        help="Input pathogen sheet",
    )

    parser.add_argument(
        "--assembly-summary",
        type=Path,
        default="/Users/tas1/Documents/APHL/taxtriage/assets/assembly_summary_refseq.txt",
        help="Location where the RefSeq assembly summary is downloaded",
    )

    parser.add_argument(
        "--output",
        type=Path,
        default="/Users/tas1/Documents/APHL/taxtriage/assets/pathogen_sheet_with_assembly_metrics.csv",
        help="Output pathogen sheet containing the added assembly metrics",
    )

    parser.add_argument(
        "--report-directory",
        type=Path,
        default="/Users/tas1/Documents/APHL/taxtriage/assets/assembly_reports",
        help="Directory used to cache downloaded assembly reports",
    )

    parser.add_argument(
        "--use-existing-assembly-summary",
        action="store_true",
        help=(
            "Use the existing local assembly summary instead of downloading "
            "the newest version from NCBI"
        ),
    )

    return parser.parse_args()
def main() -> None:
    args = parse_arguments()

    pathogen_sheet = pd.read_csv(
        args.pathogen_sheet,
        dtype={"taxid": str, "assembly_accession": str},
        low_memory=False,
    )

    if "assembly_accession" not in pathogen_sheet.columns:
        raise ValueError(
            "The pathogen sheet does not contain an assembly_accession column"
        )

    # Create the report cache directory before any report downloads.
    args.report_directory.mkdir(parents=True, exist_ok=True)

    if args.use_existing_assembly_summary:
        if not args.assembly_summary.exists():
            raise FileNotFoundError(
                "--use-existing-assembly-summary was specified, but the "
                f"summary file does not exist: {args.assembly_summary}"
            )
        print(f"Using existing assembly summary: {args.assembly_summary}")
    else:
        print(
            "Downloading newest RefSeq assembly summary to "
            f"{args.assembly_summary}"
        )
        download_file(
            url=ASSEMBLY_SUMMARY_REFSEQ_URL,
            output_path=args.assembly_summary,
            force=True,
        )

    assemblies = read_assembly_summary(args.assembly_summary)

    enriched_sheet = enrich_pathogen_sheet(
        pathogen_sheet=pathogen_sheet,
        assemblies=assemblies,
        report_directory=args.report_directory,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    enriched_sheet.to_csv(args.output, index=False)

    print(f"Enriched pathogen sheet: {args.output}")
    print(f"Assembly report cache: {args.report_directory}")


if __name__ == "__main__":
    main()
