#!/usr/bin/env python3
##############################################################################################
# Copyright 2022 The Johns Hopkins University Applied Physics Laboratory LLC
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

"""Extract top hits from a Kraken2 report with optional rank filtering,
specific taxid limits, and distribution-based z-score filtering."""

import argparse
import csv
import logging
import re
import sys
from pathlib import Path

import pandas as pd

from distributions import body_site_map, import_distributions

logger = logging.getLogger()

KRAKEN_HEADER = [
    "abundance",
    "clade_fragments_covered",
    "number_fragments_assigned",
    "rank",
    "taxid",
    "name",
    "parents",
]

OUTPUT_HEADER = [
    "abundance",
    "clade_fragments_covered",
    "number_fragments_assigned",
    "rank",
    "taxid",
    "name",
]


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Extract top-hit organisms from a Kraken2 report.",
        epilog=(
            "Example:\n"
            "  python get_top_hits.py -i report.tsv -o top_hits.tsv -t 5 "
            "--ranks S S1 -s 10239:100:S 2:100:S"
        ),
    )
    parser.add_argument(
        "-i", "--file_in",
        type=Path, required=True,
        help="Input Kraken2 report (tab-delimited).",
    )
    parser.add_argument(
        "-o", "--file_out",
        type=Path, required=True,
        help="Output TSV file for the top-hit organisms.",
    )
    parser.add_argument(
        "--ranks",
        nargs="+", default=None, metavar="RANK",
        help=(
            "Only include specific Kraken ranks in the general top-hit fill "
            "(e.g. S S1 S2 G). If omitted, all ranks are included. "
            "Force-included taxids (pathogens, distributions, -s overrides) "
            "always appear regardless of this filter."
        ),
    )
    parser.add_argument(
        "-s", "--top_hits_string",
        nargs="+", default=[], metavar="TAXID:COUNT:RANK",
        help=(
            "Per-taxid overrides as TAXID:COUNT:RANK. "
            "These always take priority over -t and --ranks. "
            "Example: -s 10239:100:S 2:100:S"
        ),
    )
    parser.add_argument(
        "-t", "--top_per_rank",
        type=int, default=5,
        help="Default max hits per rank (default: 5). Overridden by -s for specific taxids.",
    )
    parser.add_argument(
        "-d", "--distributions",
        type=Path, default=None,
        help="Optional HMP distributions file (TSV) with z-score data per taxid.",
    )
    parser.add_argument(
        "-p", "--pathogens",
        type=Path, default=None,
        help="Optional CSV with annotated pathogens (matched by taxid).",
    )
    parser.add_argument(
        "-z", "--zscore",
        type=float, default=1.5,
        help="Z-score threshold for distribution filtering (default: 1.5).",
    )
    parser.add_argument(
        "-b", "--body_site",
        type=str, default="Unknown",
        help="Body site for distribution look-up (default: Unknown).",
    )
    parser.add_argument(
        "--include_taxids",
        type=int, nargs="+", default=[],
        help="Taxids to force-include regardless of rank or distribution.",
    )
    parser.add_argument(
        "--remove_commensals",
        action="store_true", default=False,
        help="Remove commensal organisms listed in the pathogen sheet.",
    )
    parser.add_argument(
        "-l", "--log-level",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
        help="Logging verbosity (default: WARNING).",
    )
    return parser.parse_args(argv)


# ---------------------------------------------------------------------------
# Kraken report parsing
# ---------------------------------------------------------------------------

def import_file(filepath):
    """Parse a Kraken2 report into a list of entry dicts.

    Always returns the FULL unfiltered mapping. Rank filtering is applied
    later so that force-included taxids (pathogens, distributions, -s)
    are never silently dropped.
    """
    header = KRAKEN_HEADER[:6]  # parents are computed, not in file
    depth_regex = re.compile(r"^(\s+)(.+)")

    mapping = []
    last_parents = {}

    with open(filepath, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            entry = {}
            for idx, col in enumerate(header):
                if idx >= len(row):
                    entry[col] = ""
                    continue
                val = row[idx]
                if col in ("abundance", "clade_fragments_covered", "number_fragments_assigned"):
                    entry[col] = float(val)
                elif col == "taxid":
                    entry[col] = int(val)
                else:
                    entry[col] = val

            # Determine depth from leading whitespace in the name field
            match = depth_regex.search(entry["name"])
            if match:
                depth = len(match.group(1)) // 2
                entry["name"] = match.group(2)
            else:
                depth = 0

            entry["depth"] = depth

            # Build parent lineage (walk back up the tree)
            parents = []
            for d in range(depth - 1, 0, -1):
                if d in last_parents:
                    parents.append(int(last_parents[d]))
            entry["parents"] = parents

            last_parents[depth] = entry["taxid"]
            mapping.append(entry)

    return mapping


# ---------------------------------------------------------------------------
# Top-hit selection
# ---------------------------------------------------------------------------

def select_top_hits(full_mapping, specific_limits, top_per_rank,
                    extra_taxids=None, allowed_ranks=None):
    """Select top organisms from *full_mapping*.

    Priority order (each step ignores the --ranks filter):
      1. *extra_taxids*  – always included (pathogens, distributions, --include_taxids).
      2. *specific_limits* – per-taxid overrides from -s, always honoured.
      3. General fill     – up to *top_per_rank* per rank, restricted to
                            *allowed_ranks* when set.

    Returns a dict {taxid: row}.
    """
    extra_taxids = extra_taxids or []

    # Index the full mapping by taxid (for force-includes + specific limits)
    by_taxid = {row["taxid"]: row for row in full_mapping}

    # Group ALL rows by rank, sorted by abundance descending
    by_rank_full = {}
    for row in full_mapping:
        by_rank_full.setdefault(row["rank"], []).append(row)
    for rank in by_rank_full:
        by_rank_full[rank].sort(key=lambda r: r["abundance"], reverse=True)

    selected = {}       # taxid -> row
    rank_counts = {}    # rank  -> count selected so far

    # ── Step 1: Force-include extra taxids ──────────────────────────────
    # These bypass --ranks entirely. Pathogens, distribution outliers,
    # and --include_taxids all land here.
    for tid in extra_taxids:
        if tid in by_taxid and tid not in selected:
            selected[tid] = by_taxid[tid]
            logger.debug("Force-including taxid %d (%s)", tid, by_taxid[tid]["name"])

    # ── Step 2: Specific per-taxid limits (-s) ──────────────────────────
    # These also bypass --ranks. If you say -s 10239:100:S, we pull 100 S
    # descendants of 10239 from the FULL mapping regardless of --ranks.
    for parent_taxid, spec in specific_limits.items():
        rank = spec["rank"]
        limit = spec["limit"]
        if rank not in by_rank_full:
            continue
        count = 0
        for row in by_rank_full[rank]:
            if parent_taxid in row["parents"]:
                if count >= limit:
                    break
                if row["taxid"] not in selected:
                    selected[row["taxid"]] = row
                    rank_counts[rank] = rank_counts.get(rank, 0) + 1
                count += 1

    # ── Step 3: General fill (respects --ranks) ─────────────────────────
    # Build a rank-filtered view for the general fill only.
    if allowed_ranks is not None:
        rank_set = set(allowed_ranks)
        general_ranks = {r: rows for r, rows in by_rank_full.items() if r in rank_set}
    else:
        general_ranks = by_rank_full

    for rank, rows in general_ranks.items():
        current = rank_counts.get(rank, 0)
        for row in rows:
            if current >= top_per_rank:
                break
            if row["taxid"] not in selected:
                selected[row["taxid"]] = row
                current += 1
        rank_counts[rank] = current

    return selected


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_output(selected, outpath):
    """Write the selected top-hit rows to a TSV file."""
    outpath.parent.mkdir(parents=True, exist_ok=True)
    with open(outpath, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(OUTPUT_HEADER)
        for row in selected.values():
            writer.writerow([
                row[col].strip() if col == "name" else row[col]
                for col in OUTPUT_HEADER
            ])
    logger.info("Wrote top-hits report to %s", outpath)


# ---------------------------------------------------------------------------
# Distribution z-score filtering
# ---------------------------------------------------------------------------

def _get_dist_stats(row, stats_dict):
    """Look up distribution stats for a (tax_id, body_site) pair."""
    return stats_dict.get((row["tax_id"], row["body_site"]))


def apply_distribution_filter(full_mapping, args, body_sites):
    """Run the HMP distribution z-score filter on the full (unfiltered) mapping.

    The distribution file contains per-taxid abundance statistics from HMP
    (columns: rank, mean, std, tax_id, name, body_site, site_count, abundances, ...).

    For each species-level taxid present in BOTH the Kraken report and the
    distribution file, we compute:

        norm_abundance = sum(observed_abundances) / site_count
            → population mean that accounts for samples where the taxid
              was absent (implicit zeros)

        zscore = (sample_abundance - norm_abundance) / std

    Taxids with zscore > threshold are considered *out of line* relative to
    healthy baseline — these are force-included in top hits.

    Taxids with zscore <= threshold look commensal-like and may be removed
    if --remove_commensals is set.

    Returns (passing_taxids, below_threshold_taxids).
    """
    dists, _site_counts = import_distributions(
        args.distributions, "tax_id", body_sites,
    )

    # Work on species-level rows from the FULL mapping (no rank filter applied)
    df = pd.DataFrame(full_mapping)
    df = df[df["rank"].str.contains("S")]

    if df.empty:
        logger.warning("No species-level rows in Kraken report for distribution filter.")
        return [], []

    df = df.rename(columns={"taxid": "tax_id"})
    df["body_site"] = (body_sites[0] if body_sites else None)
    if df["body_site"].dtype == object:
        df["body_site"] = df["body_site"].str.lower()

    # Keep rows that either:
    #   a) have a matching entry in the distribution stats, OR
    #   b) have never been seen in distributions at all (novel taxid)
    known_taxids = {k[0] for k in dists}
    df = df[
        df.apply(
            lambda r: (r["tax_id"], r["body_site"]) in dists
            or r["tax_id"] not in known_taxids,
            axis=1,
        )
    ]

    if df.empty:
        logger.warning("No data remaining after distribution filter.")
        return [], []

    # Compute z-scores
    df["stats"] = df.apply(lambda r: _get_dist_stats(r, dists), axis=1)
    df["norm_abundance"] = df["stats"].apply(lambda s: s["norm_abundance"] if s else None)
    df["std"] = df["stats"].apply(lambda s: s["std"] if s else None)

    # zscore = (observed - population_norm) / std
    # High positive zscore → abundance is unusually elevated for this body site
    df["zscore"] = (df["abundance"] - df["norm_abundance"]) / df["std"]
    df["zscore"] = df["zscore"].fillna(-1)

    # Taxids above the threshold are outliers → force-include
    passing = df[df["zscore"] > args.zscore]["tax_id"].tolist()
    # Taxids at or below threshold look commensal → candidate for removal
    below_threshold = (
        df[df["zscore"] <= args.zscore]["tax_id"].tolist()
        if args.remove_commensals
        else []
    )

    logger.info(
        "Distribution filter: %d taxids above z=%.1f, %d below",
        len(passing), args.zscore, len(below_threshold),
    )
    return passing, below_threshold


# ---------------------------------------------------------------------------
# Pathogen sheet processing
# ---------------------------------------------------------------------------

def _translate_sites(sites_str):
    """Map a comma-separated body-site string through body_site_map and deduplicate."""
    sites = [s.strip() for s in sites_str.split(",") if s.strip()]
    translated = set()
    for site in sites:
        mapped = body_site_map(site.lower())
        if isinstance(mapped, list):
            translated.update(mapped)
        else:
            translated.add(mapped)
    return ", ".join(sorted(translated))


def process_pathogens(args, body_sites, seen_taxids):
    """Load the pathogen sheet and return (pathogen_taxids, commensal_taxids_to_remove).

    Pathogen taxids that appear in the Kraken report (*seen_taxids*) are
    returned for force-inclusion — they bypass the --ranks filter entirely.
    """
    pathogen_taxids = []
    remove_taxids = []

    if not args.pathogens:
        return pathogen_taxids, remove_taxids

    try:
        with open(args.pathogens, "r", encoding="utf-8", errors="replace") as fh:
            sheet = pd.read_csv(fh, sep=",")

        # Normalise body-site columns
        sheet["pathogenic_sites"] = (
            sheet["pathogenic_sites"].fillna("Unknown").apply(_translate_sites)
        )
        sheet["commensal_sites"] = (
            sheet["commensal_sites"].fillna("").apply(_translate_sites).str.lower()
        )

        # Filter to pathogens relevant for this body site
        if body_sites:
            sheet = sheet[
                sheet["pathogenic_sites"].str.lower().isin(body_sites)
                | ("sterile" in body_sites)
            ]

        pathogen_classes = {
            "primary", "opportunistic", "potential", "oportunistic",
        }
        pathogen_taxids = (
            sheet.loc[sheet["general_classification"].isin(pathogen_classes), "taxid"]
            .dropna().astype(int).tolist()
        )
        # Keep only those actually present in the Kraken report
        pathogen_taxids = [t for t in pathogen_taxids if t in seen_taxids]
        logger.info("Found %d pathogen taxids in Kraken report", len(pathogen_taxids))

        # Commensal removal (only if body site is not sterile)
        if args.remove_commensals and "sterile" not in body_sites:
            commensal_taxids = (
                sheet.loc[sheet["general_classification"] == "commensal", "taxid"]
                .dropna().astype(int).tolist()
            )
            remove_taxids.extend(commensal_taxids)

    except Exception as exc:
        logger.error("Error reading pathogens file: %s", exc)

    return pathogen_taxids, remove_taxids


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv=None):
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")

    if not args.file_in.is_file():
        logger.error("Input file not found: %s", args.file_in)
        sys.exit(2)

    # --- Parse -s specific limits (TAXID:COUNT:RANK) -----------------------
    specific_limits = {}
    for entry in args.top_hits_string:
        parts = entry.split(":")
        if len(parts) >= 3:
            specific_limits[int(parts[0])] = {
                "limit": int(parts[1]),
                "rank": parts[2],
            }

    # --- Import the FULL Kraken report (no rank filter yet) ----------------
    full_mapping = import_file(args.file_in)
    seen_taxids = {row["taxid"] for row in full_mapping}

    # --- Body site ---------------------------------------------------------
    body_sites = (
        [body_site_map(args.body_site.lower())]
        if args.body_site and args.body_site != "Unknown"
        else []
    )

    # --- Force-include taxids (from --include_taxids) ----------------------
    extra_orgs = [t for t in args.include_taxids if t in seen_taxids]

    # --- Pathogens (uses full mapping's seen_taxids) -----------------------
    pathogen_taxids, remove_taxids = process_pathogens(args, body_sites, seen_taxids)
    extra_orgs.extend(pathogen_taxids)

    # --- HMP Distributions (uses full mapping, species-level) --------------
    if args.distributions:
        dist_passing, dist_remove = apply_distribution_filter(
            full_mapping, args, body_sites,
        )
        extra_orgs.extend(dist_passing)
        if args.remove_commensals:
            remove_taxids.extend(dist_remove)

    # --- Apply commensal removal -------------------------------------------
    remove_set = set(remove_taxids)
    if remove_set:
        full_mapping = [m for m in full_mapping if m["taxid"] not in remove_set]
        extra_orgs = [t for t in extra_orgs if t not in remove_set]
        for tid in remove_set & seen_taxids:
            logger.info("Removing commensal taxid %d", tid)

    extra_orgs = list(set(extra_orgs))

    # --- Build the effective rank filter for general fill ------------------
    # --ranks restricts ONLY the general fill (step 3). Force-includes from
    # pathogens, distributions, --include_taxids, and -s all bypass it.
    allowed_ranks = set(args.ranks) if args.ranks is not None else None

    # --- Select top hits ---------------------------------------------------
    selected = select_top_hits(
        full_mapping, specific_limits, args.top_per_rank,
        extra_taxids=extra_orgs, allowed_ranks=allowed_ranks,
    )

    # --- Write output ------------------------------------------------------
    write_output(selected, args.file_out)


if __name__ == "__main__":
    sys.exit(main())
