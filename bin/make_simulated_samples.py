#!/usr/bin/env python3
"""
Generate abundance files and reference FASTAs for read simulation from top_report.tsv files.

This script produces the INPUT files needed by InSilicoSeq (ISS) and NanoSim:
1. Parses top_report.tsv to extract species-level taxids with abundances
2. Maps taxids to accessions via the merged.taxid.tsv mapping file
3. Extracts reference sequences from FASTA file(s)
4. Creates per-sample abundance.tsv files (accession<TAB>abundance)

The actual read simulation (ISS / NanoSim) is handled OUTSIDE this script
by the calling bash scripts (get_refs.sh, generate_batch_samples.sh).
"""

import argparse
import csv
import random
import shutil
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional

import pysam


def dirichlet_weights(n: int) -> List[float]:
    """
    Simple Dirichlet(1,...,1) sampler using exponential variables.
    Produces n positive weights that sum to 1.0.
    """
    if n <= 0:
        return []
    xs = [random.expovariate(1.0) for _ in range(n)]
    s = sum(xs)
    return [x / s for x in xs]


def parse_top_report(
    top_report_path: Path,
    ranks: List[str],
    exclude_taxids: Set[int],
    include_taxids: Optional[Set[int]] = None,
    use_abundance: bool = False,
    minreads_threshold: int = 3
) -> Dict[int, float]:
    """
    Parse a top_report.tsv file and extract taxa at specified ranks.

    Args:
        top_report_path: Path to top_report.tsv file
        ranks: List of rank codes to extract (e.g., ['S', 'S1', 'S2'])
        exclude_taxids: Set of taxids to exclude (e.g., host contamination)
        use_abundance: If True, use abundance column; otherwise use number_fragments_assigned

    Returns:
        Dictionary mapping taxid -> abundance or fragment count

    Expected format (tab-delimited with header):
        abundance	clade_fragments_covered	number_fragments_assigned	rank	taxid	name
    """
    taxid_values = {}

    with top_report_path.open('r', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')

        for row in reader:
            try:
                abundance = float(row.get('abundance', 0))
                num_frags = float(row.get('number_fragments_assigned', 0))
                num_frags_covered = float(row.get('clade_fragments_covered', 0))
                rank = row.get('rank', '').strip()
                taxid = int(row.get('taxid', 0))
            except (ValueError, KeyError):
                continue
            # if taxid == 57975 or taxid == 271848:
            #     print(abundance, num_frags, num_frags_covered, rank, taxid)
            #     exit()
            if taxid not in include_taxids:
                if num_frags_covered < minreads_threshold:
                    continue
                # Check if this is a rank we want
                if rank not in ranks:
                    continue
                # Exclude specific taxids
                if taxid in exclude_taxids:
                    continue
            # Choose value based on flag
            if use_abundance:
                value = abundance
            else:
                value = num_frags_covered

            if value > 0:
                taxid_values[taxid] = value

    return taxid_values


def parse_merged_taxid_file(merged_taxid_path: Path) -> Dict[int, List[str]]:
    """
    Parse the merged.taxid.tsv file to create taxid -> accessions mapping.

    Expected format (tab-delimited with header):
        Acc	Assembly	Organism_Name	Description	Mapped_Value

    Returns:
        Dictionary mapping taxid -> list of accessions
    """
    taxid_to_accs = defaultdict(list)

    with merged_taxid_path.open('r', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')

        for row in reader:
            try:
                # Try multiple possible column names for accession
                acc = (row.get('Acc') or row.get('Accession') or '').strip()
                # Try multiple possible column names for taxid
                taxid_str = (row.get('Mapped_Value') or row.get('TaxID') or row.get('taxid') or '').strip()
                taxid = int(taxid_str)
            except (ValueError, KeyError, AttributeError):
                continue

            if acc and taxid:
                taxid_to_accs[taxid].append(acc)

    # Deduplicate accessions per taxid
    for taxid in taxid_to_accs:
        taxid_to_accs[taxid] = list(set(taxid_to_accs[taxid]))

    return dict(taxid_to_accs)


def build_fasta_index(fasta_path: Path) -> Dict[str, int]:
    """
    Build an in-memory lookup of sequence IDs in a FASTA file using pysam.

    If an .fai index doesn't exist, pysam will create one automatically.

    Args:
        fasta_path: Path to the FASTA file

    Returns:
        Set of reference names available in the FASTA
    """
    fasta = pysam.FastaFile(str(fasta_path))
    refs = set(fasta.references)
    fasta.close()
    return refs


def discover_organisms_from_fasta(
    fasta_path: Path,
    taxid_to_accs: Dict[int, List[str]],
    taxid_values: Dict[int, float],
    default_min_count: float = 1.0,
) -> Dict[int, float]:
    """
    Discover all organisms present in the reference FASTA by cross-referencing
    accessions with the merged_taxid mapping.  This ensures every organism that
    has a downloaded reference gets included in the simulation, even if its
    taxid doesn't appear in the top_report (e.g. due to taxid level mismatch
    between species-level Kraken2 report and strain-level reference accessions).

    For organisms already in taxid_values (from top_report), the original
    abundance/count is kept.  For newly discovered organisms, a minimum count
    is assigned so they still receive simulated reads.

    Args:
        fasta_path: Path to the merged reference FASTA
        taxid_to_accs: Dictionary of taxid -> list of accessions (from merged_taxid.tsv)
        taxid_values: Existing taxid -> value dict (from top_report parsing)
        default_min_count: Minimum count assigned to newly discovered organisms

    Returns:
        Updated dictionary mapping taxid -> abundance/count (superset of taxid_values)
    """
    # Get all accessions present in the reference FASTA
    fasta_accessions = build_fasta_index(fasta_path)
    if not fasta_accessions:
        return taxid_values

    # Build reverse map: accession -> taxid
    acc_to_taxid = {}
    for taxid, accs in taxid_to_accs.items():
        for acc in accs:
            acc_to_taxid[acc] = taxid

    # Find all taxids that have at least one accession in the FASTA
    fasta_taxids = set()
    for acc in fasta_accessions:
        if acc in acc_to_taxid:
            fasta_taxids.add(acc_to_taxid[acc])
        else:
            # Try partial match (some FASTAs have "ACC.1 description" headers)
            for ref_acc, taxid in acc_to_taxid.items():
                if ref_acc.startswith(acc) or acc.startswith(ref_acc):
                    fasta_taxids.add(taxid)
                    break

    # Merge: keep existing values, add newly discovered taxids with min count
    merged = dict(taxid_values)
    newly_discovered = 0
    for taxid in fasta_taxids:
        if taxid not in merged:
            merged[taxid] = default_min_count
            newly_discovered += 1

    if newly_discovered > 0:
        print(f"  FASTA discovery: found {newly_discovered} additional taxids "
              f"in reference FASTA (not in top_report)")
        print(f"  Total organisms for simulation: {len(merged)}")

    return merged


def extract_sequences(fasta_path: Path, accessions: List[str], output_fasta: Path,
                      _fasta_handle: 'pysam.FastaFile' = None) -> None:
    """
    Extract sequences from FASTA file by accession using pysam.

    Uses the FASTA index (.fai) for fast random access. If the index
    doesn't exist, pysam creates it automatically on first open.

    Args:
        fasta_path: Input FASTA file
        accessions: List of accessions to extract
        output_fasta: Output FASTA file path
        _fasta_handle: Optional pre-opened pysam.FastaFile for reuse
    """
    close_handle = False
    fa = _fasta_handle
    if fa is None:
        fa = pysam.FastaFile(str(fasta_path))
        close_handle = True

    available = set(fa.references)
    found = 0

    try:
        with output_fasta.open('w') as out:
            for acc in accessions:
                # Try exact match first, then try matching just the first token
                if acc in available:
                    seq = fa.fetch(acc)
                    out.write(f">{acc}\n{seq}\n")
                    found += 1
                else:
                    # Try partial match: some FASTAs have "ACC.1 description" as the full header
                    for ref in available:
                        if ref.startswith(acc) or acc.startswith(ref):
                            seq = fa.fetch(ref)
                            out.write(f">{ref}\n{seq}\n")
                            found += 1
                            break
        if found == 0:
            print(f"  WARNING: No sequences found for {len(accessions)} accessions in {fasta_path.name}")
    finally:
        if close_handle:
            fa.close()


def create_abundance_file(
    taxid_values: Dict[int, float],
    taxid_to_accs: Dict[int, List[str]],
    output_path: Path,
    random_abundance: bool = False
) -> Tuple[Dict[str, float], float]:
    """
    Create an abundance file for ISS.

    Args:
        taxid_values: Dictionary of taxid -> abundance/count
        taxid_to_accs: Dictionary of taxid -> list of accessions
        output_path: Path to write abundance file
        random_abundance: If True, assign random abundances per taxid

    Returns:
        Tuple of (accession -> abundance dict, total value)
    """
    acc_abundance = {}

    if random_abundance:
        # Assign random abundance per taxid using Dirichlet
        taxids = list(taxid_values.keys())
        taxid_weights = dirichlet_weights(len(taxids))

        for taxid, weight in zip(taxids, taxid_weights):
            accs = taxid_to_accs.get(taxid, [])
            if not accs:
                continue

            # Split taxid abundance across its accessions
            acc_weights = dirichlet_weights(len(accs))
            for acc, acc_w in zip(accs, acc_weights):
                acc_abundance[acc] = weight * acc_w
    else:
        # Use values from top_report (normalize later)
        total_value = sum(taxid_values.values())

        for taxid, value in taxid_values.items():
            accs = taxid_to_accs.get(taxid, [])
            if not accs:
                continue

            taxid_proportion = value / total_value if total_value > 0 else 0

            # Split taxid proportion across accessions using Dirichlet
            acc_weights = dirichlet_weights(len(accs))
            for acc, acc_w in zip(accs, acc_weights):
                acc_abundance[acc] = taxid_proportion * acc_w

    # Normalize to sum to 1.0
    total_abu = sum(acc_abundance.values())
    if total_abu > 0:
        for acc in acc_abundance:
            acc_abundance[acc] /= total_abu

    # Write abundance file
    with output_path.open('w') as f:
        for acc, abu in acc_abundance.items():
            f.write(f"{acc}\t{abu}\n")

    return acc_abundance, sum(taxid_values.values())


def main():
    parser = argparse.ArgumentParser(
        description="Generate abundance files and reference FASTAs for ISS/NanoSim",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  # Basic usage – produces abundance.tsv + reference FASTA per sample
  %(prog)s --top_report top_report.tsv --merged_taxid all.merged.taxid.tsv \\
    --fasta all.fasta --outdir simulated_reads --nreads 100000

  # Multiple samples with random abundances
  %(prog)s --top_report top_report.tsv --merged_taxid all.merged.taxid.tsv \\
    --fasta all.fasta --outdir simulated_reads --nreads 100000 \\
    --nsamples 5 --random_abundance

Outputs per sample directory:
  abundance.tsv        – accession<TAB>fractional_abundance (sums to 1.0)
  reference.fasta      – extracted reference sequences for simulation
  metadata.tsv         – sample parameters
  taxa.tsv             – taxid-level summary
        """
    )

    # Input files
    parser.add_argument(
        "--top_report",
        type=Path,
        required=False,
        default=None,
        help="Path to top_report.tsv file (default input source)"
    )
    parser.add_argument(
        "--abundance_input",
        type=Path,
        required=False,
        default=None,
        help="Direct abundance input: 2-column TSV (taxid<TAB>abundance). "
             "Bypasses --top_report parsing. Still requires --merged_taxid and --fasta."
    )
    parser.add_argument(
        "--minreads",
        type=int,
        required=False,
        help="Threhsold for min clade covered in top_report.tsv (default: 3)",
        default=3
    )
    parser.add_argument(
        "--merged_taxid",
        type=Path,
        required=True,
        help="Path to merged.taxid.tsv mapping file"
    )
    parser.add_argument(
        "--fasta",
        type=Path,
        required=True,
        help="Path to FASTA file containing reference sequences"
    )

    # Output
    parser.add_argument(
        "--outdir",
        type=Path,
        required=True,
        help="Output directory for simulated reads"
    )

    # Parsing options
    parser.add_argument(
        "--ranks",
        nargs="+",
        default=["S", "S1", "S2"],
        help="Rank codes to extract from top_report (default: S S1 S2)"
    )
    parser.add_argument(
        "--exclude_taxids",
        nargs="+",
        type=int,
        default=[9606],
        help="TaxIDs to exclude (default: 9606 for Homo sapiens)"
    )
    parser.add_argument(
        "--include_taxids",
        nargs="+",
        type=int,
        default=[],
        help="TaxIDs to ALWAYS include if present in the report regardless of the top report filtering"
    )

    # Read simulation
    parser.add_argument(
        "--nreads",
        type=int,
        default=None,
        help="Total number of reads to simulate. If not specified, uses sum of number_fragments_assigned from top_report"
    )
    parser.add_argument(
        "--use_abundance",
        action="store_true",
        help="Use abundance column instead of number_fragments_assigned. Only valid when --nreads is specified. "
             "When set, abundances are used for proportions and --nreads determines total reads."
    )
    parser.add_argument(
        "--random_abundance",
        action="store_true",
        help="Assign random abundances per taxid instead of using top_report values"
    )

    # Batch processing
    parser.add_argument(
        "--nsamples",
        type=int,
        default=1,
        help="Number of samples to generate (creates subdirectories)"
    )

    # Other
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility (default: 42)"
    )

    args = parser.parse_args()

    # Validate inputs: must provide either --top_report or --abundance_input
    if not args.top_report and not args.abundance_input:
        raise SystemExit(
            "ERROR: Must provide either --top_report or --abundance_input"
        )

    # Validate use_abundance flag
    if args.use_abundance and args.nreads is None:
        raise SystemExit(
            "ERROR: --use_abundance can only be used when --nreads is specified.\n"
            "Either:\n"
            "  1) Add --nreads <number> to use abundance column for proportions, or\n"
            "  2) Remove --use_abundance to use number_fragments_assigned (default)"
        )

    # Set random seed
    random.seed(args.seed)

    # Create output directory
    args.outdir.mkdir(parents=True, exist_ok=True)

    # Convert exclude_taxids to set
    exclude_taxids = set(args.exclude_taxids)
    include_taxids = set(args.include_taxids)

    # Check that input files exist
    if args.top_report and not args.top_report.exists():
        raise SystemExit(f"Top report file not found: {args.top_report}")
    if args.abundance_input and not args.abundance_input.exists():
        raise SystemExit(f"Abundance input file not found: {args.abundance_input}")
    if not args.merged_taxid.exists():
        raise SystemExit(f"Merged taxid file not found: {args.merged_taxid}")
    if not args.fasta.exists():
        raise SystemExit(f"FASTA file not found: {args.fasta}")

    # Parse mapping file
    print("Parsing merged taxid mapping file...")
    taxid_to_accs = parse_merged_taxid_file(args.merged_taxid)
    print(f"Found {len(taxid_to_accs)} unique taxids in mapping file")

    # Parse taxid values from either --abundance_input or --top_report
    if args.abundance_input:
        # Custom abundance input: 2-column TSV (taxid<TAB>abundance)
        print(f"Reading custom abundance input: {args.abundance_input}")
        taxid_values = {}
        with args.abundance_input.open() as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) >= 2:
                    try:
                        taxid = int(parts[0])
                        value = float(parts[1])
                        if taxid not in exclude_taxids and value > 0:
                            taxid_values[taxid] = value
                    except ValueError:
                        continue
        print(f"Read {len(taxid_values)} taxids from custom abundance input")
    else:
        # Parse top report (original behavior)
        print(f"Parsing top_report.tsv for ranks: {', '.join(args.ranks)}...")
        if args.use_abundance:
            print("Using 'abundance' column for proportions")
            use_abundance_col = True
        else:
            print("Using 'number_fragments_assigned' column for counts")
            use_abundance_col = False

        taxid_values = parse_top_report(args.top_report, args.ranks, exclude_taxids, include_taxids, use_abundance_col, args.minreads)

    if not taxid_values:
        print("WARNING: No taxa found from top_report/abundance_input "
              f"(after excluding taxids: {exclude_taxids}). "
              "Will attempt to discover organisms from the reference FASTA.")
        taxid_values = {}

    print(f"Found {len(taxid_values)} unique taxids from top_report/abundance_input")

    # ── FASTA-first organism discovery ──────────────────────────────────────
    # Cross-reference the reference FASTA with merged_taxid.tsv to discover
    # ALL organisms that have downloaded references.  This catches organisms
    # whose taxids in merged_taxid are at a different rank level (e.g. strain)
    # than the top_report (e.g. species).
    taxid_values = discover_organisms_from_fasta(
        args.fasta, taxid_to_accs, taxid_values,
        default_min_count=max(1.0, args.minreads)
    )

    if not taxid_values:
        raise SystemExit(
            "No taxa found after parsing input and scanning reference FASTA "
            f"(after excluding taxids: {exclude_taxids})"
        )

    # Filter to taxids present in mapping file
    taxids_with_accessions = {
        tid for tid in taxid_values.keys() if tid in taxid_to_accs
    }

    if not taxids_with_accessions:
        raise SystemExit(
            "No taxids found in mapping file. "
            "Check that your mapping file covers the taxa in your report/reference."
        )

    print(f"{len(taxids_with_accessions)} taxids have accessions in mapping file")

    # Filter values to only mapped taxids
    filtered_values = {
        tid: taxid_values[tid] for tid in taxids_with_accessions
    }

    # Determine number of reads to simulate
    total_counts_from_report = sum(filtered_values.values())
    if args.nreads is None:
        # Use the sum of number_fragments_assigned from top_report
        num_reads = int(round(total_counts_from_report))
        print(f"No --nreads specified, using sum of number_fragments_assigned from top_report: {num_reads}")
    else:
        num_reads = args.nreads
        if args.use_abundance:
            print(f"Using user-specified read count ({num_reads}) with abundance-based proportions")
        else:
            print(f"Using user-specified read count: {num_reads}")

    if num_reads <= 0:
        raise SystemExit("Number of reads must be positive")

    # Get all accessions we need
    all_accessions = []
    for taxid in taxids_with_accessions:
        all_accessions.extend(taxid_to_accs[taxid])
    all_accessions = list(set(all_accessions))  # deduplicate

    print(f"Total accessions to extract: {len(all_accessions)}")

    # Extract sequences to a shared reference FASTA
    shared_reference = args.outdir / "reference_sequences.fasta"
    print(f"Extracting sequences to {shared_reference}...")
    extract_sequences(args.fasta, all_accessions, shared_reference)

    # Check if extraction was successful
    if not shared_reference.exists() or shared_reference.stat().st_size == 0:
        raise SystemExit(
            f"Failed to extract sequences. Check that accessions exist in FASTA file."
        )

    # Generate samples (abundance files + reference FASTAs only)
    manifest_path = args.outdir / "manifest.tsv"
    with manifest_path.open('w') as mf:
        mf.write("sample_id\tnum_taxids\tnum_reads\trandom_abundance\n")

        for i in range(1, args.nsamples + 1):
            sample_id = f"sample_{i:03d}"
            sample_dir = args.outdir / sample_id
            sample_dir.mkdir(parents=True, exist_ok=True)

            print(f"\n{'='*60}")
            print(f"Generating {sample_id} ({i}/{args.nsamples})")
            print(f"{'='*60}")

            # Create abundance file (used by BOTH ISS and NanoSim)
            abundance_file = sample_dir / "abundance.tsv"
            acc_abundance, total_val = create_abundance_file(
                filtered_values,
                taxid_to_accs,
                abundance_file,
                random_abundance=args.random_abundance
            )

            print(f"Created abundance file with {len(acc_abundance)} accessions")

            # Copy reference FASTA to sample directory
            sample_reference = sample_dir / "reference.fasta"
            shutil.copy2(str(shared_reference), str(sample_reference))

            # Create sample metadata
            metadata_path = sample_dir / "metadata.tsv"
            with metadata_path.open('w') as f:
                f.write("key\tvalue\n")
                f.write(f"sample_id\t{sample_id}\n")
                f.write(f"num_taxids\t{len(filtered_values)}\n")
                f.write(f"num_reads\t{num_reads}\n")
                f.write(f"num_reads_source\t{'top_report_sum' if args.nreads is None else 'user_specified'}\n")
                f.write(f"value_column\t{'abundance' if args.use_abundance else 'number_fragments_assigned'}\n")
                f.write(f"top_report_value_sum\t{total_counts_from_report}\n")
                f.write(f"random_abundance\t{args.random_abundance}\n")
                f.write(f"ranks\t{','.join(args.ranks)}\n")
                f.write(f"excluded_taxids\t{','.join(map(str, args.exclude_taxids))}\n")
                f.write(f"included_taxids\t{','.join(map(str, args.include_taxids))}\n")
                f.write(f"seed\t{args.seed}\n")

            # Create taxa report
            taxa_report = sample_dir / "taxa.tsv"
            with taxa_report.open('w') as f:
                f.write("taxid\tabundance_value\taccessions\n")
                for taxid in sorted(filtered_values.keys()):
                    accs = taxid_to_accs[taxid]
                    f.write(f"{taxid}\t{filtered_values[taxid]}\t{','.join(accs)}\n")

            # Update manifest
            mf.write(f"{sample_id}\t{len(filtered_values)}\t{num_reads}\t{args.random_abundance}\n")

            print(f"Sample {sample_id} complete!")
            print(f"  abundance.tsv : {abundance_file}")
            print(f"  reference.fasta: {sample_reference}")

    print(f"\n{'='*60}")
    print("All sample abundance files and references generated!")
    print(f"Manifest written to: {manifest_path}")
    print(f"Run ISS and/or NanoSim on the output abundance.tsv + reference.fasta files.")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
