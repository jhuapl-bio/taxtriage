#!/usr/bin/env python3
"""
Prepare NanoSim metagenome inputs from ISS-style abundance files.

Converts accession-level abundance + reference FASTA into:
1. genome_list.tsv  - organism_name<TAB>path_to_organism_fasta
2. size_file.tsv    - header "Size<TAB>N" then organism_name<TAB>abundance_pct

Requires: pysam (for indexed FASTA extraction)
"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import pysam


def extract_sequences_for_accessions(fa: pysam.FastaFile, accessions: list,
                                     output_fasta: Path) -> int:
    """
    Extract sequences from an indexed FASTA for a list of accessions.

    Args:
        fa: Pre-opened pysam.FastaFile handle
        accessions: List of accession IDs to extract
        output_fasta: Path to write the output FASTA

    Returns:
        Number of sequences successfully extracted
    """
    available = set(fa.references)
    found = 0

    with output_fasta.open('w') as out:
        for acc in accessions:
            if acc in available:
                seq = fa.fetch(acc)
                out.write(f">{acc}\n{seq}\n")
                found += 1
            else:
                # Try partial match (some FASTAs have versioned accessions like ACC.1)
                for ref in available:
                    if ref.startswith(acc) or acc.startswith(ref):
                        seq = fa.fetch(ref)
                        out.write(f">{ref}\n{seq}\n")
                        found += 1
                        break

    return found


def main():
    parser = argparse.ArgumentParser(
        description="Prepare NanoSim metagenome simulation inputs from ISS-style abundance"
    )
    parser.add_argument(
        "--abundance", type=Path, required=True,
        help="ISS-style abundance file (accession<TAB>abundance)"
    )
    parser.add_argument(
        "--reference", type=Path, required=True,
        help="Reference FASTA file containing all sequences"
    )
    parser.add_argument(
        "--merged_taxid", type=Path, required=True,
        help="Merged taxid mapping file (Acc<TAB>Assembly<TAB>Organism_Name<TAB>...)"
    )
    parser.add_argument(
        "--num_reads", type=int, required=True,
        help="Total number of ONT reads to simulate"
    )
    parser.add_argument(
        "--outdir", type=Path, default=Path("."),
        help="Output directory (default: current directory)"
    )
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    genomes_dir = args.outdir / "genomes"
    genomes_dir.mkdir(parents=True, exist_ok=True)

    # 1. Parse merged_taxid to get accession -> organism mapping
    acc_to_org = {}
    with args.merged_taxid.open('r', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            acc = (row.get('Acc') or row.get('Accession') or '').strip()
            org = (row.get('Organism_Name') or row.get('organism') or '').strip()
            if acc and org:
                acc_to_org[acc] = org

    # 2. Parse abundance file (accession -> fractional abundance)
    acc_abundance = {}
    with args.abundance.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                try:
                    acc_abundance[parts[0]] = float(parts[1])
                except ValueError:
                    continue

    # 3. Aggregate accession abundance to organism level
    org_abundance = defaultdict(float)
    org_accs = defaultdict(list)
    for acc, abu in acc_abundance.items():
        org = acc_to_org.get(acc, acc)
        org_abundance[org] += abu
        org_accs[org].append(acc)

    # 4. Open reference FASTA once (pysam auto-creates .fai if missing)
    fa = pysam.FastaFile(str(args.reference))

    # 5. Create per-organism FASTAs and build genome_list + size_file
    genome_list_path = args.outdir / "genome_list.tsv"
    size_file_path = args.outdir / "size_file.tsv"

    genome_entries = []
    size_entries = []

    for org, abu in sorted(org_abundance.items()):
        slug = org.replace(' ', '_').replace('/', '_').replace('(', '').replace(')', '')
        slug = ''.join(c for c in slug if c.isalnum() or c in ('_', '-', '.'))
        org_fasta = genomes_dir / f"{slug}.fasta"
        accs = org_accs[org]

        found = extract_sequences_for_accessions(fa, accs, org_fasta)

        if found > 0 and org_fasta.exists() and org_fasta.stat().st_size > 0:
            genome_entries.append(f"{org}\t{org_fasta}")
            abu_pct = abu * 100.0
            size_entries.append(f"{org}\t{abu_pct:.6f}")
        else:
            print(f"WARNING: No sequences extracted for {org} ({len(accs)} accessions), skipping")

    fa.close()

    if not genome_entries:
        raise SystemExit("ERROR: No organism genomes were successfully extracted. "
                         "Check that accessions in abundance.tsv match the reference FASTA headers.")

    # Write genome_list.tsv
    with genome_list_path.open('w') as f:
        for entry in genome_entries:
            f.write(entry + '\n')

    # Write size_file.tsv (NanoSim format)
    with size_file_path.open('w') as f:
        f.write(f"Size\t{args.num_reads}\n")
        for entry in size_entries:
            f.write(entry + '\n')

    print(f"Prepared NanoSim inputs:")
    print(f"  genome_list.tsv : {len(genome_entries)} organisms")
    print(f"  size_file.tsv   : {args.num_reads} total reads")
    print(f"  Genomes dir     : {genomes_dir}")


if __name__ == "__main__":
    main()
