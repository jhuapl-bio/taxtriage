#!/usr/bin/env python3
"""
Reconstruct a subsampled in-silico FASTQ from a master FASTQ + a read-index file.

The TaxTriage subsampling step (subsample_insilico_reads.py) describes every
derived dataset by a gzip'd list of 0-based read indices into the master FASTQ.
That lets you delete the (large) subsampled FASTQs after analysis and rebuild
them on demand from the much smaller master + index file.

Usage
-----
  # Single-end
  reconstruct_insilico_reads.py \\
      --index sample_ss_randomized_c1000_r1.idx.gz \\
      --r1 master.fastq.gz \\
      --out-r1 rebuilt.fastq.gz

  # Paired-end
  reconstruct_insilico_reads.py \\
      --index sample_ss_randomized_c1000_r1.idx.gz \\
      --r1 master_R1.fastq.gz --r2 master_R2.fastq.gz \\
      --out-r1 rebuilt_R1.fastq.gz --out-r2 rebuilt_R2.fastq.gz

An index refers to a read *pair* for paired data (same record position in R1/R2).
"""

import argparse
import gzip
from pathlib import Path
from typing import List


def open_maybe_gzip(path: Path, mode: str = "rt"):
    p = str(path)
    if p.endswith(".gz"):
        return gzip.open(p, mode)
    return open(p, mode)


def load_indices(index_path: Path) -> set:
    idx = set()
    with open_maybe_gzip(index_path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if line:
                idx.add(int(line))
    return idx


def extract(master: Path, out: Path, members: set) -> int:
    """Stream master, write records whose position is in `members`. Returns count."""
    written = 0
    with open_maybe_gzip(master, "rt") as fh, gzip.open(out, "wt") as oh:
        rec_idx = 0
        while True:
            l1 = fh.readline()
            if not l1:
                break
            l2 = fh.readline()
            l3 = fh.readline()
            l4 = fh.readline()
            if rec_idx in members:
                oh.write(l1 + l2 + l3 + l4)
                written += 1
            rec_idx += 1
    return written


def main():
    parser = argparse.ArgumentParser(
        description="Rebuild a subsampled FASTQ from a master FASTQ + read-index file."
    )
    parser.add_argument("--index", required=True, type=Path, help="Read-index file (.idx.gz).")
    parser.add_argument("--r1", required=True, type=Path, help="Master R1/single-end FASTQ.")
    parser.add_argument("--r2", type=Path, default=None, help="Master R2 FASTQ (paired).")
    parser.add_argument("--out-r1", required=True, type=Path, help="Output R1/single-end FASTQ.gz.")
    parser.add_argument("--out-r2", type=Path, default=None, help="Output R2 FASTQ.gz (paired).")
    args = parser.parse_args()

    members = load_indices(args.index)
    if not members:
        raise ValueError(f"Index file {args.index} contained no indices.")

    n1 = extract(args.r1, args.out_r1, members)
    if n1 != len(members):
        raise ValueError(
            f"Expected {len(members)} reads from R1 but wrote {n1} — "
            f"master may not match this index file."
        )

    if args.r2 is not None:
        if args.out_r2 is None:
            raise ValueError("--r2 given but --out-r2 missing.")
        n2 = extract(args.r2, args.out_r2, members)
        if n2 != len(members):
            raise ValueError(
                f"Expected {len(members)} reads from R2 but wrote {n2}."
            )

    print(f"Reconstructed {n1} reads into {args.out_r1}" + (f" / {args.out_r2}" if args.r2 else ""))


if __name__ == "__main__":
    main()
