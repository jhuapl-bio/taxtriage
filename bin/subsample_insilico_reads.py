#!/usr/bin/env python3
"""
Subsample a master in-silico (synthetic) FASTQ into a series of derived datasets.

This powers two related use-cases for the TaxTriage in-silico workflow:

  1. Spike-in datasets (--mode consistent | randomized)
       consistent : nested prefixes of the master dataset. The dataset of size
                    n+1 contains exactly the same reads as the dataset of size n,
                    plus one more. (First-N reads of the master, in file order.)
       randomized : each dataset is an independent random draw from the master.
                    Two datasets of the same size will not share reads except by
                    chance.

  2. Dilution-series datasets (--replicates x, randomized mode)
       For each requested read count, produce x independently sampled datasets.

Every derived dataset is defined by a list of 0-based read indices into the
master FASTQ. Those indices are written to a tiny gzip'd ``.idx.gz`` file and
recorded in a manifest. Because the selection is fully described by the index
file, the (large) subsampled FASTQs can be deleted after analysis and rebuilt
later with ``reconstruct_insilico_reads.py`` from the master FASTQ + index file.

For paired-end data (R1/R2), an "index" refers to a read *pair*: the same record
position is taken from both mates so pairs stay together.

Outputs
-------
  <outdir>/<dataset_id>.subsample_R1.fastq.gz   (and _R2 for paired)
  <outdir>/<dataset_id>.subsample.fastq.gz      (single-end)
  <index_dir>/<dataset_id>.idx.gz               (0-based read indices, one per line)
  <manifest>                                    (TSV describing every dataset)

The dataset_id encodes everything needed to regroup datasets downstream:
  <parent>_ss_<mode>_c<count>_r<replicate>
"""

import argparse
import gzip
import random
import sys
from pathlib import Path
from typing import Dict, List, Optional


def open_maybe_gzip(path: Path, mode: str = "rt"):
    """Open a file transparently whether or not it is gzip-compressed."""
    p = str(path)
    if p.endswith(".gz"):
        return gzip.open(p, mode)
    return open(p, mode)


def count_records(fastq: Path) -> int:
    """Count FASTQ records (every 4 lines = 1 record)."""
    n_lines = 0
    with open_maybe_gzip(fastq, "rt") as fh:
        for _ in fh:
            n_lines += 1
    if n_lines % 4 != 0:
        raise ValueError(
            f"{fastq} has {n_lines} lines, not a multiple of 4 — is it a valid FASTQ?"
        )
    return n_lines // 4


def parse_counts(counts_str: str, total: int) -> List[int]:
    """Parse a comma/space separated list of read counts, clamped to `total`."""
    raw = [c for c in counts_str.replace(",", " ").split() if c]
    out = []
    for c in raw:
        try:
            v = int(c)
        except ValueError:
            raise ValueError(f"Invalid read count '{c}' in series '{counts_str}'")
        if v <= 0:
            continue
        out.append(min(v, total))
    # de-duplicate while preserving order, then sort ascending for nested logic
    seen = set()
    deduped = []
    for v in out:
        if v not in seen:
            seen.add(v)
            deduped.append(v)
    return sorted(deduped)


def build_index_sets(
    mode: str,
    counts: List[int],
    replicates: int,
    total: int,
    seed: int,
) -> List[Dict]:
    """
    Return a list of dataset descriptors, each:
        {count, replicate, indices (sorted list of 0-based record indices), seed}
    """
    datasets = []

    if mode == "consistent":
        # Nested prefixes in original file order. Replicates are meaningless
        # here (the prefix is deterministic), so we always emit a single one.
        for c in counts:
            datasets.append(
                {
                    "count": c,
                    "replicate": 1,
                    "indices": list(range(c)),  # already sorted
                    "seed": "NA",  # deterministic, no RNG used
                }
            )
    elif mode == "randomized":
        for c in counts:
            for r in range(1, replicates + 1):
                # Derive a reproducible per-(count, replicate) seed so the whole
                # series is reproducible from the single top-level --seed.
                sub_seed = hash((seed, c, r)) & 0xFFFFFFFF
                rng = random.Random(sub_seed)
                idx = sorted(rng.sample(range(total), c))
                datasets.append(
                    {
                        "count": c,
                        "replicate": r,
                        "indices": idx,
                        "seed": sub_seed,
                    }
                )
    else:
        raise ValueError(f"Unknown --mode '{mode}' (expected consistent|randomized)")

    return datasets


def write_index_file(index_path: Path, indices: List[int]) -> None:
    """Write sorted 0-based indices, one per line, gzip-compressed."""
    with gzip.open(index_path, "wt") as fh:
        for i in indices:
            fh.write(f"{i}\n")


def stream_write_datasets(
    master: Path,
    datasets: List[Dict],
    fastq_paths: List[Path],
    write_fastq: bool,
) -> None:
    """
    Single pass over the master FASTQ, writing each record to every dataset
    whose index set contains that record's position.

    `datasets` and `fastq_paths` are aligned lists.
    """
    if not write_fastq:
        return

    # Membership: for each dataset, a set for O(1) lookup during streaming.
    member_sets = [set(d["indices"]) for d in datasets]
    handles = [gzip.open(p, "wt") for p in fastq_paths]

    try:
        with open_maybe_gzip(master, "rt") as fh:
            rec_idx = 0
            while True:
                l1 = fh.readline()
                if not l1:
                    break
                l2 = fh.readline()
                l3 = fh.readline()
                l4 = fh.readline()
                record = l1 + l2 + l3 + l4
                for ds_i, members in enumerate(member_sets):
                    if rec_idx in members:
                        handles[ds_i].write(record)
                rec_idx += 1
    finally:
        for h in handles:
            h.close()


def main():
    parser = argparse.ArgumentParser(
        description="Subsample a master synthetic FASTQ into spike-in / dilution-series datasets."
    )
    parser.add_argument("--r1", required=True, type=Path, help="Master FASTQ (R1 or single-end).")
    parser.add_argument("--r2", type=Path, default=None, help="Master R2 FASTQ (paired-end only).")
    parser.add_argument(
        "--mode",
        choices=["consistent", "randomized"],
        default="randomized",
        help="consistent = nested prefixes; randomized = independent random draws.",
    )
    parser.add_argument(
        "--counts",
        required=True,
        help="Comma/space separated list of target read counts, e.g. '100,500,1000'.",
    )
    parser.add_argument(
        "--replicates",
        type=int,
        default=1,
        help="Datasets produced per count (randomized mode). Forced to 1 for consistent.",
    )
    parser.add_argument("--seed", type=int, default=42, help="RNG seed for randomized mode.")
    parser.add_argument("--parent", required=True, help="Parent sample id (prefix for dataset ids).")
    parser.add_argument("--outdir", type=Path, default=Path("datasets"), help="Output dir for FASTQs.")
    parser.add_argument("--index-dir", type=Path, default=Path("indices"), help="Output dir for index files.")
    parser.add_argument("--manifest", type=Path, default=Path("subsample_manifest.tsv"))
    parser.add_argument(
        "--no-fastq",
        action="store_true",
        help="Only write index files + manifest (do not materialise FASTQs).",
    )
    args = parser.parse_args()

    paired = args.r2 is not None
    args.outdir.mkdir(parents=True, exist_ok=True)
    args.index_dir.mkdir(parents=True, exist_ok=True)

    # ── Count master records (and sanity-check paired mates match) ────────────
    total = count_records(args.r1)
    if paired:
        total_r2 = count_records(args.r2)
        if total_r2 != total:
            raise ValueError(
                f"R1 has {total} records but R2 has {total_r2} — mates are out of sync."
            )
    if total == 0:
        raise ValueError(f"Master FASTQ {args.r1} contains no reads.")

    counts = parse_counts(args.counts, total)
    if not counts:
        raise ValueError(f"No valid read counts parsed from '{args.counts}'.")

    replicates = args.replicates if args.mode == "randomized" else 1
    if args.mode == "consistent" and args.replicates > 1:
        sys.stderr.write(
            "[subsample] NOTE: --mode consistent ignores --replicates (nested prefixes are deterministic).\n"
        )

    datasets = build_index_sets(args.mode, counts, replicates, total, args.seed)

    # ── Assign dataset ids + output paths, write index files ──────────────────
    fastq_paths_flat: List[Path] = []
    per_dataset_fastqs: List[Path] = []  # first fastq of each dataset (for streaming alignment)
    stream_targets: List[Dict] = []
    stream_paths: List[List[Path]] = []

    manifest_rows = []
    for d in datasets:
        ds_id = f"{args.parent}_ss_{args.mode}_c{d['count']}_r{d['replicate']}"
        idx_path = args.index_dir / f"{ds_id}.idx.gz"
        write_index_file(idx_path, d["indices"])

        if paired:
            r1_out = args.outdir / f"{ds_id}.subsample_R1.fastq.gz"
            r2_out = args.outdir / f"{ds_id}.subsample_R2.fastq.gz"
            out_paths = [r1_out, r2_out]
        else:
            se_out = args.outdir / f"{ds_id}.subsample.fastq.gz"
            out_paths = [se_out]

        manifest_rows.append(
            {
                "dataset_id": ds_id,
                "parent_id": args.parent,
                "platform": "paired" if paired else "single",
                "mode": args.mode,
                "target_count": d["count"],
                "actual_count": len(d["indices"]),
                "replicate": d["replicate"],
                "seed": d["seed"],
                "total_master_reads": total,
                "master_r1": args.r1.name,
                "master_r2": args.r2.name if paired else "NA",
                "index_file": idx_path.name,
            }
        )

        # For paired: R1 uses R1 master, R2 uses R2 master. Handle separately.
        if paired:
            stream_targets.append(d)
            stream_paths.append([out_paths[0]])  # R1
        else:
            stream_targets.append(d)
            stream_paths.append([out_paths[0]])

    write_fastq = not args.no_fastq

    # Stream R1 (or single-end) master once, distributing to all datasets' R1.
    stream_write_datasets(
        args.r1,
        stream_targets,
        [p[0] for p in stream_paths],
        write_fastq,
    )
    # Stream R2 master for paired datasets.
    if paired:
        r2_paths = [args.outdir / f"{row['dataset_id']}.subsample_R2.fastq.gz" for row in manifest_rows]
        stream_write_datasets(args.r2, stream_targets, r2_paths, write_fastq)

    # ── Write manifest ────────────────────────────────────────────────────────
    cols = [
        "dataset_id",
        "parent_id",
        "platform",
        "mode",
        "target_count",
        "actual_count",
        "replicate",
        "seed",
        "total_master_reads",
        "master_r1",
        "master_r2",
        "index_file",
    ]
    with open(args.manifest, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for row in manifest_rows:
            fh.write("\t".join(str(row[c]) for c in cols) + "\n")

    sys.stderr.write(
        f"[subsample] mode={args.mode} total_master={total} "
        f"datasets={len(datasets)} counts={counts} replicates={replicates} "
        f"fastq={'written' if write_fastq else 'skipped (index-only)'}\n"
    )


if __name__ == "__main__":
    main()
