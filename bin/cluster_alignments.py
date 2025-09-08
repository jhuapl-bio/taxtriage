#!/usr/bin/env python3
"""
Extract "unique" reads from a BAM by clustering on reference positions and (optionally)
sequence similarity via MMseqs2, then write representatives to a FASTA.

Two stages (you can use either or both):
  1) Position clustering: group reads by (ref, start, end, strand) with an optional
     tolerance window so similar starts/ends fall into the same bin.
  2) Sequence clustering (optional): run `mmseqs easy-cluster` on the candidate reads
     and keep the representative sequence of each cluster.

Examples
--------
# 1) Position-only clustering (exact start/end):
python bam_unique_reads.py \
  --bam example.bam --out unique.fasta --method pos --min-mapq 20

# 2) Position clustering using 25bp window, then MMseqs2 at 97% id and 0.8 coverage:
python bam_unique_reads.py \
  --bam example.bam --out unique.fasta --method pos+mmseqs \
  --window 25 --min-mapq 10 --min-seq-id 0.97 --cov 0.8

# 3) Run MMseqs2 on all aligned reads directly (no position clustering):
python bam_unique_reads.py \
  --bam example.bam --out unique.fasta --method mmseqs --min-seq-id 0.99 --cov 0.9

Dependencies
------------
- pysam (required): `pip install pysam`
- MMseqs2 (optional, for --method mmseqs or pos+mmseqs): tool must be in $PATH

Notes
-----
- By default we export the *aligned* portion of each read (no soft clips). Use
  --extract full to export the full query sequence instead.
- Supplementary/secondary alignments are excluded by default; enable via flags
  if you need them.
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from math import ceil
from math import floor
from typing import Dict, List, Tuple, Iterable

try:
    import pysam
except ImportError as e:
    sys.stderr.write("ERROR: pysam is required. Install with `pip install pysam`.\n")
    raise


@dataclass
class ReadRec:
    name: str
    seq: str
    qual: str
    rname: str
    start: int
    end: int
    strand: str
    mapq: int
    aln_len: int

@dataclass
class RepRec:
    rec: ReadRec
    cluster_size: int  # number of reads in the positional cluster this rec represents


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Extract unique reads from BAM via positional and/or MMseqs2 clustering.")
    p.add_argument("--bam", required=True, help="Input BAM (indexed).")
    p.add_argument("--out", required=True, help="Output FASTA path.")
    p.add_argument("--method", choices=["pos", "mmseqs", "pos+mmseqs"], default="pos+mmseqs",
                   help="Clustering strategy: position-only, MMseqs-only, or both.")
    p.add_argument("--window", type=int, default=0,
                   help="Bin size (bp) for start/end when forming positional clusters (0 = exact).")
    p.add_argument("--min-cluster-size", type=int, default=1,
                   help="Minimum reads required to keep a positional cluster.")
    p.add_argument("--min-mapq", type=int, default=0, help="Minimum MAPQ to include a read.")
    p.add_argument("--include-secondary", action="store_true", help="Include secondary alignments.")
    p.add_argument("--include-supplementary", action="store_true", help="Include supplementary alignments.")
    p.add_argument("--extract", choices=["aligned", "full"], default="aligned",
                   help="Use only the aligned portion (no soft clips) or full query sequence.")

    # MMseqs params
    p.add_argument("--min-seq-id", type=float, default=0.97, help="MMseqs2: minimum sequence identity (0-1).")
    p.add_argument("--cov", type=float, default=0.8, help="MMseqs2: minimum coverage fraction (0-1).")
    p.add_argument("--mmseqs-threads", type=int, default=1, help="MMseqs2: number of threads.")
    p.add_argument("--keep-temp", action="store_true", help="Do not remove temp directory (debugging).")

    # Uniqueness selection
    p.add_argument("--top-unique-pct", type=float, default=0.0,
                   help="If >0, select the top X percent most unique positional clusters per reference (e.g., 5 for top 5%).")
    p.add_argument("--min-unique-per-ref", type=int, default=1,
                   help="Minimum number of clusters to keep per reference when using --top-unique-pct.")

    return p.parse_args()


def mmseqs_available() -> bool:
    return shutil.which("mmseqs") is not None


def bam_iter(bam_path: str,
             min_mapq: int = 0,
             include_secondary: bool = False,
             include_supplementary: bool = False,
             extract: str = "aligned") -> tuple[list[ReadRec], int]:
    """Yield ReadRec for each primary alignment passing filters."""
    reads = []
    total_reads = 0
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for aln in bam.fetch(until_eof=True):
            if aln.is_unmapped:
                continue
            if aln.mapping_quality < min_mapq:
                continue
            if (not include_secondary) and aln.is_secondary:
                continue
            if (not include_supplementary) and aln.is_supplementary:
                continue

            rname = bam.get_reference_name(aln.reference_id)
            start = int(aln.reference_start)
            end = int(aln.reference_end)
            strand = "-" if aln.is_reverse else "+"
            mapq = int(aln.mapping_quality)
            aln_len = end - start if end is not None and start is not None else 0

            if extract == "aligned":
                seq = aln.query_alignment_sequence
                qual = (aln.query_alignment_qualities or [])
            else:
                seq = aln.query_sequence
                qual = (aln.query_qualities or [])

            if not seq:
                continue

            total_reads += 1
            name = aln.query_name
            if aln.is_paired:
                if aln.is_read1:
                    name = f"{name}/1"
                elif aln.is_read2:
                    name = f"{name}/2"

            qual_str = "".join(chr(q + 33) for q in qual) if qual else ""
            reads.append(ReadRec(name, seq, qual_str, rname, start, end, strand, mapq, aln_len))

    return reads, total_reads

def bin_pos(v: int, window: int) -> int:
    if window <= 0:
        return v
    return int(floor(v / window) * window)


def position_clusters(reads: Iterable[ReadRec], window: int, min_cluster_size: int) -> Dict[Tuple[str, int, int, str], List[ReadRec]]:
    clusters: Dict[Tuple[str, int, int, str], List[ReadRec]] = {}
    for r in reads:
        s_bin = bin_pos(r.start, window)
        e_bin = bin_pos(r.end, window)
        key = (r.rname, s_bin, e_bin, r.strand)
        clusters.setdefault(key, []).append(r)
    # filter small clusters
    if min_cluster_size > 1:
        clusters = {k: v for k, v in clusters.items() if len(v) >= min_cluster_size}
    return clusters


def choose_representative(rs: List[ReadRec]) -> ReadRec:
    """Pick a representative read from a positional cluster: highest MAPQ, then longest alignment."""
    return sorted(rs, key=lambda r: (r.mapq, r.aln_len), reverse=True)[0]


def write_fasta(path: str, recs: Iterable[ReadRec]) -> None:
    with open(path, "w") as fh:
        for r in recs:
            header = f">{r.name}|{r.rname}:{r.start}-{r.end}({r.strand})|MAPQ={r.mapq}"
            fh.write(header + "\n")
            seq = r.seq
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i+80] + "\n")


def select_top_unique_by_ref(reps: List[RepRec], top_pct: float, min_per_ref: int) -> List[RepRec]:
    """Keep the most unique positional clusters per reference.

    Uniqueness metric: smaller positional-cluster size = more unique.
    For each reference, keep max(min_per_ref, ceil(N * top_pct/100)).
    """
    if top_pct <= 0:
        return reps

    by_ref: Dict[str, List[RepRec]] = {}
    for rr in reps:
        by_ref.setdefault(rr.rec.rname, []).append(rr)

    kept: List[RepRec] = []
    for ref, items in by_ref.items():
        n = len(items)
        k = max(min_per_ref, ceil(n * (top_pct / 100.0)))
        # smaller cluster_size is more unique
        items_sorted = sorted(items, key=lambda x: (x.cluster_size, -x.rec.mapq, -x.rec.aln_len))
        kept.extend(items_sorted[:k])
    return kept


def run_mmseqs_easy_cluster(in_fasta: str, out_dir: str, min_seq_id: float, cov: float, threads: int) -> Tuple[str, str]:
    """Run mmseqs easy-cluster and return (rep_fasta, tsv) paths.
    Raises subprocess.CalledProcessError on failure.
    """
    out_base = os.path.join(out_dir, "mmseqs_out")
    tmp_dir = os.path.join(out_dir, "mmseqs_tmp")
    os.makedirs(tmp_dir, exist_ok=True)

    cmd = [
        "mmseqs", "easy-cluster", in_fasta, out_base, tmp_dir,
        "--min-seq-id", str(min_seq_id), "-c", str(cov), "--threads", str(threads)
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    rep_fa = out_base + "_rep_seq.fasta"
    clu_tsv = out_base + "_cluster.tsv"
    if not os.path.exists(rep_fa):
        raise FileNotFoundError(f"Expected MMseqs2 output not found: {rep_fa}")
    if not os.path.exists(clu_tsv):
        # Newer MMseqs may emit *_cluster.tsv OR *_all_seqs.fasta only; tolerate missing TSV
        clu_tsv = ""
    return rep_fa, clu_tsv


def main():
    args = parse_args()

    if not os.path.exists(args.bam):
        sys.stderr.write(f"ERROR: BAM not found: {args.bam}\n")
        sys.exit(2)

    # Stage A: prepare candidate sequences (either all reads, or positional representatives)
    reads_iter, tt_reads = bam_iter(
        args.bam,
        min_mapq=args.min_mapq,
        include_secondary=args.include_secondary,
        include_supplementary=args.include_supplementary,
        extract=args.extract,
    )
    print(f"Total reads: {tt_reads}")

    # For uniqueness selection, we need positional cluster sizes. Build clusters if method includes pos.
    reps: List[RepRec] = []
    if args.method in ("pos", "pos+mmseqs"):
        clusters = position_clusters(reads_iter, window=args.window, min_cluster_size=args.min_cluster_size)
        for key, rs in clusters.items():
            rep = choose_representative(rs)
            reps.append(RepRec(rep, cluster_size=len(rs)))
        # Select top unique clusters per reference if requested
        reps = select_top_unique_by_ref(reps, args.top_unique_pct, args.min_unique_per_ref)
        candidates: List[ReadRec] = [rr.rec for rr in reps]
    else:
        # mmseqs only: take all reads (no positional uniqueness available)
        candidates = list(reads_iter)


    if not candidates:
        sys.stderr.write("No reads passed the filters; nothing to write.\n")
        # Still create an empty file for predictability
        open(args.out, "w").close()
        return

    # If method is position-only, we're done
    if args.method == "pos":
        # candidates already filtered by top-unique-pct if provided
        write_fasta(args.out, candidates)
        return

    # Otherwise, run MMseqs2 on the candidates
    if not mmseqs_available():
        sys.stderr.write("WARNING: mmseqs not found in PATH. Falling back to position-only representatives.\n")
        write_fasta(args.out, candidates)
        return

    with tempfile.TemporaryDirectory(prefix="bamuniq_") as tdir:
        if args.keep_temp:
            # Re-create persistent temp dir by copying at end (simpler than disabling context manager)
            persistent = tdir
        in_fa = os.path.join(tdir, "candidates.fasta")
        write_fasta(in_fa, candidates)

        try:
            rep_fa, _ = run_mmseqs_easy_cluster(
                in_fasta=in_fa,
                out_dir=tdir,
                min_seq_id=args.min_seq_id,
                cov=args.cov,
                threads=args.mmseqs_threads,
            )
        except subprocess.CalledProcessError as e:
            sys.stderr.write("ERROR: MMseqs2 failed. Falling back to position-only representatives.\n")
            sys.stderr.write(e.stderr.decode(errors="ignore") if e.stderr else str(e) + "\n")
            write_fasta(args.out, candidates)
            return
        except Exception as e:
            sys.stderr.write(f"ERROR: {e}. Falling back to position-only representatives.\n")
            write_fasta(args.out, candidates)
            return

        # Copy representative fasta to desired output
        shutil.copyfile(rep_fa, args.out)

        if args.keep_temp:
            keep_dir = os.path.abspath("mmseqs_debug")
            os.makedirs(keep_dir, exist_ok=True)
            for fn in os.listdir(tdir):
                src = os.path.join(tdir, fn)
                dst = os.path.join(keep_dir, fn)
                try:
                    if os.path.isdir(src):
                        shutil.copytree(src, dst, dirs_exist_ok=True)
                    else:
                        shutil.copy2(src, dst)
                except Exception:
                    pass
            sys.stderr.write(f"Kept intermediate files in: {keep_dir}\n")


if __name__ == "__main__":
    main()
