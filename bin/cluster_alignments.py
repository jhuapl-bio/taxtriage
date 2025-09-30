#!/usr/bin/env python
"""
Extract "unique" reads from a BAM by clustering on reference positions, then write representatives to a FASTA.

Position-only clustering with optional splitting of long reads.

Usage examples
--------------
# Position-only clustering (exact start/end), max segment length 2000 bp:
python bam_unique_reads.py \
  --bam example.bam --out unique.fasta --window 25 --min-mapq 20 --max-length 2000

Dependencies
------------
- pysam: install with `pip install pysam`

Notes
-----
- By default we export the *aligned* portion of each read (no soft clips). Use
  --extract full to export the full query sequence instead.
- Supplementary/secondary alignments are excluded by default; enable via flags
  if you need them.
"""

import argparse
import os
import sys
from dataclasses import dataclass
from math import floor, ceil
from typing import Dict, List, Tuple, Iterable

try:
    import pysam
except ImportError:
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
    cluster_size: int


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Extract unique reads from BAM via positional clustering and split long reads."
    )
    p.add_argument("--bam", required=True, help="Input BAM (indexed).")
    p.add_argument("--out", required=True, help="Output FASTA path.")
    p.add_argument("--window", type=int, default=0,
                   help="Bin size (bp) for start/end when clustering (0 = exact).")
    p.add_argument("--min-cluster-size", type=int, default=1,
                   help="Minimum reads required to keep a positional cluster.")
    p.add_argument("--min-mapq", type=int, default=0,
                   help="Minimum MAPQ to include a read.")
    p.add_argument("--include-secondary", action="store_true",
                   help="Include secondary alignments.")
    p.add_argument("--include-supplementary", action="store_true",
                   help="Include supplementary alignments.")
    p.add_argument("--extract", choices=["aligned", "full"], default="aligned",
                   help="Export aligned portion or full query sequence.")
    p.add_argument("--top-unique-pct", type=float, default=0.0,
                   help="If >0, select top X%% most unique clusters per reference.")
    p.add_argument("--min-unique-per-ref", type=int, default=1,
                   help="Minimum clusters to keep per reference when using --top-unique-pct.")
    p.add_argument("--max-length", type=int, default=2000,
                   help="Maximum length of each output read; splits longer reads evenly.")
    return p.parse_args()
def bam_iter(bam_path: str,
             min_mapq: int = 0,
             include_secondary: bool = False,
             include_supplementary: bool = False,
             extract: str = "aligned") -> tuple[List[ReadRec], int]:
    reads: List[ReadRec] = []
    total: int = 0
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for aln in bam.fetch(until_eof=True):
            if aln.mapping_quality < min_mapq:
                continue
            if not include_secondary and aln.is_secondary:
                continue
            if not include_supplementary and aln.is_supplementary:
                continue

            is_unmapped = aln.is_unmapped

            # Choose sequence/quality: no aligned portion for unmapped
            if extract == "aligned" and not is_unmapped:
                seq = aln.query_alignment_sequence
                qual_list = aln.query_alignment_qualities or []
            else:
                seq = aln.query_sequence
                qual_list = aln.query_qualities or []

            if not seq:
                continue

            total += 1

            name = aln.query_name
            if aln.is_paired:
                if aln.is_read1:
                    name += "/1"
                elif aln.is_read2:
                    name += "/2"

            if is_unmapped:
                rname = "*"
                start = None
                end = None
                strand = "."
                aln_len = len(seq)
            else:
                # reference_id can be -1 for unmapped; we’re already in mapped branch
                rname = bam.get_reference_name(aln.reference_id)
                start, end = aln.reference_start, aln.reference_end
                strand = '-' if aln.is_reverse else '+'
                aln_len = (end - start) if end is not None and start is not None else 0

            qual_str = "".join(chr(q + 33) for q in qual_list) if qual_list else ""

            reads.append(ReadRec(name, seq, qual_str,
                                 rname, start, end, strand, aln.mapping_quality, aln_len))
    return reads, total
def bin_pos(v: int, window: int):
    if v is None:
        return None
    return int(floor(v / window) * window) if window > 0 else v

def position_clusters(reads: Iterable[ReadRec], window: int,
                      min_cluster_size: int) -> Dict[Tuple[str,int,int,str], List[ReadRec]]:
    clusters: Dict[Tuple[str,int,int,str], List[ReadRec]] = {}
    for r in reads:
        key = (r.rname, bin_pos(r.start, window), bin_pos(r.end, window), r.strand)
        clusters.setdefault(key, []).append(r)
    if min_cluster_size > 1:
        clusters = {k: v for k, v in clusters.items() if len(v) >= min_cluster_size}
    return clusters


def choose_representative(rs: List[ReadRec]) -> ReadRec:
    """Pick representative: highest MAPQ, then longest alignment."""
    return sorted(rs, key=lambda r: (r.mapq, r.aln_len), reverse=True)[0]


def split_read(r: ReadRec, max_length: int) -> List[ReadRec]:
    """
    Split a ReadRec into multiple segments if seq length > max_length.
    Splits into n = ceil(len/max_length) segments of approximately equal size.
    """
    seq, qual = r.seq, r.qual
    L = len(seq)
    if L <= max_length:
        return [r]
    n = ceil(L / max_length)
    parts: List[ReadRec] = []
    idx = 0
    for i in range(n):
        remain = L - idx
        left = n - i
        size = ceil(remain / left)
        seg_seq = seq[idx:idx+size]
        seg_qual = qual[idx:idx+size] if qual else ""
        part_name = f"{r.name}_part{i+1}"
        parts.append(ReadRec(
            part_name, seg_seq, seg_qual,
            r.rname, r.start, r.end, r.strand, r.mapq, len(seg_seq)
        ))
        idx += size
    return parts


def write_fasta(path: str, recs: Iterable[ReadRec]) -> None:
    with open(path, "w") as fh:
        for r in recs:
            header = f">{r.name}|{r.rname}:{r.start}-{r.end}({r.strand})|MAPQ={r.mapq}"
            fh.write(header + "\n")
            for i in range(0, len(r.seq), 80):
                fh.write(r.seq[i:i+80] + "\n")


def select_top_unique_by_ref(reps: List[RepRec], top_pct: float,
                             min_per_ref: int) -> List[RepRec]:
    if top_pct <= 0:
        return reps
    by_ref: Dict[str, List[RepRec]] = {}
    for rr in reps:
        by_ref.setdefault(rr.rec.rname, []).append(rr)
    kept: List[RepRec] = []
    for ref, items in by_ref.items():
        k = max(min_per_ref, ceil(len(items) * top_pct / 100.0))
        sorted_items = sorted(items, key=lambda x: (x.cluster_size, -x.rec.mapq, -x.rec.aln_len))
        kept.extend(sorted_items[:k])
    return kept


def main():
    args = parse_args()
    if not os.path.exists(args.bam):
        sys.stderr.write(f"ERROR: BAM not found: {args.bam}\n")
        sys.exit(2)

    reads, total = bam_iter(
        args.bam,
        min_mapq=args.min_mapq,
        include_secondary=args.include_secondary,
        include_supplementary=args.include_supplementary,
        extract=args.extract
    )
    print(f"Total reads: {total}")

    clusters = position_clusters(reads, window=args.window,
                                 min_cluster_size=args.min_cluster_size)
    reps: List[RepRec] = [RepRec(choose_representative(rs), len(rs)) for rs in clusters.values()]
    reps = select_top_unique_by_ref(reps, args.top_unique_pct, args.min_unique_per_ref)
    candidates: List[ReadRec] = [rr.rec for rr in reps]

    # Split long reads
    final_recs: List[ReadRec] = []
    for r in candidates:
        final_recs.extend(split_read(r, args.max_length))

    if not final_recs:
        sys.stderr.write("No reads passed filters; output empty FASTA.\n")
        open(args.out, "w").close()
        return

    write_fasta(args.out, final_recs)

if __name__ == "__main__":
    main()
