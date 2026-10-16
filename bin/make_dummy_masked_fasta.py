#!/usr/bin/env python3
"""
make_dummy_masked_fasta.py -- write a masked copy of a reference FASTA for
testing TaxTriage region masking (--mask_fasta).

Regions either come from a TSV (-r: seq_id, start, end, kind, note; 0-based,
half-open; end may be 'end') or are drawn at random (--random), optionally
restricted to the accessions in a TaxTriage mapping file (--mapping, first
column). kind = hard (replace with N) or soft (lowercase). Sequences without
regions are copied unchanged. The regions actually applied are always written
next to the output as <out>.regions.tsv, so results can be checked.

  # fixed regions
  make_dummy_masked_fasta.py -i refs.fasta -r regions.tsv -o refs.masked.fasta.gz
  # 1-3 random regions per accession in the mapping file, reproducible
  make_dummy_masked_fasta.py -i refs.fasta --mapping run.dwnld.gcfmapping.tsv \\
      --random 3 --seed 42 -o refs.masked.fasta.gz
"""
import argparse
import gzip
import random
import sys
from collections import OrderedDict


def opener(p, mode):
    return gzip.open(p, mode + "t") if p.endswith(".gz") else open(p, mode)


def read_fasta(path):
    recs = OrderedDict()
    name, hdr, chunks = None, None, []
    with opener(path, "r") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    recs[name] = (hdr, "".join(chunks))
                hdr, name, chunks = line, line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        recs[name] = (hdr, "".join(chunks))
    return recs


def random_regions(recs, names, labels, a):
    """1..a.random non-overlapping regions per sequence.

    Region size scales with the sequence: a.min_len .. max(a.min_len,
    a.max_frac * length / a.random), capped at a.max_len, so a 30 kb virus gets
    a few hundred bp masked and a 5 Mb chromosome tens of kb. Each region is
    soft (lowercase) with probability a.soft_prob, otherwise hard (N).
    """
    rng = random.Random(a.seed)
    out = []
    for n in names:
        L = len(recs[n][1])
        if L < 2 * a.min_len:
            continue
        k = rng.randint(1, a.random)
        hi = int(min(a.max_len, max(a.min_len, a.max_frac * L / k)))
        taken = []
        for _ in range(k * 20):                    # retries to avoid overlaps
            if len(taken) == k:
                break
            size = rng.randint(a.min_len, max(a.min_len, hi))
            s = rng.randint(0, L - size)
            e = s + size
            if any(s < te and ts < e for ts, te in taken):
                continue
            taken.append((s, e))
        for s, e in sorted(taken):
            kind = "soft" if rng.random() < a.soft_prob else "hard"
            out.append([n, str(s), str(e), kind,
                        f"random (seed {a.seed}) - {labels.get(n, '')}".strip(" -")])
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--input", required=True, help="Reference FASTA(.gz)")
    ap.add_argument("-o", "--output", required=True, help="Masked FASTA(.gz)")
    ap.add_argument("-r", "--regions", help="Region TSV: seq_id start end kind note")
    ap.add_argument("--mapping", help="TaxTriage accession mapping TSV; random regions only on its accessions")
    ap.add_argument("--random", type=int, default=0, help="Up to N random regions per sequence")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--min_len", type=int, default=200, help="Smallest random region (bp)")
    ap.add_argument("--max_len", type=int, default=200_000, help="Largest random region (bp)")
    ap.add_argument("--max_frac", type=float, default=0.05,
                    help="Random regions cover at most about this fraction of each sequence")
    ap.add_argument("--soft_prob", type=float, default=0.2,
                    help="Probability a random region is soft-masked (lowercase) instead of N")
    ap.add_argument("--width", type=int, default=80)
    a = ap.parse_args()
    if not a.regions and not a.random:
        ap.error("give -r/--regions or --random N")

    recs = read_fasta(a.input)
    labels = {}
    names = list(recs)
    if a.mapping:
        with open(a.mapping) as fh:
            wanted = []
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f and f[0] and not line.startswith("#"):
                    wanted.append(f[0])
                    labels[f[0]] = f[-1] if len(f) > 1 else ""
        missing = [w for w in wanted if w not in recs]
        if missing:
            print(f"WARNING: mapping accessions not in FASTA: {', '.join(missing)}", file=sys.stderr)
        names = [w for w in wanted if w in recs]

    rows = []
    if a.regions:
        with open(a.regions) as fh:
            for line in fh:
                if line.strip() and not line.startswith("#"):
                    rows.append(line.rstrip("\n").split("\t"))
    if a.random:
        rows += random_regions(recs, names, labels, a)

    by_seq = {}
    for f in rows:
        by_seq.setdefault(f[0], []).append(f)
    missing = sorted(set(by_seq) - set(recs))
    if missing:
        print(f"WARNING: not in FASTA: {', '.join(missing)}", file=sys.stderr)

    stats = {"hard": 0, "soft": 0}
    with opener(a.output, "w") as out:
        for name, (hdr, seq) in recs.items():
            if name in by_seq:
                b = bytearray(seq, "ascii")
                for f in by_seq[name]:
                    s = max(0, min(int(f[1]), len(b)))
                    e = len(b) if f[2] == "end" else max(s, min(int(f[2]), len(b)))
                    if f[3] == "hard":
                        b[s:e] = b"N" * (e - s)
                    else:
                        b[s:e] = b[s:e].lower()
                    stats[f[3]] += e - s
                seq = b.decode()
            out.write(hdr)
            for i in range(0, len(seq), a.width):
                out.write(seq[i:i + a.width] + "\n")

    base = a.output[:-3] if a.output.endswith(".gz") else a.output
    with open(base + ".regions.tsv", "w") as fh:
        fh.write("#seq_id\tstart\tend\tkind\tnote\n")
        for f in rows:
            fh.write("\t".join(f) + "\n")
    print(f"{len(rows)} regions on {len(by_seq)} sequences: {stats['hard']:,} bp hard (N), "
          f"{stats['soft']:,} bp soft (lowercase) -> {a.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
