#!/usr/bin/env python3
"""
Keep only the named records of a FASTA (spike-in `record` column).

Ids match the header's first token exactly, or without the version suffix
(NC_003310 selects NC_003310.1). Every requested id must be found — a typo
should fail the run rather than silently spike nothing.
"""
import argparse
import re
import sys


def base(x):
    return re.sub(r"\.\d+$", "", x)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--records", required=True, help="';' or ',' separated record ids")
    ap.add_argument("--output", required=True)
    a = ap.parse_args()

    want = [r for r in re.split(r"[;,\s]+", a.records) if r]
    exact = set(want)
    loose = {base(w): w for w in want}
    hit = set()
    keep = False
    seen_ids = []
    with open(a.fasta) as fin, open(a.output, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                rid = line[1:].split(None, 1)[0] if line[1:].strip() else ""
                seen_ids.append(rid)
                m = rid if rid in exact else loose.get(base(rid)) if base(rid) in loose else None
                keep = m is not None
                if keep:
                    hit.add(m)
            if keep:
                fout.write(line)
    missing = [w for w in want if w not in hit]
    if missing:
        sys.exit(f"ERROR: record(s) not found in {a.fasta}: {', '.join(missing)}\n"
                 f"       it holds {len(seen_ids)} record(s): {', '.join(seen_ids[:12])}"
                 + (" ..." if len(seen_ids) > 12 else ""))
    print(f"[extract_records] kept {len(hit)} of {len(seen_ids)} record(s): {', '.join(want)}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
