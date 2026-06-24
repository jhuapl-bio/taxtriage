#!/usr/bin/env python3
##############################################################################################
# Copyright 2022 The Johns Hopkins University Applied Physics Laboratory LLC
# (header trimmed -- copy the full APL header from bin/novelty_score.py)
##############################################################################################
"""
Propagate a representative's taxonomy call to every member of its linclust cluster.

The novelty branch clusters predicted ORFs with mmseqs linclust and runs the (expensive)
translated-search taxonomy on the cluster REPRESENTATIVES only. This script expands those
per-representative results back to per-member rows so the downstream novelty score sees
member-level abundance instead of representative-level abundance.

mmseqs reduces FASTA headers to their first whitespace-delimited token in both the cluster
TSV and the createtsv/convertalis outputs, so the representative key matches across files.

Inputs
  --clusters   mmseqs linclust *_cluster.tsv : representative<TAB>member (rep is a member of itself)
  --lca        mmseqs createtsv output       : query, taxid, rank, name, [lineage]
  --tophit     mmseqs convertalis output     : query, target, pident, ...   (optional)

Outputs (same column layout as the inputs, one row per member)
  --out-lca, --out-tophit
"""
import argparse
import sys
from collections import defaultdict


def load_clusters(path):
    """representative -> [members]  (mmseqs lists the rep as a member of itself)."""
    rep_members = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            rep, member = parts[0], parts[1]
            rep_members[rep].append(member)
    return rep_members


def load_keyed_best(path, key_col=0):
    """First (best) row per key. mmseqs emits hits best-first, so first-seen == best."""
    d = {}
    if not path:
        return d
    try:
        with open(path) as fh:
            for line in fh:
                if not line.strip():
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) <= key_col:
                    continue
                k = cols[key_col]
                if k not in d:
                    d[k] = cols
    except FileNotFoundError:
        return {}
    return d


def main(argv=None):
    ap = argparse.ArgumentParser(description="Expand linclust representative taxonomy to members.")
    ap.add_argument("--clusters", required=True)
    ap.add_argument("--lca", required=True)
    ap.add_argument("--tophit", default=None)
    ap.add_argument("--out-lca", required=True)
    ap.add_argument("--out-tophit", required=True)
    a = ap.parse_args(argv)

    rep_members = load_clusters(a.clusters)
    lca = load_keyed_best(a.lca)
    top = load_keyed_best(a.tophit) if a.tophit else {}

    n_members = n_lca = n_top = 0
    with open(a.out_lca, "w") as ol, open(a.out_tophit, "w") as ot:
        for rep, members in rep_members.items():
            lrow = lca.get(rep)
            trow = top.get(rep)
            for m in members:
                n_members += 1
                if lrow is not None:
                    out = list(lrow)
                    out[0] = m
                    ol.write("\t".join(out) + "\n")
                    n_lca += 1
                if trow is not None:
                    out = list(trow)
                    out[0] = m
                    ot.write("\t".join(out) + "\n")
                    n_top += 1

    sys.stderr.write(
        f"[expand_clusters] members={n_members} expanded_lca_rows={n_lca} "
        f"expanded_tophit_rows={n_top} (reps with lca={len(lca)})\n"
    )


if __name__ == "__main__":
    sys.exit(main())
