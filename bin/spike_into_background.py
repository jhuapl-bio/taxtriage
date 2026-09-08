#!/usr/bin/env python3
"""
Draw exact read counts from per-organism simulated pools and write the spike
payload for every (level, replicate) dataset of a spike-in series.

This does NOT touch the background FASTQ. The background is identical in every
dataset of a spike-in series — that is the whole point of the design — so the
Nextflow module concatenates it onto these spike payloads (gzip members
concatenate cleanly). Keeping the background out of Python means memory stays
flat no matter how deep the background is.

For each level L and replicate R, and each organism A spiked at count N:
  * take N reads (or read PAIRS, when paired) from A's simulated pool,
  * chosen without replacement using a seed derived from (seed, level, rep, A),
    so replicates draw different reads from the same pool and the whole thing is
    reproducible,
  * renamed <dataset>_spike_<A>_<i> so spiked reads are identifiable in the BAM
    and can never collide with background read names.

Emits one manifest row per dataset, recording exactly what went in — the report
reads this to tell a spike-in series from a depth series and to know the true
spiked count at each level.
"""

import argparse
import gzip
import hashlib
import io
import json
import os
import random
import sys
from collections import OrderedDict


def openr(path):
    """Open a FASTQ, gzipped or not."""
    if str(path).endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def count_records(path):
    """Number of FASTQ records, counted by lines/4 without holding the file."""
    n = 0
    with openr(path) as fh:
        for _ in fh:
            n += 1
    if n % 4:
        raise SystemExit(f"ERROR: {path} is not a well-formed FASTQ ({n} lines, not a multiple of 4)")
    return n // 4


def derived_seed(base, *parts):
    """A stable per-(dataset, organism) seed, so runs are reproducible."""
    h = hashlib.sha256(("|".join([str(base)] + [str(p) for p in parts])).encode()).hexdigest()
    return int(h[:16], 16)


def pick_indices(total, want, seed):
    """`want` distinct record indices from [0, total), or all of them if want >= total."""
    rng = random.Random(seed)
    if want >= total:
        return list(range(total)), True   # exhausted: caller warns
    return sorted(rng.sample(range(total), want)), False


def stream_take(pool_path, wanted_by_dataset, writers, name_prefix_by_dataset):
    """
    One pass over a pool FASTQ, writing the selected records to each dataset that
    wants them. `wanted_by_dataset` maps dataset -> sorted index list.
    """
    # position -> [(dataset, ordinal), ...]
    plan = {}
    for ds, idxs in wanted_by_dataset.items():
        for ordinal, i in enumerate(idxs):
            plan.setdefault(i, []).append((ds, ordinal))
    if not plan:
        return
    last = max(plan)
    with openr(pool_path) as fh:
        i = 0
        while i <= last:
            head = fh.readline()
            if not head:
                break
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            if not qual:
                break
            targets = plan.get(i)
            if targets:
                for ds, ordinal in targets:
                    w = writers.get(ds)
                    if w is None:
                        continue
                    w.write("@%s_%d\n%s%s%s" % (name_prefix_by_dataset[ds], ordinal, seq, plus, qual))
            i += 1


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--spikein", required=True,
                    help="Normalised spike-in TSV from parse_spikein_sheet.py")
    ap.add_argument("--pool", action="append", default=[], metavar="ACC=R1[,R2]",
                    help="Simulated read pool for one accession. Repeat per accession.")
    ap.add_argument("--parent", required=True,
                    help="Dataset id prefix, e.g. '<background>_background'")
    ap.add_argument("--mode", default="randomized", choices=["randomized", "consistent"],
                    help="Recorded in the dataset id; 'consistent' reuses replicate 1's draw")
    ap.add_argument("--paired", action="store_true", help="Pools and output are paired-end")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--outdir", default="spike")
    ap.add_argument("--manifest", default="spikein_manifest.tsv")
    ap.add_argument("--background-reads", type=int, default=0,
                    help="Record count of the background these payloads are mixed into")
    ap.add_argument("--background-name", default="", help="Recorded in the manifest")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # ── pools ────────────────────────────────────────────────────────────────
    pools = {}
    for spec in args.pool:
        if "=" not in spec:
            raise SystemExit(f"ERROR: --pool {spec!r} must look like ACC=R1[,R2]")
        acc, paths = spec.split("=", 1)
        parts = [p for p in paths.split(",") if p]
        if args.paired and len(parts) < 2:
            raise SystemExit(f"ERROR: --pool for {acc} needs R1,R2 in paired mode (got {paths!r})")
        pools[acc.strip()] = parts

    # ── spike sheet ──────────────────────────────────────────────────────────
    levels = OrderedDict()   # level_count -> {"reps": r, "members": [(acc, count, name)]}
    with open(args.spikein) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            lc = int(f[idx["level_count"]])
            entry = levels.setdefault(lc, {"reps": 1, "key": f[idx["level_key"]], "members": []})
            entry["reps"] = max(entry["reps"], int(f[idx["replicates"]]))
            entry["members"].append((f[idx["accession"]], int(f[idx["count"]]),
                                     f[idx["name"]] if "name" in idx and len(f) > idx["name"] else ""))

    missing = sorted({acc for e in levels.values() for acc, _, _ in e["members"]} - set(pools))
    if missing:
        raise SystemExit(
            "ERROR: no simulated pool was provided for: " + ", ".join(missing) +
            ". Every accession in the spike-in sheet must have reads simulated for it."
        )

    # ── plan every dataset ───────────────────────────────────────────────────
    datasets = []
    for lc in sorted(levels):
        e = levels[lc]
        for rep in range(1, e["reps"] + 1):
            ds = f"{args.parent}_ss_{args.mode}_c{lc}_r{rep}"
            datasets.append({"id": ds, "level_count": lc, "level_key": e["key"],
                             "replicate": rep, "members": e["members"]})

    # Per accession: which records each dataset takes.
    per_acc = {}       # acc -> {dataset_id: [indices]}
    exhausted = []
    pool_sizes = {}
    for acc, paths in pools.items():
        pool_sizes[acc] = count_records(paths[0])
    for d in datasets:
        for acc, count, _name in d["members"]:
            # 'consistent' makes every replicate draw the SAME reads, so the only
            # thing changing across a consistent series is the level.
            rep_for_seed = 1 if args.mode == "consistent" else d["replicate"]
            seed = derived_seed(args.seed, d["level_count"], rep_for_seed, acc)
            idxs, ran_out = pick_indices(pool_sizes[acc], count, seed)
            if ran_out:
                exhausted.append((d["id"], acc, count, pool_sizes[acc]))
            per_acc.setdefault(acc, {})[d["id"]] = idxs

    for ds, acc, want, have in exhausted:
        print(f"[spike] WARNING: {ds}: asked for {want} reads of {acc} but its pool holds "
              f"only {have}; using all of them. Raise --spikein_pool_factor or sim_nreads.",
              file=sys.stderr)

    # ── write ────────────────────────────────────────────────────────────────
    suffixes = ["_R1", "_R2"] if args.paired else [""]
    writers = {sfx: {} for sfx in suffixes}
    try:
        for d in datasets:
            for si, sfx in enumerate(suffixes):
                path = os.path.join(args.outdir, f"{d['id']}.spike{sfx}.fastq.gz")
                writers[sfx][d["id"]] = gzip.open(path, "wt")

        for acc in sorted(per_acc):
            for si, sfx in enumerate(suffixes):
                pool_path = pools[acc][si]
                prefixes = {ds: f"{ds}_spike_{acc}" for ds in per_acc[acc]}
                stream_take(pool_path, per_acc[acc], writers[sfx], prefixes)
    finally:
        for sfx in suffixes:
            for w in writers[sfx].values():
                try:
                    w.close()
                except Exception:
                    pass

    # ── manifest ─────────────────────────────────────────────────────────────
    # target_count / actual_count keep the same meaning the depth-series manifest
    # gives them (what was asked for, what was delivered) so the report's existing
    # reader works unchanged; kind + spiked_count tell it this is a spike series.
    # level_count is the c<N> in the dataset id (the series x position); target vs
    # actual are the dataset's TOTAL requested and delivered spike, so they compare
    # like for like exactly as they do for a depth series.
    cols = ["dataset_id", "parent_id", "kind", "platform", "mode", "level_key",
            "level_count", "target_count", "actual_count", "replicate", "seed",
            "spiked_count", "background_reads", "total_master_reads",
            "background_name", "spike_detail"]
    with open(args.manifest, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for d in datasets:
            detail = []
            actual = 0
            requested = 0
            for acc, count, name in d["members"]:
                got = len(per_acc[acc][d["id"]])
                actual += got
                requested += count
                detail.append({"accession": acc, "requested": count, "spiked": got, "name": name})
            row = {
                "dataset_id": d["id"],
                "parent_id": args.parent,
                "kind": "spikein",
                "platform": "paired" if args.paired else "single",
                "mode": args.mode,
                "level_key": d["level_key"],
                "level_count": d["level_count"],
                "target_count": requested,
                "actual_count": actual,
                "replicate": d["replicate"],
                "seed": args.seed,
                "spiked_count": actual,
                "background_reads": args.background_reads,
                "total_master_reads": args.background_reads + actual,
                "background_name": args.background_name,
                "spike_detail": json.dumps(detail, separators=(",", ":")),
            }
            fh.write("\t".join(str(row[c]) for c in cols) + "\n")

    print(f"[spike] {len(datasets)} dataset(s) across {len(levels)} level(s); "
          f"pools: " + ", ".join(f"{a}={pool_sizes[a]}" for a in sorted(pool_sizes)),
          file=sys.stderr)


if __name__ == "__main__":
    main()
