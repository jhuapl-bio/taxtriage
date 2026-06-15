#!/usr/bin/env python3
import argparse, json, sys, re
from collections import defaultdict, OrderedDict
import pandas as pd
import os


# Reference accession embedded in a clustered-read sequence_id, e.g.
#   "BA000007.2_1297_1/1|NC_000913.3:780733-781131(+)|MAPQ=60"  ->  "NC_000913.3"
# The accession is the token immediately preceding a "<start>-<end>" coordinate span.
_ACC_RE = re.compile(r'([A-Za-z0-9_.|]*?)([A-Za-z][A-Za-z0-9_]*\.\d+):\d+-\d+')


def parse_accession(sequence_id):
    """Return the reference accession parsed from a sequence_id, or None."""
    if not sequence_id:
        return None
    for part in str(sequence_id).split('|'):
        m = re.match(r'^([A-Za-z][A-Za-z0-9_]*\.\d+):\d+-\d+', part.strip())
        if m:
            return m.group(1)
    # Fallback: search anywhere for "<acc>:<start>-<end>"
    m = _ACC_RE.search(str(sequence_id))
    return m.group(2) if m else None


def finest_confidence(ann, ranks_fine_to_coarse):
    """Pick the model's most specific *concrete* top-1 probability for one record.

    `ann` is {rank: (label, prob)}. Prefer the finest rank with a concrete label
    (not None / not the "-1" unknown sentinel); fall back to the finest rank with
    any numeric probability. Returns float prob or None.
    """
    fallback = None
    for rank in ranks_fine_to_coarse:
        label, prob = ann.get(rank, (None, None))
        if prob is None:
            continue
        try:
            p = float(prob)
        except (TypeError, ValueError):
            continue
        if fallback is None:
            fallback = p
        if label is not None and str(label) not in ("-1", ""):
            return p
    return fallback


def read_fasta_headers_in_order(fa_path):
    """Return a list of FASTA headers (without '>') in the order they appear."""
    headers = []
    with open(fa_path, 'r') as fh:
        for line in fh:
            if line.startswith('>'):
                headers.append(line[1:].strip())
    return headers

def read_clusters_tsv(tsv_path):
    """
    Read clusters.tsv from mmseqs easy-cluster:
      representative<TAB>member
    Returns: dict rep -> list[members]
    """
    rep_to_members = defaultdict(list)
    with open(tsv_path, 'r') as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                # be tolerant but warn
                sys.stderr.write(f"WARNING: Skipping malformed line (expected 2+ columns): {line}\n")
                continue
            rep, member = parts[0].strip(), parts[1].strip()
            rep_to_members[rep].append(member)
    return rep_to_members

def read_json_labels(json_path):
    """
    JSON structure expected:
    [
      [ { 'label_column': 'superkingdom', 'topk_predictions': [{'label': 2, 'probability': 0.99}] },
        { 'label_column': 'phylum',        'topk_predictions': [{'label': 1224, 'probability': 0.94}] }, ... ],
      [ ... for next representative in FASTA order ... ],
      ...
    ]
    Returns:
      annotations_by_index: list[dict(rank -> (label, prob))]
      all_ranks: ordered list of unique rank names in first-seen order
    """
    with open(json_path, 'r') as fh:
        data = json.load(fh)

    all_ranks_order = OrderedDict()  # preserves first-seen order
    annotations_by_index = []
    for entry_idx, entry in enumerate(data):
        rank_map = {}
        if entry is None:
            annotations_by_index.append(rank_map)
            continue
        for item in entry.get('predictions', []):
            # be robust to missing keys
            rank = item.get('label_column')
            topk = item.get('topk_predictions') or []
            if not rank or not topk:
                continue
            # take top-1
            top = topk[0]
            label = top.get('label')
            prob = top.get('probability')
            rank_map[rank] = (label, prob)
            all_ranks_order.setdefault(rank, None)
        annotations_by_index.append(rank_map)

    return annotations_by_index, list(all_ranks_order.keys())

def main():
    ap = argparse.ArgumentParser(description="Map representative sequence JSON labels to cluster members and output TSV.")
    ap.add_argument("--fasta", required=True, help="Representative sequences FASTA (headers in order).")
    ap.add_argument("--json", required=True, help="JSON with per-representative predictions (indexed to FASTA order).")
    ap.add_argument("--clusters", required=True, help="MMseqs clusters.tsv (representative<TAB>member).")
    ap.add_argument("--out", required=True, help="Output TSV path.")
    ap.add_argument("--modelname", required=False, help="Name of the MicrobeRT model.")
    ap.add_argument("--include-prob", action="store_true",
                    help="Also include *_prob columns for each rank (top-1 probability).")
    # NEW: optional taxa report
    ap.add_argument("--taxa-report", default=None,
                    help="Optional TSV path for rank × taxid probability profile (avg, median, std).")
    # NEW: optional per-source-accession report. Aggregates the model's per-sequence
    # confidence by the reference accession the sequence was aligned to (parsed from
    # the sequence_id). This lets downstream steps attach a MicrobeRT probability to
    # each DETECTED organism via its reference accession(s), instead of relying on the
    # model's predicted-taxid namespace (which need not match NCBI taxids).
    ap.add_argument("--accession-report", default=None,
                    help="Optional TSV path for per-accession probability profile "
                         "(accession, avg, median, std, n, model).")
    args = ap.parse_args()

    # 1) representatives in order from FASTA
    reps_in_order = read_fasta_headers_in_order(args.fasta)
    if not reps_in_order:
        sys.stderr.write("ERROR: No FASTA headers found.\n")
        sys.exit(2)

    # 2) JSON annotations (indexed to FASTA order)
    annotations_by_index, all_ranks = read_json_labels(args.json)
    if not annotations_by_index:
        sys.stderr.write("ERROR: JSON appears empty or malformed.\n")
        sys.exit(2)

    if len(annotations_by_index) != len(reps_in_order):
        sys.stderr.write(
            f"WARNING: JSON entries ({len(annotations_by_index)}) != FASTA headers ({len(reps_in_order)}). "
            f"Proceeding with min length.\n"
        )
    N = min(len(annotations_by_index), len(reps_in_order))

    # 3) clusters map
    rep_to_members = read_clusters_tsv(args.clusters)

    # 4) Build rows
    rows = []
    modelname = os.path.basename(args.modelname) if args.modelname else None
    for i in range(N):
        rep = reps_in_order[i]
        ann = annotations_by_index[i]  # dict rank -> (label, prob)
        # cluster members for this representative (default: just the rep itself if not in file)
        members = rep_to_members.get(rep, [rep])

        # Prepare an output row per member with rank columns
        for member in members:
            row = {
                "sequence_id": member,
                "representative_id": rep
            }
            # fill all ranks; if missing, None
            for rank in all_ranks:
                label, prob = ann.get(rank, (None, None))
                row[rank] = label
                if args.include_prob:
                    row[f"{rank}_prob"] = prob
            rows.append(row)

    # 5) DataFrame and write TSV
    # Ensure deterministic column order: id columns + ranks (+ prob columns if requested)
    cols = ["sequence_id", "representative_id"] + all_ranks
    if args.include_prob:
        cols += [f"{r}_prob" for r in all_ranks]

    df = pd.DataFrame(rows, columns=cols)
    df.to_csv(args.out, sep="\t", index=False)
    # 6) Optional taxa profile/report (rank × taxid: avg, median, std of top-1 prob)
    if args.taxa_report:
        # Build per-representative (not per-member) long table to avoid double-counting
        long_rows = []
        for i in range(N):
            ann = annotations_by_index[i]
            for rank in all_ranks:
                label, prob = ann.get(rank, (None, None))
                if label is None or prob is None:
                    continue
                long_rows.append({"rank": rank, "taxid": label, "probability": float(prob), 'model': modelname})

        if long_rows:
            ldf = pd.DataFrame(long_rows)
            prof = (
                ldf.groupby(["rank", "taxid"])["probability"]
                   .agg(avg="mean", median="median", std="std")
                   .reset_index()
            )
            prof = prof[["rank", "taxid", "avg", "median", "std"]]
            # fill prof with 0 if empty or NaN
            prof = prof.fillna(0.0)
            # add the basename of the model to the prof df
            prof['model'] = modelname
            prof.to_csv(args.taxa_report, sep="\t", index=False)
        else:
            # empty file with header
            pd.DataFrame(columns=["rank", "taxid", "avg", "median", "std", 'model']) \
              .to_csv(args.taxa_report, sep="\t", index=False)

    # 7) Optional per-accession report: aggregate the model's per-sequence confidence
    #    by the reference accession each clustered read was aligned to. One row per
    #    member sequence (so accessions are weighted by read support), keyed to the
    #    detected organism downstream via its reference accession(s).
    if args.accession_report:
        # Map representative -> finest-rank confidence (members inherit the rep's call).
        rep_conf = {}
        ranks_fine_to_coarse = list(reversed(all_ranks))
        for i in range(N):
            rep_conf[reps_in_order[i]] = finest_confidence(annotations_by_index[i], ranks_fine_to_coarse)

        acc_rows = []
        for rep, members in rep_to_members.items():
            conf = rep_conf.get(rep)
            if conf is None:
                continue
            for member in members:
                acc = parse_accession(member)
                if acc:
                    acc_rows.append({"accession": acc, "probability": conf})

        if acc_rows:
            adf = pd.DataFrame(acc_rows)
            prof = (
                adf.groupby("accession")["probability"]
                   .agg(avg="mean", median="median", std="std", n="count")
                   .reset_index()
            )
            prof = prof.fillna(0.0)
            prof["model"] = modelname
            prof = prof[["accession", "avg", "median", "std", "n", "model"]]
            prof.to_csv(args.accession_report, sep="\t", index=False)
        else:
            pd.DataFrame(columns=["accession", "avg", "median", "std", "n", "model"]) \
              .to_csv(args.accession_report, sep="\t", index=False)

if __name__ == "__main__":
    main()
