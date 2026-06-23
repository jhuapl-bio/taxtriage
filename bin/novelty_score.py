#!/usr/bin/env python3
##############################################################################################
# Copyright 2022 The Johns Hopkins University Applied Physics Laboratory LLC
# (header trimmed for the sketch -- copy the full APL header from bin/mergeConfidence.py)
##############################################################################################
"""
Per-sample NOVELTY scoring for taxtriage.

Reads things that fell through Kraken2 + reference alignment and quantifies how much of
the sample looks like it has *no near neighbor in the DB*, plus which taxa are only
assignable at genus-or-higher rank.

Two kinds of output:
  1. <sample>.novelty.summary.tsv   -- one row, the per-sample novelty signal (joins to
                                       your confidence/mqc tables by sample id)
  2. <sample>.novelty.candidates.tsv-- one row per candidate taxon assignable at genus+
                                       but not species, with read support and divergence.

The score is a transparent, weighted combination of three z-scored components so it sits
naturally next to the existing custom scoring rather than being a black box:

    novelty = w_dark * z(dark_fraction)
            + w_rank * z(highrank_only_fraction)
            + w_idnt * z(lowident_tail_mass)

  dark_fraction          fraction of input reads/bases unassigned by EVERYTHING
                         (not classified by K2, not aligned to a reference, and not
                          assigned even at protein level)
  highrank_only_fraction fraction of reads that the translated search COULD place, but
                         only at genus/family/order/class (no species). This is the
                         "we see it at the genus level" signal.
  lowident_tail_mass     fraction of best-hit identities sitting below `idnt_cut`
                         (default 50% aa) -- divergent-but-homologous content.

z() is computed against a baseline: negative/no-template controls if provided
(--baseline), otherwise the in-run distribution across samples (--baseline auto).
A sample is flagged when novelty >= --flag-threshold (default 2.0 ~ "2 sigma out").
"""

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

# ranks we treat as "high rank only" (assignable, but not to species)
# NOTE: NCBI's 2024+ taxonomy renamed "superkingdom" -> "domain"; keep BOTH so the
# top rank is counted regardless of which dump produced the LCA. Dropping "domain"
# silently discarded those hits from both candidates and highrank_only_fraction.
HIGH_RANKS = {"genus", "family", "order", "class", "phylum",
              "superkingdom", "domain", "kingdom"}
SPECIES_RANKS = {"species", "subspecies", "strain", "serotype", "no rank"}


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="Per-sample novelty scoring for taxtriage.")
    p.add_argument("-s", "--sample", required=True, help="sample id (meta.id)")
    p.add_argument("--lca", type=Path, required=True,
                   help="mmseqs/kaiju per-read LCA tsv: query, taxid, rank, name, [lineage]")
    p.add_argument("--tophit", type=Path, default=None,
                   help="best-hit aln table with pident column (mmseqs convertalis / diamond outfmt6)")
    p.add_argument("--total-reads", type=int, required=True,
                   help="total QC'd reads for the sample (denominator for fractions)")
    p.add_argument("--k2-classified", type=int, default=0,
                   help="reads classified by Kraken2 (any rank)")
    p.add_argument("--ref-aligned", type=int, default=0,
                   help="reads aligned to a reference fasta in the alignment subworkflow")
    p.add_argument("--idnt-col", type=int, default=2,
                   help="0-based column of percent identity in --tophit (default 2 = mmseqs)")
    p.add_argument("--count-col", type=int, default=None,
                   help="0-based column in --lca holding a per-row read count (count-weighted "
                        "backends like bracken). When unset, each LCA row counts as 1 "
                        "(per-query backends like mmseqs2/kaiju).")
    p.add_argument("--idnt-cut", type=float, default=50.0,
                   help="aa %% identity below which a hit counts toward the low-identity tail")
    p.add_argument("--min-reads", type=int, default=2,
                   help="absolute floor: min hits at a taxon to report it as a genus+ "
                        "candidate (default 2). This is the LOWER bound; the effective "
                        "cutoff scales up with sample depth via --min-cand-frac.")
    p.add_argument("--min-cand-frac", type=float, default=0.002,
                   help="depth-scaled cutoff: a candidate must also hold at least this "
                        "fraction of all genus+ (placeable) hits in the sample. The effective "
                        "threshold is max(--min-reads, ceil(min_cand_frac * highrank_total)), "
                        "so deep samples stay strict (~old default of 10 at ~5k placeable hits) "
                        "while shallow / long-read / ORF-level samples still surface their top "
                        "genera instead of emitting an empty table (default 0.002).")
    p.add_argument("--baseline", default="auto",
                   help="'auto' (z vs in-run samples, supply --run-summaries) or a controls "
                        "summary tsv produced by earlier runs of this script")
    p.add_argument("--run-summaries", type=Path, default=None,
                   help="optional: concatenated summary.tsv of sibling samples for 'auto' baseline")
    p.add_argument("--weights", default="0.5,0.3,0.2",
                   help="w_dark,w_rank,w_idnt")
    p.add_argument("--flag-threshold", type=float, default=2.0)
    p.add_argument("--classifier", default="",
                   help="novelty backend that produced the LCA (mmseqs2|kaiju|bracken); "
                        "passed through to the report so the Novelty tab can label its source.")
    p.add_argument("--gene-mode", action="store_true",
                   help="set when the query was Pyrodigal-predicted genes (--novelty_gene) rather "
                        "than whole contigs; the report uses this to label counts as genes/seqs "
                        "instead of contigs.")
    p.add_argument("-o", "--out-prefix", required=True)
    return p.parse_args(argv)


def read_lca(path, count_col=None):
    """Return (reads_assigned_total, by_taxon_highrank, by_taxon_species, highrank_only_reads).

    by_taxon_highrank[(taxid, rank, name)] = n — genus/family/order/class/phylum hits only.
    by_taxon_species[(taxid, rank, name)] = n — species/no-rank hits with a real taxid.

    Both dicts count contigs (kaiju/mmseqs2) or reads (bracken) per taxon.  They are kept
    separate so the novelty *score* uses only the high-rank pool (unchanged behaviour) while
    the candidates *table* can surface species-level protein hits too — those represent taxa
    that the reference alignment never saw, which is exactly the novelty signal.

    Per-query backends (mmseqs2/kaiju) emit one row per query, so each row weighs 1.
    Count-weighted backends (bracken) pass `count_col`: each row carries its own read
    count, so a single row stands in for that many reads.
    """
    assigned = 0
    highrank_only = 0
    by_taxon_highrank = defaultdict(int)
    by_taxon_species = defaultdict(int)
    with open(path) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 4:
                continue
            taxid, rank, name = cols[1], cols[2].strip().lower(), cols[3]
            if taxid in ("0", "", "no rank") and rank in ("", "no rank"):
                continue  # unassigned -> contributes to dark matter, counted elsewhere
            if count_col is not None and len(cols) > count_col:
                try:
                    w = int(float(cols[count_col]))
                except ValueError:
                    w = 1
            else:
                w = 1
            if w <= 0:
                continue
            assigned += w
            if rank in HIGH_RANKS:
                highrank_only += w
                by_taxon_highrank[(taxid, rank, name)] += w
            elif rank in SPECIES_RANKS and taxid not in ("0", ""):
                # Species-level protein hits: real taxid but not genus+ — still potentially
                # novel (protein search found them; reference alignment did not).
                by_taxon_species[(taxid, rank, name)] += w
    return assigned, by_taxon_highrank, by_taxon_species, highrank_only


def lowident_tail(path, idnt_col, cut):
    """Fraction of best hits with identity below `cut`. Returns (frac, n)."""
    if path is None or not Path(path).exists():
        return 0.0, 0
    n = below = 0
    seen = set()
    with open(path) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) <= idnt_col:
                continue
            q = cols[0]
            if q in seen:        # keep one (best) row per query
                continue
            seen.add(q)
            try:
                pid = float(cols[idnt_col])
            except ValueError:
                continue
            n += 1
            if pid < cut:
                below += 1
    return (below / n if n else 0.0), n


def zscore(x, mean, sd):
    if sd is None or sd <= 1e-9:
        # no spread to compare against -> fall back to raw magnitude so a real
        # signal in a clean run still surfaces instead of being divided away
        return x * 10.0
    return (x - mean) / sd


def baseline_stats(run_summaries, fields):
    """mean/sd per field from sibling-sample or control summaries."""
    import statistics
    cols = defaultdict(list)
    if run_summaries and Path(run_summaries).exists():
        with open(run_summaries) as fh:
            r = csv.DictReader(fh, delimiter="\t")
            for row in r:
                # exclude flagged/positive samples from the baseline if marked
                if row.get("is_control", "").lower() in ("0", "false", ""):
                    if row.get("use_as_baseline", "true").lower() in ("0", "false"):
                        continue
                for f in fields:
                    try:
                        cols[f].append(float(row[f]))
                    except (KeyError, ValueError):
                        pass
    stats = {}
    for f in fields:
        vals = cols[f]
        if len(vals) >= 2:
            stats[f] = (statistics.mean(vals), statistics.pstdev(vals))
        else:
            stats[f] = (0.0, None)
    return stats


def main(argv=None):
    a = parse_args(argv)
    w_dark, w_rank, w_idnt = (float(x) for x in a.weights.split(","))

    assigned, by_taxon, by_taxon_species, highrank_only = read_lca(a.lca, a.count_col)
    tail_frac, tail_n = lowident_tail(a.tophit, a.idnt_col, a.idnt_cut)

    total = max(a.total_reads, 1)
    # "dark matter": not classified by K2, not ref-aligned, not protein-assigned.
    # Conservative upper bound: reads explained by none of the three.
    explained = max(a.k2_classified, a.ref_aligned) + assigned
    dark_fraction = max(0.0, (total - explained) / total)
    highrank_only_fraction = highrank_only / total

    fields = ["dark_fraction", "highrank_only_fraction", "lowident_tail_mass"]
    raw = {
        "dark_fraction": dark_fraction,
        "highrank_only_fraction": highrank_only_fraction,
        "lowident_tail_mass": tail_frac,
    }
    stats = baseline_stats(a.run_summaries, fields)
    z = {f: zscore(raw[f], *stats[f]) for f in fields}

    novelty = w_dark * z["dark_fraction"] + w_rank * z["highrank_only_fraction"] + w_idnt * z["lowident_tail_mass"]
    flagged = novelty >= a.flag_threshold

    # Method-level read accounting, exported so the report can render a per-sample
    # "where did every read go" breakdown (closed-set alignment / kraken2 / mmseqs rescue
    # / dark) without re-deriving it client-side. `assigned` is everything mmseqs placed;
    # split into species-or-finer vs genus+ (the "rescued above species" bucket).
    assigned_highrank = highrank_only
    assigned_species  = max(0, assigned - highrank_only)
    ref_aligned_frac  = min(1.0, a.ref_aligned / total)
    k2_frac           = min(1.0, a.k2_classified / total)
    mmseqs_frac       = min(1.0, assigned / total)

    # ---- summary row (one line; joins to mqc/confidence by sample) ----
    summary_path = f"{a.out_prefix}.novelty.summary.tsv"
    cols = ["sample", "classifier", "gene_mode", "total_reads", "dark_fraction", "highrank_only_fraction",
            "lowident_tail_mass", "z_dark", "z_highrank", "z_lowident",
            "novelty_score", "novelty_flag",
            # extended method-accounting columns (additive; older readers ignore them)
            "ref_aligned", "k2_classified", "mmseqs_assigned",
            "mmseqs_assigned_species", "mmseqs_assigned_highrank",
            "ref_aligned_frac", "k2_frac", "mmseqs_frac"]
    with open(summary_path, "w", newline="") as fh:
        wri = csv.writer(fh, delimiter="\t")
        wri.writerow(cols)
        wri.writerow([a.sample, a.classifier, int(a.gene_mode), total,
                      f"{dark_fraction:.4f}", f"{highrank_only_fraction:.4f}", f"{tail_frac:.4f}",
                      f"{z['dark_fraction']:.3f}", f"{z['highrank_only_fraction']:.3f}", f"{z['lowident_tail_mass']:.3f}",
                      f"{novelty:.3f}", int(flagged),
                      a.ref_aligned, a.k2_classified, assigned,
                      assigned_species, assigned_highrank,
                      f"{ref_aligned_frac:.4f}", f"{k2_frac:.4f}", f"{mmseqs_frac:.4f}"])

    # ---- candidate taxa (sorted by support) ----
    # Two pools are written:
    #   1. Genus+ (high-rank-only) hits — taxa the protein search could place at genus or above
    #      but not to species.  These are the classic "novel genus" signal.
    #   2. Species-level hits — taxa the protein search placed at species (or no rank with a
    #      real taxid) that the reference alignment never saw.  For kaiju/mmseqs2 these are
    #      protein-search species hits on assembled contigs; they represent real biology that
    #      just wasn't in the reference DB.
    #
    # Depth-scaled cutoff for genus+ pool: `highrank_only` is the count of all placeable
    # (genus+) hits; scaling against it keeps deep runs strict while still surfacing the top
    # genera in shallow / long-read / ORF-level samples.
    import math
    eff_min_highrank = max(a.min_reads, math.ceil(a.min_cand_frac * highrank_only))
    # Species pool uses a simpler absolute floor; these hits are already species-resolved so
    # depth-scaling is less critical.
    eff_min_species = a.min_reads
    cand_path = f"{a.out_prefix}.novelty.candidates.tsv"
    n_written = 0
    with open(cand_path, "w", newline="") as fh:
        wri = csv.writer(fh, delimiter="\t")
        # frac_of_sample: contig-count / total-residual-reads for kaiju/mmseqs2 (note: mixed
        #   units; kept for backward compat but the report labels it as "ctgs / residual reads"
        #   when the classifier is contig-based).
        # frac_of_highrank: share of all genus+-placeable hits (genus+ pool only; empty for
        #   species entries where the concept doesn't apply).
        wri.writerow(["sample", "taxid", "rank", "name", "reads",
                      "frac_of_sample", "frac_of_highrank"])
        for (taxid, rank, name), n in sorted(by_taxon.items(), key=lambda kv: -kv[1]):
            if n >= eff_min_highrank:
                foh = (n / highrank_only) if highrank_only else 0.0
                wri.writerow([a.sample, taxid, rank, name, n,
                              f"{n/total:.4f}", f"{foh:.4f}"])
                n_written += 1
        # Species-level protein hits (not in reference alignment, resolved to species by the
        # protein classifier).  frac_of_highrank left empty — these aren't in the genus+ pool.
        for (taxid, rank, name), n in sorted(by_taxon_species.items(), key=lambda kv: -kv[1]):
            if n >= eff_min_species:
                wri.writerow([a.sample, taxid, rank, name, n,
                              f"{n/total:.4f}", ""])
                n_written += 1

    sys.stderr.write(
        f"[novelty] {a.sample}: dark={dark_fraction:.3f} "
        f"highrank_only={highrank_only_fraction:.3f} idtail={tail_frac:.3f} "
        f"score={novelty:.2f} flag={int(flagged)} "
        f"| candidates={n_written} (genus+ eff_min={eff_min_highrank} of {highrank_only} placeable; "
        f"species eff_min={eff_min_species} of {assigned - highrank_only} species hits)\n"
    )


if __name__ == "__main__":
    sys.exit(main())
