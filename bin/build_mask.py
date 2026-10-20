#!/usr/bin/env python3
# ##############################################################################################
# # Copyright 2022 The Johns Hopkins University Applied Physics Laboratory LLC
# # All rights reserved.
# # Permission is hereby granted, free of charge, to any person obtaining a copy of this
# # software and associated documentation files (the "Software"), to deal in the Software
# # without restriction, including without limitation the rights to use, copy, modify,
# # merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
# # permit persons to whom the Software is furnished to do so.
# #
# # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# # INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
# # PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# # LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# # TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
# # OR OTHER DEALINGS IN THE SOFTWARE.
# #
"""
build_mask.py -- build a TaxTriage region mask (BED) for TASS scoring.

Sources (any combination; everything is merged into one BED):

  NCBI FCS (contamination) -- regions to EXCLUDE
    --fcs_details   NCBI aggregated details, e.g. fcs_details_genbank.txt[.gz],
                    fcs_details_refseq.txt[.gz] from
                    https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/
                    columns: seq_id start_pos end_pos seq_len action
                             contam_type coverage contam_details assembly_accession
    --fcs_gx_report per-assembly FCS-GX reports (<asm>_fcs_report.txt from the
                    assembly's FTP folder, or *.fcs_gx_report.txt from running
                    FCS-GX yourself). Columns: seq_id start_pos end_pos seq_len
                    action div agg_cont_cov top_tax_name
    --fcs_adaptor_report  FCS-adaptor fcs_adaptor_report.txt
                    columns: accession length action range name
                    (range like "1..45,1200..1260")
    --fcs_summary   NCBI fcs_summary_{genbank,refseq}.txt[.gz]: assembly-level
                    contamination stats, written to --assembly_out for the
                    assemblies in --mapping (informational / reporting).

  Annotation -- regions to TARGET (use with --mask_mode target)
    --feature_table NCBI *_feature_table.txt[.gz]
    --gff           GFF3[.gz]
    --feature_select  case-insensitive regex matched against gene/locus_tag/
                    product/name/Name/ID attributes, e.g. "^(stx[12][ab]?|eae)$"
    --feature_types   feature types to consider (default CDS,gene)

  Anything else
    --bed           your own BED files, merged in as-is

Restricting to the references you actually aligned to (strongly recommended
for the 4+ GB FCS details file -- it is streamed once and only matching rows
are kept):
    --refs          FASTA(.gz) / .fai / BAM / plain list of accessions
    --mapping       TaxTriage accession->assembly map (col1 accession,
                    col2 GCF_/GCA_ assembly). Enables assembly-level matching
                    (GCF_x <-> GCA_x) and the --assembly_out table.

Accession matching is version- and RefSeq-prefix tolerant: an FCS row for
CP012345.1 masks NZ_CP012345.1 in the BAM and vice-versa. Output BED column 1
always uses the reference name as it appears in --refs when --refs is given.

FCS actions (https://github.com/ncbi/fcs/wiki/FCS-GX-taxonomy-and-actions):
  EXCLUDE  whole sequence is contaminant          -> masked by default
  TRIM     contaminant at a sequence end         -> masked by default
  FIX      contaminant inside a sequence         -> masked by default
  REVIEW   possible contaminant, low confidence  -> add with --actions
  INFO     informational (e.g. organelle)        -> add with --actions
  ACTION_TRIM / ACTION_EXCLUDE (adaptor reports) -> treated as TRIM/EXCLUDE

Outputs
  --out           merged mask BED (chrom, start, end, label, score, strand)
  --summary_out   per-reference TSV (length, masked bp, %, intervals, sources)
  --assembly_out  per-assembly FCS summary TSV (needs --fcs_summary + --mapping)
  --subset_out    the matching raw FCS details rows (a small cache you can pass
                  back as --fcs_details next time instead of the full file)

Examples
  # FCS exclude-mask for the references of a TaxTriage run
  build_mask.py --fcs_details fcs_details_genbank.txt fcs_details_refseq.txt.gz \\
      --fcs_summary fcs_summary_refseq.txt --mapping sample.combined.gcfmap.tsv \\
      --refs sample.dwnld.references.fasta --out fcs_mask.bed

  # Targets: Shiga toxin + intimin genes from NCBI feature tables
  build_mask.py --feature_table GCF_000008865.2_ASM886v2_feature_table.txt.gz \\
      --feature_select '^(stx1a|stx1b|stx2a|stx2b|eae)$' --out stec_targets.bed
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mask_regions import accession_aliases, merge_intervals, open_text  # noqa: E402

DEFAULT_ACTIONS = "EXCLUDE,TRIM,FIX"


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Build a TaxTriage region-mask BED from NCBI FCS reports, annotation, or BEDs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split("Examples", 1)[-1] if "Examples" in __doc__ else None,
    )
    src = p.add_argument_group("FCS sources (exclude)")
    src.add_argument("--fcs_details", nargs="+", default=[], help="NCBI aggregated fcs_details_*.txt[.gz]")
    src.add_argument("--fcs_gx_report", nargs="+", default=[], help="FCS-GX per-assembly report(s)")
    src.add_argument("--fcs_adaptor_report", nargs="+", default=[], help="FCS-adaptor report(s)")
    src.add_argument("--fcs_summary", nargs="+", default=[], help="NCBI fcs_summary_*.txt[.gz] (assembly level)")
    src.add_argument("--actions", default=DEFAULT_ACTIONS,
                     help=f"Comma-separated FCS actions to mask (default {DEFAULT_ACTIONS}); 'all' for every action")
    src.add_argument("--divisions", default=None,
                     help="Optional regex on the FCS contaminant division (e.g. 'prok|anml:primates|synt'); "
                          "default: every division")
    src.add_argument("--min_contam_len", type=int, default=0, help="Skip FCS spans shorter than this (bp)")

    ann = p.add_argument_group("Annotation sources (targets)")
    ann.add_argument("--feature_table", nargs="+", default=[], help="NCBI *_feature_table.txt[.gz]")
    ann.add_argument("--gff", nargs="+", default=[], help="GFF3[.gz]")
    ann.add_argument("--feature_select", default=None,
                     help="Case-insensitive regex on gene / locus_tag / product / Name / ID")
    ann.add_argument("--feature_types", default="CDS,gene", help="Feature types to keep (default CDS,gene)")

    other = p.add_argument_group("Other")
    other.add_argument("--bed", nargs="+", default=[], help="Extra BED files merged in as-is")
    other.add_argument("--masked_fasta", nargs="+", default=[],
                       help="Masked FASTA(.gz) of the references: masked bases become mask regions. Use the "
                            "same sequence IDs as the BAM (NZ_/version differences are tolerated).")
    other.add_argument("--masked_fasta_case", choices=["hard", "soft", "both"], default="hard",
                       help="Which masking to read from --masked_fasta: hard = runs of N/n (e.g. FCS-GX "
                            "'clean' output, dustmasker -outfmt fasta with hard masking), soft = lowercase "
                            "bases (RepeatMasker / dustmasker soft masking), both. Default hard.")
    other.add_argument("--masked_fasta_min_run", type=int, default=1,
                       help="Ignore masked runs shorter than this many bp (e.g. 10 to skip isolated "
                            "ambiguity Ns). Default 1.")
    other.add_argument("--refs", nargs="+", default=[],
                       help="Restrict to these references: FASTA(.gz), .fai, BAM/CRAM, or text list of accessions")
    other.add_argument("--mapping", nargs="+", default=[],
                       help="accession<TAB>assembly map(s), e.g. TaxTriage *.gcfmap.tsv")
    other.add_argument("--pad", type=int, default=0, help="Pad every interval by N bp each side")

    out = p.add_argument_group("Outputs")
    out.add_argument("-o", "--out", required=True, help="Output BED")
    out.add_argument("--summary_out", default=None, help="Per-reference summary TSV")
    out.add_argument("--assembly_out", default=None, help="Per-assembly FCS summary TSV")
    out.add_argument("--subset_out", default=None, help="Write matching FCS details rows here (cache)")
    out.add_argument("--json_out", default=None, help="Machine-readable run summary JSON")
    return p.parse_args(argv)


# ─────────────────────────────────────────────────────────────────────────────
# Reference universe
# ─────────────────────────────────────────────────────────────────────────────
def load_refs(paths):
    """Return {accession: length or None}."""
    refs = {}
    for path in paths:
        low = path.lower()
        if low.endswith((".bam", ".cram", ".sam")):
            try:
                import pysam
            except ImportError:
                sys.exit("pysam is required to read references from a BAM/CRAM")
            with pysam.AlignmentFile(path, "rb" if not low.endswith(".sam") else "r") as fh:
                for name, ln in zip(fh.references, fh.lengths):
                    refs[name] = int(ln)
            continue
        if low.endswith(".fai"):
            with open(path) as fh:
                for line in fh:
                    f = line.split("\t")
                    if len(f) >= 2:
                        refs[f[0]] = int(f[1])
            continue
        with open_text(path) as fh:
            first = fh.readline()
        if first.startswith(">"):
            name, ln = None, 0
            with open_text(path) as fh:
                for line in fh:
                    if line.startswith(">"):
                        if name is not None:
                            refs[name] = ln
                        name = line[1:].split()[0] if len(line) > 1 else None
                        ln = 0
                    else:
                        ln += len(line.strip())
            if name is not None:
                refs[name] = ln
        else:
            with open_text(path) as fh:
                for line in fh:
                    t = line.strip().split()
                    if t and not t[0].startswith("#"):
                        refs.setdefault(t[0], None)
    return refs


def load_mapping(paths):
    """accession -> assembly accession (GCF_/GCA_)."""
    m = {}
    for path in paths:
        with open_text(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 2 or line.startswith("#"):
                    continue
                if re.match(r"^GC[AF]_\d+(\.\d+)?$", f[1]):
                    m[f[0]] = f[1]
    return m


def asm_aliases(asm):
    """GCF_000001405.40 <-> GCA_000001405.40 (+ unversioned)."""
    asm = asm.strip()
    out = [asm]
    if asm.startswith("GCF_"):
        out.append("GCA_" + asm[4:])
    elif asm.startswith("GCA_"):
        out.append("GCF_" + asm[4:])
    out += [a.rsplit(".", 1)[0] for a in list(out) if "." in a]
    return out


class Universe:
    """Which reference names we care about, with alias resolution."""

    def __init__(self, refs: dict):
        self.refs = refs
        self.alias = {}
        for r in refs:
            for a in accession_aliases(r):
                self.alias.setdefault(a, r)

    @property
    def open(self):
        return not self.refs

    def resolve(self, seq_id):
        if self.open:
            return seq_id
        for a in accession_aliases(seq_id):
            r = self.alias.get(a)
            if r is not None:
                return r
        return None


# ─────────────────────────────────────────────────────────────────────────────
# Parsers -> rows of (ref, start0, end, label, source)
# ─────────────────────────────────────────────────────────────────────────────
def _action_ok(action, actions):
    if actions is None:
        return True
    a = action.upper().replace("ACTION_", "")
    return a in actions


def parse_fcs_details(path, uni, actions, div_re, min_len, subset_fh, stats):
    """Stream an NCBI aggregated details file (or a GX report with the same lead columns)."""
    n = kept = 0
    with open_text(path) as fh:
        header = None
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                if subset_fh is not None and not stats.get("_subset_header"):
                    subset_fh.write(line)
                    stats["_subset_header"] = True
                continue
            n += 1
            sid, rest = line.split("\t", 1)
            ref = uni.resolve(sid)
            if ref is None:
                continue
            f = line.rstrip("\n").split("\t")
            try:
                s1, e1 = int(f[1]), int(f[2])
            except (ValueError, IndexError):
                continue
            action = f[4] if len(f) > 4 else "EXCLUDE"
            div = f[5] if len(f) > 5 else ""
            if not _action_ok(action, actions):
                stats["skipped_action"][action] += 1
                continue
            if div_re is not None and not div_re.search(div):
                stats["skipped_division"] += 1
                continue
            if e1 - s1 + 1 < min_len:
                continue
            if subset_fh is not None:
                subset_fh.write(line)
            asm = f[8] if len(f) > 8 else ""
            detail = f[7] if len(f) > 7 else ""
            label = f"FCS:{action}:{div}" + (f":{detail}" if detail else "")
            label = label.replace("\t", " ").replace(" ", "_")[:120]
            stats["by_action"][action] += 1
            if asm:
                stats["assemblies_hit"].add(asm)
            kept += 1
            yield ref, s1 - 1, e1, label, os.path.basename(path)
    stats["rows_scanned"] += n
    print(f"[build_mask] {os.path.basename(path)}: scanned {n:,} rows, kept {kept:,}", file=sys.stderr)


def parse_fcs_adaptor(path, uni, actions, stats):
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            ref = uni.resolve(f[0])
            if ref is None:
                continue
            action = f[2]
            if not _action_ok(action, actions):
                stats["skipped_action"][action] += 1
                continue
            name = f[4] if len(f) > 4 else "adaptor"
            try:
                length = int(f[1])
            except ValueError:
                length = None
            rng = f[3].strip()
            spans = []
            if rng and rng != "NA":
                for part in rng.split(","):
                    if ".." in part:
                        a, b = part.split("..", 1)
                        spans.append((int(a) - 1, int(b)))
            elif "EXCLUDE" in action.upper() and length:
                spans.append((0, length))
            for a, b in spans:
                stats["by_action"]["ADAPTOR:" + action] += 1
                yield ref, a, b, f"FCS-adaptor:{action}:{name}"[:120].replace(" ", "_"), os.path.basename(path)


def _attrs_gff(s):
    d = {}
    for kv in s.strip().split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k.strip()] = v.strip()
    return d


def parse_feature_table(path, uni, sel_re, ftypes, stats):
    with open_text(path) as fh:
        header = None
        for line in fh:
            if line.startswith("#"):
                header = line.lstrip("# ").rstrip("\n").split("\t")
                continue
            if header is None:
                continue
            f = dict(zip(header, line.rstrip("\n").split("\t")))
            if ftypes and f.get("feature") not in ftypes:
                continue
            acc = f.get("genomic_accession") or ""
            ref = uni.resolve(acc)
            if ref is None:
                continue
            fields = " ".join(f.get(k, "") for k in ("symbol", "locus_tag", "name", "product_accession", "attributes"))
            names = [f.get("symbol", ""), f.get("locus_tag", ""), f.get("name", ""), f.get("product_accession", "")]
            if sel_re is not None and not any(sel_re.search(x) for x in names if x) and not sel_re.search(fields):
                continue
            try:
                a, b = int(f["start"]) - 1, int(f["end"])
            except (KeyError, ValueError):
                continue
            label = f.get("symbol") or f.get("locus_tag") or f.get("name") or "feature"
            stats["features"] += 1
            yield ref, a, b, f"target:{f.get('feature')}:{label}"[:120].replace(" ", "_"), os.path.basename(path)


def parse_gff(path, uni, sel_re, ftypes, stats):
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            if ftypes and f[2] not in ftypes:
                continue
            ref = uni.resolve(f[0])
            if ref is None:
                continue
            at = _attrs_gff(f[8])
            names = [at.get(k, "") for k in ("gene", "Name", "locus_tag", "ID", "product", "gene_synonym")]
            if sel_re is not None and not any(sel_re.search(x) for x in names if x):
                continue
            label = at.get("gene") or at.get("Name") or at.get("locus_tag") or at.get("ID") or "feature"
            stats["features"] += 1
            yield ref, int(f[3]) - 1, int(f[4]), f"target:{f[2]}:{label}"[:120].replace(" ", "_"), os.path.basename(path)


def parse_masked_fasta(path, uni, case, min_run, stats):
    """Turn masked bases of a FASTA into intervals.

    hard: N/n runs. soft: lowercase runs (a lowercase 'n' counts for both).
    Streams the file line by line, so multi-GB FASTAs are fine; sequences not
    in --refs are skipped without being scanned in detail.
    """
    want_hard = case in ("hard", "both")
    want_soft = case in ("soft", "both")
    hard_re = re.compile(r"[Nn]+")
    soft_re = re.compile(r"[a-z]+")
    label_h = f"masked_fasta:hard:{os.path.basename(path)}"
    label_s = f"masked_fasta:soft:{os.path.basename(path)}"
    ref, pos = None, 0
    open_h = open_s = None          # run start carried across line breaks

    def flush(kind, start, end):
        if start is None or ref is None or end - start < min_run:
            return None
        stats["masked_fasta_runs"] += 1
        return (ref, start, end, label_h if kind == "h" else label_s, os.path.basename(path))

    with open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                for kind, st in (("h", open_h), ("s", open_s)):
                    r = flush(kind, st, pos)
                    if r:
                        yield r
                ref = uni.resolve(line[1:].split()[0]) if len(line) > 1 else None
                stats["masked_fasta_seqs"] += 1 if ref is not None else 0
                pos, open_h, open_s = 0, None, None
                continue
            seq = line.strip()
            if ref is None or not seq:
                pos += len(seq)
                continue
            L = len(seq)
            for kind, want, rx in (("h", want_hard, hard_re), ("s", want_soft, soft_re)):
                if not want:
                    continue
                cur = open_h if kind == "h" else open_s
                for m in rx.finditer(seq):
                    a, b = m.start(), m.end()
                    if cur is not None and a == 0:
                        pass                      # run continues from previous line
                    else:
                        r = flush(kind, cur, pos) if cur is not None else None
                        if r:
                            yield r
                        cur = pos + a
                    if b < L:                     # run ends inside this line
                        r = flush(kind, cur, pos + b)
                        if r:
                            yield r
                        cur = None
                if cur is not None and not (seq[-1:] and rx.match(seq[-1])):
                    r = flush(kind, cur, pos)
                    if r:
                        yield r
                    cur = None
                if kind == "h":
                    open_h = cur
                else:
                    open_s = cur
            pos += L
    for kind, st in (("h", open_h), ("s", open_s)):
        r = flush(kind, st, pos)
        if r:
            yield r


def parse_bed(path, uni, stats):
    with open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                f = line.split()
            if len(f) < 3:
                continue
            ref = uni.resolve(f[0])
            if ref is None:
                stats["bed_unmatched"] += 1
                continue
            yield ref, int(f[1]), int(f[2]), (f[3] if len(f) > 3 else "bed"), os.path.basename(path)


def load_fcs_summary(paths, wanted_asms):
    """Return {assembly: row-dict} for assemblies of interest (or all if none)."""
    out = {}
    want = set()
    for a in wanted_asms:
        want.update(asm_aliases(a))
    for path in paths:
        with open_text(path) as fh:
            header = None
            for line in fh:
                if line.startswith("#"):
                    header = line.lstrip("#").rstrip("\n").split("\t")
                    continue
                aid = line.split("\t", 1)[0]
                if want and aid not in want and aid.rsplit(".", 1)[0] not in want:
                    continue
                row = dict(zip(header or [], line.rstrip("\n").split("\t")))
                row["_source"] = os.path.basename(path)
                out[aid] = row
    return out


# ─────────────────────────────────────────────────────────────────────────────
def main(argv=None):
    args = parse_args(argv)
    actions = None if args.actions.strip().lower() == "all" else {
        a.strip().upper().replace("ACTION_", "") for a in args.actions.split(",") if a.strip()}
    div_re = re.compile(args.divisions) if args.divisions else None
    sel_re = re.compile(args.feature_select, re.I) if args.feature_select else None
    ftypes = {t.strip() for t in args.feature_types.split(",") if t.strip()}

    refs = load_refs(args.refs)
    uni = Universe(refs)
    mapping = load_mapping(args.mapping)
    if refs:
        print(f"[build_mask] restricting to {len(refs):,} reference sequences", file=sys.stderr)
    elif args.fcs_details:
        print("[build_mask] WARNING: no --refs given; EVERY FCS row will be written "
              "(the GenBank details file alone yields millions of intervals)", file=sys.stderr)

    stats = {
        "rows_scanned": 0, "by_action": Counter(), "skipped_action": Counter(),
        "skipped_division": 0, "assemblies_hit": set(), "features": 0, "bed_unmatched": 0,
        "masked_fasta_runs": 0, "masked_fasta_seqs": 0,
    }
    subset_fh = open(args.subset_out, "w") if args.subset_out else None

    raw = defaultdict(list)       # ref -> [(s, e)]
    labels = defaultdict(list)    # ref -> [(s, e, label)]
    sources = defaultdict(Counter)

    def add(gen):
        for ref, s, e, label, src in gen:
            if e <= s:
                continue
            s = max(0, s - args.pad)
            e = e + args.pad
            ln = refs.get(ref)
            if ln:
                e = min(e, ln)
            raw[ref].append((s, e))
            labels[ref].append((s, e, label))
            sources[ref][src] += 1

    for p in args.fcs_details + args.fcs_gx_report:
        add(parse_fcs_details(p, uni, actions, div_re, args.min_contam_len, subset_fh, stats))
    for p in args.fcs_adaptor_report:
        add(parse_fcs_adaptor(p, uni, actions, stats))
    for p in args.feature_table:
        add(parse_feature_table(p, uni, sel_re, ftypes, stats))
    for p in args.gff:
        add(parse_gff(p, uni, sel_re, ftypes, stats))
    for p in args.masked_fasta:
        add(parse_masked_fasta(p, uni, args.masked_fasta_case, args.masked_fasta_min_run, stats))
    for p in args.bed:
        add(parse_bed(p, uni, stats))
    if subset_fh is not None:
        subset_fh.close()

    # Write BED: one line per ORIGINAL interval (labels preserved) -- the
    # scorer merges overlaps itself. Sorted for tabix/IGV friendliness.
    n_lines = 0
    with open(args.out, "w") as out:
        _srcs = [k for k, v in (("fcs_details", args.fcs_details), ("fcs_gx_report", args.fcs_gx_report),
                                 ("fcs_adaptor_report", args.fcs_adaptor_report), ("feature_table", args.feature_table),
                                 ("gff", args.gff), ("bed", args.bed),
                                 ("masked_fasta", args.masked_fasta)) if v]
        _fcs_note = f"; fcs_actions={args.actions}" if (args.fcs_details or args.fcs_gx_report or args.fcs_adaptor_report) else ""
        _sel_note = f"; feature_select={args.feature_select}" if args.feature_select else ""
        out.write(f"#TaxTriage region mask -- build_mask.py; sources={','.join(_srcs) or 'none'}{_fcs_note}{_sel_note}\n")
        for ref in sorted(labels):
            for s, e, lab in sorted(labels[ref]):
                out.write(f"{ref}\t{s}\t{e}\t{lab}\t0\t.\n")
                n_lines += 1

    per_ref = []
    for ref in sorted(raw):
        merged = merge_intervals(raw[ref])
        bp = sum(b - a for a, b in merged)
        ln = refs.get(ref)
        per_ref.append({
            "reference": ref,
            "assembly": mapping.get(ref, ""),
            "length": ln if ln is not None else "",
            "region_bp": bp,
            "region_pct": round(100.0 * bp / ln, 4) if ln else "",
            "intervals": len(merged),
            "labels": ";".join(sorted({l.split(":")[0] + ":" + l.split(":")[1] if l.count(":") >= 1 else l
                                       for _, _, l in labels[ref]})),
            "sources": ";".join(f"{k}({v})" for k, v in sources[ref].items()),
        })
    if args.summary_out:
        with open(args.summary_out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["reference", "assembly", "length", "region_bp", "region_pct",
                                               "intervals", "labels", "sources"], delimiter="\t")
            w.writeheader()
            for r in per_ref:
                w.writerow(r)

    asm_rows = {}
    if args.fcs_summary:
        wanted = set(mapping.values()) if mapping else set()
        asm_rows = load_fcs_summary(args.fcs_summary, wanted)
        if args.assembly_out:
            # one row per assembly in the run (RefSeq row preferred over GenBank)
            asms = sorted(wanted) if wanted else sorted(asm_rows)
            cols = ["assembly", "fcs_source", "fraction_contaminated", "total_contaminant_length",
                    "total_contaminant_accessions", "exclude_ct", "trim_ct", "fix_ct", "review_ct",
                    "asserted_division", "top_inferred_division", "fcs_gx_db_rev", "masked_refs", "masked_bp"]
            by_asm_mask = defaultdict(lambda: [0, 0])
            for r in per_ref:
                if r["assembly"]:
                    by_asm_mask[r["assembly"]][0] += 1
                    by_asm_mask[r["assembly"]][1] += r["region_bp"]
            with open(args.assembly_out, "w") as fh:
                fh.write("\t".join(cols) + "\n")
                for a in asms:
                    row = None
                    for al in asm_aliases(a):
                        if al in asm_rows:
                            row = asm_rows[al]
                            break
                    row = row or {}
                    vals = [a, row.get("_source", "not_screened"), row.get("fraction_contaminated", ""),
                            row.get("total_contaminant_length", ""), row.get("total_contaminant_accessions", ""),
                            row.get("exclude_ct", ""), row.get("trim_ct", ""), row.get("fix_ct", ""),
                            row.get("review_ct", ""), row.get("asserted_division", ""),
                            row.get("top_inferred_division", ""), row.get("fcs_gx_db_rev", ""),
                            str(by_asm_mask[a][0]), str(by_asm_mask[a][1])]
                    fh.write("\t".join(str(v) for v in vals) + "\n")

    tot_bp = sum(r["region_bp"] for r in per_ref)
    summary = {
        "bed": os.path.basename(args.out),
        "references_considered": len(refs) if refs else None,
        "references_with_intervals": len(per_ref),
        "bed_lines": n_lines,
        "region_bp": tot_bp,
        "fcs_rows_scanned": stats["rows_scanned"],
        "fcs_rows_by_action": dict(stats["by_action"]),
        "fcs_rows_skipped_by_action": dict(stats["skipped_action"]),
        "fcs_rows_skipped_by_division": stats["skipped_division"],
        "fcs_assemblies_hit": sorted(stats["assemblies_hit"]),
        "annotation_features": stats["features"],
        "bed_rows_unmatched": stats["bed_unmatched"],
        "masked_fasta_sequences_matched": stats["masked_fasta_seqs"],
        "masked_fasta_runs": stats["masked_fasta_runs"],
        "actions": args.actions,
        "divisions": args.divisions,
        "assemblies_with_fcs_summary": len(asm_rows),
        "assemblies_contaminated": len({
            k.split("_", 1)[-1] for k, r in asm_rows.items()
            if float((r.get("total_contaminant_length") or 0) or 0) > 0}),
    }
    if args.json_out:
        with open(args.json_out, "w") as fh:
            json.dump(summary, fh, indent=2)
    print(f"[build_mask] {len(per_ref):,} references with mask intervals, {tot_bp:,} bp, "
          f"{n_lines:,} BED lines -> {args.out}", file=sys.stderr)
    if refs and not per_ref:
        print("[build_mask] NOTE: none of the references matched any source row. For FCS this is "
              "normal for RefSeq (already screened/cleaned) references -- see the --assembly_out table.",
              file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
