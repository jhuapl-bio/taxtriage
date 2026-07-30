#!/usr/bin/env python3
"""
Rebuild per-sample ``<sample>.paths.json`` inputs from a generated report.

bin/make_report.py flattens each ``*.paths.json`` (metadata + nested
organisms) into the report's ``window.HEATMAP_BOOT`` payload. This script walks
that back: it extracts the bootstrap, then re-nests the flat ``records`` into
the genus → species → strain hierarchy and re-emits files in the same shape and
style as the originals (``{"metadata": {...}, "organisms": [...]}``, 4-space
indent, original key order).

Fidelity
--------
* ``metadata``  — lossless: the report stores each sample's metadata verbatim.
* ``organisms`` — partial: flattening keeps ~45 of the ~85 organism fields and
  rounds several of them (Coverage → 1 dp, Mean Depth → 2 dp, Gini → 3 dp, …).
  Unrecoverable fields (hmp_*, minhash gates, control_comparison, accessions,
  rpm_norm_scale, …) are omitted unless --include-nulls is given.
* ``contigs`` / ``depth_histogram`` / ``breadth_histogram`` are restored from
  the payload's contig_data when present.
* ``protein_annotations_genus`` is rebuilt on genus nodes from the report's
  protein tables (--no-protein to skip).

Run with --report to see exactly which fields were recovered vs. dropped.

Usage
-----
    python scripts/reconstruct_paths_json.py all.odr.html -d rebuilt/
    python scripts/reconstruct_paths_json.py all.odr.html -d rebuilt/ --sample Miseq_Run_A
    python scripts/reconstruct_paths_json.py all.odr.html -d rebuilt/ --include-nulls --report
    python scripts/reconstruct_paths_json.py all.odr.html --combined -o all.samples.json
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from collections import defaultdict
from pathlib import Path

# reuse the extractor that lives next to this script
_EXTRACTOR = Path(__file__).resolve().parent / "extract_heatmap_boot.py"
_spec = importlib.util.spec_from_file_location("extract_heatmap_boot", _EXTRACTOR)
_ehb = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_ehb)

# ── field maps ─────────────────────────────────────────────────────────────

# Canonical organism key order, taken from a real ALIGNMENT_PER_SAMPLE
# paths.json, so rebuilt files read like the originals.
CANONICAL_ORDER = [
    "toplevelkey", "name", "subkey", "subkeyname", "toplevelname", "numreads",
    "length", "covered_bases", "k2_reads", "numreads_strain_sum", "accessions",
    "coverage", "covered_bases_union", "covered_bases_strain_sum",
    "coverage_pooled", "rpm_confidence_weight", "read_fraction", "rpkm",
    "meandepth", "minhash_reduction", "diamond_identity", "minhash_score",
    "hmp_percentile", "plasmid_score", "mapq_score", "meanmapq",
    "k2_disparity_score", "meanbaseq", "rpm", "highmapq_fraction",
    "abundance_confidence", "has_plasmid", "highmapq_fraction_effective",
    "breadth_log_score", "gini_coefficient", "gini_contig_frac",
    "gini_contig_penalty", "gini_mapq_scale", "gini_depth_conc_penalty",
    "n_contigs_total", "n_contigs_covered", "minhash_confidence",
    "minhash_reduction_pre_gate", "minhash_coverage_gate",
    "minhash_breadth_gate", "minhash_gate_floor", "minhash_reduction_raw",
    "rpm_share_in_group", "rpm_norm_scale", "solo_exclusivity_boost",
    "score_power_scale", "members", "key", "sample_name", "mmbert",
    "mmbert_model", "zscore", "hmp_norm_abundance", "hmp_mean", "hmp_std",
    "observed_abundance", "hmp_site_count", "hmp_num_samples", "disparity",
    "tass_score", "high_cons", "mol_type", "is_pathogen", "microbial_category",
    "annClass", "is_annotated", "status", "ref", "matched_taxid",
    "matched_rank", "normalized_sample_site", "commensal_sites",
    "pathogenic_sites", "sampletype", "total_reads", "control_comparison",
    "protein_annotations_genus", "taxonomy", "high_ani_matches", "contigs",
    "depth_histogram", "breadth_histogram",
]

# report column -> (organism field, transform)
_ident = lambda v: v
COLUMN_MAP = {
    "Detected Organism":     ("name", _ident),
    "Taxonomic ID #":        ("key", str),
    "Subkey":                ("subkey", str),
    "TASS Score":            ("tass_score", float),
    "# Reads Aligned":       ("numreads", float),
    "Covered Bases":         ("covered_bases", float),
    "Genome Length (bp)":    ("length", float),
    "Coverage":              ("coverage", lambda v: float(v) / 100.0),
    "Mean Depth":            ("meandepth", float),
    "Gini Coefficient":      ("gini_coefficient", float),
    "Mean MapQ":             ("meanmapq", float),
    "Mean BaseQ":            ("meanbaseq", float),
    "Minhash Score":         ("minhash_reduction", float),
    "Breadth Score":         ("breadth_log_score", float),
    "MapQ Score":            ("mapq_score", float),
    "Disparity Score":       ("disparity", float),
    "Diamond Identity":      ("diamond_identity", float),
    "K2 Reads":              ("k2_reads", float),
    "RPM":                   ("rpm", float),
    "RPKM":                  ("rpkm", float),
    "MicrobeRT Probability": ("mmbert", lambda v: None if v in (None, "") else float(v) / 100.0),
    "MicrobeRT Model":       ("mmbert_model", lambda v: v or None),
    "Microbial Category":    ("microbial_category", _ident),
    "Ann Class":             ("annClass", _ident),
    "IsAnnotated":           ("is_annotated", _ident),
    "High Consequence":      ("high_cons", bool),
    "Mol Type":              ("mol_type", _ident),
    "Status":                ("status", _ident),
    "Passes Threshold":      ("passes_threshold", bool),
    "High ANI Matches":      ("high_ani_matches", _ident),
}

_TAX_COLS = {
    "Domain": "domain", "Kingdom": "kingdom", "Superkingdom": "superkingdom",
    "Phylum": "phylum", "Class": "class", "Order": "order",
    "Family": "family", "Genus": "genus",
}

# Organism fields the flattening step throws away entirely.
UNRECOVERABLE = [
    "numreads_strain_sum", "accessions", "covered_bases_union",
    "covered_bases_strain_sum", "coverage_pooled", "rpm_confidence_weight",
    "read_fraction", "minhash_score", "hmp_percentile", "plasmid_score",
    "k2_disparity_score", "highmapq_fraction", "abundance_confidence",
    "has_plasmid", "highmapq_fraction_effective", "gini_contig_frac",
    "gini_contig_penalty", "gini_mapq_scale", "gini_depth_conc_penalty",
    "n_contigs_total", "n_contigs_covered", "minhash_confidence",
    "minhash_reduction_pre_gate", "minhash_coverage_gate",
    "minhash_breadth_gate", "minhash_gate_floor", "minhash_reduction_raw",
    "rpm_share_in_group", "rpm_norm_scale", "solo_exclusivity_boost",
    "score_power_scale", "zscore", "hmp_norm_abundance", "hmp_mean", "hmp_std",
    "observed_abundance", "hmp_site_count", "hmp_num_samples", "is_pathogen",
    "ref", "matched_taxid", "matched_rank", "normalized_sample_site",
    "commensal_sites", "pathogenic_sites", "control_comparison",
]

# Rounded during flattening — recovered value is close, not exact.
LOSSY = {
    "coverage": "1 dp on the 0-100 % scale",
    "meandepth": "2 dp", "gini_coefficient": "3 dp", "meanmapq": "1 dp",
    "meanbaseq": "1 dp", "minhash_reduction": "3 dp", "breadth_log_score": "3 dp",
    "mapq_score": "3 dp", "disparity": "3 dp", "diamond_identity": "1 dp",
    "rpm": "2 dp", "rpkm": "4 dp", "numreads": "truncated to int",
    "covered_bases": "truncated to int", "length": "truncated to int",
    "k2_reads": "truncated to int", "mmbert": "2 dp on the 0-100 % scale",
}


def _order(d: dict) -> dict:
    """Reorder an organism dict into the canonical paths.json field order."""
    rank = {k: i for i, k in enumerate(CANONICAL_ORDER)}
    return {k: d[k] for k in sorted(d, key=lambda k: (rank.get(k, 10_000), k))}


# ── node construction ──────────────────────────────────────────────────────

def _node_from_row(row: dict, meta: dict, *, include_nulls: bool) -> dict:
    node: dict = {}
    for col, (field, fn) in COLUMN_MAP.items():
        if col not in row:
            continue
        val = row[col]
        if val is None and field not in ("mmbert", "mmbert_model"):
            continue
        node[field] = fn(val)

    # High ANI Matches is only meaningful when the run was ANI-annotated.
    if not row.get("ANI Annotated") and not node.get("high_ani_matches"):
        node.pop("high_ani_matches", None)

    tax = {std: row[col] for col, std in _TAX_COLS.items() if row.get(col)}
    if tax:
        node["taxonomy"] = tax

    node["sample_name"] = row.get("Specimen ID")
    # On a genus row make_report sets "Species Name" to the genus group itself,
    # so it carries no real subkey name — leave it out rather than write a wrong
    # value (the originals name a representative species there).
    if row.get("Level") != "Genus":
        node["subkeyname"] = row.get("Species Name") or node.get("name")
    node["toplevelname"] = row.get("Genus Name") or tax.get("genus", "")

    if include_nulls:
        for field in UNRECOVERABLE:
            node.setdefault(field, None)
    return node


def _decorate_genus(node: dict, row: dict, meta: dict) -> None:
    """Fields the originals carry only on genus (toplevelkey) nodes."""
    node["sampletype"] = row.get("Sample Type") or meta.get("sample_type")
    if meta.get("total_reads") is not None:
        node["total_reads"] = meta["total_reads"]


def _protein_index(payload: dict) -> dict:
    """(sample, genus) -> [protein_annotations_genus entries] from the report tables."""
    idx: dict[tuple[str, str], list] = defaultdict(list)
    prot = payload.get("prot_data") or {}

    def _num(v):
        try:
            return float(v)
        except (TypeError, ValueError):
            return None

    # The report's protein tables list one row per organism the hit was seen on,
    # so the same annotation repeats across strains of a genus. The originals
    # store it once per genus — dedupe on the identity of the hit.
    seen: set = set()

    for table, class_col in (("per_gene_hits", "Description"), ("amr_genes", "Classification")):
        for r in prot.get(table) or []:
            sample = r.get("Specimen ID") or r.get("Sample")
            genus = r.get("Genus")
            if not sample or not genus:
                continue
            sig = (sample, genus, r.get("Source ID"), r.get("Gene"), r.get("Species"),
                   _num(r.get("%id")), _num(r.get("E-value")))
            if sig in seen:
                continue
            seen.add(sig)
            idx[(sample, genus)].append({
                "sseqid":            r.get("Source ID"),
                "pident":            _num(r.get("%id")),
                "evalue":            _num(r.get("E-value")),
                "bitscore":          _num(r.get("Bitscore")),
                "source_id":         r.get("Source ID"),
                "gene_name":         r.get("Gene"),
                "product":           r.get("Product"),
                "classification":    r.get(class_col),
                "antibiotics_class": r.get("Antibiotics Class"),
                "antibiotics":       r.get("Antibiotics"),
                "organism":          r.get("Organism"),
                "genus":             genus,
                "species":           r.get("Species"),
                "property":          r.get("Property") or ("Antibiotic Resistance" if table == "amr_genes" else None),
                "source":            r.get("Source"),
                "level":             r.get("Level"),
                "host_name":         r.get("Host"),
            })
    return idx


def _contig_index(payload: dict) -> dict:
    idx = {}
    for cd in payload.get("contig_data") or []:
        idx[(cd.get("sample"), cd.get("organism"), str(cd.get("taxon_id", "")))] = cd
    return idx


# ── hierarchy rebuild ──────────────────────────────────────────────────────

def rebuild_sample(sample: str, rows: list, payload: dict, *,
                   include_nulls: bool = False, with_protein: bool = True) -> dict:
    meta = (payload.get("sample_meta") or {}).get(sample, {})
    contigs = _contig_index(payload)
    proteins = _protein_index(payload) if with_protein else {}

    by_level = defaultdict(list)
    for r in rows:
        by_level[r.get("Level", "Strain")].append(r)

    genus_nodes: dict[str, dict] = {}      # taxid -> node
    genus_by_name: dict[str, dict] = {}
    for r in by_level["Genus"]:
        node = _node_from_row(r, meta, include_nulls=include_nulls)
        _decorate_genus(node, r, meta)
        node["toplevelkey"] = str(r.get("Taxonomic ID #", ""))
        node["members"] = []
        genus_nodes[node["toplevelkey"]] = node
        genus_by_name.setdefault(node.get("name", ""), node)

    def _genus_for(row: dict) -> dict:
        """Find (or synthesize) the genus group a row belongs under."""
        gname = row.get("Genus Name") or row.get("Genus") or ""
        node = genus_by_name.get(gname)
        if node is None:
            node = {
                "toplevelkey": str(row.get("Taxonomic ID #", "")),
                "name": gname or row.get("Detected Organism", "Unknown"),
                "toplevelname": gname,
                "tass_score": row.get("Genus TASS"),
                "sample_name": sample,
                "members": [],
            }
            _decorate_genus(node, row, meta)
            genus_by_name[gname] = node
            genus_nodes[node["toplevelkey"] + "|synth|" + gname] = node
        return node

    species_nodes: dict[str, dict] = {}    # subkey/taxid -> node
    for r in by_level["Species"]:
        node = _node_from_row(r, meta, include_nulls=include_nulls)
        node["members"] = []
        gnode = _genus_for(r)
        node["toplevelkey"] = gnode["toplevelkey"]
        node["toplevelname"] = gnode.get("name", node.get("toplevelname"))
        gnode["members"].append(node)
        species_nodes[str(r.get("Taxonomic ID #", ""))] = node

    for r in by_level["Strain"]:
        node = _node_from_row(r, meta, include_nulls=include_nulls)
        subkey = str(r.get("Subkey", "") or "")
        parent = species_nodes.get(subkey)
        if parent is None:
            # Leaf subkey node: make_report emits a childless species group as a
            # "Strain" row, so it belongs directly under its genus group.
            gnode = _genus_for(r)
            node["toplevelkey"] = gnode["toplevelkey"]
            node["toplevelname"] = gnode.get("name", node.get("toplevelname"))
            node["members"] = []
            gnode["members"].append(node)
            species_nodes.setdefault(str(r.get("Taxonomic ID #", "")), node)
            target = node
        else:
            node["toplevelkey"] = parent["toplevelkey"]
            node["toplevelname"] = parent.get("toplevelname")
            parent["members"].append(node)
            target = node

        cd = contigs.get((sample, r.get("Detected Organism"), str(r.get("Taxonomic ID #", ""))))
        if cd:
            if cd.get("contigs"):
                target["contigs"] = cd["contigs"]
            if cd.get("depth_histogram"):
                target["depth_histogram"] = cd["depth_histogram"]
            if cd.get("breadth_histogram"):
                target["breadth_histogram"] = cd["breadth_histogram"]

    organisms = []
    for gnode in genus_nodes.values():
        if with_protein:
            hits = proteins.get((sample, gnode.get("name", "")))
            if hits:
                gnode["protein_annotations_genus"] = hits
        gnode["members"] = [
            _order({**sp, "members": [_order(st) for st in sp.get("members", [])]})
            for sp in gnode["members"]
        ]
        organisms.append(_order(gnode))

    organisms.sort(key=lambda o: -(o.get("tass_score") or 0))
    return {"metadata": meta, "organisms": organisms}


# ── cli ────────────────────────────────────────────────────────────────────

def _fidelity_report(payload: dict) -> str:
    lines = ["Field fidelity per organism node:", "",
             "  exact (verbatim from the report):"]
    exact = [f for f, _ in COLUMN_MAP.values() if f not in LOSSY]
    lines.append("    " + ", ".join(sorted(exact)))
    lines.append("  approximate (rounded when the report was built):")
    lines += [f"    {f:<22} {why}" for f, why in sorted(LOSSY.items())]
    lines.append("  dropped by flattening (null with --include-nulls):")
    lines.append("    " + ", ".join(sorted(UNRECOVERABLE)))
    lines.append("")
    lines.append("  metadata block is byte-for-byte the original.")
    return "\n".join(lines)


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("html", type=Path, help="report .html (or heatmap_boot.js / pages.js)")
    ap.add_argument("-d", "--outdir", type=Path, default=Path("."),
                    help="directory for the rebuilt <sample>.paths.json files")
    ap.add_argument("--sample", action="append",
                    help="only rebuild this sample (repeatable)")
    ap.add_argument("--combined", action="store_true",
                    help="write one combined all.samples.json instead of per-sample files")
    ap.add_argument("-o", "--output", type=Path, help="output path when --combined is used")
    ap.add_argument("--include-nulls", action="store_true",
                    help="emit unrecoverable fields as null so the key set matches the original")
    ap.add_argument("--no-protein", action="store_true",
                    help="skip rebuilding protein_annotations_genus")
    ap.add_argument("--indent", type=int, default=4, help="JSON indent (default 4, as produced upstream)")
    ap.add_argument("--report", action="store_true", help="print a field-fidelity summary")
    ap.add_argument("--stdout", action="store_true", help="print the JSON instead of writing files")
    args = ap.parse_args()

    if not args.html.exists():
        sys.exit(f"ERROR: not found: {args.html}")

    payload, kind = _ehb.extract(args.html)
    records = payload.get("records") or []
    if not records:
        sys.exit("ERROR: payload has no 'records' — nothing to rebuild")

    by_sample = defaultdict(list)
    for r in records:
        by_sample[r.get("Specimen ID")].append(r)

    wanted = args.sample or list(by_sample)
    missing = [s for s in wanted if s not in by_sample]
    if missing:
        sys.exit(f"ERROR: no records for sample(s): {', '.join(missing)}. "
                 f"Available: {', '.join(by_sample)}")

    built = {
        s: rebuild_sample(s, by_sample[s], payload,
                          include_nulls=args.include_nulls,
                          with_protein=not args.no_protein)
        for s in wanted
    }

    def _dump(obj) -> str:
        return json.dumps(obj, indent=args.indent, ensure_ascii=False, allow_nan=False)

    if args.combined:
        combined = {"taxtriage_combined": True,
                    "samples": [built[s] for s in wanted]}
        if args.stdout:
            print(_dump(combined))
        else:
            out = args.output or (args.outdir / "all.samples.json")
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_text(_dump(combined) + "\n", encoding="utf-8")
            print(f"✓ {out}  ({len(wanted)} samples)")
    elif args.stdout:
        for s in wanted:
            print(_dump(built[s]))
    else:
        args.outdir.mkdir(parents=True, exist_ok=True)
        for s in wanted:
            out = args.outdir / f"{s}.paths.json"
            out.write_text(_dump(built[s]) + "\n", encoding="utf-8")
            n_org = len(built[s]["organisms"])
            n_leaf = sum(len(m.get("members") or []) or 1
                         for g in built[s]["organisms"] for m in g.get("members") or [])
            print(f"✓ {out}  ({n_org} genus groups, {n_leaf} leaf organisms)")

    print(f"  source: {args.html.name} ({kind})", file=sys.stderr)
    if args.report:
        print("", file=sys.stderr)
        print(_fidelity_report(payload), file=sys.stderr)


if __name__ == "__main__":
    main()
