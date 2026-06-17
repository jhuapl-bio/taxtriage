#!/usr/bin/env python3
##############################################################################################
# Copyright 2025 The Johns Hopkins University Applied Physics Laboratory LLC
# All rights reserved.
# Permission is hereby granted, free of charge, to any person obtaining a copy of this
# software and associated documentation files (the "Software"), to deal in the Software
# without restriction, including without limitation the rights to use, copy, modify,
# merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
# PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
# OR OTHER DEALINGS IN THE SOFTWARE.
#

"""
make_report.py
==============
Build the all.comparison.report.html from either:
  a) one or more *.paths.json files produced by ALIGNMENT_PER_SAMPLE  (preferred)
  b) a single TSV/XLSX tabular file produced by ORGANISM_MERGE_REPORT (fallback)

Optionally accepts one or more protein-annotation XLSX files produced by the
--annotate_proteins / --annotate_meta steps (NOT use_diamond / get_features).
"""

import argparse
import glob
import json
import math
import os
import re
import sys

import pandas as pd


# ──────────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────────

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Build TaxTriage multi-run comparison HTML report.",
    )
    parser.add_argument(
        "-i", "--input", required=True, nargs="+",
        metavar="FILE",
        help="Input file(s): one or more .paths.json files  OR  a single "
             "TSV/XLSX tabular report (auto-detected by extension).",
    )
    parser.add_argument(
        "-p", "--protein_annotations", nargs="*", default=[],
        metavar="XLSX",
        help="Optional: protein-annotation XLSX file(s) produced by "
             "--annotate_proteins / --annotate_meta. "
             "Do NOT pass files from use_diamond or get_features here.",
    )
    parser.add_argument(
        "--mintass", default=1.0, type=float,
        help="Minimum TASS score for inclusion in the report (default: 1.0, i.e. include all). "
             "This is a hard filter; organisms below this threshold will be excluded entirely. "
             "Note that the UI filter slider is pre-populated from best_cutoffs in the input data, "
             "so you can use that to set a more conservative default while still allowing users to see all organisms if they wish."
    )
    parser.add_argument(
        "-t", "--template",
        metavar="TEMPLATE", default="heatmap.html",
        help="Input HTML template file (default: heatmap.html).",
    )
    parser.add_argument(
        "-pident", '--pident', default=0.0, type=float,
        help="Minimum percent identity (0-100) for inclusion in the report (default: 0.0, i.e. include all). Only used if you have VF/AMR results"
    )
    parser.add_argument(
        "-o", "--output",
        metavar="OUTPUT", default="all.comparison.report.html",
        help="Output HTML file.",
    )
    parser.add_argument(
        "-mc", "--microbial_category", nargs="+",
        default=["all"],
        metavar="CAT",
        help="Microbial category filter (default: Primary). "
             "Accepted values: Primary, Commensal, Opportunistic, Potential, Unknown. "
             "Use 'all' to include every category. "
             "Multiple values are allowed, e.g. -mc Primary Potential",
    )
    parser.add_argument(
        "-n", "--novelty", default=None, metavar="JSON",
        help="Optional: combined all.novelty.json produced by NOVELTY_COLLECT. When given, the "
             "report shows a Novelty Detection panel (per-sample score/flag + candidate taxa).",
    )
    parser.add_argument(
        "--novelty-downloads", nargs="*", default=[], metavar="FILE",
        help="Optional: per-sample + combined novelty JSON/XLSX files to expose as download links. "
             "Only their basenames are used (the files are published alongside the report).",
    )
    parser.add_argument(
        "--embedding", nargs="*", default=[], metavar="FILE",
        help="Optional: per-sample *.umap.tsv, *.cluster_summary.tsv, and *.clusters.tsv files "
             "from NOVEL_HOMOLOGS. When present, the report shows a Cluster Explorer subtab "
             "in the Novelty panel. Pass all files together; they are matched by sample name "
             "(the stem before the first dot).",
    )
    return parser.parse_args(argv)


# ──────────────────────────────────────────────────────────────────────────────
# JSON ingestion
# ──────────────────────────────────────────────────────────────────────────────

_TAX_RANKS = ["domain", "kingdom", "phylum", "class", "order", "family", "genus"]


def _collect_best_cutoffs(sample_meta):
    """Aggregate best_cutoffs across all loaded samples.

    For each granularity level (key/subkey/toplevelkey), takes the minimum
    best_threshold so the UI starts at the most conservative recommended cutoff.
    Returns a dict in the same shape as a single sample's best_cutoffs, or None.
    """
    levels = ("key", "subkey", "toplevelkey")
    aggregated = {}
    for level in levels:
        thresholds = []
        for meta in sample_meta.values():
            bc = meta.get("best_cutoffs") or {}
            t = (bc.get(level) or {}).get("best_threshold")
            if t is not None:
                thresholds.append(float(t))
        if thresholds:
            # Use the sample whose best_threshold is the minimum as the representative
            best_meta = min(
                (m for m in sample_meta.values() if (m.get("best_cutoffs") or {}).get(level)),
                key=lambda m: float((m.get("best_cutoffs", {}).get(level) or {}).get("best_threshold") or 9999),
            )
            aggregated[level] = dict((best_meta.get("best_cutoffs") or {}).get(level) or {})
            aggregated[level]["best_threshold"] = min(thresholds)
    return aggregated or None

_VALID_MICROBIAL_CATS = {"Primary", "Commensal", "Opportunistic", "Potential", "Unknown"}


def _resolve_microbial_cats(cat_args):
    """Return a set of accepted microbial categories, or None to mean 'all'."""
    if not cat_args:
        return None
    lowered = [c.strip().lower() for c in cat_args]
    if "all" in lowered:
        return None  # no filtering
    result = set()
    for raw in cat_args:
        canon = raw.strip().capitalize()
        # handle "opportunistic" -> "Opportunistic", title-case normalisation
        # do a case-insensitive lookup against valid set
        matched = next(
            (v for v in _VALID_MICROBIAL_CATS if v.lower() == raw.strip().lower()),
            None,
        )
        if matched:
            result.add(matched)
        else:
            print(
                f"[make_report] WARNING: --microbial_category value {raw!r} is not recognised "
                f"(valid: {sorted(_VALID_MICROBIAL_CATS)} or 'all'); ignoring.",
                file=sys.stderr,
            )
    return result or None  # if everything was invalid, fall back to no filter


def _flatten_organism(org, sample_name, sample_type, total_reads,
                      species_parent=None, genus_parent=None, level="Strain"):
    """
    Flatten one organism entry (any hierarchy level) into a flat dict
    suitable for the tabular view and all plots.

    species_parent / genus_parent: the subkey (species) and toplevelkey (genus)
    group objects this organism rolls up into. Their tass_score is attached as
    Species TASS / Genus TASS so the UI can show that a strain failing its own
    threshold may still be detected at the species or genus level (LCA-aware
    rollup). For a row that IS already the species/genus level, the parent of
    that same level points to itself.
    """
    strain_reads = float(org.get("numreads", 0) or 0)
    pct = strain_reads / max(1, total_reads) * 100.0
    covered = int(org.get("covered_bases", 0) or 0)
    genome_len = int(org.get("length", 0) or 0)
    breadth_pct = round(min(100.0, covered / genome_len * 100), 2) if genome_len > 0 else 0.0
    tass = float(org.get("tass_score", 0) or 0)

    tax = org.get("taxonomy", {})

    # ── ANI annotation (set by match_paths.py when the ANI matrix is enabled) ──
    # Presence of the 'high_ani_matches' key — even as an empty list — signals
    # that ANI was computed for this run. Its absence means the data predates
    # ANI support, so ANI-dependent views fall back to an "out of date /
    # unsupported" state for this sample.
    _ani_annotated = 'high_ani_matches' in org
    _ani_list = [
        {"key": str(m.get("key", "")), "ani_pct": round(float(m.get("ani_pct", 0) or 0), 2)}
        for m in (org.get('high_ani_matches') or []) if isinstance(m, dict)
    ]

    # Parent-level TASS (species = subkey, genus = toplevelkey). Fall back to the
    # organism's own TASS when a parent level is absent so the rollup never
    # under-reports the row itself.
    _species_src = species_parent if species_parent is not None else org
    _genus_src = genus_parent if genus_parent is not None else (species_parent if species_parent is not None else org)
    species_tass = float(_species_src.get("tass_score", tass) or 0)
    genus_tass = float(_genus_src.get("tass_score", tass) or 0)
    species_name = str(_species_src.get("name", org.get("name", "")) or "")
    genus_name = str(_genus_src.get("name", "") or tax.get("genus", "") or "")

    return {
        "Specimen ID":         sample_name,
        "Sample Type":         sample_type,
        "Detected Organism":   org.get("name", "Unknown"),
        "TASS Score":          round(tass, 100),  # 4th col — 0–100 scale for display
        "Taxonomic ID #":      str(org.get("key", "")),
        "Subkey":              str(org.get("subkey", org.get("key", ""))),
        "Microbial Category":  org.get("microbial_category", "Unknown"),
        "Ann Class":           org.get("annClass", ""),
        "IsAnnotated":         "Yes" if org.get("is_annotated", "No") == "Yes" else "No",
        "High Consequence":    bool(org.get("high_cons", False)),
        "Mol Type":            org.get("mol_type", ""),
        "Status":              org.get("status", ""),
        "# Reads Aligned":     int(strain_reads),
        "% Reads":             round(pct, 4),
        "Coverage":            round(min(100.0, (org.get("coverage", 0) or 0) * 100), 1),
        "Covered Bases":       covered,
        "Genome Length (bp)":  genome_len,
        "Breadth %":           breadth_pct,
        "Mean Depth":          round(float(org.get("meandepth", 0) or 0), 2),
        "Gini Coefficient":    round(float(org.get("gini_coefficient", 0) or 0), 3),
        "Mean MapQ":           round(float(org.get("meanmapq", 0) or 0), 1),
        "Mean BaseQ":          round(float(org.get("meanbaseq", 0) or 0), 1),
        "Minhash Score":       round(float(org.get("minhash_reduction", 0) or 0), 3),
        "Breadth Score":       round(float(org.get("breadth_log_score", 0) or 0), 3),
        "MapQ Score":          round(float(org.get("mapq_score", 0) or 0), 3),
        "Disparity Score":     round(float(org.get("disparity", 0) or 0), 3),
        "Diamond Identity":    round(float(org.get("diamond_identity", 0) or 0), 1),
        # MicrobeRT (mmbert) classifier probability + model name. mmbert is a
        # 0–1 probability; surfaced here as a 0–100 % to match the PDF report.
        # Kept None when absent so the column renders blank rather than 0.
        "MicrobeRT Probability": (round(float(org.get("mmbert")) * 100, 2)
                                  if org.get("mmbert") not in (None, "")
                                  else None),
        "MicrobeRT Model":     org.get("mmbert_model") or "",
        "K2 Reads":            int(org.get("k2_reads", 0) or 0),
        "RPM":                 round(float(org.get("rpm", 0) or 0), 2),
        "RPKM":                round(float(org.get("rpkm", 0) or 0), 4),
        "Passes Threshold":    bool(org.get("passes_threshold", False)),
        # ANI annotation: capability flag + list of high-ANI partner taxa
        # ({key, ani_pct}). Consumed by the cross-sample Feature Compare view and
        # by client-side capability detection (absence ⇒ ANI unsupported).
        "ANI Annotated":       _ani_annotated,
        "High ANI Matches":    _ani_list,
        # Taxonomic rollup level of this row: "Strain" (key), "Species" (subkey),
        # or "Genus" (toplevelkey). Lets the UI switch the view granularity and
        # surface a species/genus summary row when its children fail their own
        # threshold. Strain rows are the default view.
        "Level":               level,
        # Parent-level rollup TASS — lets the UI show a strain that fails its own
        # threshold but is still detected at the species (subkey) or genus
        # (toplevelkey) level. Pass/fail itself is computed client-side against
        # the active threshold, so only the scores are emitted here.
        "Species TASS":        round(species_tass, 4),
        "Genus TASS":          round(genus_tass, 4),
        "Species Name":        species_name,
        "Genus Name":          genus_name,
        # taxonomy
        "Kingdom":             tax.get("kingdom", ""),
        "Domain":             tax.get("domain", ""),
        "Superkingdom":        tax.get("superkingdom", ""),
        "Phylum":              tax.get("phylum", ""),
        "Class":               tax.get("class", ""),
        "Order":               tax.get("order", ""),
        "Family":              tax.get("family", ""),
        "Genus":               tax.get("genus", ""),
    }


def _iter_organisms(json_data, sample_name, mintass=0, microbial_cats=None):
    """Yield flat organism dicts from a parsed paths JSON.

    microbial_cats: set of accepted microbial_category strings, or None to allow all.
    """
    meta = json_data.get("metadata", {})
    sample_type = meta.get("sample_type", "unknown")
    total_reads = int(meta.get("total_reads", 1) or 1)

    def _cat_ok(org):
        if microbial_cats is None:
            return True
        return org.get("microbial_category", "Unknown") in microbial_cats

    # grp = toplevelkey (genus) group; sk_m = subkey (species) group; strain = key.
    #
    # Leaf nodes (the organisms actually shown in the default view) are emitted as
    # Level="Strain". For every group that has children we ALSO emit a summary row
    # one level up — Level="Species" for a subkey group with strain members, and
    # Level="Genus" for a genus group with members — carrying that group's own
    # aggregate TASS/coverage. These summary rows let the UI roll up to the
    # species/genus level and surface a species hit even when every child strain
    # falls below the active cutoff. They are hidden in the default Strain view
    # unless promoted by the rollup.
    for grp in json_data.get("organisms", []):
        grp_has_members = bool(grp.get("members"))
        for sk_m in grp.get("members", []):
            sk_has_members = bool(sk_m.get("members"))
            for strain in sk_m.get("members", []):
                if float(strain.get("tass_score", 0) or 0) >= mintass and _cat_ok(strain):
                    yield _flatten_organism(strain, sample_name, sample_type, total_reads,
                                            species_parent=sk_m, genus_parent=grp,
                                            level="Strain")
            if sk_has_members:
                # species summary row over its strain members
                if float(sk_m.get("tass_score", 0) or 0) >= mintass and _cat_ok(sk_m):
                    yield _flatten_organism(sk_m, sample_name, sample_type, total_reads,
                                            species_parent=sk_m, genus_parent=grp,
                                            level="Species")
            else:
                # no nested strains: the subkey node is itself the leaf
                if float(sk_m.get("tass_score", 0) or 0) >= mintass and _cat_ok(sk_m):
                    yield _flatten_organism(sk_m, sample_name, sample_type, total_reads,
                                            species_parent=sk_m, genus_parent=grp,
                                            level="Strain")
        if grp_has_members:
            # genus summary row over its species/strain members
            if float(grp.get("tass_score", 0) or 0) >= mintass and _cat_ok(grp):
                yield _flatten_organism(grp, sample_name, sample_type, total_reads,
                                        species_parent=grp, genus_parent=grp,
                                        level="Genus")
        else:
            # no members at all: the genus node is itself the leaf
            if float(grp.get("tass_score", 0) or 0) >= mintass and _cat_ok(grp):
                yield _flatten_organism(grp, sample_name, sample_type, total_reads,
                                        species_parent=grp, genus_parent=grp,
                                        level="Strain")


def load_json_inputs(paths, mintass=0, microbial_cats=None):
    """Return flat organism rows, per-sample metadata, and per-organism contig data."""
    rows = []
    sample_meta = {}
    # contig_data: dict keyed by "<sample>||<organism_name>||<taxon_id>"
    # value: {contigs: [...], depth_histogram: {...}}
    contig_data = {}

    for path in paths:
        path = path.strip()
        if not os.path.isfile(path):
            expanded = glob.glob(path)
            if not expanded:
                print(f"[make_report] WARNING: cannot find {path!r}, skipping", file=sys.stderr)
                continue
            for p in expanded:
                rows_, sm_, cd_ = load_json_inputs([p], mintass=mintass, microbial_cats=microbial_cats)
                rows.extend(rows_)
                sample_meta.update(sm_)
                contig_data.update(cd_)
            continue

        try:
            with open(path) as fh:
                data = json.load(fh)
        except Exception as exc:
            print(f"[make_report] WARNING: failed to parse {path}: {exc}", file=sys.stderr)
            continue

        # ── Combined all.samples.json: expand into per-sample loads ──────────
        if data.get("taxtriage_combined") and isinstance(data.get("samples"), list):
            print(f"[make_report] Detected combined JSON {path!r} with "
                  f"{len(data['samples'])} sample(s); expanding...")
            for sample_data in data["samples"]:
                s_meta = sample_data.get("metadata", {})
                s_name = s_meta.get("sample_name",
                                     os.path.basename(path).split(".")[0])
                sample_meta[s_name] = s_meta
                for row in _iter_organisms(sample_data, s_name,
                                           mintass=mintass,
                                           microbial_cats=microbial_cats):
                    rows.append(row)
                _STRIP = {'members', 'subkey', 'key', 'toplevelkey'}
                for grp in sample_data.get("organisms", []):
                    for sk_m in grp.get("members", []):
                        if float(sk_m.get('tass_score', 0) or 0) < mintass:
                            continue
                        if microbial_cats is not None and sk_m.get('microbial_category', 'Unknown') not in microbial_cats:
                            continue
                        for strain in sk_m.get("members", []):
                            if float(strain.get('tass_score', 0) or 0) < mintass:
                                continue
                            if microbial_cats is not None and strain.get('microbial_category', 'Unknown') not in microbial_cats:
                                continue
                            _contigs = strain.get("contigs")
                            _dhist   = strain.get("depth_histogram")
                            _bhist   = strain.get("breadth_histogram")
                            if _contigs or _dhist or _bhist:
                                _key = f"{s_name}||{strain.get('name','')}||{strain.get('key','')}"
                                _cd_entry = {
                                    "sample":          s_name,
                                    "organism":        strain.get("name", "Unknown"),
                                    "taxon_id":        str(strain.get("key", "")),
                                    "contigs":         [{k: v for k, v in c.items() if k not in _STRIP} for c in (_contigs or [])],
                                    "depth_histogram": _dhist or {},
                                }
                                if _bhist:
                                    _cd_entry["breadth_histogram"] = _bhist
                                contig_data[_key] = _cd_entry
            continue

        meta = data.get("metadata", {})
        sample_name = meta.get("sample_name", os.path.basename(path).split(".")[0])
        sample_meta[sample_name] = meta

        for row in _iter_organisms(data, sample_name, mintass=mintass, microbial_cats=microbial_cats):
            rows.append(row)

        # Extract per-contig and depth-histogram data from each strain
        _STRIP = {'members', 'subkey', 'key', 'toplevelkey'}
        for grp in data.get("organisms", []):
            for sk_m in grp.get("members", []):
                if float(sk_m.get('tass_score', 0) or 0) < mintass:
                    continue
                if microbial_cats is not None and sk_m.get('microbial_category', 'Unknown') not in microbial_cats:
                    continue
                for strain in sk_m.get("members", []):
                    if float(strain.get('tass_score', 0) or 0) < mintass:
                        continue
                    if microbial_cats is not None and strain.get('microbial_category', 'Unknown') not in microbial_cats:
                        continue
                    _contigs = strain.get("contigs")
                    _dhist   = strain.get("depth_histogram")
                    _bhist   = strain.get("breadth_histogram")
                    if _contigs or _dhist or _bhist:
                        _key = f"{sample_name}||{strain.get('name','')}||{strain.get('key','')}"
                        _cd_entry = {
                            "sample":          sample_name,
                            "organism":        strain.get("name", "Unknown"),
                            "taxon_id":        str(strain.get("key", "")),
                            "contigs":         [{k: v for k, v in c.items() if k not in _STRIP} for c in (_contigs or [])],
                            "depth_histogram": _dhist or {},
                        }
                        if _bhist:
                            _cd_entry["breadth_histogram"] = _bhist
                        contig_data[_key] = _cd_entry

    return rows, sample_meta, contig_data


# ──────────────────────────────────────────────────────────────────────────────
# TSV / XLSX fallback ingestion
# ──────────────────────────────────────────────────────────────────────────────

def load_tabular_input(path, mintass=0, microbial_cats=None):
    ext = os.path.splitext(path)[1].lower()
    if ext in (".xlsx", ".xls"):
        df = pd.read_excel(path, dtype=str)
    else:
        df = pd.read_csv(path, sep="\t", dtype=str)
    df = df[df['TASS Score'].astype(float) >= mintass].copy()
    if microbial_cats is not None and 'Microbial Category' in df.columns:
        df = df[df['Microbial Category'].isin(microbial_cats)].copy()
    df.columns = df.columns.str.strip()
    if "Detected Organism" in df.columns:
        df["Detected Organism"] = df["Detected Organism"].str.replace("°", "", regex=False).str.strip()
    if "Index" in df.columns:
        df = df.drop(columns=["Index"], errors="ignore")
    if "index" in df.columns:
        df = df.drop(columns=["index"], errors="ignore")
    df = df.where(pd.notnull(df), None)
    return df.to_dict(orient="records"), {}


# ──────────────────────────────────────────────────────────────────────────────
# Protein annotation ingestion  (only --annotate_proteins / --annotate_meta)
# ──────────────────────────────────────────────────────────────────────────────

def load_novelty(novelty_json, download_paths):
    """
    Read the combined all.novelty.json (NOVELTY_COLLECT) and return
      (novelty_payload, downloads)
    where novelty_payload = {"samples": {<sample>: {"summary": {...}, "candidates": [...]}}}
    and downloads = [{"label": <sample or "All samples">, "kind": "json"|"xlsx",
                      "filename": <basename>}, ...] for the report's download links.
    Returns ({"samples": {}}, []) when nothing usable is provided.
    """
    payload = {"samples": {}}
    if novelty_json:
        p = novelty_json.strip()
        if p and os.path.isfile(p):
            try:
                with open(p, encoding="utf-8") as fh:
                    data = json.load(fh)
                payload["samples"] = data.get("samples", {}) or {}
            except Exception as exc:  # noqa: BLE001
                print(f"[make_report] WARNING: cannot read novelty json {p}: {exc}",
                      file=sys.stderr)

    # Build download link descriptors from the staged files (basenames only — the files are
    # published next to all.odr.html, so relative links resolve).
    downloads = []
    seen = set()
    for path in (download_paths or []):
        name = os.path.basename((path or "").strip())
        if not name or name in seen or name.startswith("NO_FILE") or name.startswith("~"):
            continue
        if not (name.endswith(".json") or name.endswith(".xlsx")):
            continue
        seen.add(name)
        kind = "xlsx" if name.endswith(".xlsx") else "json"
        if name.startswith("all.novelty"):
            label = "All samples (combined)"
        else:
            # strip the .novelty.json / .novelty.xlsx suffix to recover the sample name
            label = name.split(".novelty.")[0]
        downloads.append({"label": label, "kind": kind, "filename": name})

    # Stable order: combined first, then per-sample alphabetical, xlsx before json within a group.
    downloads.sort(key=lambda d: (not d["filename"].startswith("all.novelty"),
                                  d["label"].lower(), d["kind"]))
    return payload, downloads


def load_embedding(embed_files, max_scatter_pts=8000, max_singleton_pts=1500):
    """
    Read per-sample *.umap.tsv, *.cluster_summary.tsv, and *.clusters.tsv files
    produced by the NOVEL_HOMOLOGS process and return a compact embedding payload:

      {
        "samples": {
          "<sample>": {
            "umap":    [{"id","set","cluster_id","x","y"}, ...],   # sampled
            "summary": [{"cluster_id","size","novel_family","dominant_ref_class",
                         "mean_outlier_score","cx","cy","r"}, ...],  # with centroid
            "clusters": [{"query_id","cluster_id","novel_flag","outlier_score",
                          "cluster_size"}, ...],               # for hover/colour
          }
        }
      }

    Performance strategy:
      - Pre-compute cluster centroid (cx, cy) and enclosure radius (r) in Python so
        the JS never has to iterate all points to draw circles.
      - Keep all non-singleton clustered points (usually <20k) for the scatter.
      - Randomly sample singletons down to max_singleton_pts.
      - Total scatter points capped at max_scatter_pts (excess singletons trimmed first).
    Returns {"samples": {}} when nothing usable is provided.
    """
    import math as _math
    import random as _random

    payload = {"samples": {}}
    if not embed_files:
        return payload

    # Group files by sample name (stem before first dot that isn't a known suffix)
    umap_files, summary_files, cluster_files = {}, {}, {}
    for p in (embed_files or []):
        p = (p or "").strip()
        if not p or not os.path.isfile(p):
            continue
        base = os.path.basename(p)
        if base.startswith("NO_FILE") or base.startswith("~"):
            continue
        # Derive sample name: everything before the first recognised suffix
        for suffix in (".umap.tsv", ".cluster_summary.tsv", ".clusters.tsv"):
            if base.endswith(suffix):
                sname = base[: -len(suffix)]
                if suffix == ".umap.tsv":
                    umap_files[sname] = p
                elif suffix == ".cluster_summary.tsv":
                    summary_files[sname] = p
                elif suffix == ".clusters.tsv":
                    cluster_files[sname] = p
                break

    all_samples = set(umap_files) | set(summary_files) | set(cluster_files)
    for sname in sorted(all_samples):
        s_data = {}

        # ── cluster_summary ──────────────────────────────────────────────────
        summary_rows = []
        if sname in summary_files:
            try:
                df = pd.read_csv(summary_files[sname], sep="\t", dtype=str)
                summary_rows = df.where(pd.notnull(df), None).to_dict(orient="records")
            except Exception as exc:
                print(f"[make_report] WARNING: cannot read {summary_files[sname]}: {exc}",
                      file=sys.stderr)

        # ── clusters (per-ORF novelty/cluster assignment) ────────────────────
        cluster_rows = []
        if sname in cluster_files:
            try:
                df = pd.read_csv(cluster_files[sname], sep="\t", dtype=str)
                keep = [c for c in ("query_id", "cluster_id", "novel_flag",
                                    "outlier_score", "cluster_size") if c in df.columns]
                cluster_rows = df[keep].where(pd.notnull(df[keep]), None).to_dict(orient="records")
            except Exception as exc:
                print(f"[make_report] WARNING: cannot read {cluster_files[sname]}: {exc}",
                      file=sys.stderr)

        # ── umap (x/y coordinates) + pre-compute centroids ──────────────────
        umap_rows = []
        if sname in umap_files:
            try:
                df = pd.read_csv(umap_files[sname], sep="\t", dtype=str)
                df["_x"] = pd.to_numeric(df["x"], errors="coerce")
                df["_y"] = pd.to_numeric(df["y"], errors="coerce")
                df = df.dropna(subset=["_x", "_y"])

                # Separate singletons (cluster_id == -1) from clustered points
                df["_cid"] = df["cluster_id"].astype(str)
                clustered  = df[df["_cid"] != "-1"]
                singletons = df[df["_cid"] == "-1"]

                # Pre-compute centroid + 90th-pct enclosure radius per cluster
                centroid_map = {}
                for cid, grp in clustered.groupby("_cid"):
                    cx = grp["_x"].mean()
                    cy = grp["_y"].mean()
                    dists = sorted(
                        ((grp["_x"] - cx)**2 + (grp["_y"] - cy)**2).pow(0.5).tolist()
                    )
                    idx90 = min(len(dists) - 1, int(0.9 * len(dists)))
                    r = round(dists[idx90] * 1.25 + 0.05, 4)
                    centroid_map[cid] = {"cx": round(cx, 4), "cy": round(cy, 4), "r": r}

                # Merge centroids into summary rows
                for row in summary_rows:
                    cid = str(row.get("cluster_id", ""))
                    if cid in centroid_map:
                        row.update(centroid_map[cid])

                # Sample singletons
                sing_sample = singletons
                if len(singletons) > max_singleton_pts:
                    sing_sample = singletons.sample(max_singleton_pts, random_state=42)

                # Budget for clustered points: max_scatter_pts minus singletons
                clustered_budget = max(max_scatter_pts - len(sing_sample), 500)
                clust_sample = clustered
                if len(clustered) > clustered_budget:
                    # Proportional stratified sample: at least 1 pt per cluster
                    n_clusters = clustered["_cid"].nunique()
                    min_per = max(1, min(3, clustered_budget // max(1, n_clusters)))
                    remain = clustered_budget - n_clusters * min_per
                    if remain < 0:
                        remain = 0
                    sampled_parts = []
                    for cid, grp in clustered.groupby("_cid"):
                        quota = min_per
                        if remain > 0 and len(grp) > min_per:
                            extra = max(0, round(remain * len(grp) / max(1, len(clustered))))
                            quota = min_per + extra
                        sampled_parts.append(grp.sample(min(quota, len(grp)), random_state=42))
                    clust_sample = pd.concat(sampled_parts)

                pts = pd.concat([clust_sample, sing_sample]).reset_index(drop=True)

                cols = [c for c in ("id", "set", "cluster_id", "x", "y") if c in pts.columns]
                umap_rows = pts[cols].where(pd.notnull(pts[cols]), None).to_dict(orient="records")

            except Exception as exc:
                print(f"[make_report] WARNING: cannot read {umap_files[sname]}: {exc}",
                      file=sys.stderr)

        if umap_rows or summary_rows or cluster_rows:
            s_data["umap"]    = umap_rows
            s_data["summary"] = summary_rows
            s_data["clusters"] = cluster_rows
            payload["samples"][sname] = s_data
            n_pts = len(umap_rows)
            n_cl  = len(summary_rows)
            print(f"[make_report] Embedding [{sname}]: {n_pts} scatter pts, {n_cl} clusters")

    return payload


def load_protein_annotations(paths, pident=0):
    """
    Read one or more protein-annotation XLSX files (sheets: Genus Summary,
    Per-Gene Hits, Sample Overview, AMR Genes) and return a dict:
      {
        "genus_summary":   [row, ...],
        "per_gene_hits":   [row, ...],
        "sample_overview": [row, ...],
        "amr_genes":       [row, ...],
      }
    Rows from multiple files are concatenated.
    """
    out = {
        "genus_summary": [],
        "per_gene_hits": [],
        "sample_overview": [],
        "amr_genes": [],
    }
    sheet_map = {
        "Genus Summary":   "genus_summary",
        "Per-Gene Hits":   "per_gene_hits",
        "Sample Overview": "sample_overview",
        "AMR Genes":       "amr_genes",
    }

    for path in (paths or []):
        path = path.strip()
        if not path or not os.path.isfile(path):
            continue
        try:
            wb = pd.read_excel(path, sheet_name=None, dtype=str)
        except Exception as exc:
            print(f"[make_report] WARNING: cannot read {path}: {exc}", file=sys.stderr)
            continue
        for sheet_name, key in sheet_map.items():
            if sheet_name in wb:
                df = wb[sheet_name].where(pd.notnull(wb[sheet_name]), None)
                # filter where %id is >= pident (if that column exists)
                if '%id' in df.columns:
                    df = df[df['%id'].astype(float) >= pident].copy()
                out[key].extend(df.to_dict(orient="records"))
    return out


# ──────────────────────────────────────────────────────────────────────────────
# JSON sanitisation helpers
# ──────────────────────────────────────────────────────────────────────────────

def _sanitize(obj):
    if isinstance(obj, dict):
        return {k: _sanitize(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_sanitize(v) for v in obj]
    try:
        if pd.isna(obj):
            return None
    except Exception:
        pass
    if isinstance(obj, float):
        return obj if math.isfinite(obj) else None
    if hasattr(obj, "item"):       # numpy scalar
        return obj.item()
    return obj


# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    # ── detect input mode ─────────────────────────────────────────────────────
    is_json_mode = all(
        f.strip().endswith(".json")
        for f in args.input
        if f.strip()
    )

    microbial_cats = _resolve_microbial_cats(args.microbial_category)
    if microbial_cats is None:
        print("[make_report] Microbial category filter: all")
    else:
        print(f"[make_report] Microbial category filter: {sorted(microbial_cats)}")

    # All organisms with TASS > 0 are always included; the UI pre-populates its
    # filter slider from best_cutoffs baked into the BOOT payload at load time.
    mintass = args.mintass
    print("[make_report] TASS threshold: 0.0 (all organisms included; UI filter set from best_cutoffs)")

    if is_json_mode:
        rows, sample_meta, contig_data = load_json_inputs(
            args.input, mintass=mintass, microbial_cats=microbial_cats
        )
        print(f"[make_report] Loaded {len(rows)} organism rows from "
              f"{len(args.input)} JSON file(s); {len(contig_data)} organisms have contig data")
    else:
        contig_data = {}
        if len(args.input) > 1:
            print("[make_report] WARNING: multiple non-JSON inputs given; "
                  "using only the first.", file=sys.stderr)
        rows, sample_meta = load_tabular_input(args.input[0], mintass, microbial_cats=microbial_cats)
        print(f"[make_report] Loaded {len(rows)} rows from tabular file "
              f"{args.input[0]!r}")

    # ── derive column lists ────────────────────────────────────────────────────
    # These fields are carried on each record for client-side analysis (the
    # Feature Compare view + capability detection) but are NOT human-displayable
    # table columns — 'High ANI Matches' is a nested list — so keep them out of
    # the column picker / detections table.
    _NON_DISPLAY_COLS = {"High ANI Matches", "ANI Annotated"}
    all_cols = [c for c in (rows[0].keys() if rows else []) if c not in _NON_DISPLAY_COLS]
    numeric_cols = []
    if rows:
        for col in all_cols:
            vals = [r[col] for r in rows if r.get(col) is not None]
            if vals and all(
                isinstance(v, (int, float)) or
                (isinstance(v, str) and _is_numeric_str(v))
                for v in vals[:50]
            ):
                numeric_cols.append(col)

    # ── protein annotations ───────────────────────────────────────────────────
    prot_data = load_protein_annotations(args.protein_annotations, pident=args.pident)
    has_prot = any(len(v) > 0 for v in prot_data.values())
    print(f"[make_report] Protein annotations loaded: {has_prot} "
          f"({sum(len(v) for v in prot_data.values())} total rows)")

    # ── novelty detection (reference-free LCA) ────────────────────────────────
    novelty_data, novelty_downloads = load_novelty(args.novelty, args.novelty_downloads)
    has_novelty = bool(novelty_data.get("samples"))
    print(f"[make_report] Novelty loaded: {has_novelty} "
          f"({len(novelty_data.get('samples', {}))} sample(s), "
          f"{len(novelty_downloads)} download link(s))")

    # ── embedding / UMAP cluster data (NOVEL_HOMOLOGS) ───────────────────────
    embedding_data = load_embedding(args.embedding)
    has_embedding = bool(embedding_data.get("samples"))
    print(f"[make_report] Embedding loaded: {has_embedding} "
          f"({len(embedding_data.get('samples', {}))} sample(s))")

    # ── extract run-level metadata for map / metadata panel ──────────────────
    # These fields are part of the fixed pipeline/sample identity — not run metadata.
    _META_PIPELINE_KEYS = {
        "sample_name", "sample_type", "platform", "workflow_revision", "commit_id",
        "total_reads", "aligned_reads", "total_organism_reads", "num_species_groups",
        "num_keys", "num_subkeys", "num_toplevelkeys", "minmapq", "mapq_breadth_power",
        "weights", "min_conf_applied", "best_cutoffs", "best_cutoffs_source",
        "preferred_granularity", "control_type", "negative_controls_used",
        "positive_controls_used", "control_fold_threshold", "missing_positive_controls",
        "insilico_controls_used", "insilico_simulator_types", "missing_insilico_controls",
        "missing_insilico_by_type",
    }
    run_metadata_records = []
    for sname, smeta in sample_meta.items():
        rec = {"sample_name": sname}
        # Include ALL metadata fields that are not fixed pipeline/identity fields
        for k, v in smeta.items():
            if k not in _META_PIPELINE_KEYS and k != "sample_name":
                rec[k] = v
        # Only add to records if at least one field has a non-null, non-empty value
        if any(v is not None and v != "" for k, v in rec.items() if k != "sample_name"):
            run_metadata_records.append(rec)

    has_geo = any(
        r.get("latitude") is not None and r.get("longitude") is not None
        for r in run_metadata_records
    )

    # ── collect best_cutoffs for UI pre-population ────────────────────────────
    best_cutoffs_payload = _collect_best_cutoffs(sample_meta)
    if best_cutoffs_payload:
        thresh = (best_cutoffs_payload.get("subkey") or {}).get("best_threshold")
        print(f"[make_report] UI will pre-set TASS filter to: {thresh} (from best_cutoffs.subkey)")
    else:
        print("[make_report] No best_cutoffs found; UI TASS filter will default to 0")

    # ── build bootstrap payload ───────────────────────────────────────────────
    payload = _sanitize({
        "records":               rows,
        "all_cols":              all_cols,
        "numeric_cols":          numeric_cols,
        "sample_meta":           sample_meta,
        "prot_data":             prot_data,
        "has_prot":              has_prot,
        "contig_data":           list(contig_data.values()),   # list of organism contig objects
        "best_cutoffs":          best_cutoffs_payload,         # for UI filter pre-population
        "run_metadata_records":  run_metadata_records,         # per-sample run metadata
        "has_geo":               has_geo,                      # true if any sample has lat/lon
        "novelty":               novelty_data,                 # {samples: {<s>: {summary, candidates}}}
        "has_novelty":           has_novelty,                  # true if any novelty sample present
        "novelty_downloads":     novelty_downloads,            # [{label, kind, filename}] for links
        "embedding":             embedding_data,               # {samples: {<s>: {umap, summary, clusters}}}
        "has_embedding":         has_embedding,                # true if any embedding data present
    })

    bootstrap_json = json.dumps(payload, ensure_ascii=False, allow_nan=False, separators=(',', ':'))

    # Build inline JS instead of writing a separate heatmap_boot.js file
    bootstrap_script = f"<script>\nwindow.HEATMAP_BOOT = {bootstrap_json};\n</script>"

    # ── render template ───────────────────────────────────────────────────────
    with open(args.template, "r", encoding="utf-8") as fh:
        tpl = fh.read()

    # Use a function replacement (not a plain string) so backslashes / `\g`-like
    # sequences inside the JSON payload are inserted verbatim and never treated
    # as regex backreferences.
    _repl = lambda m: bootstrap_script

    # ── Strip the pages.js demo loader (active or commented) ──────────────────
    # The real dataset is injected inline below. pages.js is not copied next to
    # the report, so a leftover tag 404s in the browser console.
    html = re.sub(
        r'[ \t]*<!--\s*<script[^>]+src=["\']pages\.js["\'][^>]*>\s*</script>\s*-->[ \t]*\n?',
        '', tpl, flags=re.IGNORECASE,
    )
    html = re.sub(
        r'[ \t]*<script[^>]+src=["\']pages\.js["\'][^>]*>\s*</script>[ \t]*\n?',
        '', html, flags=re.IGNORECASE,
    )

    # ── Inject the inline bootstrap at the heatmap_boot.js anchor ─────────────
    # Handle the commented placeholder, an active tag, or an existing BOOTSTRAP
    # block. A *commented* anchor MUST be matched together with its <!-- … -->
    # so the boot script is not left commented out — that would blank the report.
    anchor_patterns = [
        r'<!--\s*<script[^>]+src=["\']heatmap_boot\.js["\'][^>]*>\s*</script>\s*-->',
        r'<script[^>]+src=["\']heatmap_boot\.js["\'][^>]*>\s*</script>',
        r'<!--\s*<script id=["\']BOOTSTRAP["\'][^>]*>.*?</script>\s*-->',
        r'<script id=["\']BOOTSTRAP["\'][^>]*>.*?</script>',
    ]
    n_replaced = 0
    for pat in anchor_patterns:
        html, n_replaced = re.subn(pat, _repl, html, count=1, flags=re.DOTALL | re.IGNORECASE)
        if n_replaced:
            break

    if not n_replaced and "__BOOTSTRAP_SCRIPT__" in html:
        # Fallback: explicit placeholder token.
        html = html.replace("__BOOTSTRAP_SCRIPT__", bootstrap_script)
        n_replaced = 1

    if not n_replaced:
        raise SystemExit(
            "[make_report] ERROR: no heatmap_boot.js / BOOTSTRAP anchor found in template"
        )

    # ── Guard: the boot script must not have landed inside an HTML comment ────
    _idx = html.find("window.HEATMAP_BOOT =")
    if _idx != -1:
        _open = html.rfind("<!--", 0, _idx)
        _close = html.rfind("-->", 0, _idx)
        if _open != -1 and _open > _close:
            raise SystemExit(
                "[make_report] ERROR: bootstrap script was inserted inside an HTML "
                "comment — the report would be blank. Check the template anchor."
            )

    with open(args.output, "w", encoding="utf-8") as fh:
        fh.write(html)

    print(f"[make_report] Written: {args.output}")


def _is_numeric_str(s):
    try:
        float(s)
        return True
    except (ValueError, TypeError):
        return False


if __name__ == "__main__":
    main()
