#!/usr/bin/env python3
"""
HMP "healthy abundance" outlier annotation for the interactive report.

The ODR PDF (create_report.py) marks each detection whose abundance sits inside
the range expected for a healthy subject at that body site: a faded row plus a
diamond showing how many healthy reference samples carried the organism
(``◆ n (p%)◆``). Those numbers come from the HMP abundance distributions in
``assets/taxid_abundance_stats.hmp.tsv.gz`` (one row per taxid x body site,
with the per-sample relative abundances the mean / std were computed from).

This module computes the same figures for the interactive report
(make_report.py) at report-build time, independently of TASS scoring, so turning
it on never changes a score:

    observed  = % of the sample's reads assigned to the row (the report's "% Reads")
    z         = (observed - healthy mean) / healthy std        (summed across sites)
    elevated  = z >= threshold   (default 2.0, same as create_report --zscore_threshold)

Healthy means/stds are taken over the reference samples in which the organism
was present (that is how the HMP table is built), so prevalence — the share of
ALL reference samples for the site that carried it — is reported alongside. An
organism the reference never saw at that site is "Absent in healthy" (the PDF's
z = 3 fallback, i.e. not faded). Normally sterile sites (blood / plasma / WB,
CSF, ...) have no HMP distribution because nothing is expected there in health,
so every detection is reported as "Sterile site" (never faded).

The HMP table stores abundances in PERCENT (0-100), so observed abundance is
compared in percent as well.
"""
from __future__ import annotations

import csv
import gzip
import math
from bisect import bisect_right

csv.field_size_limit(1 << 30)

# Body sites that exist in the HMP reference table.
HMP_SITES = ("nasal", "oral", "skin", "stool", "throat", "vaginal")

# Sample-type -> HMP body site(s). Mirrors distributions.body_site_map (the
# mapping match_paths/create_report apply to the same table) plus the common
# clinical shorthands. Matched on the whole string first, then on its words
# ("NP swab", "whole_blood", "stool-2" ...). Falls back to
# body_site_normalization.normalize_body_site for anything else.
_SITE_ALIASES = {
    # stool / gut
    "gut": ["stool"], "feces": ["stool"], "fecal": ["stool"], "faeces": ["stool"],
    "faecal": ["stool"], "rectal": ["stool"], "rectum": ["stool"], "fece": ["stool"],
    "intestinal": ["stool"], "gi": ["stool"],
    # nasal / nasopharyngeal
    "nose": ["nasal"], "nasal": ["nasal"], "nares": ["nasal"], "naris": ["nasal"],
    "nostril": ["nasal"], "sinus": ["nasal"], "np": ["nasal"], "nps": ["nasal"],
    "naso": ["nasal"], "nasopharyngeal": ["nasal"], "nasopharynx": ["nasal"],
    "nasopharangyeal": ["nasal"], "nasopharangeal": ["nasal"], "nasophyngeal": ["nasal"],
    "mid turbinate": ["nasal"], "midturbinate": ["nasal"], "anterior nasal": ["nasal"], "resp": ["nasal"], "respiratory": ["nasal"],
    "rhinitis": ["nasal"],
    # vaginal
    "vagina": ["vaginal"], "vag": ["vaginal"], "cervical": ["vaginal"], "cervix": ["vaginal"],
    # oral / throat
    "teeth": ["oral"], "sputum": ["oral"], "mouth": ["oral"], "saliva": ["oral"],
    "buccal": ["oral"], "dental": ["oral"], "gingival": ["oral"], "tongue": ["oral"],
    "oropharyngeal": ["throat"], "op": ["throat"], "pharyngeal": ["throat"],
    "pharynx": ["throat"], "tonsil": ["throat"], "tonsillar": ["throat"],
    # skin
    "abscess": ["skin"], "absscess": ["skin"], "ear": ["skin"], "urogenital": ["skin"],
    "wound": ["skin"], "cutaneous": ["skin"], "dermal": ["skin"],
    # lower respiratory -> both upper-airway references
    "lung": ["oral", "nasal"], "bal": ["oral", "nasal"], "bronchoalveolar": ["oral", "nasal"],
}

# Normally sterile sites. HMP has no healthy distribution for these because a
# healthy subject carries (essentially) no microbes there: ANY detection is
# outside the healthy expectation. Value = label shown in the report.
_STERILE_ALIASES = {
    "blood": "blood", "wb": "blood", "whole blood": "blood", "wholeblood": "blood",
    "bld": "blood", "plasma": "blood", "serum": "blood", "peripheral blood": "blood", "bloodstream": "blood", "bacteremia": "blood",
    "sepsis": "blood", "septicemia": "blood", "blood culture": "blood",
    "dbs": "blood", "buffy coat": "blood", "pbmc": "blood", "cfdna": "blood",
    "csf": "csf", "cerebrospinal": "csf", "cerebrospinal fluid": "csf", "spinal fluid": "csf",
    "brain": "brain", "sterile": "sterile site", "clean": "sterile site",
    "peritoneal": "sterile site", "pleural fluid": "sterile site",
    "pericardial": "sterile site", "synovial": "sterile site", "synovial fluid": "sterile site",
    "joint fluid": "sterile site", "bone marrow": "sterile site", "tissue": "sterile site",
}
# normalize_body_site() categories that count as sterile
_STERILE_NORMALIZED = {"blood": "blood", "csf": "csf", "brain": "brain", "sterile": "sterile site"}

# Histogram of log10(% abundance) shipped to the report for the distribution
# plot: bins of 0.5 decades from 1e-6 % to 100 %.
HIST_LO, HIST_HI, HIST_STEP = -6.0, 2.0, 0.5
HIST_NBINS = int(round((HIST_HI - HIST_LO) / HIST_STEP))

# Report columns added to every detection row (display order).
HMP_COLUMNS = [
    "HMP Status",
    "HMP Z-Score",
    "HMP Healthy Percentile",
    "HMP Healthy Mean %",
    "HMP Prevalence %",
    "HMP Healthy Samples",
    "HMP Reference Samples",
    "HMP Body Site",
]
# Carried per row for the report's lookups but never shown as a table column.
HMP_HIDDEN_COLUMNS = {"HMP Ref Key"}

STATUS_ELEVATED = "Elevated"
STATUS_WITHIN = "Within healthy range"
STATUS_ABSENT = "Absent in healthy"
STATUS_STERILE = "Sterile site"

try:  # same normaliser create_report.py uses for pathogen/commensal site calls
    from body_site_normalization import normalize_body_site as _normalize_body_site
except Exception:  # pragma: no cover
    _normalize_body_site = None


def resolve_site(sample_type):
    """Classify a sample type / body site string.

    Returns (kind, value):
      ("hmp", [sites])      compare against these HMP healthy distribution(s)
      ("sterile", label)    normally sterile site: nothing is expected in health
      (None, "")            no healthy reference (environmental, unknown, ...)
    """
    s = str(sample_type or "").strip().lower()
    if not s or s in ("unknown", "none", "nan", "na", "n/a", "other"):
        return None, ""
    flat = " ".join(s.replace("_", " ").replace("-", " ").replace("/", " ").split())

    def _lookup(t):
        if t in HMP_SITES:
            return "hmp", [t]
        if t in _SITE_ALIASES:
            return "hmp", list(_SITE_ALIASES[t])
        if t in _STERILE_ALIASES:
            return "sterile", _STERILE_ALIASES[t]
        return None

    # whole string, then 2-word phrases, then single words ("NP swab", "whole_blood")
    for cand in (s, flat, flat.replace(" ", "")):
        hit = _lookup(cand)
        if hit:
            return hit
    toks = flat.split()
    for i in range(len(toks) - 1):
        hit = _lookup(toks[i] + " " + toks[i + 1])
        if hit:
            return hit
    for t in toks:
        hit = _lookup(t)
        if hit:
            return hit
    if _normalize_body_site is not None:
        n = _normalize_body_site(s)
        if n in _STERILE_NORMALIZED:
            return "sterile", _STERILE_NORMALIZED[n]
        hit = _lookup(n)
        if hit:
            return hit
    return None, ""


def hmp_sites_for(sample_type) -> list:
    """HMP site(s) a sample type compares to ([] for sterile / unmapped)."""
    kind, val = resolve_site(sample_type)
    return list(val) if kind == "hmp" else []


def _hist(abund):
    counts = [0] * HIST_NBINS
    for a in abund:
        if a <= 0:
            continue
        i = int(math.floor((math.log10(a) - HIST_LO) / HIST_STEP))
        counts[min(max(i, 0), HIST_NBINS - 1)] += 1
    return counts


def _quantiles(sorted_abund, qs=(0.05, 0.25, 0.5, 0.75, 0.95)):
    if not sorted_abund:
        return []
    n = len(sorted_abund)
    out = []
    for q in qs:
        pos = q * (n - 1)
        lo = int(math.floor(pos))
        hi = min(lo + 1, n - 1)
        out.append(sorted_abund[lo] + (sorted_abund[hi] - sorted_abund[lo]) * (pos - lo))
    return out


def load_hmp_reference(path, taxids, sites):
    """Read only the (taxid, site) rows the run needs from the HMP table.

    Returns (ref, site_totals):
      ref         {(site, taxid_str): {mean, std, n, N, abund(sorted), name, rank}}
      site_totals {site: number of reference samples for that site}
    """
    taxids = {str(t) for t in taxids if str(t).strip()}
    sites = set(sites)
    ref, site_totals = {}, {}
    if not path or not sites:
        return ref, site_totals
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        header = next(rdr)
        ix = {h: i for i, h in enumerate(header)}
        need = ("mean", "std", "tax_id", "body_site", "site_count", "abundances")
        missing = [c for c in need if c not in ix]
        if missing:
            raise ValueError(f"HMP table {path} is missing column(s): {', '.join(missing)}")
        i_tax, i_site, i_sc = ix["tax_id"], ix["body_site"], ix["site_count"]
        for row in rdr:
            if len(row) <= max(ix.values()):
                continue
            site = row[i_site].strip().lower()
            if site not in sites:
                continue
            try:
                sc = int(float(row[i_sc]))
            except ValueError:
                sc = 0
            if sc > site_totals.get(site, 0):
                site_totals[site] = sc
            tid = row[i_tax].strip()
            if tid.endswith(".0"):
                tid = tid[:-2]
            if tid not in taxids:
                continue
            try:
                abund = sorted(float(x) for x in row[ix["abundances"]].split(",") if x.strip())
                mean = float(row[ix["mean"]])
                std = float(row[ix["std"]])
            except ValueError:
                continue
            ref[(site, tid)] = {
                "mean": mean, "std": std, "n": len(abund), "N": sc,
                "abund": abund,
                "name": row[ix["name"]] if "name" in ix else "",
                "rank": row[ix["rank"]] if "rank" in ix else "",
            }
    return ref, site_totals


def compute_hmp_fields(observed_pct, sites, ref_taxid, ref, site_totals, threshold):
    """HMP outlier fields for one detection.

    observed_pct : abundance of the detection in percent of the sample's reads
    sites        : HMP sites the sample type maps to (may be empty)
    ref_taxid    : taxid to look up (already resolved for strain->species fallback)
    """
    blank = {c: None for c in HMP_COLUMNS}
    blank["HMP Status"] = ""
    blank["HMP Body Site"] = "+".join(sites)
    blank["HMP Ref Key"] = ""
    if not sites:
        return blank
    N_total = sum(site_totals.get(s, 0) for s in sites)
    if N_total <= 0:
        return blank
    entries = [(s, ref[(s, ref_taxid)]) for s in sites if (s, ref_taxid) in ref]
    obs = max(0.0, float(observed_pct or 0))
    out = dict(blank)
    out["HMP Reference Samples"] = N_total
    if not entries:
        out.update({
            "HMP Status": STATUS_ABSENT,
            "HMP Healthy Samples": 0,
            "HMP Prevalence %": 0.0,
            "HMP Healthy Percentile": 100.0 if obs > 0 else 0.0,
        })
        return out
    mean = sum(e["mean"] for _, e in entries)
    std = sum(e["std"] for _, e in entries)
    n = sum(e["n"] for _, e in entries)
    if std > 0:
        z = (obs - mean) / std
    else:
        z = 3.0 if obs > mean else 0.0
    # empirical percentile across ALL reference samples (absent samples = 0 %)
    le = sum((site_totals.get(s, 0) - e["n"]) + bisect_right(e["abund"], obs) for s, e in entries)
    pct = 100.0 * le / max(1, sum(site_totals.get(s, 0) for s, _ in entries))
    out.update({
        "HMP Status": STATUS_ELEVATED if z >= threshold else STATUS_WITHIN,
        "HMP Z-Score": round(z, 3),
        "HMP Healthy Percentile": round(min(100.0, pct), 2),
        "HMP Healthy Mean %": round(mean, 6),
        "HMP Prevalence %": round(100.0 * n / N_total, 3),
        "HMP Healthy Samples": n,
        "HMP Ref Key": "|".join(f"{s}:{ref_taxid}" for s, _ in entries),
    })
    return out


def reference_payload(ref, site_totals, used_keys, threshold, source):
    """Compact reference block shipped to the report (only the keys rows use)."""
    out = {}
    for key in sorted(used_keys):
        site, tid = key.split(":", 1)
        e = ref.get((site, tid))
        if not e:
            continue
        out[key] = {
            "m": e["mean"], "s": e["std"], "n": e["n"], "N": site_totals.get(site, e["N"]),
            "h": _hist(e["abund"]),
            "q": [round(v, 8) for v in _quantiles(e["abund"])],
            "name": e["name"], "rank": e["rank"],
        }
    return {
        "threshold": threshold,
        "source": source,
        "sites": site_totals,
        "bins": {"lo": HIST_LO, "hi": HIST_HI, "step": HIST_STEP},
        "ref": out,
    }


def annotate_rows(rows, hmp_path, threshold=2.0):
    """Stamp HMP outlier columns onto report rows in place.

    Rows are make_report.py flat records. A row whose own taxid has no healthy
    reference falls back to its species (Subkey), compared against the
    species' abundance in that sample. Returns the BOOT.hmp payload, or None
    when nothing could be compared (no table, or no sample maps to an HMP site).
    """
    sites_needed, taxids = set(), set()
    any_sterile = False
    for r in rows:
        kind, val = resolve_site(r.get("Sample Type"))
        if kind == "hmp":
            sites_needed.update(val)
            taxids.add(str(r.get("Taxonomic ID #", "")).strip())
            taxids.add(str(r.get("Subkey", "")).strip())
        elif kind == "sterile":
            any_sterile = True
    if not sites_needed and not any_sterile:
        for r in rows:
            r.update(compute_hmp_fields(0, [], "", {}, {}, threshold))
        return None
    ref, site_totals = (load_hmp_reference(hmp_path, taxids, sites_needed)
                        if sites_needed else ({}, {}))

    # species-level abundance per (sample, taxid) for the strain -> species fallback
    pct_by = {}
    for r in rows:
        k = (r.get("Specimen ID"), str(r.get("Taxonomic ID #", "")).strip())
        try:
            v = float(r.get("% Reads") or 0)
        except (TypeError, ValueError):
            v = 0.0
        if v > pct_by.get(k, -1):
            pct_by[k] = v

    used = set()
    for r in rows:
        kind, val = resolve_site(r.get("Sample Type"))
        if kind == "sterile":
            f = compute_hmp_fields(0, [], "", {}, {}, threshold)
            f["HMP Status"] = STATUS_STERILE
            f["HMP Body Site"] = val
            f["HMP Healthy Samples"] = 0
            f["HMP Prevalence %"] = 0.0
            r.update(f)
            continue
        sites = list(val) if kind == "hmp" else []
        tid = str(r.get("Taxonomic ID #", "")).strip()
        sub = str(r.get("Subkey", "")).strip()
        try:
            obs = float(r.get("% Reads") or 0)
        except (TypeError, ValueError):
            obs = 0.0
        look = tid
        if sites and sub and sub != tid and not any((s, tid) in ref for s in sites) \
                and any((s, sub) in ref for s in sites):
            look = sub
            obs = pct_by.get((r.get("Specimen ID"), sub), obs)
        f = compute_hmp_fields(obs, sites, look, ref, site_totals, threshold)
        r.update(f)
        if f["HMP Ref Key"]:
            used.update(f["HMP Ref Key"].split("|"))
    return reference_payload(ref, site_totals, used, threshold, str(hmp_path))


def annotate_rows_from_json(rows, hmp_src, threshold=2.0):
    """Fallback when no HMP table is given: reuse the per-organism HMP fields
    match_paths.py wrote into the paths JSON (only present when it ran with
    --hmp). ``hmp_src`` maps id(row) -> the organism's raw JSON fields."""
    any_ref = False
    for r in rows:
        src = hmp_src.get(id(r)) or {}
        N = int(src.get("hmp_site_count") or 0)
        f = compute_hmp_fields(0, [], "", {}, {}, threshold)
        if N > 0:
            any_ref = True
            n = int(src.get("hmp_num_samples") or 0)
            z = float(src.get("zscore") or 0)
            f.update({
                "HMP Status": STATUS_ELEVATED if z >= threshold else STATUS_WITHIN,
                "HMP Z-Score": round(z, 3),
                "HMP Healthy Mean %": src.get("hmp_mean"),
                "HMP Prevalence %": round(100.0 * n / N, 3),
                "HMP Healthy Samples": n,
                "HMP Reference Samples": N,
                "HMP Body Site": str(src.get("normalized_sample_site") or ""),
            })
        r.update(f)
    if not any_ref:
        return None
    return {"threshold": threshold, "source": "paths.json (match_paths --hmp)",
            "sites": {}, "bins": {"lo": HIST_LO, "hi": HIST_HI, "step": HIST_STEP}, "ref": {}}
