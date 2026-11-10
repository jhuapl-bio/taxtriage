#!/usr/bin/env python3
# ##############################################################################################
# # Copyright 2025 The Johns Hopkins University Applied Physics Laboratory LLC
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
index_hopping.py -- flag detections that Illumina index hopping could explain.

On patterned Illumina flow cells a small fraction of reads swap an index with
another library in the same lane, so a sample can pick up a few reads of an
organism that is abundant in a neighbour. This script reads the demultiplexing
index-hopping table (bcl-convert ``Index_Hopping_Counts.csv`` or the same table
exported as TSV / XLSX), groups the run's samples by lane, and for every
organism in every sample asks: could the reads seen here have hopped in from
the other samples of the lane that carry the same organism?

    expected = sum over donors D of  reads_D(organism) x rate(D, R)
    p        = P(X >= observed | X ~ Poisson(factor x expected))
    suspect  = expected > 0 and p >= alpha     (hopping cannot be ruled out)

rate(D, R), first that applies, per shared lane (largest over shared lanes):
    pair     the file lists hopped index combinations of D and R:
             (hopped reads i7_D+i5_R + i7_R+i5_D + 1) / (reads_D + reads_R)
    lane     the lane has hopped rows but none for this pair: the lane-wide
             rate, hopped / (hopped + correctly paired)  -- an upper bound
    default  the file has no hopped rows (or no index2): --rate

Flag only: TASS, calls and every existing field are left exactly as they are.
Each organism node (genus, species and strain) gains an ``index_hop`` block and
each sample's metadata an ``index_hopping`` block; samples outside the file get
``metadata.index_hopping.applied = false``.

Stdlib only (the XLSX reader is a small zipfile/XML parser), so it runs in any
of the pipeline's Python containers.
"""

import argparse
import csv
import gzip
import io
import json
import math
import os
import re
import sys
import zipfile
import xml.etree.ElementTree as ET
from collections import defaultdict

ALL_LANES = "all"
# Donor-side notes only for a donor expected to have sent at least this many
# reads into a flagged detection (keeps negligible shares off the report).
# Same value in assets/src/js/54_index_hopping.js.
MIN_SOURCE_READS = 0.5

# ──────────────────────────────────────────────────────────────────────────────
# Table reading (CSV / TSV / XLSX)
# ──────────────────────────────────────────────────────────────────────────────
_NS_MAIN = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
_NS_REL = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
_NS_PKG = "http://schemas.openxmlformats.org/package/2006/relationships"


def _col_index(ref):
    """'C12' -> 2 (0-based column)."""
    letters = re.match(r"[A-Za-z]+", ref or "")
    if not letters:
        return None
    n = 0
    for ch in letters.group(0).upper():
        n = n * 26 + (ord(ch) - 64)
    return n - 1


def read_xlsx_sheets(path):
    """Return [(sheet_name, rows)] with rows as lists of strings (stdlib only)."""
    out = []
    with zipfile.ZipFile(path) as z:
        names = set(z.namelist())
        shared = []
        if "xl/sharedStrings.xml" in names:
            root = ET.fromstring(z.read("xl/sharedStrings.xml"))
            for si in root.findall(f"{{{_NS_MAIN}}}si"):
                shared.append("".join(t.text or "" for t in si.iter(f"{{{_NS_MAIN}}}t")))
        wb = ET.fromstring(z.read("xl/workbook.xml"))
        relmap = {}
        if "xl/_rels/workbook.xml.rels" in names:
            rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
            for r in rels.findall(f"{{{_NS_PKG}}}Relationship"):
                relmap[r.get("Id")] = r.get("Target", "")
        sheets_el = wb.find(f"{{{_NS_MAIN}}}sheets")
        for i, s in enumerate(list(sheets_el) if sheets_el is not None else []):
            target = relmap.get(s.get(f"{{{_NS_REL}}}id"), f"worksheets/sheet{i + 1}.xml")
            target = target.lstrip("/")
            if not target.startswith("xl/"):
                target = "xl/" + target
            if target not in names:
                continue
            root = ET.fromstring(z.read(target))
            rows = []
            for row in root.iter(f"{{{_NS_MAIN}}}row"):
                cells = {}
                for c in row.findall(f"{{{_NS_MAIN}}}c"):
                    idx = _col_index(c.get("r")) if c.get("r") else len(cells)
                    t = c.get("t")
                    v = c.find(f"{{{_NS_MAIN}}}v")
                    if t == "s" and v is not None:
                        val = shared[int(v.text)] if v.text and int(v.text) < len(shared) else ""
                    elif t == "inlineStr":
                        val = "".join(x.text or "" for x in c.iter(f"{{{_NS_MAIN}}}t"))
                    else:
                        val = v.text if v is not None and v.text is not None else ""
                    cells[idx] = val
                if cells:
                    width = max(cells) + 1
                    rows.append([cells.get(j, "") for j in range(width)])
                else:
                    rows.append([])
            out.append((s.get("name", f"sheet{i + 1}"), rows))
    return out


def _open_text(path):
    if path.lower().endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8-sig", newline="")
    return open(path, encoding="utf-8-sig", newline="")


def read_delimited(path):
    """CSV / TSV (optionally .gz). The delimiter is taken from the header-ish
    lines: tab when they hold tabs, otherwise comma (then ';')."""
    with _open_text(path) as fh:
        text = fh.read()
    lines = [ln for ln in text.splitlines() if ln.strip()][:20]
    tabs = sum(ln.count("\t") for ln in lines)
    commas = sum(ln.count(",") for ln in lines)
    semis = sum(ln.count(";") for ln in lines)
    lower = path.lower()
    if tabs and (tabs >= commas or lower.endswith((".tsv", ".tsv.gz", ".txt", ".txt.gz"))):
        delim = "\t"
    elif commas:
        delim = ","
    elif semis:
        delim = ";"
    else:
        delim = "\t"
    return list(csv.reader(io.StringIO(text), delimiter=delim))


def read_table_rows(path):
    """All candidate row lists in a file (one per XLSX sheet, one otherwise)."""
    with open(path, "rb") as fh:
        magic = fh.read(4)
    if magic == b"PK\x03\x04" or path.lower().endswith((".xlsx", ".xlsm")):
        return [rows for _, rows in read_xlsx_sheets(path)]
    if path.lower().endswith(".xls"):
        raise ValueError(f"{path}: legacy .xls is not supported; save it as .xlsx, .csv or .tsv")
    return [read_delimited(path)]


# ──────────────────────────────────────────────────────────────────────────────
# Header recognition + parsing
# ──────────────────────────────────────────────────────────────────────────────
def _norm_header(h):
    h = str(h or "").strip().lower().replace("#", " ").replace("%", " pct ")
    return re.sub(r"[^a-z0-9]+", "_", h).strip("_")


_HEADER_ALIASES = {
    "lane": {"lane"},
    "sample": {"sampleid", "sample_id", "sample", "sample_name", "samplename", "sample_id_s"},
    "project": {"sample_project", "project", "sampleproject"},
    "i7": {"index", "index1", "index_1", "i7", "i7_index", "index_i7", "i7_index_id_sequence"},
    "i5": {"index2", "index_2", "i5", "i5_index", "index_i5"},
    "reads": {"reads", "read_count", "readcount", "num_reads", "clusters", "count",
              "number_of_reads", "reads_count"},
    "pct_hopped": {"pct_of_hopped_reads", "pct_hopped_reads", "pct_hopped", "of_hopped_reads"},
    "pct_all": {"pct_of_all_reads", "pct_all_reads", "of_all_reads"},
}


def _match_header(row):
    """{field: column index} when `row` looks like the table header, else None."""
    cols = {}
    for i, raw in enumerate(row):
        h = _norm_header(raw)
        for field, names in _HEADER_ALIASES.items():
            if h in names and field not in cols:
                cols[field] = i
    if "reads" in cols and ("i7" in cols or "sample" in cols):
        return cols
    return None


def _num(v):
    s = str(v or "").strip().replace(",", "").replace(" ", "")
    if not s:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def _seq(v):
    return re.sub(r"[^ACGTN]", "", str(v or "").upper())


def parse_hopping_file(path):
    """Read one index-hopping table.

    Returns {"path", "lanes": {lane: {"samples": {sid: {"i7","i5","reads","project"}},
                                      "hopped": [(i7, i5, reads)]}},
             "has_index2": bool}.

    A row with a SampleID whose index pair is that sample's own (the pair with
    the most reads for the sample in that lane) is a correctly paired count;
    every other row (no SampleID, or a SampleID with another pair) is a hopped
    combination.
    """
    tables = read_table_rows(path)
    header, rows = None, None
    for t in tables:
        for k, r in enumerate(t):
            cols = _match_header(r)
            if cols:
                header, rows = cols, t[k + 1:]
                break
        if header:
            break
    if not header:
        raise ValueError(
            f"{path}: no index-hopping header found (need a '# Reads' column plus "
            f"'index' or 'SampleID'; see docs/index-hopping.md)")

    def cell(r, field):
        i = header.get(field)
        return r[i].strip() if i is not None and i < len(r) and r[i] is not None else ""

    records = []
    for r in rows:
        if not any(str(x).strip() for x in r):
            continue
        reads = _num(cell(r, "reads"))
        if reads is None:
            # a second table / footer after the data: stop at the first row
            # that is not a data row
            if _match_header(r):
                continue
            continue
        records.append({
            "lane": cell(r, "lane") or ALL_LANES,
            "sample": cell(r, "sample"),
            "project": cell(r, "project"),
            "i7": _seq(cell(r, "i7")),
            "i5": _seq(cell(r, "i5")),
            "reads": reads,
        })

    lanes = {}
    # pass 1: each sample's own pair per lane = its row with the most reads
    own = {}
    for rec in records:
        if not rec["sample"]:
            continue
        k = (rec["lane"], rec["sample"])
        if k not in own or rec["reads"] > own[k]["reads"]:
            own[k] = rec
    for (lane, sid), rec in own.items():
        L = lanes.setdefault(lane, {"samples": {}, "hopped": []})
        L["samples"][sid] = {"i7": rec["i7"], "i5": rec["i5"], "reads": rec["reads"],
                             "project": rec["project"]}
    # pass 2: everything else is a hopped combination
    for rec in records:
        if rec["sample"] and own.get((rec["lane"], rec["sample"])) is rec:
            continue
        L = lanes.setdefault(rec["lane"], {"samples": {}, "hopped": []})
        L["hopped"].append((rec["i7"], rec["i5"], rec["reads"]))
    has_index2 = any(s["i5"] for L in lanes.values() for s in L["samples"].values())
    return {"path": path, "lanes": lanes, "has_index2": has_index2}


# ──────────────────────────────────────────────────────────────────────────────
# Hop rates
# ──────────────────────────────────────────────────────────────────────────────
class LaneRates:
    """Hop-rate lookup for one lane of one file."""

    def __init__(self, lane, data):
        self.lane = lane
        self.samples = data["samples"]
        by_i7, by_i5 = defaultdict(list), defaultdict(list)
        for sid, s in self.samples.items():
            if s["i7"]:
                by_i7[s["i7"]].append(sid)
            if s["i5"]:
                by_i5[s["i5"]].append(sid)
        self.pair_hops = defaultdict(float)  # frozenset({A, B}) -> hopped reads
        self.hopped_total = 0.0               # hops between two known samples
        self.unattributed = 0.0               # hopped rows whose indexes match no sample
        self.n_hopped_rows = len(data["hopped"])
        for i7, i5, reads in data["hopped"]:
            src, dst = by_i7.get(i7, []), by_i5.get(i5, [])
            combos = [(a, b) for a in src for b in dst if a != b]
            if not combos:
                self.unattributed += reads
                continue
            self.hopped_total += reads
            share = reads / len(combos)
            for a, b in combos:
                self.pair_hops[frozenset((a, b))] += share
        self.proper_total = sum(s["reads"] for s in self.samples.values())
        denom = self.hopped_total + self.proper_total
        self.lane_rate = (self.hopped_total / denom) if (self.n_hopped_rows and denom > 0) else None

    def rate(self, a, b, default_rate):
        """(rate, source) for reads of sample a landing in sample b."""
        if self.lane_rate is None:
            return default_rate, "default"
        key = frozenset((a, b))
        if key in self.pair_hops:
            na = self.samples[a]["reads"]
            nb = self.samples[b]["reads"]
            if na + nb > 0:
                return (self.pair_hops[key] + 1.0) / (na + nb), "pair"
        return self.lane_rate, "lane"


# ──────────────────────────────────────────────────────────────────────────────
# Sample matching
# ──────────────────────────────────────────────────────────────────────────────
def norm_id(s):
    return re.sub(r"[^a-z0-9]+", "_", str(s or "").strip().lower()).strip("_")


_BCL_SUFFIX = re.compile(r"_s\d+(_l\d{3})?(_r[12](_001)?)?$")


def id_candidates(sample_name, override=None):
    """Normalised names a TaxTriage sample may appear under in the hop file."""
    if override:
        return [norm_id(override)]
    base = norm_id(sample_name)
    cands = [base]
    stripped = _BCL_SUFFIX.sub("", base)
    if stripped and stripped != base:
        cands.append(stripped)
    return cands


# ──────────────────────────────────────────────────────────────────────────────
# Organism tree helpers
# ──────────────────────────────────────────────────────────────────────────────
def iter_nodes(organisms):
    """Yield (depth, node, is_leaf) for every genus / species / strain node."""
    stack = [(0, n) for n in reversed(organisms or [])]
    while stack:
        depth, node = stack.pop()
        if not isinstance(node, dict):
            continue
        members = node.get("members") or []
        yield depth, node, not members
        for m in reversed(members):
            stack.append((depth + 1, m))


def node_key(node):
    return str(node.get("key", node.get("subkey", "")) or "")


def node_reads(node):
    """Whole reads, truncated like make_report's "# Reads Aligned"."""
    try:
        return int(float(node.get("numreads", 0) or 0))
    except (TypeError, ValueError):
        return 0


def level_name(depth, is_leaf):
    if is_leaf:
        return "strain"
    return {0: "genus", 1: "species"}.get(depth, "strain")


def poisson_sf(k, lam):
    """P(X >= k) for X ~ Poisson(lam)."""
    if k <= 0:
        return 1.0
    if lam <= 0:
        return 0.0
    # far in the upper tail: effectively zero (avoids long loops)
    if k > lam + 12.0 * math.sqrt(lam) + 30.0:
        return 0.0
    term = math.exp(-lam)
    cdf = term
    for i in range(1, int(k)):
        term *= lam / i
        cdf += term
    return max(0.0, min(1.0, 1.0 - cdf))


# ──────────────────────────────────────────────────────────────────────────────
# Main analysis
# ──────────────────────────────────────────────────────────────────────────────
def analyse(samples, hop_files, sample_files, overrides, default_rate, factor, alpha,
            max_donors=5, log=print):
    """Annotate samples in place.

    samples       {sample_id: parsed paths.json}
    hop_files     [parsed hopping files]
    sample_files  {sample_id: [file indexes]} (None -> every file for every sample)
    overrides     {sample_id: SampleID in the file}
    Returns summary rows (one per suspect node).
    """
    # 1. resolve each TaxTriage sample to a SampleID in each of its files
    lookup = []  # per file: {normalised SampleID: SampleID}
    for hf in hop_files:
        ids = {}
        for L in hf["lanes"].values():
            for sid in L["samples"]:
                ids.setdefault(norm_id(sid), sid)
        lookup.append(ids)

    membership = defaultdict(dict)  # sample -> {file index: SampleID}
    for sname, data in samples.items():
        files = sample_files.get(sname) if sample_files is not None else None
        if files is None:
            files = list(range(len(hop_files)))
        meta_name = (data.get("metadata") or {}).get("sample_name") or sname
        cands = id_candidates(sname, overrides.get(sname))
        if not overrides.get(sname):
            for c in id_candidates(meta_name):
                if c not in cands:
                    cands.append(c)
        for fi in files:
            for c in cands:
                if c in lookup[fi]:
                    membership[sname][fi] = lookup[fi][c]
                    break
        if files and not membership[sname]:
            want = overrides.get(sname) or sname
            log(f"[index_hopping] WARNING: sample {sname!r} not found (as {want!r}) in "
                f"{', '.join(os.path.basename(hop_files[i]['path']) for i in files)}")

    lane_rates = {}
    for fi, hf in enumerate(hop_files):
        for lane, L in hf["lanes"].items():
            lane_rates[(fi, lane)] = LaneRates(lane, L)

    # 2. donors: samples that share a lane of the same file
    def shared_lanes(a, b):
        out = []
        for fi, sid_a in membership[a].items():
            sid_b = membership[b].get(fi)
            if sid_b is None:
                continue
            for lane, L in hop_files[fi]["lanes"].items():
                if sid_a in L["samples"] and sid_b in L["samples"]:
                    out.append((fi, lane, sid_a, sid_b))
        return out

    names = sorted(samples)
    pair_rate = {}  # (donor, recipient) -> (rate, source, lane)
    for d in names:
        for r in names:
            if d == r or not membership[d] or not membership[r]:
                continue
            best = None
            for fi, lane, sid_d, sid_r in shared_lanes(d, r):
                rate, src = lane_rates[(fi, lane)].rate(sid_d, sid_r, default_rate)
                if best is None or rate > best[0]:
                    best = (rate, src, lane)
            if best:
                pair_rate[(d, r)] = best

    # 3. per-sample organism read counts at (level, key). Whole reads, as the
    #    report rows carry them ("# Reads Aligned"), so the in-browser version
    #    (assets/src/js/54_index_hopping.js) gives the same answer.
    reads_at = {}
    node_at = {}  # the node each reads_at entry came from (donor-side notes)
    for s, data in samples.items():
        idx, nodes = {}, {}
        for depth, node, leaf in iter_nodes(data.get("organisms")):
            key = node_key(node)
            if not key:
                continue
            k = (level_name(depth, leaf), key)
            if k not in idx or node_reads(node) > idx[k]:
                idx[k] = node_reads(node)
                nodes[k] = node
        reads_at[s] = idx
        node_at[s] = nodes
    # donor side: (donor, level, key) -> [{sample, expected, observed}] for every
    # flagged detection elsewhere that this organism's reads could explain
    sent = defaultdict(list)

    summary = []
    for r in names:
        data = samples[r]
        meta = data.setdefault("metadata", {})
        if not membership[r]:
            files = sample_files.get(r) if sample_files is not None else None
            files = list(range(len(hop_files))) if files is None else files
            meta["index_hopping"] = {
                "applied": False,
                "reason": ("sample not found in the index hopping file(s)" if files
                           else "no index hopping file for this sample"),
                "files": [os.path.basename(hop_files[i]["path"]) for i in files],
            }
            continue
        donors = sorted({d for (d, rr) in pair_rate if rr == r})
        counts = defaultdict(int)
        n_nodes = defaultdict(int)
        for depth, node, leaf in iter_nodes(data.get("organisms")):
            key = node_key(node)
            obs = node_reads(node)
            lvl = level_name(depth, leaf)
            n_nodes[lvl] += 1
            contrib = []
            for d in donors:
                nd = reads_at[d].get((lvl, key), 0)
                if nd <= obs or nd <= 0:
                    continue  # a donor with no more reads cannot explain these
                rate, src, lane = pair_rate[(d, r)]
                contrib.append({"sample": d, "reads": nd, "rate": rate,
                                "rate_source": src, "lane": lane, "expected": nd * rate})
            expected = sum(c["expected"] for c in contrib)
            obs_i = obs
            p = poisson_sf(obs_i, factor * expected) if expected > 0 else 0.0
            suspect = bool(expected > 0 and obs_i > 0 and p >= alpha)
            contrib.sort(key=lambda c: -c["expected"])
            block = {
                "suspect": suspect,
                "observed_reads": obs_i,
                "expected_reads": round(expected, 4),
                "p_value": round(p, 6),
                "explained_fraction": round(min(1.0, expected / obs), 4) if obs > 0 else 0.0,
                "top_donor": contrib[0]["sample"] if contrib else None,
                "top_donor_reads": contrib[0]["reads"] if contrib else None,
                "top_donor_rate": contrib[0]["rate"] if contrib else None,
                "rate_source": contrib[0]["rate_source"] if contrib else None,
                "donors": [{k: (round(v, 8) if isinstance(v, float) else v) for k, v in c.items()}
                           for c in contrib[:max_donors]],
            }
            node["index_hop"] = block
            if suspect:
                counts[lvl] += 1
                for c in contrib:
                    if c["expected"] < MIN_SOURCE_READS:
                        continue
                    sent[(c["sample"], lvl, key)].append(
                        {"sample": r, "expected": c["expected"], "observed": obs_i})
                summary.append({
                    "sample": r, "level": lvl, "key": key, "name": node.get("name", ""),
                    "observed_reads": obs_i, "expected_reads": block["expected_reads"],
                    "p_value": block["p_value"], "explained_fraction": block["explained_fraction"],
                    "top_donor": block["top_donor"], "top_donor_reads": block["top_donor_reads"],
                    "top_donor_rate": block["top_donor_rate"], "rate_source": block["rate_source"],
                    "tass_score": node.get("tass_score"),
                    "passes_threshold": node.get("passes_threshold"),
                })
        file_info = []
        for fi, sid in sorted(membership[r].items()):
            hf = hop_files[fi]
            lanes = sorted(l for l, L in hf["lanes"].items() if sid in L["samples"])
            file_info.append({
                "file": os.path.basename(hf["path"]),
                "sample_id": sid,
                "lanes": lanes,
                "lane_rates": {l: lane_rates[(fi, l)].lane_rate for l in lanes},
                "reads_in_file": sum(hf["lanes"][l]["samples"][sid]["reads"] for l in lanes),
            })
        meta["index_hopping"] = {
            "applied": True,
            "files": file_info,
            "donor_samples": donors,
            "default_rate": default_rate,
            "factor": factor,
            "alpha": alpha,
            "rate_sources": sorted({pair_rate[(d, r)][1] for d in donors}),
            "n_suspect": dict(counts),
            "n_nodes": dict(n_nodes),
        }
    # Donor-side note: a light annotation on the organism the hopped reads
    # likely came FROM. Informational only (never a flag).
    n_source = defaultdict(lambda: defaultdict(int))
    for (d, lvl, key), recips in sent.items():
        node = node_at.get(d, {}).get((lvl, key))
        if node is None:
            continue
        recips = sorted(recips, key=lambda x: -x["expected"])
        node["index_hop_source"] = {
            "n_recipients": len(recips),
            "expected_out": round(sum(x["expected"] for x in recips), 4),
            "recipients": [{"sample": x["sample"], "expected_reads": round(x["expected"], 4),
                            "observed_reads": x["observed"]} for x in recips[:max_donors]],
        }
        n_source[d][lvl] += 1
    for d, by in n_source.items():
        ih = samples[d].get("metadata", {}).get("index_hopping")
        if isinstance(ih, dict) and ih.get("applied"):
            ih["n_source"] = dict(by)
    return summary


# ──────────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────────
SUMMARY_COLS = ["sample", "level", "key", "name", "observed_reads", "expected_reads",
                "p_value", "explained_fraction", "top_donor", "top_donor_reads",
                "top_donor_rate", "rate_source", "tass_score", "passes_threshold"]


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--jsons", nargs="+", required=True,
                   help="per-sample *.paths.json files (match_paths.py output)")
    p.add_argument("--hopping", nargs="+", required=True,
                   help="index hopping table(s): CSV, TSV or XLSX (bcl-convert Index_Hopping_Counts layout)")
    p.add_argument("--sample-map", default=None,
                   help="JSON {\"samples\": {sample_id: {\"files\": [index into --hopping], "
                        "\"index_sample_id\": str|null}}}. Without it every file applies to every sample.")
    p.add_argument("--rate", type=float, default=0.001,
                   help="hop rate used when a file has no hopped rows (default 0.001 = 0.1%%)")
    p.add_argument("--factor", type=float, default=2.0,
                   help="multiplier on the expected hopped reads, for rate uncertainty (default 2)")
    p.add_argument("--alpha", type=float, default=0.01,
                   help="flag when P(reads >= observed | hopping) >= alpha (default 0.01)")
    p.add_argument("-o", "--outdir", default="annotated", help="annotated JSONs are written here")
    p.add_argument("--summary", default="index_hopping.tsv", help="TSV of suspect detections")
    p.add_argument("--report", default="index_hopping.json", help="JSON run summary")
    a = p.parse_args(argv)
    if not 0 <= a.rate <= 1:
        p.error("--rate must be between 0 and 1")
    if a.factor <= 0:
        p.error("--factor must be > 0")
    if not 0 < a.alpha <= 1:
        p.error("--alpha must be in (0, 1]")
    return a


def main(argv=None):
    a = parse_args(argv)
    hop_files = []
    for path in a.hopping:
        hf = parse_hopping_file(path)
        n_s = len({s for L in hf["lanes"].values() for s in L["samples"]})
        n_h = sum(len(L["hopped"]) for L in hf["lanes"].values())
        print(f"[index_hopping] {os.path.basename(path)}: {len(hf['lanes'])} lane(s), "
              f"{n_s} sample(s), {n_h} hopped combination row(s)"
              + ("" if hf["has_index2"] else "; no index2 (default rate used)"))
        hop_files.append(hf)

    sample_files, overrides = None, {}
    if a.sample_map:
        with open(a.sample_map) as fh:
            smap = json.load(fh).get("samples", {})
        sample_files = {k: [int(i) for i in (v.get("files") or [])] for k, v in smap.items()}
        overrides = {k: v.get("index_sample_id") for k, v in smap.items() if v.get("index_sample_id")}

    samples, paths, passthrough = {}, {}, []
    for path in a.jsons:
        with open(path) as fh:
            data = json.load(fh)
        if not isinstance(data, dict) or "organisms" not in data:
            passthrough.append(path)
            continue
        base = os.path.basename(path)
        sname = re.sub(r"\.paths\.json$|\.json$", "", base)
        if sample_files is not None and sname not in sample_files:
            # map keyed by meta.id; the file stem normally equals it
            alt = (data.get("metadata") or {}).get("sample_name")
            if alt in sample_files:
                sname = alt
        samples[sname] = data
        paths[sname] = base

    if sample_files is not None:
        for s in samples:
            sample_files.setdefault(s, [])

    summary = analyse(samples, hop_files, sample_files, overrides, a.rate, a.factor, a.alpha)

    os.makedirs(a.outdir, exist_ok=True)
    for s, data in samples.items():
        with open(os.path.join(a.outdir, paths[s]), "w") as fh:
            json.dump(data, fh, indent=4)
    for path in passthrough:
        with open(path) as src, open(os.path.join(a.outdir, os.path.basename(path)), "w") as dst:
            dst.write(src.read())

    with open(a.summary, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=SUMMARY_COLS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for row in sorted(summary, key=lambda r: (r["sample"], r["level"], -r["explained_fraction"])):
            w.writerow(row)

    report = {
        "files": [os.path.basename(h["path"]) for h in hop_files],
        "default_rate": a.rate, "factor": a.factor, "alpha": a.alpha,
        "samples": {s: d["metadata"].get("index_hopping") for s, d in sorted(samples.items())},
    }
    with open(a.report, "w") as fh:
        json.dump(report, fh, indent=2)

    n_app = sum(1 for d in samples.values() if d["metadata"]["index_hopping"].get("applied"))
    n_strain = sum(1 for r in summary if r["level"] == "strain")
    print(f"[index_hopping] {n_app}/{len(samples)} sample(s) matched; "
          f"{n_strain} strain-level detection(s) possibly explained by index hopping")
    return 0


if __name__ == "__main__":
    sys.exit(main())
