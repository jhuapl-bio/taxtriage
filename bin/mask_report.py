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
mask_report.py -- visual summary of a masked (FCS / targeted) TaxTriage run.

Reads the per-sample *.paths.json files written by match_paths.py with
--mask_bed, and optionally the matching UNMASKED JSONs (--baseline) to show how
each organism's TASS moved. Writes:

  <prefix>.html   self-contained page (no external libraries):
                    run settings + per-sample read accounting
                    organism table (TASS masked vs unmasked, reads / bp masked,
                    penalty, NCBI FCS assembly flag), sortable and filterable
                    coverage tracks for the most-affected organisms, with the
                    masked (exclude) or targeted (target) spans drawn over the
                    positional breadth histogram
                    NCBI FCS assembly table (from build_mask.py --assembly_out)
  <prefix>.tsv    the organism table

Example
  mask_report.py --jsons masked/*.paths.json --baseline unmasked/*.paths.json \\
      --assembly_table fcs_assemblies.tsv --mask_summary fcs_mask.summary.tsv \\
      --prefix all.masking
"""
from __future__ import annotations

import argparse
import base64
import csv
import html
import json
import os
import sys


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--jsons", nargs="+", required=True, help="Masked per-sample *.paths.json")
    p.add_argument("--baseline", nargs="*", default=[], help="Unmasked per-sample *.paths.json (optional)")
    p.add_argument("--assembly_table", default=None, help="build_mask.py --assembly_out TSV")
    p.add_argument("--mask_summary", default=None, help="build_mask.py --summary_out TSV")
    p.add_argument("--mask_json", default=None, help="build_mask.py --json_out")
    p.add_argument("--top", type=int, default=24, help="Coverage tracks to draw (most-affected first)")
    p.add_argument("--min_tass", type=float, default=0.0, help="Hide organisms below this TASS in both runs")
    p.add_argument("--prefix", default="all.masking", help="Output prefix")
    return p.parse_args(argv)


def _f(x, d=0.0):
    try:
        return float(x)
    except (TypeError, ValueError):
        return d


def _sample_name(d, path):
    return (d.get("metadata") or {}).get("sample_name") or os.path.basename(path).split(".")[0]


def _iter_strains(d):
    for g in d.get("organisms", []) or []:
        for sk in g.get("members", []) or []:
            for st in sk.get("members", []) or []:
                yield g, sk, st


def load(paths):
    out = {}
    for p in paths:
        if not p or os.path.basename(p).startswith("NO_FILE"):
            continue
        try:
            with open(p) as fh:
                d = json.load(fh)
        except Exception as e:  # noqa: BLE001
            print(f"[mask_report] skip {p}: {e}", file=sys.stderr)
            continue
        if isinstance(d, dict) and "samples" in d:
            for s in d["samples"]:
                out[_sample_name(s, p)] = s
        else:
            out[_sample_name(d, p)] = d
    return out


def read_tsv(path):
    if not path or not os.path.exists(path):
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def b64_bins(s):
    if not s:
        return []
    raw = base64.b64decode(s)
    return list(raw)


# ─────────────────────────────────────────────────────────────────────────────
# SVG coverage track
# ─────────────────────────────────────────────────────────────────────────────
def track_svg(bh, mode, width=760, height=92):
    bins = b64_bins(bh.get("b64"))
    mbins = b64_bins(bh.get("mask_b64"))
    n = len(bins)
    if not n:
        return ""
    ml, mr, mt, mb = 34, 6, 6, 18
    iw, ih = width - ml - mr, height - mt - mb
    bw = iw / n
    parts = [f'<svg viewBox="0 0 {width} {height}" width="100%" role="img" '
             f'aria-label="positional breadth with {"targets" if mode == "target" else "masked regions"}">']
    parts.append(f'<rect x="{ml}" y="{mt}" width="{iw}" height="{ih}" class="trk-bg"/>')
    # mask layer first (behind bars). exclude: shade the masked share of each
    # bin; target: mask_bins holds the share OUTSIDE targets, so shade 100-m.
    # Consecutive bins with the same (5%-rounded) value merge into one rect.
    layer = []
    for i in range(n):
        m = mbins[i] if i < len(mbins) else 0
        v = (100 - m) if mode == "target" else m
        layer.append(int(round(v / 5.0)) * 5)
    cls = "trk-target" if mode == "target" else "trk-mask"
    i = 0
    while i < n:
        v = layer[i]
        j = i
        while j + 1 < n and layer[j + 1] == v:
            j += 1
        if v > 0:
            x0 = ml + i * bw
            w = max((j - i + 1) * bw, 1.0)
            parts.append(f'<rect x="{x0:.2f}" y="{mt}" width="{w:.2f}" height="{ih}" class="{cls}" '
                         f'opacity="{0.25 + 0.75 * v / 100:.2f}" shape-rendering="crispEdges"/>')
        i = j + 1
    for i, v in enumerate(bins):
        if v <= 0:
            continue
        h = ih * min(100, v) / 100.0
        parts.append(f'<rect x="{ml + i * bw:.2f}" y="{mt + ih - h:.2f}" width="{max(bw, 0.6):.2f}" '
                     f'height="{h:.2f}" class="trk-bar"/>')
    for b in (bh.get("breaks") or [])[1:]:
        x = ml + b * bw
        parts.append(f'<line x1="{x:.2f}" x2="{x:.2f}" y1="{mt}" y2="{mt + ih}" class="trk-brk"/>')
    tot = int(bh.get("total_len") or n * (bh.get("bin_size") or 1))
    for frac in (0, 0.25, 0.5, 0.75, 1.0):
        x = ml + iw * frac
        parts.append(f'<text x="{x:.1f}" y="{height - 4}" class="trk-ax" text-anchor="'
                     f'{"start" if frac == 0 else "end" if frac == 1 else "middle"}">{_fmt_bp(tot * frac)}</text>')
    parts.append(f'<text x="{ml - 4}" y="{mt + 8}" class="trk-ax" text-anchor="end">100%</text>')
    parts.append(f'<text x="{ml - 4}" y="{mt + ih}" class="trk-ax" text-anchor="end">0</text>')
    parts.append("</svg>")
    return "".join(parts)


def _fmt_bp(v):
    v = float(v)
    if v >= 1e6:
        return f"{v / 1e6:.2f} Mb"
    if v >= 1e3:
        return f"{v / 1e3:.1f} kb"
    return f"{int(v)} bp"


# ─────────────────────────────────────────────────────────────────────────────
def main(argv=None):
    a = parse_args(argv)
    masked = load(a.jsons)
    base = load(a.baseline)
    asm_rows = read_tsv(a.assembly_table)
    mask_sum = read_tsv(a.mask_summary)
    build_info = {}
    if a.mask_json and os.path.exists(a.mask_json):
        with open(a.mask_json) as fh:
            build_info = json.load(fh)

    # baseline lookup: sample -> strain key -> tass
    base_idx = {}
    for s, d in base.items():
        for _, _, st in _iter_strains(d):
            base_idx[(s, str(st.get("key")))] = st
    rows = []
    tracks = []
    sample_rows = []
    mode = None
    for s, d in sorted(masked.items()):
        m = (d.get("metadata") or {}).get("masking") or {}
        mode = mode or m.get("mode")
        sample_rows.append({
            "sample": s,
            "mode": m.get("mode", "—"),
            "references_masked": m.get("references_masked", 0),
            "references_in_bam": m.get("references_in_bam", 0),
            "masked_bp": m.get("masked_bp", 0),
            "kept_bp": m.get("kept_bp", 0),
            "reads_considered": m.get("reads_considered", 0),
            "reads_masked": m.get("reads_masked", 0),
            "penalty_power": m.get("penalty_power", 0),
            "fcs_contaminated": m.get("fcs_assemblies_contaminated", ""),
        })
        seen = set()
        for g, sk, st in _iter_strains(d):
            key = str(st.get("key"))
            seen.add(key)
            b = base_idx.get((s, key))
            t_m = _f(st.get("tass_score"))
            t_b = _f(b.get("tass_score")) if b else None
            if max(t_m, t_b or 0) < a.min_tass:
                continue
            tot = _f(st.get("mask_reads_total"))
            msk = _f(st.get("mask_reads_masked"))
            fcs = st.get("fcs_assembly") or {}
            row = {
                "sample": s,
                "genus": g.get("name", ""),
                "species": sk.get("name", ""),
                "organism": st.get("name", ""),
                "taxid": key,
                "tass_masked": round(t_m, 2),
                "tass_unmasked": round(t_b, 2) if t_b is not None else "",
                "tass_delta": round(t_m - t_b, 2) if t_b is not None else "",
                "reads_kept": int(_f(st.get("numreads"))),
                "reads_masked": int(msk),
                "reads_masked_pct": round(100 * msk / tot, 2) if tot else 0.0,
                "scored_bp": int(_f(st.get("length"))),
                "masked_bp": int(_f(st.get("masked_bp"))),
                "masked_bp_pct": round(100 * _f(st.get("mask_bp_fraction")), 3),
                "breadth_pct": round(100 * _f(st.get("coverage")), 4),
                "penalty": round(_f(st.get("mask_penalty"), 1.0), 4),
                "fcs_assembly": fcs.get("assembly", ""),
                "fcs_contaminant_bp": fcs.get("contaminant_bp") if fcs.get("contaminant_bp") is not None else "",
                "fcs_fraction": fcs.get("fraction_contaminated") if fcs.get("fraction_contaminated") is not None else "",
            }
            rows.append(row)
            bh = st.get("breadth_histogram") or {}
            if bh.get("b64") and (bh.get("mask_b64") or st.get("mask_regions")):
                impact = row["reads_masked"] + abs(_f(row["tass_delta"])) * 10 + row["tass_masked"] / 10
                tracks.append((impact, s, row, bh, bh.get("mask_mode") or m.get("mode"), st.get("mask_regions") or []))
        # organisms that vanished entirely after masking
        for (bs, key), b in base_idx.items():
            if bs != s or key in seen:
                continue
            t_b = _f(b.get("tass_score"))
            if t_b < a.min_tass:
                continue
            rows.append({
                "sample": s, "genus": "", "species": "", "organism": b.get("name", ""), "taxid": key,
                "tass_masked": "", "tass_unmasked": round(t_b, 2), "tass_delta": "removed",
                "reads_kept": 0, "reads_masked": "", "reads_masked_pct": "", "scored_bp": "", "masked_bp": "",
                "masked_bp_pct": "", "breadth_pct": "", "penalty": "", "fcs_assembly": "",
                "fcs_contaminant_bp": "", "fcs_fraction": "",
            })

    cols = ["sample", "organism", "taxid", "tass_masked", "tass_unmasked", "tass_delta", "reads_kept",
            "reads_masked", "reads_masked_pct", "scored_bp", "masked_bp", "masked_bp_pct", "breadth_pct",
            "penalty", "fcs_assembly", "fcs_contaminant_bp", "fcs_fraction", "species", "genus"]
    with open(a.prefix + ".tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    tracks.sort(key=lambda t: -t[0])
    html_out = render(a, mode, sample_rows, rows, cols, tracks[: a.top], asm_rows, mask_sum, build_info, bool(base))
    with open(a.prefix + ".html", "w") as fh:
        fh.write(html_out)
    print(f"[mask_report] {len(masked)} samples, {len(rows)} organism rows, "
          f"{min(len(tracks), a.top)} tracks -> {a.prefix}.html / .tsv", file=sys.stderr)


CSS = """
:root{--bg:#fbfbfa;--fg:#1d1d1f;--muted:#6b6b70;--card:#fff;--line:#e3e3e0;--accent:#2f6fde;
--bar:#2f6fde;--mask:#9a9aa1;--target:#f2b84b;--up:#1a7f4b;--down:#c2412d;--hdr:#f3f3f1}
@media (prefers-color-scheme: dark){:root:not([data-theme=light]){--bg:#151517;--fg:#ececee;--muted:#9b9ba3;
--card:#1d1d20;--line:#323238;--accent:#6c9bf0;--bar:#6c9bf0;--mask:#6e6e76;--target:#c9942f;--up:#4cc38a;
--down:#f07a64;--hdr:#232327}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);
font:14px/1.45 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}
main{max-width:1180px;margin:0 auto;padding:28px 16px 60px}
h1{font-size:22px;margin:0 0 4px}h2{font-size:16px;margin:30px 0 10px}
.sub{color:var(--muted);margin:0 0 18px}
.kpis{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:10px}
.kpi{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:10px 12px}
.kpi b{display:block;font-size:20px;font-variant-numeric:tabular-nums}.kpi span{color:var(--muted);font-size:12px}
.wrap{overflow-x:auto;background:var(--card);border:1px solid var(--line);border-radius:10px}
table{border-collapse:collapse;width:100%;font-size:12.5px;font-variant-numeric:tabular-nums}
th,td{padding:6px 9px;border-bottom:1px solid var(--line);text-align:right;white-space:nowrap}
th{background:var(--hdr);position:sticky;top:0;cursor:pointer;user-select:none;font-weight:600}
td.l,th.l{text-align:left}tr:last-child td{border-bottom:0}
.up{color:var(--up)}.down{color:var(--down)}.flag{color:var(--down);font-weight:600}
.tools{display:flex;gap:10px;align-items:center;margin:0 0 8px;flex-wrap:wrap}
input[type=search]{padding:6px 9px;border:1px solid var(--line);border-radius:7px;background:var(--card);color:var(--fg);min-width:240px}
.tracks{display:grid;gap:12px}
.trk{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:10px 12px}
.trk h3{font-size:13px;margin:0 0 2px}.trk p{margin:0 0 6px;color:var(--muted);font-size:12px}
.trk-bg{fill:var(--bg);stroke:var(--line)}.trk-bar{fill:var(--bar)}.trk-mask{fill:var(--mask)}
.trk-target{fill:var(--target)}.trk-brk{stroke:var(--muted);stroke-dasharray:2 2}
.trk-ax{fill:var(--muted);font-size:10px}
.legend{display:flex;gap:14px;color:var(--muted);font-size:12px;margin:0 0 10px}
.sw{display:inline-block;width:11px;height:11px;border-radius:2px;margin-right:5px;vertical-align:-1px}
code{background:var(--hdr);padding:1px 5px;border-radius:4px}
"""

JS = """
document.querySelectorAll('table.sortable').forEach(t=>{
  t.querySelectorAll('th').forEach((th,i)=>th.addEventListener('click',()=>{
    const tb=t.tBodies[0],rows=[...tb.rows],asc=th.dataset.dir!=='asc';
    t.querySelectorAll('th').forEach(h=>delete h.dataset.dir);th.dataset.dir=asc?'asc':'desc';
    rows.sort((a,b)=>{const x=a.cells[i].dataset.v??a.cells[i].textContent,y=b.cells[i].dataset.v??b.cells[i].textContent;
      const nx=parseFloat(x),ny=parseFloat(y);
      if(!isNaN(nx)&&!isNaN(ny))return asc?nx-ny:ny-nx;return asc?x.localeCompare(y):y.localeCompare(x)});
    rows.forEach(r=>tb.appendChild(r));}));});
const q=document.getElementById('q');if(q)q.addEventListener('input',()=>{const v=q.value.toLowerCase();
  document.querySelectorAll('#orgs tbody tr').forEach(r=>r.style.display=r.textContent.toLowerCase().includes(v)?'':'none')});
"""


def _td(v, cls=""):
    if isinstance(v, int) and not isinstance(v, bool) and abs(v) >= 10000:
        e = f"{v:,}"
    else:
        e = html.escape(str(v))
    dv = ""
    if isinstance(v, (int, float)):
        dv = f' data-v="{v}"'
    return f'<td class="{cls}"{dv}>{e}</td>'


def render(a, mode, sample_rows, rows, cols, tracks, asm_rows, mask_sum, build_info, has_base):
    E = html.escape
    tot_reads = sum(_f(r["reads_considered"]) for r in sample_rows)
    tot_masked = sum(_f(r["reads_masked"]) for r in sample_rows)
    moved = [r for r in rows if isinstance(r.get("tass_delta"), (int, float)) and abs(r["tass_delta"]) >= 1]
    removed = [r for r in rows if r.get("tass_delta") == "removed"]
    asm_flag = [r for r in asm_rows if _f(r.get("total_contaminant_length")) > 0]
    mode_txt = {"exclude": "Exclude — listed regions are ignored", "target": "Target — only listed regions are scored"}.get(mode or "", "—")
    o = [f"<!doctype html><html lang=en><head><meta charset=utf-8>"
         f"<meta name=viewport content='width=device-width,initial-scale=1'><title>TaxTriage masking report</title>"
         f"<style>{CSS}</style></head><body><main>"]
    o.append("<h1>TaxTriage masking report</h1>")
    o.append(f"<p class=sub>Mode: <b>{E(mode_txt)}</b>. Scores below were computed only on unmasked bases"
             f"{' and compared against an unmasked run of the same samples' if has_base else ''}.</p>")
    o.append("<div class=kpis>")
    for val, lab in [(len(sample_rows), "samples"),
                     (f"{int(tot_masked):,} / {int(tot_reads):,}", "reads removed by mask"),
                     (f"{(100 * tot_masked / tot_reads if tot_reads else 0):.2f}%", "of considered reads"),
                     (len({r['organism'] for r in rows if _f(r.get('reads_masked')) > 0}), "organisms with masked reads"),
                     (len(moved) if has_base else "—", "TASS moved ≥ 1 point"),
                     (len(removed) if has_base else "—", "detections removed"),
                     (len(asm_flag) if asm_rows else "—", "assemblies NCBI FCS flags")]:
        o.append(f"<div class=kpi><b>{E(str(val))}</b><span>{E(lab)}</span></div>")
    o.append("</div>")

    o.append("<h2>Per-sample read accounting</h2><div class=wrap><table class=sortable><thead><tr>")
    sc = ["sample", "references_masked", "references_in_bam", "masked_bp", "kept_bp", "reads_considered",
          "reads_masked", "penalty_power", "fcs_contaminated"]
    o.append("".join(f"<th class='{'l' if c == 'sample' else ''}'>{E(c.replace('_', ' '))}</th>" for c in sc))
    o.append("</tr></thead><tbody>")
    for r in sample_rows:
        o.append("<tr>" + "".join(_td(r[c], "l" if c == "sample" else "") for c in sc) + "</tr>")
    o.append("</tbody></table></div>")

    o.append("<h2>Organisms</h2><div class=tools><input id=q type=search placeholder='Filter organisms, samples, assemblies…'>"
             f"<span class=sub style='margin:0'>{len(rows)} rows · click a column to sort · also in <code>{E(os.path.basename(a.prefix))}.tsv</code></span></div>")
    o.append("<div class=wrap style='max-height:560px'><table id=orgs class=sortable><thead><tr>")
    show = [c for c in cols if c not in ("species", "genus") and (has_base or c not in ("tass_unmasked", "tass_delta"))]
    o.append("".join(f"<th class='{'l' if c in ('sample', 'organism', 'fcs_assembly') else ''}'>{E(c.replace('_', ' '))}</th>" for c in show))
    o.append("</tr></thead><tbody>")
    for r in sorted(rows, key=lambda r: (-_f(r.get("reads_masked")), -_f(r.get("tass_masked")))):
        cells = []
        for c in show:
            v = r.get(c, "")
            cls = "l" if c in ("sample", "organism", "fcs_assembly") else ""
            if c == "tass_delta" and isinstance(v, (int, float)):
                cls += " up" if v > 0 else " down" if v < 0 else ""
            if c == "tass_delta" and v == "removed":
                cls += " flag"
            if c == "fcs_contaminant_bp" and _f(v) > 0:
                cls += " flag"
            cells.append(_td(v, cls))
        o.append("<tr>" + "".join(cells) + "</tr>")
    o.append("</tbody></table></div>")

    if tracks:
        o.append(f"<h2>Coverage tracks — {len(tracks)} most-affected organisms</h2>")
        if mode == "target":
            o.append("<div class=legend><span><i class=sw style='background:var(--target)'></i>targeted (scored)</span>"
                     "<span><i class=sw style='background:var(--bar)'></i>% of bin covered</span>"
                     "<span>dashed = contig break</span></div>")
        else:
            o.append("<div class=legend><span><i class=sw style='background:var(--mask)'></i>masked (not scored)</span>"
                     "<span><i class=sw style='background:var(--bar)'></i>% of bin covered</span>"
                     "<span>dashed = contig break</span></div>")
        o.append("<div class=tracks>")
        for _, s, r, bh, tmode, regs in tracks:
            lab = []
            for reg in regs[:6]:
                n_iv = len(reg.get("intervals") or [])
                lab.append(f"{reg['accession']}: {n_iv} {'target' if tmode == 'target' else 'masked'} interval{'s' if n_iv != 1 else ''}")
            meta = (f"{E(s)} · TASS {r['tass_masked']}"
                    + (f" (unmasked {r['tass_unmasked']})" if has_base and r['tass_unmasked'] != '' else "")
                    + f" · {r['reads_kept']} reads kept, {r['reads_masked']} masked · {E('; '.join(lab))}")
            o.append(f"<div class=trk><h3>{E(r['organism'])}</h3><p>{meta}</p>{track_svg(bh, tmode)}</div>")
        o.append("</div>")

    if asm_rows:
        o.append("<h2>NCBI FCS — assembly screen of the references used</h2>")
        o.append("<p class=sub>From NCBI's precomputed FCS-GX summaries. A contaminated assembly can still be "
                 "the right call; its flagged spans are masked when per-sequence details are available "
                 "(<code>masked_bp</code> &gt; 0).</p>")
        ac = list(asm_rows[0].keys())
        o.append("<div class=wrap style='max-height:420px'><table class=sortable><thead><tr>"
                 + "".join(f"<th class='{'l' if c in ('assembly', 'fcs_source', 'asserted_division', 'top_inferred_division') else ''}'>{E(c.replace('_', ' '))}</th>" for c in ac)
                 + "</tr></thead><tbody>")
        for r in sorted(asm_rows, key=lambda r: -_f(r.get("total_contaminant_length"))):
            o.append("<tr>" + "".join(
                _td(r.get(c, ""), ("l" if c in ("assembly", "fcs_source", "asserted_division", "top_inferred_division") else "")
                    + (" flag" if c == "total_contaminant_length" and _f(r.get(c)) > 0 else ""))
                for c in ac) + "</tr>")
        o.append("</tbody></table></div>")

    if mask_sum:
        o.append("<h2>Mask contents by reference</h2><div class=wrap style='max-height:420px'><table class=sortable><thead><tr>")
        mc = list(mask_sum[0].keys())
        o.append("".join(f"<th class='{'l' if c in ('reference', 'assembly', 'labels', 'sources') else ''}'>{E(c.replace('_', ' '))}</th>" for c in mc))
        o.append("</tr></thead><tbody>")
        for r in mask_sum:
            r = {k: (int(v) if isinstance(v, str) and v.isdigit() else v) for k, v in r.items()}
            o.append("<tr>" + "".join(_td(r.get(c, ""), "l" if c in ("reference", "assembly", "labels", "sources") else "") for c in mc) + "</tr>")
        o.append("</tbody></table></div>")
    if build_info:
        o.append("<h2>Mask build</h2><div class=wrap><table><tbody>")
        for k, v in build_info.items():
            if isinstance(v, list) and len(v) > 12:
                v = ", ".join(map(str, v[:12])) + f" … (+{len(v) - 12})"
            o.append(f"<tr><td class=l>{E(k)}</td><td class=l>{E(json.dumps(v) if isinstance(v, (dict, list)) else str(v))}</td></tr>")
        o.append("</tbody></table></div>")
    o.append(f"<script>{JS}</script></main></body></html>")
    return "".join(o)


if __name__ == "__main__":
    main()
