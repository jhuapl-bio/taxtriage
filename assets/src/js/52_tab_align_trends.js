/* ═══════════════════════════════════════════════════════════════════════════
       -  §  TRENDS SUB-TAB: ALIGNMENT TRENDS   (data-metasub="align")
       -     Cross-sample frequency of zero / low / high depth windows along
       -     each reference, from the per-strain `depth_profile` that
       -     match_paths.py embeds (bin/depth_profile.py documents the format).
       -     Mirrors bin/alignment_trends.py, so the numbers here match the
       -     pipeline's --alignment_trends tables for the same settings.
       -
       -     Per sample, window depth is normalised by that sample's mean depth
       -     on the reference, so deep and shallow libraries compare. A sample
       -     only counts toward a window when it EXPECTS enough reads there
       -     (sample reads x window / reference length >= "min reads / window"); with
       -     fewer, an empty window is Poisson noise rather than a real gap.
       -     Windows where a state recurs in >= the recurrence fraction of the
       -     counted samples merge into regions: recurrent zero / low regions
       -     mark stretches of the reference the sample set does not carry
       -     (deletions, divergent or novel loci, assembly artefacts); recurrent
       -     high regions mark repeats, rRNA operons, mobile elements, plasmid
       -     copy number or contamination-prone loci.
       -
       -     Two views: "Live analysis" (this file's analysis over the visible
       -     samples, every cutoff adjustable) and, when the run used
       -     --alignment_trends, "Pipeline results" — the embedded
       -     all.alignment_trends.json reported as the pipeline computed it.
═══════════════════════════════════════════════════════════════════════════ */
(function () {
  const AT = {
    cache: new Map(), // `${refKey}|${sig}` -> analysis
    zoom: null, // [w0, w1) window-index domain, null = full
    selected: "",
    wired: false,
  };
  const MAX_COLS = 900; // display columns for the heatmap / tracks
  const esc = (s) =>
    String(s == null ? "" : s)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;");
  const fmtBp = (v) =>
    v >= 1e6 ? (v / 1e6).toFixed(v >= 1e7 ? 1 : 2) + " Mb" : v >= 1e3 ? (v / 1e3).toFixed(v >= 1e4 ? 0 : 1) + " kb" : v + " bp";
  const fmtInt = (v) => Number(v || 0).toLocaleString();
  const $ = (id) => document.getElementById(id);

  /* ── profile decoding (bin/depth_profile.py) ─────────────────────────── */
  function _b64(s) {
    const raw = atob(s || "");
    const out = new Uint8Array(raw.length);
    for (let i = 0; i < raw.length; i++) out[i] = raw.charCodeAt(i);
    return out;
  }
  const _DEP_LUT = (() => {
    const t = new Float32Array(256);
    for (let b = 1; b < 256; b++) t[b] = Math.pow(2, b / 16) - 1;
    return t;
  })();
  // Decoded profiles live in a WeakMap, NOT on the profile object, so session
  // save (29_session_state.js serialises CONTIG_DATA) never sees typed arrays.
  const _DEC = new WeakMap();
  function _decode(prof) {
    if (_DEC.has(prof)) return _DEC.get(prof);
    const w = +prof.window || 100;
    const contigs = new Map();
    (prof.contigs || []).forEach((e) => {
      const acc = String(e[0]);
      const len = +e[1] || 0;
      const n = Math.max(1, Math.ceil(len / w));
      const dep = new Float32Array(n);
      const br = new Float32Array(n);
      if (e.length >= 4 && e[2]) {
        const db = _b64(e[2]);
        const bb = _b64(e[3]);
        for (let i = 0; i < n && i < db.length; i++) {
          dep[i] = _DEP_LUT[db[i]];
          br[i] = bb[i] || 0;
        }
      }
      contigs.set(acc, { len, dep, br });
    });
    const dec = { w, contigs };
    _DEC.set(prof, dec);
    return dec;
  }
  function _rebin(vals, f, len, w, n) {
    const out = new Float32Array(n);
    if (f <= 1) {
      out.set(vals.subarray(0, Math.min(n, vals.length)));
      return out;
    }
    for (let o = 0; o < n; o++) {
      let tot = 0,
        ws = 0;
      for (let j = o * f; j < Math.min(vals.length, (o + 1) * f); j++) {
        const lo = j * w;
        const bl = Math.max(1, Math.min(lo + w, len) - lo);
        tot += vals[j] * bl;
        ws += bl;
      }
      out[o] = ws ? tot / ws : 0;
    }
    return out;
  }

  /* ── settings ────────────────────────────────────────────────────────── */
  function _num(id, dflt) {
    const el = $(id);
    const v = el ? parseFloat(el.value) : NaN;
    return isFinite(v) ? v : dflt;
  }
  function _opts() {
    return {
      lowFrac: _num("at-low", 0.2),
      highFrac: _num("at-high", 3),
      minRpw: _num("at-minrpw", 1),
      minFreq: _num("at-recur", 50) / 100,
      minReads: _num("at-minreads", 3),
      minSamples: Math.max(1, Math.round(_num("at-minsamples", 2))),
      coarsen: Math.max(1, Math.round(_num("at-coarsen", 4))),
      useFilters: !!($("at-use-filters") || { checked: true }).checked,
      sort: ($("at-sort") || { value: "order" }).value,
    };
  }
  const _optSig = (o) =>
    [o.lowFrac, o.highFrac, o.minRpw, o.minFreq, o.minReads, o.minSamples, o.coarsen, o.useFilters].join(",");

  /* ── which (sample, reference) entries are in play ──────────────────── */
  function _passingTaxa() {
    const set = new Set();
    filteredData().forEach((r) => {
      const pi = typeof rowPassInfo === "function" ? rowPassInfo(r) : null;
      if (pi && !isNaN(pi.strain) && !pi.strainPass) return;
      const tax = r["Taxonomic ID #"];
      set.add(`${r["Specimen ID"]}||${tax}`);
      if (Array.isArray(r.__mergedFrom)) r.__mergedFrom.forEach((s) => set.add(`${s}||${tax}`));
    });
    return set;
  }
  function _refIndex(o) {
    const pass = o.useFilters ? _passingTaxa() : null;
    const refs = new Map();
    (typeof CONTIG_DATA !== "undefined" ? CONTIG_DATA : []).forEach((cd) => {
      if (!cd || !cd.depth_profile || !cd.depth_profile.contigs) return;
      if (sampleHidden[cd.sample]) return;
      if (pass && !pass.has(`${cd.sample}||${cd.taxon_id}`)) return;
      const k = String(cd.taxon_id || cd.organism);
      if (!refs.has(k)) refs.set(k, { key: k, name: cd.organism || k, entries: [] });
      refs.get(k).entries.push(cd);
    });
    return refs;
  }
  function _readsOf(cd) {
    const c = cd.contigs || [];
    let n = 0;
    c.forEach((x) => (n += +x.reads || 0));
    return n;
  }

  /* ── analysis (mirror of bin/alignment_trends.py) ───────────────────── */
  function _analyse(ref, o) {
    const sig = _optSig(o) + "|" + ref.entries.map((e) => e.sample).sort().join(",");
    const ck = ref.key + "|" + sig;
    if (AT.cache.has(ck)) return AT.cache.get(ck);
    const dec = ref.entries.map((cd) => ({ cd, d: _decode(cd.depth_profile) }));
    let W = 0;
    dec.forEach((x) => (W = Math.max(W, x.d.w)));
    for (let i = 1; i < o.coarsen; i *= 2) W *= 2;
    const lens = new Map();
    dec.forEach((x) => x.d.contigs.forEach((c, acc) => lens.set(acc, Math.max(lens.get(acc) || 0, c.len))));
    const layout = [...lens.keys()].sort().map((acc) => {
      const len = lens.get(acc);
      return { acc, len, n: Math.max(1, Math.ceil(len / W)) };
    });
    let nwin = 0;
    layout.forEach((c) => {
      c.off = nwin;
      nwin += c.n;
    });
    const wlen = new Float32Array(nwin);
    const wStart = new Float64Array(nwin); // concatenated bp coordinate of each window start
    let bpOff = 0;
    layout.forEach((c) => {
      for (let i = 0; i < c.n; i++) {
        wlen[c.off + i] = Math.max(1, Math.min((i + 1) * W, c.len) - i * W);
        wStart[c.off + i] = bpOff + i * W;
      }
      c.bpOff = bpOff;
      bpOff += c.len;
    });
    // Full reference size: fragmented assemblies list zero-read contigs only up
    // to a cap (bin/depth_profile.py), so contigs no sample covered can be
    // absent from the layout — they are reported as never-covered bp.
    let totalLen = bpOff,
      nContigsFull = layout.length;
    dec.forEach((x) => {
      totalLen = Math.max(totalLen, +x.cd.depth_profile.total_len || 0);
      nContigsFull = Math.max(nContigsFull, +x.cd.depth_profile.n_contigs || 0);
    });
    const neverBp = Math.max(0, totalLen - bpOff);
    const neverContigs = Math.max(0, nContigsFull - layout.length);

    const tracks = dec.map(({ cd, d }) => {
      const f = Math.max(1, Math.round(W / d.w));
      const dep = new Float32Array(nwin);
      const br = new Float32Array(nwin);
      layout.forEach((c) => {
        const src = d.contigs.get(c.acc);
        if (!src) return;
        dep.set(_rebin(src.dep, f, src.len, d.w, c.n), c.off);
        br.set(_rebin(src.br, f, src.len, d.w, c.n), c.off);
      });
      let sd = 0,
        sb = 0;
      for (let i = 0; i < nwin; i++) {
        sd += dep[i] * wlen[i];
        sb += br[i] * wlen[i];
      }
      const mean = totalLen ? sd / totalLen : 0;
      const breadth = totalLen ? sb / totalLen : 0;
      // Reads expected per window: the profile's read count spread evenly over
      // the reference (falls back to the contig table, then to depth / read length).
      const reads = +cd.depth_profile.reads || _readsOf(cd);
      const rl = Math.max(1, +cd.depth_profile.read_len || 150);
      const perBp = reads > 0 && totalLen ? reads / totalLen : mean / rl;
      const expRpw = perBp * W;
      const informative = reads >= o.minReads && expRpw >= o.minRpw;
      const lo = o.lowFrac * mean,
        hi = o.highFrac * mean;
      const cls = new Uint8Array(nwin); // 0 normal, 1 zero, 2 low, 3 high
      const ok = new Uint8Array(nwin);
      const cnt = [0, 0, 0, 0];
      for (let i = 0; i < nwin; i++) {
        const dd = dep[i];
        cls[i] = br[i] <= 0 || dd <= 0 ? 1 : dd < lo ? 2 : dd > hi ? 3 : 0;
        cnt[cls[i]]++;
        ok[i] = informative && perBp * wlen[i] >= o.minRpw ? 1 : 0;
      }
      return {
        sample: cd.sample,
        dep,
        br,
        cls,
        ok,
        mean,
        breadth,
        reads,
        expRpw,
        informative,
        pct: cnt.map((c) => (100 * c) / nwin),
      };
    });

    const nInf = tracks.filter((t) => t.informative).length;
    const fz = new Float32Array(nwin).fill(NaN);
    const fl = new Float32Array(nwin).fill(NaN);
    const fh = new Float32Array(nwin).fill(NaN);
    const nS = new Uint16Array(nwin);
    const mNorm = new Float32Array(nwin).fill(NaN);
    if (nInf >= o.minSamples) {
      for (let i = 0; i < nwin; i++) {
        let n = 0,
          z = 0,
          l = 0,
          h = 0,
          sn = 0;
        for (const t of tracks) {
          if (!t.ok[i]) continue;
          n++;
          const c = t.cls[i];
          if (c === 1) z++;
          else if (c === 2) l++;
          else if (c === 3) h++;
          if (t.mean > 0) sn += t.dep[i] / t.mean;
        }
        nS[i] = n;
        if (n < o.minSamples) continue;
        fz[i] = z / n;
        fl[i] = (z + l) / n; // low includes zero
        fh[i] = h / n;
        mNorm[i] = sn / n;
      }
    }
    // recurrent regions — runs within one contig
    const regions = [];
    const defs = [
      ["zero", fz, [1]],
      ["low", fl, [1, 2]],
      ["high", fh, [3]],
    ];
    defs.forEach(([type, fr, member]) => {
      layout.forEach((c) => {
        let s = -1;
        const flush = (e) => {
          if (s < 0) return;
          let sf = 0,
            mx = 0,
            sn = 0,
            ns = 0;
          const aff = new Set();
          for (let i = s; i < e; i++) {
            sf += fr[i];
            mx = Math.max(mx, fr[i]);
            sn += mNorm[i];
            ns = Math.max(ns, nS[i]);
            tracks.forEach((t) => {
              if (t.ok[i] && member.includes(t.cls[i])) aff.add(t.sample);
            });
          }
          const st = (s - c.off) * W;
          const en = Math.min((e - c.off) * W, c.len);
          regions.push({
            type,
            contig: c.acc,
            start: st,
            end: en,
            length: en - st,
            w0: s,
            w1: e,
            nWin: e - s,
            meanFreq: sf / (e - s),
            maxFreq: mx,
            nSamples: ns,
            meanNorm: sn / (e - s),
            wholeContig: st === 0 && en === c.len,
            affected: [...aff].sort(),
          });
          s = -1;
        };
        for (let i = c.off; i < c.off + c.n; i++) {
          if (fr[i] >= o.minFreq) {
            if (s < 0) s = i;
          } else flush(i);
        }
        flush(c.off + c.n);
      });
    });
    const sum = { zero: 0, low: 0, high: 0 };
    const nReg = { zero: 0, low: 0, high: 0 };
    regions.forEach((r) => {
      sum[r.type] += r.length;
      nReg[r.type]++;
    });
    const res = {
      ref,
      W,
      layout,
      nwin,
      wlen,
      wStart,
      totalLen,
      nContigsFull,
      neverBp,
      neverContigs,
      tracks,
      nInf,
      fz,
      fl,
      fh,
      nS,
      mNorm,
      regions,
      pct: {
        zero: totalLen ? (100 * sum.zero) / totalLen : 0,
        low: totalLen ? (100 * sum.low) / totalLen : 0,
        high: totalLen ? (100 * sum.high) / totalLen : 0,
      },
      nReg,
      analysed: nInf >= o.minSamples,
    };
    if (AT.cache.size > 300) AT.cache.clear();
    AT.cache.set(ck, res);
    return res;
  }

  /* ── controls ─────────────────────────────────────────────────────────── */
  function _wire() {
    if (AT.wired) return;
    AT.wired = true;
    const rerender = () => {
      AT.zoom = null;
      drawAlignTrends();
    };
    ["at-low", "at-high", "at-minrpw", "at-recur", "at-minreads", "at-minsamples", "at-coarsen", "at-use-filters", "at-sort"].forEach(
      (id) => {
        const el = $(id);
        if (el) el.addEventListener("change", id === "at-sort" ? () => drawAlignTrends() : rerender);
      },
    );
    const sel = $("at-ref-sel");
    if (sel)
      sel.addEventListener("change", () => {
        AT.selected = sel.value;
        AT.zoom = null;
        drawAlignTrends();
      });
    const srch = $("at-ref-search");
    if (srch)
      srch.addEventListener("input", () => {
        const prev = AT.selected;
        _fillRefSelect(_lastRefs, _lastOverview);
        if (AT.selected !== prev) {
          AT.zoom = null;
          drawAlignTrends();
        }
      });
    const rz = $("at-reset-zoom");
    if (rz)
      rz.addEventListener("click", (e) => {
        e.preventDefault();
        AT.zoom = null;
        drawAlignTrends();
      });
    const dl = (which) => (e) => {
      e.preventDefault();
      _download(which);
    };
    if ($("at-dl-regions")) $("at-dl-regions").addEventListener("click", dl("regions"));
    if ($("at-dl-windows")) $("at-dl-windows").addEventListener("click", dl("windows"));
    if ($("at-dl-overview")) $("at-dl-overview").addEventListener("click", dl("overview"));
  }

  let _lastRefs = new Map();
  let _lastOverview = [];
  let _lastRes = null;

  function _fillRefSelect(refs, overview) {
    const sel = $("at-ref-sel");
    if (!sel) return;
    const q = (($("at-ref-search") || {}).value || "").trim().toLowerCase();
    const rows = overview.filter((r) => !q || r.name.toLowerCase().includes(q) || r.key.includes(q));
    sel.innerHTML = "";
    rows.forEach((r) => {
      const opt = document.createElement("option");
      opt.value = r.key;
      opt.textContent = `${r.name} — ${r.nInf}/${r.nSamples} sample${r.nSamples === 1 ? "" : "s"}`;
      sel.appendChild(opt);
    });
    if (AT.selected && rows.some((r) => r.key === AT.selected)) sel.value = AT.selected;
    else if (rows.length) AT.selected = sel.value = rows[0].key;
  }

  /* ── main draw ────────────────────────────────────────────────────────── */
  function drawAlignTrends() {
    _wire();
    _atpInit();
    if (ATP.view === "pipeline") return drawAlignTrendsPipeline();
    const o = _opts();
    const refs = _refIndex(o);
    _lastRefs = refs;
    const empty = $("at-empty");
    const body = $("at-body");
    if (!refs.size) {
      if (empty) {
        empty.style.display = "";
        empty.innerHTML =
          "No depth profiles in this report. Per-sample JSONs from TaxTriage runs that include " +
          "<code>depth_profile</code> (default; disabled by <code>--depth_profile_windows 0</code>) light this tab up. " +
          "Hidden samples and detections removed by the active filters are excluded.";
      }
      if (body) body.style.display = "none";
      return;
    }
    if (empty) empty.style.display = "none";
    if (body) body.style.display = "";

    // Overview across references (analysis is cached per settings + sample set)
    const overview = [];
    refs.forEach((ref) => {
      const nSamples = ref.entries.length;
      let res = null;
      if (nSamples >= 1) res = _analyse(ref, o);
      overview.push({
        key: ref.key,
        name: ref.name,
        nSamples,
        nInf: res ? res.nInf : 0,
        W: res ? res.W : 0,
        totalLen: res ? res.totalLen : 0,
        nContigs: res ? res.nContigsFull : 0,
        analysed: res ? res.analysed : false,
        pct: res ? res.pct : { zero: 0, low: 0, high: 0 },
        nReg: res ? res.nReg : { zero: 0, low: 0, high: 0 },
      });
    });
    overview.sort((a, b) => b.analysed - a.analysed || b.nInf - a.nInf || b.nSamples - a.nSamples || a.name.localeCompare(b.name));
    _lastOverview = overview;
    _fillRefSelect(refs, overview);
    _drawOverview(overview, o);

    const ref = refs.get(AT.selected);
    if (!ref) return;
    const res = _analyse(ref, o);
    _lastRes = res;
    _drawKpis(res, o);
    _drawHeat(res, o);
    _drawFreq(res, o);
    _drawRegions(res, o);
    _drawSampleTable(res, o);
  }

  function _drawKpis(res, o) {
    const el = $("at-kpis");
    if (!el) return;
    const card = (label, value, sub, color) =>
      `<div class="kpi-card" style="border-left-color:${color}"><div class="kpi-label">${label}</div>` +
      `<div class="kpi-value">${value}</div><div class="kpi-sub">${sub}</div></div>`;
    el.innerHTML =
      card(
        "Informative samples",
        `${res.nInf} / ${res.tracks.length}`,
        res.analysed ? `need ≥ ${o.minSamples}` : `<span style="color:#c62828">below ${o.minSamples} — no trends</span>`,
        "#1565c0",
      ) +
      card(
        "Reference",
        fmtBp(res.totalLen),
        `${fmtInt(res.nContigsFull)} contig(s), ${fmtInt(res.nwin)} × ${fmtBp(res.W)} windows` +
          (res.neverBp ? `<br>${fmtInt(res.neverContigs)} contig(s) · ${fmtBp(res.neverBp)} with no reads in any sample (not drawn)` : ""),
        "#455a64",
      ) +
      card("Recurrent zero", res.pct.zero.toFixed(1) + "%", `${res.nReg.zero} region(s) · ${fmtBp(Math.round((res.pct.zero * res.totalLen) / 100))}`, "#212121") +
      card("Recurrent low / zero", res.pct.low.toFixed(1) + "%", `${res.nReg.low} region(s) · < ${o.lowFrac}× mean`, "#1e88e5") +
      card("Recurrent high", res.pct.high.toFixed(1) + "%", `${res.nReg.high} region(s) · > ${o.highFrac}× mean`, "#e53935");
  }

  // Display columns over the zoom domain: each column aggregates `step` windows.
  function _cols(res) {
    const [a, b] = AT.zoom || [0, res.nwin];
    const span = Math.max(1, b - a);
    const step = Math.max(1, Math.ceil(span / MAX_COLS));
    const cols = [];
    for (let s = a; s < b; s += step) cols.push([s, Math.min(b, s + step)]);
    return { a, b, step, cols };
  }
  function _winLabel(res, i) {
    const c = res.layout.find((x) => i >= x.off && i < x.off + x.n) || res.layout[0];
    const st = (i - c.off) * res.W;
    return { contig: c.acc, start: st, end: Math.min(st + res.W, c.len) };
  }
  function _sortedTracks(res, o) {
    const t = [...res.tracks];
    if (o.sort === "depth") t.sort((x, y) => y.mean - x.mean);
    else if (o.sort === "zero") t.sort((x, y) => y.pct[1] - x.pct[1]);
    else if (o.sort === "name") t.sort((x, y) => x.sample.localeCompare(y.sample));
    else {
      const ord = typeof _orderedSamples === "function" ? _orderedSamples(t.map((x) => x.sample)) : t.map((x) => x.sample);
      const idx = new Map(ord.map((s, i) => [s, i]));
      t.sort((x, y) => idx.get(x.sample) - idx.get(y.sample));
    }
    // informative first
    return t.sort((x, y) => y.informative - x.informative);
  }

  const HEAT_ML = 190;
  function _drawHeat(res, o) {
    const host = $("at-heat");
    if (!host || typeof d3 === "undefined") return;
    host.innerHTML = "";
    const tracks = _sortedTracks(res, o);
    const { cols } = _cols(res);
    const W = Math.max(600, host.clientWidth || 900);
    const rowH = tracks.length > 40 ? 10 : 16;
    const H = 18 + rowH * tracks.length + 6;
    const iw = W - HEAT_ML - 14;
    const cw = iw / cols.length;
    const color = d3.scaleDiverging(d3.interpolateRdBu).domain([3, 0, -3]);
    const svg = d3.select(host).append("svg").attr("width", W).attr("height", H).attr("font-size", 10);
    const g = svg.append("g").attr("transform", `translate(${HEAT_ML},14)`);
    const cap = typeof _sampleNameCap === "function" ? _sampleNameCap(tracks.length) : 28;
    tracks.forEach((t, ri) => {
      const y = ri * rowH;
      const label = (t.sample.length > cap ? t.sample.slice(0, cap - 1) + "…" : t.sample) + (t.informative ? "" : " *");
      svg
        .append("text")
        .attr("x", HEAT_ML - 6)
        .attr("y", 14 + y + rowH / 2 + 3)
        .attr("text-anchor", "end")
        .attr("fill", t.informative ? "#263238" : "#90a4ae")
        .text(label)
        .append("title")
        .text(
          `${t.sample}\nmean depth ${t.mean.toFixed(2)}x · breadth ${t.breadth.toFixed(1)}% · ${fmtInt(t.reads)} reads\n` +
            `expected reads / window ${t.expRpw.toFixed(1)}${t.informative ? "" : " — below cutoff, not counted"}`,
        );
      cols.forEach(([s, e], ci) => {
        // bp-weighted depth over the windows this column aggregates (zoomed
        // out), so a column of mostly-empty windows reads as low, not normal.
        let dsum = 0,
          L = 0,
          ok = 0;
        for (let i = s; i < e; i++) {
          dsum += (t.br[i] > 0 ? t.dep[i] : 0) * res.wlen[i];
          L += res.wlen[i];
          ok += t.ok[i];
        }
        const cd = L ? dsum / L : 0;
        const fill = cd <= 0 ? "#212121" : color(Math.log2(cd / (t.mean || 1) + 0.01));
        g.append("rect")
          .attr("x", ci * cw)
          .attr("y", y)
          .attr("width", Math.max(0.6, cw + 0.15))
          .attr("height", rowH - 1)
          .attr("fill", fill)
          .attr("opacity", ok ? 1 : 0.35)
          .attr("data-c", ci)
          .attr("data-r", ri);
      });
    });
    // contig dividers (only when few enough to read)
    _contigDividers(g, res, cols, cw, rowH * tracks.length);
    g.on("mousemove", (ev) => {
      const tgt = ev.target;
      if (!tgt || tgt.tagName !== "rect" || tgt.dataset.c == null) return;
      const [s, e] = cols[+tgt.dataset.c];
      const t = tracks[+tgt.dataset.r];
      const a = _winLabel(res, s),
        b = _winLabel(res, e - 1);
      let dsum = 0,
        bsum = 0,
        L = 0;
      for (let i = s; i < e; i++) {
        dsum += t.dep[i] * res.wlen[i];
        bsum += t.br[i] * res.wlen[i];
        L += res.wlen[i];
      }
      const d = L ? dsum / L : 0;
      const cls = e - s === 1 ? ["normal", "zero", "low", "high"][t.cls[s]] : "";
      showTip(
        `<b>${esc(t.sample)}</b><br>${esc(a.contig)}:${fmtInt(a.start)}–${a.contig === b.contig ? fmtInt(b.end) : esc(b.contig) + ":" + fmtInt(b.end)}` +
          `<br>depth ${d.toFixed(2)}x (${t.mean ? (d / t.mean).toFixed(2) : "–"}× sample mean)<br>breadth ${(L ? bsum / L : 0).toFixed(0)}%` +
          (cls ? `<br>class <b>${cls}</b>` : `<br>${e - s} windows aggregated`) +
          (t.ok[s] ? "" : `<br><i>not counted (too few expected reads)</i>`),
        ev,
      );
    }).on("mouseleave", () => hideTip());
    // legend
    const lg = $("at-heat-legend");
    if (lg)
      lg.innerHTML =
        `<span class="at-swatch" style="background:linear-gradient(90deg,${color(-3)},${color(0)},${color(3)})"></span> log2(depth / sample mean), −3 … +3` +
        ` &nbsp; <span class="at-swatch" style="background:#212121;width:14px"></span> zero` +
        ` &nbsp; <span class="at-swatch" style="background:#b0bec5;opacity:.5;width:14px"></span> faded = not counted` +
        ` &nbsp; * = sample below the expected-reads cutoff`;
  }

  function _contigDividers(g, res, cols, cw, h) {
    const { a, b } = _cols(res);
    const inView = res.layout.filter((c) => c.off > a && c.off < b);
    if (inView.length > 150) return;
    const colOf = (wi) => {
      for (let ci = 0; ci < cols.length; ci++) if (wi < cols[ci][1]) return ci;
      return cols.length;
    };
    inView.forEach((c) => {
      const x = colOf(c.off) * cw;
      g.append("line").attr("x1", x).attr("x2", x).attr("y1", -4).attr("y2", h).attr("stroke", "#78909c").attr("stroke-dasharray", "2,2").attr("stroke-width", 0.7).append("title").text(c.acc);
    });
  }

  function _drawFreq(res, o) {
    const host = $("at-freq");
    if (!host || typeof d3 === "undefined") return;
    host.innerHTML = "";
    const { cols } = _cols(res);
    const W = Math.max(600, host.clientWidth || 900);
    const H = 190;
    const iw = W - HEAT_ML - 14;
    const ih = H - 40;
    const cw = iw / cols.length;
    const svg = d3.select(host).append("svg").attr("width", W).attr("height", H).attr("font-size", 10);
    const g = svg.append("g").attr("transform", `translate(${HEAT_ML},10)`);
    const y = d3.scaleLinear().domain([0, 1]).range([ih, 0]);
    g.append("g").call(d3.axisLeft(y).ticks(4).tickFormat(d3.format(".0%")));
    svg.append("text").attr("x", HEAT_ML - 40).attr("y", 10 + ih / 2).attr("text-anchor", "end").attr("fill", "#546e7a").text("of samples");
    // recurrent region bands
    const { a, b } = _cols(res);
    const colOf = (wi) => {
      for (let ci = 0; ci < cols.length; ci++) if (wi < cols[ci][1]) return ci;
      return cols.length;
    };
    res.regions
      .filter((r) => r.type !== "zero" && r.w1 > a && r.w0 < b)
      .forEach((r) => {
        const x0 = colOf(Math.max(a, r.w0)) * cw;
        const x1 = (colOf(Math.min(b, r.w1) - 1) + 1) * cw;
        g.append("rect").attr("x", x0).attr("y", 0).attr("width", Math.max(1, x1 - x0)).attr("height", ih).attr("fill", r.type === "high" ? "#ffcdd2" : "#bbdefb").attr("opacity", 0.55);
      });
    // grey = no window in the column has enough counted samples
    cols.forEach(([s, e], ci) => {
      for (let i = s; i < e; i++) if (!isNaN(res.fl[i])) return;
      g.append("rect").attr("x", ci * cw).attr("y", 0).attr("width", Math.max(0.6, cw + 0.15)).attr("height", ih).attr("fill", "#eceff1");
    });
    g.append("line").attr("x1", 0).attr("x2", iw).attr("y1", y(o.minFreq)).attr("y2", y(o.minFreq)).attr("stroke", "#9e9e9e").attr("stroke-dasharray", "4,3");
    // series: max over each display column (a recurrent window must not vanish when zoomed out)
    const series = (arr) =>
      cols.map(([s, e], ci) => {
        let m = NaN;
        for (let i = s; i < e; i++) if (!isNaN(arr[i])) m = isNaN(m) ? arr[i] : Math.max(m, arr[i]);
        return [ci * cw + cw / 2, m];
      });
    const line = d3
      .line()
      .defined((d) => !isNaN(d[1]))
      .x((d) => d[0])
      .y((d) => y(d[1]));
    [
      [res.fl, "#1e88e5", 1.6, null, "low or zero"],
      [res.fz, "#212121", 1, "4,2", "zero"],
      [res.fh, "#e53935", 1.4, null, "high"],
    ].forEach(([arr, c, sw, da, lbl]) => {
      g.append("path").datum(series(arr)).attr("fill", "none").attr("stroke", c).attr("stroke-width", sw).attr("stroke-dasharray", da).attr("d", line).append("title").text(lbl);
    });
    // x axis in bp of the concatenated reference
    const sx = d3.scaleLinear().domain([res.wStart[a] || 0, b >= res.nwin ? res.totalLen : res.wStart[b]]).range([0, iw]);
    g.append("g").attr("transform", `translate(0,${ih})`).call(d3.axisBottom(sx).ticks(8).tickFormat((v) => fmtBp(Math.round(v))));
    _contigDividers(g, res, cols, cw, ih);
    if (!res.analysed) {
      g.append("text").attr("x", iw / 2).attr("y", ih / 2).attr("text-anchor", "middle").attr("fill", "#c62828").attr("font-size", 12).text(`Needs ≥ ${o.minSamples} informative samples — lower "min reads / window" or add samples`);
    }
    // brush to zoom
    const brush = d3
      .brushX()
      .extent([
        [0, 0],
        [iw, ih],
      ])
      .on("end", (ev) => {
        if (!ev.selection) return;
        const [x0, x1] = ev.selection;
        const c0 = Math.max(0, Math.floor(x0 / cw));
        const c1 = Math.min(cols.length - 1, Math.floor(x1 / cw));
        if (c1 <= c0 && cols[c0][1] - cols[c0][0] <= 1) return;
        AT.zoom = [cols[c0][0], cols[c1][1]];
        _drawHeat(res, o);
        _drawFreq(res, o);
        _zoomNote(res);
      });
    g.append("g").attr("class", "at-brush").call(brush);
    _zoomNote(res);
  }
  function _zoomNote(res) {
    const el = $("at-zoom-note");
    if (!el) return;
    if (!AT.zoom) {
      el.innerHTML = "Drag across the track to zoom; click a region row to jump to it.";
      if ($("at-reset-zoom")) $("at-reset-zoom").style.display = "none";
      return;
    }
    const a = _winLabel(res, AT.zoom[0]),
      b = _winLabel(res, AT.zoom[1] - 1);
    el.innerHTML = `Zoomed: <b>${esc(a.contig)}:${fmtInt(a.start)}</b> – <b>${esc(b.contig)}:${fmtInt(b.end)}</b>`;
    if ($("at-reset-zoom")) $("at-reset-zoom").style.display = "";
  }

  function _drawRegions(res, o) {
    const host = $("at-regions");
    if (!host) return;
    const only = ($("at-region-type") || { value: "" }).value;
    const rows = res.regions
      .filter((r) => !only || r.type === only)
      .sort((x, y) => y.length * y.meanFreq - x.length * x.meanFreq);
    const tcol = { zero: "#212121", low: "#1e88e5", high: "#e53935" };
    const shown = rows.slice(0, 500);
    host.innerHTML =
      `<table class="at-table"><thead><tr><th>Type</th><th>Contig</th><th>Start</th><th>End</th><th>Length</th>` +
      `<th>Mean freq</th><th>Max freq</th><th>Samples</th><th>Mean depth / sample mean</th><th>Affected samples</th></tr></thead><tbody>` +
      (shown.length
        ? shown
            .map(
              (r, i) =>
                `<tr data-i="${res.regions.indexOf(r)}"><td><span class="at-pill" style="background:${tcol[r.type]}">${r.type}</span>${r.wholeContig ? ' <span class="at-tag" title="Spans the whole contig">whole contig</span>' : ""}</td>` +
                `<td>${esc(r.contig)}</td><td>${fmtInt(r.start)}</td><td>${fmtInt(r.end)}</td><td>${fmtBp(r.length)}</td>` +
                `<td>${(100 * r.meanFreq).toFixed(0)}%</td><td>${(100 * r.maxFreq).toFixed(0)}%</td><td>${r.nSamples}</td>` +
                `<td>${isNaN(r.meanNorm) ? "–" : r.meanNorm.toFixed(2)}</td><td class="at-aff" title="${esc(r.affected.join(", "))}">${esc(r.affected.join(", "))}</td></tr>`,
            )
            .join("")
        : `<tr><td colspan="10" style="color:#78909c">No recurrent regions at these settings.</td></tr>`) +
      `</tbody></table>` +
      (rows.length > shown.length ? `<div class="at-note">Showing the 500 largest of ${fmtInt(rows.length)} regions — download for all.</div>` : "");
    host.querySelectorAll("tbody tr[data-i]").forEach((tr) =>
      tr.addEventListener("click", () => {
        const r = res.regions[+tr.dataset.i];
        const pad = Math.max(2, Math.round((r.w1 - r.w0) * 0.5));
        AT.zoom = [Math.max(0, r.w0 - pad), Math.min(res.nwin, r.w1 + pad)];
        _drawHeat(res, o);
        _drawFreq(res, o);
        const hw = $("at-heat");
        if (hw && hw.scrollIntoView) hw.scrollIntoView({ behavior: "smooth", block: "center" });
      }),
    );
    const sel = $("at-region-type");
    if (sel && !sel.__atWired) {
      sel.__atWired = true;
      sel.addEventListener("change", () => _lastRes && _drawRegions(_lastRes, _opts()));
    }
  }

  function _drawSampleTable(res, o) {
    const host = $("at-samples");
    if (!host) return;
    const t = _sortedTracks(res, o);
    host.innerHTML =
      `<table class="at-table"><thead><tr><th>Sample</th><th>Reads</th><th>Mean depth</th><th>Breadth</th><th>Exp. reads / window</th>` +
      `<th>Counted</th><th>% zero windows</th><th>% low</th><th>% high</th></tr></thead><tbody>` +
      t
        .map(
          (x) =>
            `<tr><td>${esc(x.sample)}</td><td>${fmtInt(x.reads)}</td><td>${x.mean.toFixed(2)}x</td><td>${x.breadth.toFixed(1)}%</td>` +
            `<td>${x.expRpw.toFixed(1)}</td><td>${x.informative ? "yes" : '<span style="color:#90a4ae">no</span>'}</td>` +
            `<td>${x.pct[1].toFixed(1)}</td><td>${x.pct[2].toFixed(1)}</td><td>${x.pct[3].toFixed(1)}</td></tr>`,
        )
        .join("") +
      `</tbody></table>`;
  }

  function _drawOverview(overview, o) {
    const host = $("at-overview");
    if (!host) return;
    host.innerHTML =
      `<table class="at-table at-clickable"><thead><tr><th>Reference</th><th>Samples (counted / total)</th><th>Size</th><th>Window</th>` +
      `<th>Recurrent zero</th><th>Recurrent low/zero</th><th>Recurrent high</th></tr></thead><tbody>` +
      overview
        .map(
          (r) =>
            `<tr data-k="${esc(r.key)}" class="${r.key === AT.selected ? "at-sel" : ""}${r.analysed ? "" : " at-dim"}"><td>${esc(r.name)}</td><td>${r.nInf} / ${r.nSamples}</td>` +
            `<td>${fmtBp(r.totalLen)}${r.nContigs > 1 ? ` · ${fmtInt(r.nContigs)} contigs` : ""}</td><td>${fmtBp(r.W)}</td>` +
            (r.analysed
              ? `<td>${r.pct.zero.toFixed(1)}% <span class="at-n">(${r.nReg.zero})</span></td><td>${r.pct.low.toFixed(1)}% <span class="at-n">(${r.nReg.low})</span></td><td>${r.pct.high.toFixed(1)}% <span class="at-n">(${r.nReg.high})</span></td>`
              : `<td colspan="3" style="color:#90a4ae">needs ≥ ${o.minSamples} counted samples</td>`) +
            `</tr>`,
        )
        .join("") +
      `</tbody></table>`;
    host.querySelectorAll("tbody tr[data-k]").forEach((tr) =>
      tr.addEventListener("click", () => {
        AT.selected = tr.dataset.k;
        AT.zoom = null;
        drawAlignTrends();
        const h = $("at-kpis");
        if (h && h.scrollIntoView) h.scrollIntoView({ behavior: "smooth", block: "start" });
      }),
    );
  }

  /* ── downloads (same columns as bin/alignment_trends.py) ──────────────── */
  function _download(which) {
    const res = _lastRes;
    const tsv = (cols, rows) => [cols.join("\t"), ...rows.map((r) => cols.map((c) => (r[c] == null ? "" : r[c])).join("\t"))].join("\n");
    const slug = (s) => String(s).replace(/[^A-Za-z0-9_.-]+/g, "_").slice(0, 60);
    if (which === "overview") {
      const rows = _lastOverview.map((r) => ({
        key: r.key,
        organism: r.name,
        n_samples: r.nSamples,
        n_informative: r.nInf,
        total_len: r.totalLen,
        n_contigs: r.nContigs,
        window: r.W,
        zero_pct: r.analysed ? r.pct.zero.toFixed(3) : "",
        n_zero_regions: r.analysed ? r.nReg.zero : "",
        low_pct: r.analysed ? r.pct.low.toFixed(3) : "",
        n_low_regions: r.analysed ? r.nReg.low : "",
        high_pct: r.analysed ? r.pct.high.toFixed(3) : "",
        n_high_regions: r.analysed ? r.nReg.high : "",
      }));
      _downloadText(tsv(Object.keys(rows[0] || { key: 1 }), rows), "alignment_trends.summary.tsv", "text/tab-separated-values");
      return;
    }
    if (!res) return;
    if (which === "regions") {
      const rows = res.regions.map((r) => ({
        key: res.ref.key,
        organism: res.ref.name,
        type: r.type,
        contig: r.contig,
        start: r.start,
        end: r.end,
        length: r.length,
        n_windows: r.nWin,
        mean_freq: r.meanFreq.toFixed(4),
        max_freq: r.maxFreq.toFixed(4),
        n_samples: r.nSamples,
        mean_norm_depth: isNaN(r.meanNorm) ? "" : r.meanNorm.toFixed(4),
        whole_contig: r.wholeContig,
        samples_affected: r.affected.join(","),
      }));
      _downloadText(tsv(Object.keys(rows[0] || { type: 1 }), rows), `alignment_trends.${slug(res.ref.name)}.regions.tsv`, "text/tab-separated-values");
    } else {
      const rows = [];
      for (let i = 0; i < res.nwin; i++) {
        const w = _winLabel(res, i);
        const row = {
          contig: w.contig,
          start: w.start,
          end: w.end,
          n_samples: res.nS[i],
          freq_zero: isNaN(res.fz[i]) ? "" : res.fz[i].toFixed(4),
          freq_low: isNaN(res.fl[i]) ? "" : res.fl[i].toFixed(4),
          freq_high: isNaN(res.fh[i]) ? "" : res.fh[i].toFixed(4),
          mean_norm_depth: isNaN(res.mNorm[i]) ? "" : res.mNorm[i].toFixed(4),
        };
        res.tracks.forEach((t) => (row[`depth:${t.sample}`] = t.dep[i].toFixed(3)));
        rows.push(row);
      }
      _downloadText(tsv(Object.keys(rows[0] || { contig: 1 }), rows), `alignment_trends.${slug(res.ref.name)}.windows.tsv`, "text/tab-separated-values");
    }
  }

  /* ══ Pipeline results view (BOOT.alignment_trends, from --alignment_trends) ══
     make_report.py embeds all.alignment_trends.json (minus the per-window
     table). This view reports the pipeline's comparison as run; clicking a
     reference opens it in the live view with the pipeline's settings applied. */
  const ATP = { view: "live", wired: false };
  const _atpData = () => (typeof BOOT !== "undefined" && BOOT && BOOT.alignment_trends) || null;

  function _atpInit() {
    const data = _atpData();
    const tog = $("at-view-toggle");
    if (tog) tog.style.display = data ? "" : "none";
    if (!data || ATP.wired) return;
    ATP.wired = true;
    document.querySelectorAll("#at-view-toggle .at-view-btn").forEach((b) =>
      b.addEventListener("click", () => _atpSetView(b.dataset.atview)),
    );
    ["atp-region-type", "atp-region-ref"].forEach((id) => {
      const el = $(id);
      if (el) el.addEventListener("change", () => (_atpRegions(), _atpSamples()));
    });
    const srch = $("atp-ref-search");
    if (srch) srch.addEventListener("input", () => _atpRefs());
    if ($("atp-dl-regions"))
      $("atp-dl-regions").addEventListener("click", (e) => {
        e.preventDefault();
        _atpDownload("regions");
      });
    if ($("atp-dl-samples"))
      $("atp-dl-samples").addEventListener("click", (e) => {
        e.preventDefault();
        _atpDownload("samples");
      });
    // Reference picker for the regions / samples tables
    const rs = $("atp-region-ref");
    if (rs)
      (data.references || []).forEach((r) => {
        const o = document.createElement("option");
        o.value = String(r.key);
        o.textContent = r.organism;
        rs.appendChild(o);
      });
  }

  function _atpSetView(v) {
    ATP.view = v === "pipeline" && _atpData() ? "pipeline" : "live";
    document
      .querySelectorAll("#at-view-toggle .at-view-btn")
      .forEach((b) => b.classList.toggle("active", b.dataset.atview === ATP.view));
    if ($("at-view-live")) $("at-view-live").style.display = ATP.view === "live" ? "" : "none";
    if ($("at-view-pipeline")) $("at-view-pipeline").style.display = ATP.view === "pipeline" ? "" : "none";
    if (ATP.view === "pipeline") drawAlignTrendsPipeline();
    else drawAlignTrends();
  }

  function drawAlignTrendsPipeline() {
    const data = _atpData();
    if (!data) return;
    const p = data.params || {};
    const refs = data.references || [];
    const regs = data.regions || [];
    const meta = $("atp-meta");
    if (meta) {
      const cut = (abs, frac, sym) => (abs != null ? `${sym} ${abs}x (absolute)` : `${sym} ${frac}× sample mean`);
      meta.innerHTML =
        `Computed by <code>bin/alignment_trends.py</code> during the run over <b>${fmtInt(data.n_samples)}</b> real sample(s) ` +
        `(controls and simulated datasets excluded; hidden samples and report filters do <i>not</i> apply here). ` +
        `Settings: low ${cut(p.low_abs, p.low_frac, "<")}, high ${cut(p.high_abs, p.high_frac, ">")}, ` +
        `≥ ${p.min_reads_per_window} expected reads / window, ≥ ${p.min_reads} reads, ≥ ${p.min_samples} counted samples, ` +
        `recurrence ≥ ${Math.round(100 * (p.min_freq || 0))}%` +
        (p.min_region_windows > 1 ? `, ≥ ${p.min_region_windows} windows per region` : "") +
        `. Source: <code>${esc(data.source || "all.alignment_trends.json")}</code>.` +
        (data.n_regions > regs.length ? ` <b>${fmtInt(regs.length)}</b> of ${fmtInt(data.n_regions)} regions embedded — the TSV on disk has all.` : "");
    }
    const k = $("atp-kpis");
    if (k) {
      const n = (t) => regs.filter((r) => r.type === t).length;
      const card = (label, value, sub, color) =>
        `<div class="kpi-card" style="border-left-color:${color}"><div class="kpi-label">${label}</div>` +
        `<div class="kpi-value">${value}</div><div class="kpi-sub">${sub}</div></div>`;
      k.innerHTML =
        card("References compared", fmtInt(refs.length), `≥ ${p.min_samples} counted samples each`, "#1565c0") +
        card("Samples", fmtInt(data.n_samples), "in the run", "#455a64") +
        card("Zero regions", fmtInt(n("zero")), "recurrent gaps", "#212121") +
        card("Low / zero regions", fmtInt(n("low")), "recurrent thin coverage", "#1e88e5") +
        card("High regions", fmtInt(n("high")), "recurrent pile-ups", "#e53935");
    }
    _atpRefs();
    _atpRegions();
    _atpSamples();
  }

  function _atpRefs() {
    const data = _atpData();
    const host = $("atp-refs");
    if (!data || !host) return;
    const q = (($("atp-ref-search") || {}).value || "").trim().toLowerCase();
    const rows = (data.references || []).filter((r) => !q || String(r.organism).toLowerCase().includes(q) || String(r.key).includes(q));
    const pc = (v, n) => `${(+v || 0).toFixed(1)}% <span class="at-n">(${fmtInt(n)})</span>`;
    host.innerHTML =
      `<table class="at-table at-clickable"><thead><tr><th>Reference</th><th>Samples (counted / total)</th><th>Size</th><th>Window</th>` +
      `<th>Recurrent zero</th><th>Recurrent low/zero</th><th>Recurrent high</th><th>Counted samples</th></tr></thead><tbody>` +
      (rows.length
        ? rows
            .map(
              (r) =>
                `<tr data-k="${esc(r.key)}"><td>${esc(r.organism)}</td><td>${r.n_informative} / ${r.n_samples}</td>` +
                `<td>${fmtBp(r.total_len)}${r.n_contigs > 1 ? ` · ${fmtInt(r.n_contigs)} contigs` : ""}</td><td>${fmtBp(r.window)}</td>` +
                `<td>${pc(r.zero_pct, r.n_zero_regions)}</td><td>${pc(r.low_pct, r.n_low_regions)}</td><td>${pc(r.high_pct, r.n_high_regions)}</td>` +
                `<td class="at-aff" title="${esc((r.informative_samples || []).join(", "))}">${esc((r.informative_samples || []).join(", "))}</td></tr>`,
            )
            .join("")
        : `<tr><td colspan="8" style="color:#78909c">No reference had enough counted samples in the pipeline run.</td></tr>`) +
      `</tbody></table>`;
    host.querySelectorAll("tbody tr[data-k]").forEach((tr) => tr.addEventListener("click", () => _atpOpenLive(tr.dataset.k)));
  }

  function _atpFilteredRegions() {
    const data = _atpData();
    const t = ($("atp-region-type") || {}).value || "";
    const rk = ($("atp-region-ref") || {}).value || "";
    return (data.regions || []).filter((r) => (!t || r.type === t) && (!rk || String(r.key) === rk));
  }
  function _atpRegions() {
    const host = $("atp-regions");
    if (!host || !_atpData()) return;
    const rows = _atpFilteredRegions()
      .slice()
      .sort((a, b) => b.length * b.mean_freq - a.length * a.mean_freq);
    const shown = rows.slice(0, 1000);
    const tcol = { zero: "#212121", low: "#1e88e5", high: "#e53935" };
    host.innerHTML =
      `<table class="at-table at-clickable"><thead><tr><th>Reference</th><th>Type</th><th>Contig</th><th>Start</th><th>End</th><th>Length</th>` +
      `<th>Mean freq</th><th>Max freq</th><th>Samples</th><th>Depth / mean</th><th>Affected samples</th></tr></thead><tbody>` +
      (shown.length
        ? shown
            .map(
              (r) =>
                `<tr data-k="${esc(r.key)}" data-c="${esc(r.contig)}" data-s="${r.start}" data-e="${r.end}"><td>${esc(r.organism)}</td>` +
                `<td><span class="at-pill" style="background:${tcol[r.type] || "#607d8b"}">${esc(r.type)}</span>${r.whole_contig ? ' <span class="at-tag">whole contig</span>' : ""}</td>` +
                `<td>${esc(r.contig)}</td><td>${fmtInt(r.start)}</td><td>${fmtInt(r.end)}</td><td>${fmtBp(r.length)}</td>` +
                `<td>${(100 * r.mean_freq).toFixed(0)}%</td><td>${(100 * r.max_freq).toFixed(0)}%</td><td>${r.n_samples}</td>` +
                `<td>${r.mean_norm_depth == null ? "–" : (+r.mean_norm_depth).toFixed(2)}</td>` +
                `<td class="at-aff" title="${esc(r.samples_affected)}">${esc(String(r.samples_affected || "").replace(/,/g, ", "))}</td></tr>`,
            )
            .join("")
        : `<tr><td colspan="11" style="color:#78909c">No recurrent regions for this selection.</td></tr>`) +
      `</tbody></table>` +
      (rows.length > shown.length ? `<div class="at-note">Showing the 1,000 largest of ${fmtInt(rows.length)} — download for all.</div>` : "");
    host.querySelectorAll("tbody tr[data-k]").forEach((tr) =>
      tr.addEventListener("click", () => _atpOpenLive(tr.dataset.k, { contig: tr.dataset.c, start: +tr.dataset.s, end: +tr.dataset.e })),
    );
  }
  function _atpSamples() {
    const data = _atpData();
    const host = $("atp-samples");
    if (!host || !data) return;
    const rk = ($("atp-region-ref") || {}).value || "";
    const analysed = new Set((data.references || []).map((r) => String(r.key)));
    const rows = (data.samples || []).filter((r) => (rk ? String(r.key) === rk : analysed.has(String(r.key))));
    host.innerHTML =
      `<table class="at-table"><thead><tr><th>Reference</th><th>Sample</th><th>Reads</th><th>Mean depth</th><th>Breadth</th>` +
      `<th>Exp. reads / window</th><th>Counted</th><th>% zero</th><th>% low</th><th>% high</th></tr></thead><tbody>` +
      rows
        .slice(0, 2000)
        .map(
          (r) =>
            `<tr><td>${esc(r.organism)}</td><td>${esc(r.sample)}</td><td>${fmtInt(r.numreads)}</td><td>${(+r.mean_depth || 0).toFixed(2)}x</td>` +
            `<td>${(+r.breadth_pct || 0).toFixed(1)}%</td><td>${r.exp_reads_per_window == null ? "–" : (+r.exp_reads_per_window).toFixed(1)}</td>` +
            `<td>${r.informative ? "yes" : '<span style="color:#90a4ae">no</span>'}</td><td>${(+r.pct_zero || 0).toFixed(1)}</td>` +
            `<td>${(+r.pct_low || 0).toFixed(1)}</td><td>${(+r.pct_high || 0).toFixed(1)}</td></tr>`,
        )
        .join("") +
      `</tbody></table>`;
  }

  // Open a pipeline reference in the live view with the pipeline's settings
  // (filters off, so the sample set matches the run as closely as the report allows).
  function _atpOpenLive(key, region) {
    const p = (_atpData() || {}).params || {};
    const set = (id, v) => {
      const el = $(id);
      if (el && v != null && isFinite(v)) el.value = v;
    };
    set("at-low", p.low_frac);
    set("at-high", p.high_frac);
    set("at-minrpw", p.min_reads_per_window);
    set("at-recur", Math.round(100 * (p.min_freq || 0.5)));
    set("at-minsamples", p.min_samples);
    set("at-minreads", p.min_reads);
    if ($("at-use-filters")) $("at-use-filters").checked = false;
    if ($("at-coarsen")) $("at-coarsen").value = "1";
    AT.selected = String(key);
    AT.zoom = null;
    _atpSetView("live");
    const res = _lastRes;
    if (!res || String(res.ref.key) !== String(key)) {
      const el = $("at-zoom-note");
      if (el) el.innerHTML = '<span style="color:#c62828">That reference is not among the report’s visible samples.</span>';
      return;
    }
    if (region) {
      const c = res.layout.find((x) => x.acc === region.contig);
      if (c) {
        const w0 = c.off + Math.floor(region.start / res.W);
        const w1 = c.off + Math.ceil(region.end / res.W);
        const pad = Math.max(2, Math.round((w1 - w0) * 0.5));
        AT.zoom = [Math.max(0, w0 - pad), Math.min(res.nwin, w1 + pad)];
        _drawHeat(res, _opts());
        _drawFreq(res, _opts());
      }
    }
    const h = $("at-kpis");
    if (h && h.scrollIntoView) h.scrollIntoView({ behavior: "smooth", block: "start" });
  }

  function _atpDownload(which) {
    const data = _atpData();
    if (!data) return;
    const rows = which === "regions" ? _atpFilteredRegions() : data.samples || [];
    if (!rows.length) return;
    const cols = Object.keys(rows[0]);
    const txt = [cols.join("\t"), ...rows.map((r) => cols.map((c) => (r[c] == null ? "" : r[c])).join("\t"))].join("\n");
    _downloadText(txt, `all.alignment_trends.${which}.tsv`, "text/tab-separated-values");
  }

  /* ══ Hover help for column headers, KPI cards and controls ════════════════
     One table of explanations keyed by the visible label, applied by a
     MutationObserver so every (re)rendered table picks them up. Numbers in the
     text follow the settings of the view the element sits in: the live
     controls, or the pipeline run's parameters under "Pipeline results". */
  function _tipSettings(el) {
    const inPipe = el && el.closest && el.closest("#at-view-pipeline");
    if (inPipe) {
      const p = (_atpData() || {}).params || {};
      return {
        low: p.low_abs != null ? `${p.low_abs}x (absolute)` : `${p.low_frac}× the sample's mean depth`,
        high: p.high_abs != null ? `${p.high_abs}x (absolute)` : `${p.high_frac}× the sample's mean depth`,
        rpw: p.min_reads_per_window,
        freq: Math.round(100 * (p.min_freq || 0)),
        minS: p.min_samples,
        minR: p.min_reads,
      };
    }
    const o = _opts();
    return {
      low: `${o.lowFrac}× the sample's mean depth`,
      high: `${o.highFrac}× the sample's mean depth`,
      rpw: o.minRpw,
      freq: Math.round(100 * o.minFreq),
      minS: o.minSamples,
      minR: o.minReads,
    };
  }
  const _COUNTED = (t) =>
    `A sample is "counted" on a reference when it has ≥ ${t.minR} reads there AND expects ≥ ${t.rpw} reads per window ` +
    `(its reads × window ÷ reference length). Below that, empty or thin windows are mostly sampling noise.`;
  const _TIPS = {
    Reference: () =>
      "The reference (strain-level assembly) the reads were aligned to. Multi-contig assemblies are analysed contig by contig and drawn end to end.",
    "Samples (counted / total)": (t) => `Counted samples / samples with a depth profile for this reference. ${_COUNTED(t)}`,
    Size: () =>
      "Total length of every contig / accession in the reference, including contigs no sample covered.",
    Window: () =>
      "Window (bin) size. Each sample's depth and breadth are averaged over windows of this many bp. Chosen as the smallest 100 × 2^k bp " +
      "that gives at most --depth_profile_windows (default 400) windows over the reference, so every sample lines up on the same grid. " +
      '"Window ×" in the live view merges windows into coarser ones.',
    "Recurrent zero": (t) =>
      `Share of the reference (and, in brackets, the number of regions) where at least ${t.freq}% of counted samples have NO reads at all. ` +
      "Points at sequence the samples do not carry: deletions, divergent or novel loci, or assembly artefacts.",
    "Recurrent low/zero": (t) =>
      `Share of the reference (and number of regions) where at least ${t.freq}% of counted samples are below ${t.low} or have no reads. ` +
      "Includes the recurrent-zero regions plus thinly covered stretches (partial deletions, divergent sequence, GC / mappability dips).",
    "Recurrent high": (t) =>
      `Share of the reference (and number of regions) where at least ${t.freq}% of counted samples are above ${t.high}. ` +
      "Typical causes: repeats, rRNA operons, mobile elements / prophage, plasmid copy number, or reads from a related organism piling onto a conserved locus.",
    "Counted samples": (t) => `The samples that counted toward this reference's frequencies. ${_COUNTED(t)}`,
    Type: (t) =>
      `zero = no reads; low = below ${t.low} (or zero); high = above ${t.high}. ` +
      `A region is consecutive windows of one contig where that state recurs in ≥ ${t.freq}% of counted samples. ` +
      '"whole contig" = the region spans the entire contig.',
    Contig: () => "Contig / accession the region lies on.",
    Start: () => "Region start on the contig (bp, 0-based). Regions are window-aligned, so edges are accurate to one window.",
    End: () => "Region end on the contig (bp, exclusive).",
    Length: () => "Region length (bp).",
    "Mean freq": () =>
      "Average, over the region's windows, of the fraction of counted samples in this state. 100% = every counted sample, in every window.",
    "Max freq": () => "Highest fraction of counted samples in this state in any single window of the region.",
    Samples: () => "Counted samples in the region's windows (the most in any one window). Short contigs and contig ends can count fewer.",
    "Mean depth / sample mean": () =>
      "Average normalised depth across the region: window depth ÷ that sample's mean depth on the reference, averaged over counted samples. 1 = typical, 0 = empty, 3 = triple.",
    "Depth / mean": () =>
      "Average normalised depth across the region: window depth ÷ that sample's mean depth on the reference, averaged over counted samples. 1 = typical, 0 = empty, 3 = triple.",
    "Affected samples": () => "Counted samples that are in this state in at least one window of the region.",
    Sample: () => "Sample (library). Hidden samples are left out of the live view.",
    Reads: () => "Reads aligned to this reference in the sample.",
    "Mean depth": () => "Mean depth over the whole reference (all contigs), the baseline that window depths are divided by.",
    Breadth: () => "Percent of the reference covered by at least one read.",
    "Exp. reads / window": (t) =>
      `Reads this sample is expected to place in one window if they were spread evenly (reads × window ÷ reference length). Must be ≥ ${t.rpw} for the sample to count.`,
    Counted: (t) => _COUNTED(t),
    "% zero windows": () => "Share of this sample's windows with no reads.",
    "% zero": () => "Share of this sample's windows with no reads.",
    "% low": (t) => `Share of this sample's windows with reads but depth below ${t.low}.`,
    "% high": (t) => `Share of this sample's windows with depth above ${t.high}.`,
    // KPI cards
    "Informative samples": (t) => `Samples that count toward this reference's frequencies, out of all with a depth profile. ${_COUNTED(t)} Trends need ≥ ${t.minS}.`,
    "Recurrent low / zero": (t) => `Share of the reference in recurrent low-or-zero regions: ≥ ${t.freq}% of counted samples below ${t.low} or empty.`,
    "References compared": (t) => `References with at least ${t.minS} counted samples in the pipeline run.`,
    "Zero regions": (t) => `Recurrent regions where ≥ ${t.freq}% of counted samples have no reads.`,
    "Low / zero regions": (t) => `Recurrent regions where ≥ ${t.freq}% of counted samples are below ${t.low} or empty.`,
    "High regions": (t) => `Recurrent regions where ≥ ${t.freq}% of counted samples are above ${t.high}.`,
  };
  // KPI cards whose label means something different from the same-named column
  const _KPI_TIPS = {
    Samples: () => "Real samples in the pipeline run (controls and simulated datasets excluded).",
    Reference: () => "Total reference length, number of contigs and the window grid used for this reference.",
  };
  // Live-control labels (matched by the input they precede)
  const _CTRL_TIPS = {
    "at-ref-sel": () =>
      "Reference to show. Listed with counted / total samples, references with enough counted samples first. The table at the bottom lists them all.",
    "at-low": (t) => `LOW cutoff: a window is low when its depth is below this multiple of the sample's own mean depth on the reference (now ${t.low}).`,
    "at-high": (t) => `HIGH cutoff: a window is high when its depth is above this multiple of the sample's mean depth (now ${t.high}).`,
    "at-minrpw": () =>
      "Minimum expected reads per window (reads × window ÷ reference length) for a sample to count — per window, so short contigs and contig ends need the same evidence.",
    "at-recur": () =>
      "Recurrence ≥: the fraction of counted samples that must share a state (zero / low / high) in a window for it to join a recurrent region.",
    "at-minsamples": () => "Counted samples needed before a reference — and each window — is analysed.",
    "at-minreads": () => "Minimum reads on the reference for a sample to count at all.",
    "at-coarsen": () => "Merge 2, 4, 8 or 16 neighbouring windows before analysis: smoother, and lets shallow samples count.",
    "at-sort": () => "Row order in the heatmap and sample table. Counted samples always come first.",
    "at-use-filters": () =>
      "Only use sample × reference pairs whose detection passes the report's active filters and TASS cutoff. Off = every visible sample with a depth profile (closest to the pipeline run).",
  };
  function _applyTips(root) {
    if (!root) return;
    root.querySelectorAll("th, .kpi-label").forEach((el) => {
      const key = el.textContent.replace(/\s+/g, " ").trim();
      const kpi = el.classList.contains("kpi-label");
      const name = kpi && _KPI_TIPS[key] ? null : Object.keys(_TIPS).find((k) => k.toLowerCase() === key.toLowerCase());
      const fn = kpi && _KPI_TIPS[key] ? _KPI_TIPS[key] : name ? _TIPS[name] : null;
      if (!fn) return;
      const txt = fn(_tipSettings(el));
      _setTip(el, txt);
    });
    Object.entries(_CTRL_TIPS).forEach(([id, fn]) => {
      const input = $(id);
      if (!input) return;
      const txt = fn(_tipSettings(input));
      _setTip(input, txt);
      // the label just before the input (or wrapping it)
      let lab = input.closest("label");
      if (!lab) {
        let p = input.previousElementSibling;
        while (p && p.tagName !== "LABEL") p = p.previousElementSibling;
        lab = p;
      }
      if (lab) _setTip(lab, txt);
    });
  }
  // Native title tooltips are unreliable here (slow to appear, and some
  // embedded viewers never show them), so the text is kept in data-attip and
  // shown by our own hover box. The title attribute is removed so the two
  // never stack.
  function _setTip(el, txt) {
    if (el.dataset.attip !== txt) el.dataset.attip = txt;
    if (el.hasAttribute("title")) el.removeAttribute("title");
    el.setAttribute("aria-label", (el.textContent || "").trim() ? `${el.textContent.trim()}: ${txt}` : txt);
    el.classList.add("at-tip");
  }
  function _tipBox() {
    let b = document.getElementById("at-tipbox");
    if (!b) {
      b = document.createElement("div");
      b.id = "at-tipbox";
      b.setAttribute("role", "tooltip");
      document.body.appendChild(b);
    }
    return b;
  }
  function _placeTip(b, ev) {
    const pad = 14,
      vw = window.innerWidth,
      vh = window.innerHeight;
    const r = b.getBoundingClientRect();
    let x = ev.clientX + pad,
      y = ev.clientY + pad;
    if (x + r.width > vw - 8) x = Math.max(8, ev.clientX - r.width - pad);
    if (y + r.height > vh - 8) y = Math.max(8, ev.clientY - r.height - pad);
    b.style.left = x + "px";
    b.style.top = y + "px";
  }
  function _wireTipHover(root) {
    let cur = null;
    root.addEventListener("mouseover", (ev) => {
      const el = ev.target.closest && ev.target.closest("[data-attip]");
      if (!el || !root.contains(el)) return;
      cur = el;
      const b = _tipBox();
      b.textContent = el.dataset.attip;
      b.style.display = "block";
      _placeTip(b, ev);
    });
    root.addEventListener("mousemove", (ev) => {
      if (!cur) return;
      _placeTip(_tipBox(), ev);
    });
    root.addEventListener("mouseout", (ev) => {
      if (!cur) return;
      const to = ev.relatedTarget;
      if (to && cur.contains(to)) return;
      cur = null;
      _tipBox().style.display = "none";
    });
    // keyboard users: show on focus of a control that carries help
    root.addEventListener("focusin", (ev) => {
      const el = ev.target.closest && ev.target.closest("[data-attip]");
      if (!el) return;
      const b = _tipBox();
      b.textContent = el.dataset.attip;
      b.style.display = "block";
      const r = el.getBoundingClientRect();
      _placeTip(b, { clientX: r.left, clientY: r.bottom });
    });
    root.addEventListener("focusout", () => (_tipBox().style.display = "none"));
  }
  (function _watchTips() {
    const root = document.getElementById("meta-subpane-align");
    if (!root || typeof MutationObserver === "undefined") return;
    let pending = false;
    new MutationObserver(() => {
      if (pending) return;
      pending = true;
      requestAnimationFrame(() => {
        pending = false;
        _applyTips(root);
      });
    }).observe(root, { childList: true, subtree: true });
    _wireTipHover(root);
    _applyTips(root);
  })();

  window.drawAlignTrends = drawAlignTrends;
  window.drawAlignTrendsPipeline = drawAlignTrendsPipeline;
  window._alignTrendsAnalyse = _analyse; // exposed for tests / console use
})();
