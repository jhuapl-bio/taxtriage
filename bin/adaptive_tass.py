#!/usr/bin/env python3
"""
adaptive_tass.py -- opt-in, per-organism adaptive TASS weighting.

Used by match_paths.py when ``--adaptive_tass`` is given. With the flag off this
module is never imported and scoring is exactly what it was before.

The standard TASS (``tass_score``) is ALWAYS computed and stays the score that
drives cutoffs, filtering and the report. The adaptive score is an additional,
informational score (``tass_score_adaptive``) written next to it.

What it does
------------
The static scorer applies ONE weight vector (resolved per sample type /
platform from sampletype_best_thresholds.json, or from the CLI) to every
organism in a sample. The adaptive scorer starts from that same vector and,
for each organism, shifts it according to the context the organism sits in:

  term           context                                    shift
  -------------  -----------------------------------------  ------------------------------
  complexity     sample richness x evenness (species level)  + minhash, + gini
  ani            a DIFFERENT species with ANI >= cutoff      + minhash, - breadth, + gini
                 is also detected                            and an ambiguity penalty
  sibling        other species of the same genus detected    + minhash, - breadth (strain)
                                                             and (smaller) ambiguity penalty
  low_abundance  < 3 % of the reads at its level             + minhash, + gini
  rollup         species / genus built from >= 2 children    + breadth, + gini, - minhash

The shifted weights are re-normalised so the core pool keeps its original
total, then the organism is scored with the normal compute_tass_score() and
multiplied by an ambiguity penalty

    penalty = 1 - max_penalty * strength * ambiguity * (1 - minhash_reduction)

so a high-ANI / sibling call is only marked down when its own unique (MinHash)
support is weak.

Formulas and magnitudes are ported from the adaptive_tass_weighting prototype
(adaptive_weighting.py). Two deliberate differences:

* "Competitors" are high-ANI organisms of a DIFFERENT species, and "siblings"
  are other species of the same genus. Strains of the same species are never
  counted against each other: which strain is present is an assignment
  question, already answered (without touching TASS) by
  compute_assignment_confidence. Counting them here would mark down a real
  species whose reads are split over several reference assemblies.
* The prototype's rank penalty is not applied for the same reason.

The prototype's adaptive cutoff is computed but only REPORTED
(``cutoff_shift``, in TASS points); the report's per-sample cutoff is not
moved. It is there so the cutoff behaviour can be calibrated before anyone
relies on it.

Safety
------
* Configuration is validated up front (AdaptiveConfig.validate) and a bad
  value stops the run instead of being silently ignored.
* strength=0 reproduces the static score exactly (tested).
* The whole sample falls back to static weights when the MinHash comparison
  signal is missing (--fast / no reference comparison), because every shift
  and the penalty lean on it.
* Any error while adapting one organism falls back to that organism's static
  score; the error is counted in metadata.adaptive_tass.errors and logged.
* Every organism carries an ``adaptive_tass`` audit block (standard and
  adaptive score, the weights used, the context and the reasons) and
  ``tass_score_adaptive`` next to the unchanged standard ``tass_score``.

CLI
---
    python adaptive_tass.py summarize <sample>.paths.json [--level strain|species|genus|all]

prints standard vs adaptive TASS, the change, and the reasons for one output JSON.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

TERMS = ("complexity", "ani", "sibling", "low_abundance", "rollup")
BASES = ("sampletype", "platform", "global")
LEVELS = ("strain", "species", "genus")

# The normalised weight pool that the adaptive shifts act on (keys of the
# match_paths `weights` dict). k2 / diamond / siblings share the pool but are
# left untouched; additive terms (plasmid bonus, abundance confidence) and
# score_power are outside the pool.
CORE_KEYS = ("breadth_weight", "minhash_weight", "gini_weight",
             "disparity_weight", "hmp_weight", "mapq_score")

# Context thresholds (prototype values).
COMPLEXITY_ONSET = 0.45        # complexity above this shifts weights
RICHNESS_SATURATION = 12       # richness at which the richness scale reaches 1
LOW_ABUNDANCE = 0.03           # relative abundance below this is "low"
SIBLING_ONSET = 0.25           # sibling pressure above this shifts strain weights
AMBIGUITY_ANI_SHARE = 0.7      # ambiguity = 0.7*ANI part + 0.3*sibling part
AMBIGUITY_SIBLING_SHARE = 0.3

# Advisory cutoff shift (prototype adaptive_threshold, minus its base).
CUT_COMPLEXITY = 0.08
CUT_AMBIGUITY = 0.10
CUT_UNIQUE_CREDIT = 0.04
CUT_ABUNDANCE_CREDIT = 0.04
CUT_ROLLUP_CREDIT = 0.05

_EPS = 1e-12


def _clamp(v, lo=0.0, hi=1.0):
    return min(hi, max(lo, v))


def _num(v, default=0.0):
    """float(v), or default for None / '' / NaN / inf / unparsable."""
    try:
        f = float(v)
    except (TypeError, ValueError):
        return default
    return f if math.isfinite(f) else default


@dataclass(frozen=True)
class AdaptiveConfig:
    strength: float = 1.0          # 0 = static scores, 1 = full prototype shifts
    terms: tuple = TERMS           # which context terms are active
    high_ani: float = 0.95         # ANI (0-1) at/above which two organisms compete
    max_penalty: float = 0.35      # largest multiplicative mark-down (strength 1)
    base: str = "sampletype"       # starting weights: sampletype | platform | global

    @classmethod
    def from_args(cls, args):
        terms = parse_terms(getattr(args, "adaptive_tass_terms", None))
        hi = getattr(args, "adaptive_tass_high_ani", None)
        if hi is None:
            hi = getattr(args, "ani_threshold", None)
        cfg = cls(
            strength=_num(getattr(args, "adaptive_tass_strength", 1.0), float("nan")),
            terms=terms,
            high_ani=_num(hi if hi is not None else 0.95, float("nan")),
            max_penalty=_num(getattr(args, "adaptive_tass_max_penalty", 0.35), float("nan")),
            base=str(getattr(args, "adaptive_tass_base", "sampletype") or "sampletype").strip().lower(),
        )
        cfg.validate()
        return cfg

    def validate(self):
        errs = []
        if not (isinstance(self.strength, (int, float)) and 0.0 <= self.strength <= 1.0):
            errs.append(f"strength must be within [0, 1] (got {self.strength!r})")
        if not (isinstance(self.high_ani, (int, float)) and 0.5 <= self.high_ani < 1.0):
            errs.append(f"high-ANI cutoff must be within [0.5, 1) as a fraction (got {self.high_ani!r})")
        if not (isinstance(self.max_penalty, (int, float)) and 0.0 <= self.max_penalty <= 0.9):
            errs.append(f"max_penalty must be within [0, 0.9] (got {self.max_penalty!r})")
        if self.base not in BASES:
            errs.append(f"base must be one of {', '.join(BASES)} (got {self.base!r})")
        bad = [t for t in self.terms if t not in TERMS]
        if bad:
            errs.append(f"unknown term(s) {bad}; valid: {', '.join(TERMS)}")
        if errs:
            raise ValueError("--adaptive_tass: " + "; ".join(errs))
        return self

    def on(self, term):
        return term in self.terms


def parse_terms(value):
    """'complexity,ani' -> ('complexity', 'ani'). None / '' / 'all' -> every
    term; 'none' -> no term (penalty and shifts off, audit fields still written).
    Unknown names raise ValueError."""
    if value is None:
        return TERMS
    if isinstance(value, (list, tuple)):
        items = [str(v).strip().lower() for v in value]
    else:
        items = [t.strip().lower() for t in str(value).replace(";", ",").split(",")]
    items = [t for t in items if t]
    if not items or items == ["all"]:
        return TERMS
    if items == ["none"]:
        return ()
    bad = [t for t in items if t not in TERMS]
    if bad:
        raise ValueError(f"--adaptive_tass_terms: unknown term(s) {bad}; valid: {', '.join(TERMS)}, all, none")
    return tuple(t for t in TERMS if t in items)  # canonical order, de-duplicated


def sample_complexity(abundances):
    """Richness / Shannon / evenness / dominance and the bounded complexity
    score (evenness x saturating richness scale) of one sample."""
    vals = [max(0.0, _num(a)) for a in abundances]
    total = sum(vals)
    rels = [v / total for v in vals if total > 0 and v > 0]
    richness = len(rels)
    shannon = -sum(r * math.log(r) for r in rels) if rels else 0.0
    evenness = shannon / math.log(richness) if richness > 1 else 0.0
    dominance = max(rels) if rels else 0.0
    richness_scale = _clamp(math.log1p(richness) / math.log1p(RICHNESS_SATURATION))
    return {
        "richness": richness,
        "shannon": shannon,
        "evenness": evenness,
        "dominance": dominance,
        "complexity": _clamp(evenness * richness_scale),
    }


def adapt_weights(base, ctx, cfg, level):
    """Shift the core weights of `base` for one organism.

    Returns (weights, reasons). `weights` is a new dict: every non-core key is
    copied unchanged; core keys are shifted (scaled by cfg.strength and by the
    pool's own total so the magnitudes match the prototype's unit-sum vector),
    clamped at 0 and re-normalised to the original core total."""
    s = cfg.strength
    w = dict(base)
    core_total = sum(max(0.0, _num(base.get(k))) for k in CORE_KEYS)
    reasons = []
    if s <= 0 or core_total <= _EPS:
        return w, reasons

    delta = {k: 0.0 for k in CORE_KEYS}
    acw_cut = 0.0
    c = ctx["complexity"]
    if cfg.on("complexity") and c > COMPLEXITY_ONSET:
        sc = (c - COMPLEXITY_ONSET) / (1.0 - COMPLEXITY_ONSET)
        delta["minhash_weight"] += 0.14 * sc
        delta["gini_weight"] += 0.08 * sc
        acw_cut += 0.03 * sc
        reasons.append(f"complex sample ({c:.2f}): more MinHash/Gini")
    if cfg.on("ani") and ctx["competitors"]:
        delta["minhash_weight"] += 0.40 * ctx["ambiguity"]
        delta["breadth_weight"] -= 0.14 * ctx["competitor_pressure"]
        delta["gini_weight"] += 0.06 * ctx["ambiguity"]
        reasons.append(f"{len(ctx['competitors'])} high-ANI competitor(s) "
                       f"(max ANI {ctx['max_ani']:.3f}): more MinHash, less breadth")
    if cfg.on("sibling") and level == "strain" and ctx["sibling_pressure"] > SIBLING_ONSET:
        sp = ctx["sibling_pressure"]
        delta["minhash_weight"] += 0.10 * sp
        delta["breadth_weight"] -= 0.05 * sp
        reasons.append(f"same-genus species share {sp:.0%} of the genus reads: more MinHash")
    if cfg.on("low_abundance") and ctx["relative_abundance"] < LOW_ABUNDANCE:
        delta["minhash_weight"] += 0.06
        delta["gini_weight"] += 0.04
        acw_cut += 0.04
        reasons.append(f"low relative abundance ({ctx['relative_abundance']:.2%}): more MinHash/Gini")
    if cfg.on("rollup") and level != "strain" and ctx["n_children"] >= 2:
        delta["breadth_weight"] += 0.08
        delta["gini_weight"] += 0.04
        delta["minhash_weight"] -= 0.05
        reasons.append(f"{level} roll-up of {ctx['n_children']} children: more breadth/Gini")

    if not any(abs(v) > _EPS for v in delta.values()) and acw_cut <= _EPS:
        return w, reasons

    shifted = {k: max(0.0, _num(base.get(k)) + s * core_total * delta[k]) for k in CORE_KEYS}
    tot = sum(shifted.values())
    if tot <= _EPS:  # every core weight clamped to 0: keep the static vector
        return dict(base), reasons + ["shift would zero every weight: static weights kept"]
    for k in CORE_KEYS:
        w[k] = shifted[k] * core_total / tot
    acw = _num(base.get("abundance_confidence_weight"))
    if acw > 0 and acw_cut > 0:
        w["abundance_confidence_weight"] = max(0.0, acw - s * acw_cut)
    return w, reasons


def check_weights(base, adapted):
    """Invariants every adapted vector must satisfy; raises AssertionError."""
    for k, v in adapted.items():
        if isinstance(v, bool) or not isinstance(v, (int, float)):
            assert adapted[k] == base.get(k), f"non-numeric weight {k} changed"
            continue
        assert math.isfinite(v), f"weight {k} is not finite"
        assert v >= 0.0 or base.get(k, 0) < 0, f"weight {k} is negative"
    for k in base:
        if k not in CORE_KEYS and k != "abundance_confidence_weight":
            assert adapted.get(k) == base.get(k), f"non-core weight {k} changed"
    b = sum(max(0.0, _num(base.get(k))) for k in CORE_KEYS)
    a = sum(_num(adapted.get(k)) for k in CORE_KEYS)
    assert abs(a - b) <= 1e-9 * max(1.0, b), f"core weight total changed ({b} -> {a})"


class AdaptiveTass:
    """Per-sample adaptive scorer. One instance per match_paths run.

    Usage (see match_paths.py):
        at = AdaptiveTass(cfg, ani_data)
        at.prepare('strain', strain_summary, strain_summary)
        for data in strain_summary.values():
            raw = at.score(data, 'strain', weights, compute_tass_score)
    """

    def __init__(self, config, ani_data=None, log=print,
                 species_field="subkey", genus_field="toplevelkey"):
        """species_field / genus_field: which organism fields hold the SPECIES
        and the GENUS (sibling group) taxid. With match_paths' default
        grouping (--rank genus, subrank species) that is subkey / toplevelkey;
        see species_fields_for_rank()."""
        self.cfg = config
        self.species_field = species_field
        self.genus_field = genus_field
        self.log = log
        self.ani = self._normalise_ani(ani_data or {})
        self._ctx = {}
        self.sample = None
        self.minhash_available = None
        self.disabled_reason = None
        self.errors = 0
        self._error_msgs = []
        self.stats = {lvl: {"n": 0, "changed": 0, "sum_delta": 0.0, "max_abs_delta": 0.0,
                            "raised": 0, "lowered": 0}
                      for lvl in LEVELS}

    # ── setup ───────────────────────────────────────────────────────────
    @staticmethod
    def _normalise_ani(ani_data):
        """Copy ani_data to {str: {str: fraction}}; percent matrices (>1) are
        divided by 100 so the cutoff is always a fraction."""
        out = {}
        pct = any(_num(v) > 1.0 for row in ani_data.values() for v in (row or {}).values())
        for a, row in ani_data.items():
            for b, v in (row or {}).items():
                f = _num(v, None)
                a_, b_ = str(a), str(b)
                if f is None or a_ == b_:
                    continue
                f = f / 100.0 if pct else f
                # symmetric: keep the larger of the two directions
                for x, y in ((a_, b_), (b_, a_)):
                    row_x = out.setdefault(x, {})
                    if f > row_x.get(y, -1.0):
                        row_x[y] = f
        return out

    def ani_value(self, a, b):
        return self.ani.get(str(a), {}).get(str(b), 0.0)

    def _ids(self, d, level):
        """(own id, species id, genus id) of an entry at `level`. The levels are
        match_paths' strain (key) -> subkey -> toplevelkey groupings; which of
        them is the species / genus depends on --rank / --subrank."""
        sf, gf = self.species_field, self.genus_field
        if level == "strain":
            own = str(d.get("key"))
            sp = str(d.get(sf, d.get("key")))
            g = str(d.get(gf, "")) if gf else ""
        elif level == "species":
            own = str(d.get("subkey", d.get("key")))
            sp = own if sf == "subkey" else str(d.get(sf, own))
            g = str(d.get(gf, "")) if gf else ""
        else:
            own = str(d.get("toplevelkey", d.get("key")))
            sp = own if sf == "toplevelkey" else None
            g = own if gf == "toplevelkey" else ""
        return own, sp, g

    @staticmethod
    def _name(d):
        return d.get("name") or d.get("subkeyname") or d.get("toplevelname") or str(d.get("key"))

    def prepare(self, level, entries, strains, minhash_available=None):
        """Compute the context of every entry at `level`.

        entries: mapping or iterable of the level's organism dicts (the same
                 objects that are scored afterwards).
        strains: the strain-level dicts (strain_summary values) -- used for the
                 sample's complexity and to map species/genus entries to the
                 reference taxids the ANI matrix is keyed by.
        minhash_available: whether the MinHash reference comparison ran
                 (match_paths passes it from the per-accession hits; the
                 aggregated strain dicts no longer carry the marker). None =
                 look for the marker on `strains`."""
        assert level in LEVELS, level
        ents = list(entries.values()) if isinstance(entries, dict) else list(entries)
        strains = list(strains.values()) if isinstance(strains, dict) else list(strains)
        self._ctx = {}

        if level == "strain" or self.sample is None:
            by_species = {}
            for st in strains:
                sp = str(st.get(self.species_field, st.get("key")))
                by_species[sp] = by_species.get(sp, 0.0) + max(0.0, _num(st.get("numreads")))
            self.sample = sample_complexity(by_species.values())
            self.sample["n_species"] = len(by_species)
            # every shift and the penalty lean on the MinHash comparison; without
            # it (--fast / no reference comparison) minhash_reduction is a
            # placeholder, so the sample is scored with static weights.
            if minhash_available is None:
                minhash_available = any("minhash_read_retention" in st for st in strains)
            self.minhash_available = bool(minhash_available)
            if not self.minhash_available and self.disabled_reason is None:
                self.disabled_reason = ("MinHash reference comparison unavailable "
                                        "(--fast or no --compare_references): static weights used")
                self.log(f"[adaptive_tass] {self.disabled_reason}")

        # reference taxids per entry (for ANI) and child counts (for roll-ups)
        refs, children = {}, {}
        for st in strains:
            k = str(st.get("key"))
            sub = str(st.get("subkey", k))
            tgt = {"strain": k, "species": sub,
                   "genus": str(st.get("toplevelkey", st.get("key")))}[level]
            refs.setdefault(tgt, set()).add(k)
            children.setdefault(tgt, set()).add(sub if level == "genus" else k)

        rows = []
        for d in ents:
            own, sp, g = self._ids(d, level)
            rows.append((d, own, sp, g, max(0.0, _num(d.get("numreads")))))
        total = sum(r[4] for r in rows)
        top = max((r[4] for r in rows), default=0.0)
        hi = self.cfg.high_ani
        sev_den = max(_EPS, 1.0 - hi)
        # reference taxid -> entry id, so each entry's ANI neighbours are found
        # by walking its own ANI rows instead of comparing every pair of entries
        ref_owner = {}
        for r in rows:
            for x in refs.get(r[1], {r[1]}):
                ref_owner.setdefault(x, set()).add(r[1])

        for d, own, sp, g, reads in rows:
            best_ani = {}
            if self.ani:
                for x in refs.get(own, {own}):
                    for y, a in self.ani.get(x, {}).items():
                        if a < hi:
                            continue
                        for own2 in ref_owner.get(y, ()):
                            if own2 != own and a > best_ani.get(own2, -1.0):
                                best_ani[own2] = a
            comps, comp_reads, max_ani = [], 0.0, 0.0
            sib_reads = 0.0
            for d2, own2, sp2, g2, reads2 in rows:
                if d2 is d or own2 == own:
                    continue
                if level != "genus" and sp2 == sp:
                    continue  # same-species strains are an assignment question, not a competitor
                a = best_ani.get(own2)
                if a is not None:
                    comps.append({"key": own2, "name": self._name(d2), "ani": round(a, 4),
                                  "reads": reads2})
                    comp_reads += reads2
                    max_ani = max(max_ani, a)
                if level != "genus" and g and g2 == g:
                    sib_reads += reads2
            cp = comp_reads / (comp_reads + reads) if (comp_reads + reads) > 0 else 0.0
            spr = sib_reads / (sib_reads + reads) if (sib_reads + reads) > 0 else 0.0
            sev = _clamp((max_ani - hi) / sev_den) if comps else 0.0
            amb_ani = cp * sev if self.cfg.on("ani") else 0.0
            amb_sib = spr if self.cfg.on("sibling") else 0.0
            ambiguity = _clamp(AMBIGUITY_ANI_SHARE * amb_ani + AMBIGUITY_SIBLING_SHARE * amb_sib)
            comps.sort(key=lambda c: (-c["ani"], -c["reads"], c["key"]))
            self._ctx[id(d)] = {
                "complexity": self.sample["complexity"],
                "relative_abundance": reads / total if total > 0 else 0.0,
                "abundance_score": (math.log1p(reads) / math.log1p(top)) if top > 0 else 0.0,
                "competitors": comps,
                "competitor_pressure": cp,
                "max_ani": max_ani,
                "ani_severity": sev,
                "sibling_pressure": spr,
                "ambiguity": ambiguity,
                "n_children": len(children.get(own, ())),
            }

    # ── scoring ─────────────────────────────────────────────────────────
    def _cutoff_shift(self, ctx, mh, level):
        """Prototype adaptive_threshold minus its base, in TASS points
        (advisory: reported, never applied)."""
        c = self.cfg
        t = CUT_COMPLEXITY * ctx["complexity"] if c.on("complexity") else 0.0
        t += CUT_AMBIGUITY * ctx["ambiguity"] * (1.0 - mh)
        t -= CUT_UNIQUE_CREDIT * max(0.0, mh - 0.75)
        t -= CUT_ABUNDANCE_CREDIT * max(0.0, ctx["abundance_score"] - 0.65)
        if c.on("rollup") and level != "strain" and ctx["n_children"] >= 2:
            t -= CUT_ROLLUP_CREDIT
        return round(100.0 * c.strength * t, 3)

    def score(self, data, level, weights, compute_fn, standard=None):
        """Raw (0-1) adaptive TASS for one organism (the static score under
        ``weights`` when it cannot be adapted). Stamps data['adaptive_tass'].

        ``weights`` is the adaptive starting point (--adaptive_tass_base);
        ``standard`` is the raw standard TASS the report keeps as its score.
        When given, the audit and the per-level change statistics compare
        against it; otherwise against the static score under ``weights``."""
        static = compute_fn(data=data, weights=weights)
        ref = static if standard is None else standard
        audit = {"level": level, "static": round(_num(static), 6),
                 "standard": round(_num(ref), 6)}
        try:
            if self.disabled_reason:
                raise _Static(self.disabled_reason)
            ctx = self._ctx.get(id(data))
            if ctx is None:
                raise _Static("no context prepared for this organism")
            mh = _clamp(_num(data.get("minhash_reduction")))
            adapted, reasons = adapt_weights(weights, ctx, self.cfg, level)
            check_weights(weights, adapted)
            score = _num(compute_fn(data=data, weights=adapted), None)
            if score is None:
                raise ValueError("compute_tass_score returned a non-number")
            penalty = 1.0 - self.cfg.max_penalty * self.cfg.strength * ctx["ambiguity"] * (1.0 - mh)
            penalty = _clamp(penalty, 1.0 - self.cfg.max_penalty, 1.0)
            if penalty < 1.0:
                reasons.append(f"ambiguity {ctx['ambiguity']:.2f} with MinHash {mh:.2f}: x{penalty:.3f}")
            final = _clamp(score * penalty)
            audit.update({
                "status": "applied",
                "adaptive": round(final, 6),
                "penalty": round(penalty, 6),
                "weights": {k: round(_num(adapted.get(k)), 6) for k in CORE_KEYS},
                "weight_shift": round(sum(abs(_num(adapted.get(k)) - _num(weights.get(k)))
                                          for k in CORE_KEYS), 6),
                "complexity": round(ctx["complexity"], 4),
                "relative_abundance": round(ctx["relative_abundance"], 6),
                "max_ani": round(ctx["max_ani"], 4),
                "competitor_pressure": round(ctx["competitor_pressure"], 4),
                "sibling_pressure": round(ctx["sibling_pressure"], 4),
                "ambiguity": round(ctx["ambiguity"], 4),
                "competitors": ctx["competitors"][:5],
                "n_children": ctx["n_children"],
                "cutoff_shift": self._cutoff_shift(ctx, mh, level),
                "reasons": reasons or ["no context term active: static weights"],
            })
            self._tally(level, ref, final)
            data["adaptive_tass"] = audit
            return final
        except _Static as why:
            audit.update({"status": "static", "adaptive": audit["static"], "reasons": [str(why)]})
        except Exception as e:  # noqa: BLE001 -- fail safe: never lose the static score
            self.errors += 1
            msg = f"{type(e).__name__}: {e}"
            if len(self._error_msgs) < 5:
                self._error_msgs.append(f"{level} {data.get('key')}: {msg}")
                self.log(f"[adaptive_tass] WARNING {level} {self._name(data)}: {msg} "
                         f"-- static score kept")
            audit.update({"status": "error", "adaptive": audit["static"], "reasons": [msg]})
        self._tally(level, ref, static)
        data["adaptive_tass"] = audit
        return static

    def _tally(self, level, static, final):
        st = self.stats[level]
        d = 100.0 * (_num(final) - _num(static))
        st["n"] += 1
        if abs(d) > 1e-9:
            st["changed"] += 1
            st["raised" if d > 0 else "lowered"] += 1
        st["sum_delta"] += d
        st["max_abs_delta"] = max(st["max_abs_delta"], abs(d))

    @staticmethod
    def floor_adaptive(parent):
        """match_paths' monotone roll-up floor, applied to tass_score_adaptive:
        a parent's adaptive score is never below its best child's."""
        own = parent.get("tass_score_adaptive")
        if own is None:
            return
        kids = [c.get("tass_score_adaptive") for c in parent.get("members", []) or []]
        best = max((_num(k) for k in kids if k is not None), default=None)
        if best is not None and best > _num(own):
            parent["tass_score_adaptive_computed"] = own
            parent["tass_score_adaptive"] = best

    def summary(self):
        """metadata.adaptive_tass block."""
        lv = {}
        for level, st in self.stats.items():
            if st["n"]:
                lv[level] = {
                    "organisms": st["n"], "changed": st["changed"],
                    "raised": st["raised"], "lowered": st["lowered"],
                    "mean_delta": round(st["sum_delta"] / st["n"], 4),
                    "max_abs_delta": round(st["max_abs_delta"], 4),
                }
        smp = {k: (round(v, 4) if isinstance(v, float) else v) for k, v in (self.sample or {}).items()}
        return {
            "enabled": True,
            "applied": self.disabled_reason is None,
            "reason": self.disabled_reason,
            "strength": self.cfg.strength,
            "terms": list(self.cfg.terms),
            "base": self.cfg.base,
            "high_ani": self.cfg.high_ani,
            "max_penalty": self.cfg.max_penalty,
            "ani_available": bool(self.ani),
            "species_field": self.species_field,
            "genus_field": self.genus_field,
            "minhash_available": self.minhash_available,
            "cutoff_shift_applied": False,
            # tass_score stays the standard score; adaptive is reported beside it
            "primary_score": "standard",
            "sample": smp,
            "levels": lv,
            "errors": self.errors,
            "error_examples": self._error_msgs or None,
        }


def species_fields_for_rank(rank, subrank="species"):
    """(species_field, genus_field, note) for match_paths' --rank / --subrank.

    match_paths groups strains (key) into subkey (--subrank, default species)
    and toplevelkey (--rank; the pipeline passes genus). The adaptive layer
    must know which field is the species so strains of one species never
    count against each other."""
    r = (rank or "").strip().lower()
    sr = (subrank or "").strip().lower()
    if r == "species":
        return "toplevelkey", None, "rank=species: species = toplevelkey; sibling term inactive"
    if sr == "species":
        note = None if r == "genus" else f"rank={r or 'none'}: siblings = same {r or 'top-level'} group"
        return "subkey", "toplevelkey", note
    return "key", "toplevelkey", (f"subrank={sr or 'none'}: species of a strain unknown, so "
                                  f"strains of one species may count as competitors")


class _Static(Exception):
    """Internal: score this organism with static weights (expected, not an error)."""


# ── CLI: inspect an output JSON ────────────────────────────────────────────
def _iter_levels(organisms):
    for g in organisms:
        yield "genus", g
        for s in g.get("members", []) or []:
            yield "species", s
            for st in s.get("members", []) or []:
                yield "strain", st


def summarize(path, level="strain", out=print):
    import json
    with open(path) as fh:
        j = json.load(fh)
    meta = (j.get("metadata") or {}).get("adaptive_tass")
    if not meta:
        out(f"{path}: no adaptive_tass block (run without --adaptive_tass)")
        return 1
    out(f"{path}: applied={meta.get('applied')} strength={meta.get('strength')} "
        f"terms={','.join(meta.get('terms') or [])} base={meta.get('base')} "
        f"ani={meta.get('ani_available')} complexity={(meta.get('sample') or {}).get('complexity')} "
        f"errors={meta.get('errors')}")
    if meta.get("reason"):
        out(f"  reason: {meta['reason']}")
    out(f"  {'level':<8}{'organism':<44}{'standard':>9}{'adaptive':>9}{'change':>8}{'cut+':>7}  reasons")
    for lvl, o in _iter_levels(j.get("organisms", [])):
        if level != "all" and lvl != level:
            continue
        a = o.get("adaptive_tass") or {}
        st, tt = _num(o.get("tass_score")), _num(o.get("tass_score_adaptive", o.get("tass_score")))
        out(f"  {lvl:<8}{str(o.get('name'))[:43]:<44}{st:>9.2f}{tt:>9.2f}{tt - st:>+8.2f}"
            f"{_num(a.get('cutoff_shift')):>+7.1f}  {'; '.join(a.get('reasons') or [])}")
    return 0


def main(argv=None):
    import argparse
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("summarize", help="standard vs adaptive TASS for one paths.json")
    s.add_argument("json")
    s.add_argument("--level", choices=list(LEVELS) + ["all"], default="strain")
    a = p.parse_args(argv)
    return summarize(a.json, a.level)


if __name__ == "__main__":
    raise SystemExit(main())
