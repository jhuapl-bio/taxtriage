#!/usr/bin/env python3
"""Tests for the opt-in adaptive TASS layer (bin/adaptive_tass.py).

The guarantees that make --adaptive_tass safe to ship as an option:

  * strength 0 and terms 'none' reproduce the static score exactly;
  * adapted weights are finite, non-negative, keep the core pool total and
    never touch a non-core weight;
  * the ambiguity penalty stays within [1 - max_penalty, 1];
  * strains of the SAME species never count against each other;
  * any failure falls back to the static score (and is counted);
  * a sample without the MinHash comparison is scored with static weights;
  * the prototype's scenarios keep their direction (the high-ANI decoy with
    weak unique support drops; the supported call does not).

Run:  python3 bin/test_adaptive_tass.py      (or: pytest bin/test_adaptive_tass.py)
"""
import json
import math
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import adaptive_tass as at  # noqa: E402

try:  # the real scorer when the pipeline deps are installed
    from optimize_weights import compute_tass_score  # noqa: E402
except Exception:  # pragma: no cover - no pandas/scipy: minimal stand-in (additive + power)
    def compute_tass_score(data={}, weights={}):
        s = sum(at._num(data.get(m, 0)) * at._num(weights.get(w, 0)) for m, w in (
            ("breadth_log_score", "breadth_weight"), ("minhash_reduction", "minhash_weight"),
            ("gini_coefficient", "gini_weight"), ("disparity", "disparity_weight"),
            ("hmp_percentile", "hmp_weight"), ("mapq_score", "mapq_score")))
        p = at._num(weights.get("score_power", 1.0), 1.0)
        s = max(0.0, min(1.0, s))
        return s ** p if p > 0 else s

HERE = os.path.dirname(os.path.abspath(__file__))
THRESHOLDS = os.path.join(HERE, "..", "assets", "sampletype_best_thresholds.json")


def base_weights():
    """match_paths' weights dict built from the all|all entry."""
    with open(THRESHOLDS) as fh:
        bw = json.load(fh)["all|all"]["best_weights"]
    w = {"mapq_score": 0.0, "disparity_weight": bw["disparity_weight"], "hmp_weight": bw["hmp_weight"],
         "gini_weight": bw["gini_weight"], "breadth_weight": bw["breadth_weight"],
         "minhash_weight": bw["minhash_weight"], "siblings_score": 0, "diamond_identity": 0.0,
         "k2_disparity_score_weight": 0.0}
    tot = sum(w.values())
    w = {k: v / tot for k, v in w.items()}
    w.update(plasmid_bonus_weight=bw["plasmid_bonus_weight"], abundance_confidence_weight=0.0,
             abundance_gate=False, score_power=bw["score_power"])
    return w


def strain(key, name, species, genus, reads, breadth, gini, minhash, hmp=0.0, disparity=0.5):
    """A strain dict carrying the fields compute_tass_score reads."""
    return {"key": key, "name": name, "subkey": species, "toplevelkey": genus, "numreads": reads,
            "breadth_log_score": breadth, "gini_coefficient": gini, "minhash_reduction": minhash,
            "hmp_percentile": hmp, "disparity": disparity, "minhash_read_retention": 0.5}


def run_strains(strains, ani=None, **cfg):
    scorer = at.AdaptiveTass(at.AdaptiveConfig(**cfg).validate(), ani_data=ani or {}, log=lambda *_: None)
    scorer.prepare("strain", strains, strains)
    w = base_weights()
    out = {}
    for s in strains:
        out[s["key"]] = (100 * compute_tass_score(data=s, weights=w),
                         100 * scorer.score(s, "strain", w, compute_tass_score))
    return out, scorer


# Prototype demo (adaptive_weighting.py build_demo), mapped onto TaxTriage fields.
def pox_high_ani():
    return [strain("10245", "Vaccinia virus", "10245", "10242", 4000, 0.80, 0.86, 0.92, 0.0, 0.62),
            strain("10243", "Cowpox virus", "10243", "10242", 3300, 0.78, 0.84, 0.18, 0.0, 0.40)]


POX_ANI = {"10245": {"10243": 0.982}}


def complex_stool():
    return [strain("562", "Escherichia coli", "562", "561", 6000, 0.80, 0.82, 0.86, 0.95, 0.45),
            strain("623", "Shigella flexneri", "623", "620", 900, 0.34, 0.68, 0.17, 0.08, 0.22),
            strain("28901", "Salmonella enterica", "28901", "590", 1400, 0.82, 0.84, 0.95, 0.12, 0.58)]


STOOL_ANI = {"562": {"623": 0.965, "28901": 0.82}, "623": {"28901": 0.81}}


# ── configuration ──────────────────────────────────────────────────────────
def test_parse_terms():
    assert at.parse_terms(None) == at.TERMS
    assert at.parse_terms("all") == at.TERMS
    assert at.parse_terms("") == at.TERMS
    assert at.parse_terms("none") == ()
    assert at.parse_terms("rollup, ani,ani") == ("ani", "rollup")  # canonical order, de-duplicated
    for bad in ("anni", "ani,foo"):
        try:
            at.parse_terms(bad)
        except ValueError:
            continue
        raise AssertionError(f"parse_terms accepted {bad!r}")


def test_config_validation():
    at.AdaptiveConfig().validate()
    for kw in ({"strength": -0.1}, {"strength": 1.5}, {"strength": float("nan")},
               {"high_ani": 95.0}, {"high_ani": 1.0}, {"high_ani": 0.2},
               {"max_penalty": 0.95}, {"max_penalty": -1}, {"base": "blood"}, {"terms": ("x",)}):
        try:
            at.AdaptiveConfig(**kw).validate()
        except ValueError:
            continue
        raise AssertionError(f"validate accepted {kw}")


def test_from_args_defaults_to_ani_threshold():
    class A:  # argparse.Namespace stand-in
        adaptive_tass_strength = 0.5
        adaptive_tass_terms = "ani"
        adaptive_tass_high_ani = None
        ani_threshold = 0.96
        adaptive_tass_max_penalty = 0.2
        adaptive_tass_base = "global"
    cfg = at.AdaptiveConfig.from_args(A())
    assert (cfg.strength, cfg.terms, cfg.high_ani, cfg.max_penalty, cfg.base) == \
        (0.5, ("ani",), 0.96, 0.2, "global")


# ── exactness of the "off" settings ───────────────────────────────────────
def test_strength_zero_is_static():
    for strains, ani in ((pox_high_ani(), POX_ANI), (complex_stool(), STOOL_ANI)):
        out, _ = run_strains(strains, ani, strength=0.0)
        for k, (s, a) in out.items():
            assert s == a, (k, s, a)


def test_terms_none_is_static():
    out, _ = run_strains(complex_stool(), STOOL_ANI, terms=())
    for k, (s, a) in out.items():
        assert s == a, (k, s, a)


# ── weight invariants (randomised) ────────────────────────────────────────
def test_weight_invariants_random():
    rng = random.Random(1234)
    base = base_weights()
    for i in range(3000):
        ctx = {"complexity": rng.random(), "relative_abundance": rng.random() ** 3,
               "competitors": [{"key": "x"}] * rng.randint(0, 3),
               "competitor_pressure": rng.random(), "max_ani": 0.95 + 0.05 * rng.random(),
               "sibling_pressure": rng.random(), "ambiguity": rng.random(),
               "n_children": rng.randint(0, 4), "abundance_score": rng.random()}
        cfg = at.AdaptiveConfig(strength=rng.random(), terms=tuple(t for t in at.TERMS if rng.random() < 0.7))
        level = rng.choice(at.LEVELS)
        w, _ = at.adapt_weights(base, ctx, cfg, level)
        at.check_weights(base, w)
        assert all(math.isfinite(w[k]) and w[k] >= 0 for k in at.CORE_KEYS)


def test_shift_directions():
    base = base_weights()
    cfg = at.AdaptiveConfig()
    plain = {"complexity": 0.0, "relative_abundance": 0.5, "competitors": [], "competitor_pressure": 0,
             "max_ani": 0, "sibling_pressure": 0, "ambiguity": 0, "n_children": 1, "abundance_score": 1}
    w, r = at.adapt_weights(base, plain, cfg, "strain")
    assert w == base and r == []
    hi = dict(plain, competitors=[{"key": "x"}], competitor_pressure=0.5, max_ani=0.98, ambiguity=0.4)
    w, _ = at.adapt_weights(base, hi, cfg, "strain")
    assert w["minhash_weight"] > base["minhash_weight"] and w["breadth_weight"] < base["breadth_weight"]
    roll = dict(plain, n_children=3)
    w, _ = at.adapt_weights(base, roll, cfg, "species")
    assert w["breadth_weight"] > base["breadth_weight"] and w["minhash_weight"] < base["minhash_weight"]
    w, _ = at.adapt_weights(base, roll, cfg, "strain")  # strains never roll up
    assert w == base


def test_penalty_bounds():
    rng = random.Random(7)
    for _ in range(300):
        strains = [strain(str(i), f"o{i}", str(i), "g", rng.randint(1, 5000), rng.random(), rng.random(),
                          rng.random(), rng.random(), rng.random()) for i in range(6)]
        ani = {str(i): {str(j): 0.95 + 0.05 * rng.random() for j in range(6) if j != i} for i in range(6)}
        mp = rng.choice((0.0, 0.1, 0.35, 0.9))
        scorer = at.AdaptiveTass(at.AdaptiveConfig(max_penalty=mp, strength=rng.random()).validate(),
                                 ani_data=ani, log=lambda *_: None)
        scorer.prepare("strain", strains, strains)
        for s in strains:
            scorer.score(s, "strain", base_weights(), compute_tass_score)
            p = s["adaptive_tass"]["penalty"]
            assert 1.0 - mp - 1e-9 <= p <= 1.0, p


# ── who counts as a competitor / sibling ──────────────────────────────────
def test_same_species_strains_never_compete():
    strains = [strain("83334", "E. coli O157:H7", "562", "561", 900, 0.9, 0.9, 0.3),
               strain("511145", "E. coli K-12", "562", "561", 1100, 0.9, 0.9, 0.3)]
    ani = {"83334": {"511145": 0.999}}
    _, sc = run_strains(strains, ani)
    for s in strains:
        a = s["adaptive_tass"]
        assert a["competitors"] == [] and a["sibling_pressure"] == 0 and a["ambiguity"] == 0, a
        assert a["penalty"] == 1.0


def test_species_fields_follow_rank():
    assert at.species_fields_for_rank("genus", "species")[:2] == ("subkey", "toplevelkey")
    assert at.species_fields_for_rank("species", "species")[:2] == ("toplevelkey", None)
    assert at.species_fields_for_rank("family", "species")[:2] == ("subkey", "toplevelkey")
    assert at.species_fields_for_rank("genus", "none")[0] == "key"


def test_rank_species_mode_same_species_never_compete():
    # --rank species: toplevelkey is the species and subkey is the strain itself
    strains = [dict(strain("83334", "E. coli O157:H7", "83334", "562", 900, 0.9, 0.9, 0.3)),
               dict(strain("511145", "E. coli K-12", "511145", "562", 1100, 0.9, 0.9, 0.3)),
               dict(strain("623", "Shigella flexneri", "623", "623", 300, 0.5, 0.6, 0.1))]
    ani = {"83334": {"511145": 0.999, "623": 0.97}, "511145": {"623": 0.97}}
    sc = at.AdaptiveTass(at.AdaptiveConfig(), ani_data=ani, log=lambda *_: None,
                         species_field="toplevelkey", genus_field=None)
    sc.prepare("strain", strains, strains)
    for s in strains:
        sc.score(s, "strain", base_weights(), compute_tass_score)
    by = {s["key"]: s["adaptive_tass"] for s in strains}
    assert [c["key"] for c in by["83334"]["competitors"]] == ["623"], by["83334"]["competitors"]
    assert {c["key"] for c in by["623"]["competitors"]} == {"83334", "511145"}
    assert all(a["sibling_pressure"] == 0 for a in by.values())


def test_same_genus_species_are_siblings():
    strains = [strain("1313", "S. pneumoniae", "1313", "1301", 600, 0.9, 0.9, 0.9),
               strain("28037", "S. mitis", "28037", "1301", 400, 0.6, 0.6, 0.4)]
    run_strains(strains)
    assert abs(strains[0]["adaptive_tass"]["sibling_pressure"] - 0.4) < 1e-9
    assert abs(strains[1]["adaptive_tass"]["sibling_pressure"] - 0.6) < 1e-9


def test_species_and_genus_levels_use_member_ani():
    strains = complex_stool()
    species = [{"key": s["subkey"], "subkey": s["subkey"], "toplevelkey": s["toplevelkey"],
                "numreads": s["numreads"], "name": s["name"], **{k: s[k] for k in (
                    "breadth_log_score", "gini_coefficient", "minhash_reduction", "hmp_percentile",
                    "disparity")}} for s in strains]
    genera = [dict(sp, key=sp["toplevelkey"]) for sp in species]
    sc = at.AdaptiveTass(at.AdaptiveConfig(), ani_data=STOOL_ANI, log=lambda *_: None)
    sc.prepare("strain", strains, strains)
    for lvl, ents in (("species", species), ("genus", genera)):
        sc.prepare(lvl, ents, strains)
        for e in ents:
            sc.score(e, lvl, base_weights(), compute_tass_score)
        by = {e["name"]: e["adaptive_tass"] for e in ents}
        assert by["Shigella flexneri"]["max_ani"] == 0.965, (lvl, by["Shigella flexneri"])
        assert by["Salmonella enterica"]["competitors"] == [], lvl


def test_percent_ani_matrix_normalised():
    sc = at.AdaptiveTass(at.AdaptiveConfig(), ani_data={"a": {"b": 97.5, "a": 100.0}})
    assert sc.ani_value("a", "b") == 0.975 and sc.ani_value("b", "a") == 0.975
    assert sc.ani_value("a", "a") == 0.0  # diagonal ignored


# ── prototype scenarios keep their direction ──────────────────────────────
def test_pox_high_ani_scenario():
    out, _ = run_strains(pox_high_ani(), POX_ANI)
    (vs, va), (cs, ca) = out["10245"], out["10243"]
    assert ca < cs, "Cowpox decoy (weak MinHash) must drop"
    assert (cs - ca) > (vs - va), "the decoy must drop more than the supported call"
    assert va >= 75.0 > ca, (va, ca)


def test_complex_stool_scenario():
    out, _ = run_strains(complex_stool(), STOOL_ANI)
    drop = {k: s - a for k, (s, a) in out.items()}
    assert drop["623"] == max(drop.values()), drop  # Shigella decoy drops most
    assert out["623"][1] < 75.0, out["623"]
    assert out["562"][1] >= 75.0 and out["28901"][1] >= 75.0, out


# ── fail-safe paths ───────────────────────────────────────────────────────
def test_error_falls_back_to_static():
    strains = pox_high_ani()
    sc = at.AdaptiveTass(at.AdaptiveConfig(), ani_data=POX_ANI, log=lambda *_: None)
    sc.prepare("strain", strains, strains)
    w = base_weights()

    def flaky(data={}, weights={}):
        if weights is not w:
            raise RuntimeError("boom")
        return compute_tass_score(data=data, weights=weights)
    for s in strains:
        assert sc.score(s, "strain", w, flaky) == compute_tass_score(data=s, weights=w)
        assert s["adaptive_tass"]["status"] == "error"
    assert sc.errors == 2 and sc.summary()["errors"] == 2


def test_missing_minhash_disables_sample():
    strains = pox_high_ani()
    sc = at.AdaptiveTass(at.AdaptiveConfig(), ani_data=POX_ANI, log=lambda *_: None)
    sc.prepare("strain", strains, strains, minhash_available=False)
    w = base_weights()
    for s in strains:
        assert sc.score(s, "strain", w, compute_tass_score) == compute_tass_score(data=s, weights=w)
        assert s["adaptive_tass"]["status"] == "static"
    summ = sc.summary()
    assert summ["applied"] is False and "MinHash" in summ["reason"]


def test_unprepared_entry_is_static():
    sc = at.AdaptiveTass(at.AdaptiveConfig(), log=lambda *_: None)
    s = pox_high_ani()[0]
    w = base_weights()
    assert sc.score(s, "strain", w, compute_tass_score) == compute_tass_score(data=s, weights=w)
    assert s["adaptive_tass"]["status"] == "static"


def test_bad_inputs_do_not_crash():
    strains = [strain("1", "a", "1", "g", None, None, "x", float("nan")),
               strain("2", "b", "2", "g", -5, 0.5, 0.5, 0.5)]
    out, sc = run_strains(strains, {"1": {"2": "n/a"}})
    for s in strains:
        assert s["adaptive_tass"]["status"] in ("applied", "static", "error")


# ── audit / summary ───────────────────────────────────────────────────────
def test_floor_static():
    parent = {"tass_score_static": 60.0, "members": [{"tass_score_static": 70.0}, {"tass_score": 90.0}]}
    at.AdaptiveTass.floor_static(parent)
    assert parent["tass_score_static"] == 90.0
    parent = {"tass_score_static": 95.0, "members": [{"tass_score_static": 70.0}]}
    at.AdaptiveTass.floor_static(parent)
    assert parent["tass_score_static"] == 95.0


def test_summary_and_audit_serialisable_and_deterministic():
    dumps = []
    for _ in range(2):
        strains = complex_stool()
        _, sc = run_strains(strains, STOOL_ANI)
        dumps.append(json.dumps({"summary": sc.summary(), "audit": [s["adaptive_tass"] for s in strains]},
                                sort_keys=True))
    assert dumps[0] == dumps[1]
    summ = json.loads(dumps[0])["summary"]
    assert summ["levels"]["strain"]["organisms"] == 3 and summ["cutoff_shift_applied"] is False


if __name__ == "__main__":
    tests = [(n, f) for n, f in sorted(globals().items()) if n.startswith("test_") and callable(f)]
    failed = 0
    for name, fn in tests:
        try:
            fn()
            print(f"PASS {name}")
        except Exception as e:  # noqa: BLE001
            failed += 1
            print(f"FAIL {name}: {type(e).__name__}: {e}")
    print(f"{len(tests) - failed}/{len(tests)} passed")
    sys.exit(1 if failed else 0)
