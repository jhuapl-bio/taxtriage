from collections import defaultdict
import math
from utils import logarithmic_weight, apply_weight, normalize_mapq
from map_taxid import get_lineage
import numpy as np
from scipy.stats import norm
from typing import Dict, Any, List, Iterable
from body_site_normalization import normalize_body_site, get_pathogen_classification


def _truth_accession_from_qname(qname: str) -> str:
    """
    Ground truth rule you described:
      NC_003310.1_3323_4  -> NC_003310.1
    """
    return qname.split("_", 1)[0]


def compute_tp_fp_counts_by_reference(bam_path: str) -> dict[str, tuple[int, int]]:
    """
    Returns: { reference_accession: (tp_reads, fp_reads) }
    TP if RNAME == truth_accession_from_qname(QNAME), else FP.
    """
    import pysam  # import here so normal runs don't require pysam unless BAM is used

    counts = defaultdict(lambda: [0, 0])  # ref -> [tp, fp]
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for aln in bam.fetch(until_eof=True):
            if aln.is_unmapped:
                continue
            # Skip non-primary to avoid double counting unless your pipeline expects otherwise
            if aln.is_secondary or aln.is_supplementary:
                continue

            ref = bam.get_reference_name(aln.reference_id)
            truth = _truth_accession_from_qname(aln.query_name)

            if ref == truth:
                counts[ref][0] += 1
            else:
                counts[ref][1] += 1

    return {k: (v[0], v[1]) for k, v in counts.items()}



import pandas as pd
CATEGORY_PRIORITY = {
    "Primary": 5,
    "Potential": 4,
    "Opportunistic": 3,
    "Commensal": 2,
    "Unknown": 1,
}

def transform_func(d):
    return math.log10(1 + d)


def build_transformed_coverage_hist(regions, genome_length):
    """
    Build a histogram of 'transformed' coverage for the entire genome.
    transform_func(depth) should return the transformed coverage value.
    """
    coverage_hist = defaultdict(int)
    total_covered_bases = 0
    for (start, end, depth) in regions:
        length = end - start
        transformed_depth = transform_func(depth)
        coverage_hist[transformed_depth] += length
        total_covered_bases += length

    # Add zeros (transformed) for uncovered portion
    uncovered = genome_length - total_covered_bases
    if uncovered > 0:
        coverage_hist[transform_func(0)] += uncovered

    return coverage_hist


def gini_coefficient_from_hist(coverage_hist):
    """
    Calculate the Gini coefficient from a histogram mapping coverage->number_of_bases.
    """
    if not coverage_hist:
        return 0.0  # no data, gini is 0 by convention

    # Sort coverage values in ascending order
    coverage_values = sorted(coverage_hist.keys())

    # N = total number of bases
    N = sum(coverage_hist[c] for c in coverage_values)
    # If somehow N is 0, return 0 to avoid division by zero
    if N == 0:
        return 0.0

    # total coverage across all bases
    total_coverage = sum(c * coverage_hist[c] for c in coverage_values)
    if total_coverage == 0:
        # All depths are zero => distribution is uniform(=0) => Gini = 0
        return 0.0

    # Build the Lorenz points (x_i, y_i)
    # x_i = cumulative fraction of bases
    # y_i = cumulative fraction of coverage
    lorenz_points = []
    pop_cum = 0
    coverage_cum = 0

    # Start from (0, 0) on the Lorenz curve
    lorenz_points.append((0.0, 0.0))

    for c in coverage_values:
        freq = coverage_hist[c]
        pop_cum_prev = pop_cum
        coverage_cum_prev = coverage_cum

        pop_cum += freq
        coverage_cum += c * freq

        x1 = pop_cum_prev / N
        y1 = coverage_cum_prev / total_coverage
        x2 = pop_cum / N
        y2 = coverage_cum / total_coverage

        # Append the point (x2, y2)
        lorenz_points.append((x2, y2))

    # Now, approximate the area under the Lorenz curve via trapezoids
    area_under_lorenz = 0.0
    for i in range(len(lorenz_points) - 1):
        x1, y1 = lorenz_points[i]
        x2, y2 = lorenz_points[i + 1]
        base = x2 - x1           # horizontal distance
        avg_height = (y1 + y2) / 2
        area_under_lorenz += base * avg_height

    # Finally, Gini = 1 - 2 * area under the Lorenz curve
    gini = 1 - 2 * area_under_lorenz
    return gini

def get_dynamic_reward_factor(genome_length, baseline=5e4, max_length=1e7, max_reward=2):
    """
    Computes a dynamic reward factor:
      - Returns 1 when genome_length is at or below baseline.
      - Returns max_reward when genome_length is at or above max_length.
      - Otherwise, linearly interpolates between 1 and max_reward.
    """
    if genome_length <= baseline:
        return 1.0
    elif genome_length >= max_length:
        return max_reward
    else:
        # Linear interpolation between 1 and max_reward
        return 1.0 + (max_reward - 1.0) * ((genome_length - baseline) / (max_length - baseline))


def getGiniCoeff(regions, genome_length, alpha=1.8, baseline=5e5, max_length=1e9, reward_factor=2, beta=0.5):
    """
    Calculate an adjusted 'Gini-based' score for the fair distribution of coverage,
    and then penalize (or boost) it according to the disparity in positions of the regions.

    Parameters:
      - regions: list of tuples (start, end, depth)
      - genome_length: total genome length
      - alpha: parameter for transforming the raw Gini
      - baseline, max_length, reward_factor: parameters for length-based scaling
      - beta: weight for the positional dispersion factor

    The final score is a product of the (transformed) Gini measure, a scaling factor
    based on genome length, and a term (1 + beta * dispersion) where dispersion is higher
    when the regions are more spread out.
    """
    # 1) Build Histograms
    coverage_hist_transformed = build_transformed_coverage_hist(regions, genome_length)

    # 2) Compute raw Gini from the transformed histogram
    gini = gini_coefficient_from_hist(coverage_hist_transformed)

    # 3) Transform the raw Gini (ensuring the result is in [0,1])
    if 0.0 <= gini <= 1.0:
        gini_log = alpha * math.sqrt(1 - gini)
        gini_log = max(0.0, min(1.0, gini_log))
    else:
        gini_log = 0.0

    # 4) Compute length-based scaling (using a log scale)
    gl_capped = min(genome_length, max_length)
    ratio = gl_capped / baseline
    ratio = max(ratio, 1.0)
    scaling_factor = 1.0 + reward_factor * math.log10(ratio)

    # 5) Compute the positional dispersion factor
    dispersion = position_dispersion_factor(regions, genome_length)

    # 6) Combine the measures.
    # The idea is to boost the score if the covered regions are spread out.
    final_score = gini_log * scaling_factor * (1 + beta * dispersion)
    final_score = min(1.0, final_score)
    return final_score

def position_dispersion_factor(regions, genome_length):
    """
    Compute a dispersion factor based on the positions of the regions.
    We use the midpoints of each region and compute their variance.
    For a uniformly distributed set of midpoints on [0, genome_length],
    the maximum variance is (genome_length^2)/12.
    We then take the square root of the normalized variance to obtain a
    factor between 0 and 1.
    """
    if not regions:
        return 0.0

    # Compute midpoints for each region
    midpoints = [(start + end) / 2.0 for (start, end, depth) in regions]
    mean_mid = sum(midpoints) / len(midpoints)
    variance = sum((m - mean_mid)**2 for m in midpoints) / len(midpoints)

    # Maximum variance for a uniform distribution in [0, genome_length]
    max_variance = (genome_length**2) / 12.0
    normalized_variance = variance / max_variance  # in [0,1]

    # Taking square root to keep the metric in a similar scale as a coefficient of variation
    dispersion = math.sqrt(normalized_variance)
    return dispersion

def gini_coefficient(values):
    """
    Compute the Gini coefficient using a sorted approach.
    """
    n = len(values)
    if n == 0:
        return 0.0

    # Sort the values in ascending order
    sorted_values = sorted(values)
    total = sum(sorted_values)

    if total == 0:
        return 0.0

    # Compute the weighted sum: sum(i * x_i) for i=1..n (using 1-indexing)
    weighted_sum = sum((i + 1) * x for i, x in enumerate(sorted_values))

    # Apply the formula:
    gini = (2 * weighted_sum) / (n * total) - (n + 1) / n
    return gini

def calculate_k2_reads_disparity(
    data = {},
    k2_mapping = {}
):
    taxid = data.get('key', None)
    # Define the rank code you're looking for, for example "G" for Genus
    dispari_k2 = calculate_k2_disparity_score(
        taxid = taxid,
        k2_mapping=k2_mapping,
        parent_rank_match = "G",
    )
    return dispari_k2.get('disparity', 0.0) if dispari_k2 else 0.0

def calculate_disparity(numreads, total_reads, variance_reads, k=1000):
        """
        Dynamically dampens the variance effect based on the proportion of reads.
        numreads: Total number of reads aligned to the organism (sum of reads)
        total_reads: Total number of reads aligned in the sample
        variance_reads: Variance of the aligned reads across all organisms
        k: Damping factor to control the influence of the proportion on the penalty
        """
        if total_reads == 0:
            return 0  # Avoid division by zero

        # Calculate the proportion of aligned reads
        proportion = numreads / total_reads

        # Dynamically adjust the variance penalty based on the proportion of reads
        dampened_variance = variance_reads / (1 + k * proportion)

        # Calculate disparity based on the proportion and the dynamically dampened variance
        disparity = proportion * (1 + dampened_variance)

        return disparity
def calculate_k2_disparity_score(
    taxid,
    k2_mapping,
    sibling_rank=None,
    parent_rank_match="F",
    value_field="clades_covered",
):
    """
    Compute disparity and proportion of `taxid` among siblings of `sibling_rank`
    that share the same ancestor at `parent_rank_match`.

    Returns a dict with keys:
      - taxid, my_rank, parent_taxid, parent_rank_match,
      - siblings_count, sibling_reads (sorted list of (value, taxid)),
      - index (0-based position of taxid in sorted sibling_reads),
      - disparity (1 - idx/(n-1) if n>1 else 0),
      - proportion (own_value / total_value, 0 if total_value == 0),
      - own_value, total_value

    If taxid not found or parent not found, returns None.
    """
    # normalize taxid to string for dict lookup
    taxid = str(taxid)

    if taxid not in k2_mapping:
        return None
    node = k2_mapping[taxid]
    my_rank = node.get("rank", None)

    # decide which rank we want to compare among; default to the node's rank
    if sibling_rank is None:
        sibling_rank = my_rank

    # find the parent taxid at the requested parent_rank_match in this node's parents
    parent_taxid = None
    for p in node.get("parents", []):
        if len(p) >= 2 and p[1] == parent_rank_match:
            parent_taxid = str(p[0])
            break

    if parent_taxid is None:
        # couldn't find an ancestor at that rank for this taxid
        return None
    # gather siblings: entries that:
    #  - have rank == sibling_rank
    #  - have the parent_taxid (with the matching rank) present in their parents list
    siblings = []
    for candidate_taxid, candidate in k2_mapping.items():
        cand_rank = candidate.get("rank", None)
        if cand_rank != sibling_rank:
            continue
        # check ancestry
        candidate_parents = candidate.get("parents", [])
        found_parent = False
        for pp in candidate_parents:
            if len(pp) >= 2 and str(pp[0]) == parent_taxid and pp[1] == parent_rank_match:
                found_parent = True
                break
        if found_parent:
            val = candidate.get(value_field, 0) or 0
            # normalize numeric types to float
            try:
                val = float(val)
            except Exception:
                val = 0.0
            siblings.append((val, str(candidate.get("taxid", candidate_taxid))))

    if not siblings:
        # no siblings found (maybe because ranks mismatch)
        return {
            "taxid": taxid,
            "my_rank": my_rank,
            "parent_taxid": parent_taxid,
            "parent_rank_match": parent_rank_match,
            "siblings_count": 0,
            "sibling_reads": [],
            "index": -1,
            "disparity": 0.0,
            "proportion": 0.0,
            "own_value": 0.0,
            "total_value": 0.0,
        }

    # sort descending by value
    siblings_sorted = sorted(siblings, key=lambda x: x[0], reverse=True)

    # total reads among siblings
    total_value = sum(v for v, _ in siblings_sorted)

    # find our taxid position
    taxid_list = [t for _, t in siblings_sorted]
    try:
        idx = taxid_list.index(taxid)
    except ValueError:
        # maybe taxid stored as int in mapping; try loose match
        try:
            idx = taxid_list.index(str(int(taxid)))
        except Exception:
            idx = -1

    own_value = 0.0
    for v, t in siblings_sorted:
        if t == taxid:
            own_value = v
            break

    # disparity: if n==1 -> 0. otherwise 1 - (index / (n-1)), so top (index=0) -> 1.0, bottom -> 0.0
    n = len(siblings_sorted)
    if n <= 1 or idx < 0:
        disparity = 0.0
    else:
        disparity = 1.0 - (idx / (n - 1))

    proportion = (own_value / total_value) if total_value > 0 else 0.0

    return {
        "taxid": taxid,
        "my_rank": my_rank,
        "parent_taxid": parent_taxid,
        "parent_rank_match": parent_rank_match,
        "siblings_count": n,
        "sibling_reads": siblings_sorted,
        "index": idx,
        "disparity": disparity,
        "proportion": proportion,
        "own_value": own_value,
        "total_value": total_value,
    }
def calculate_normalized_groups(
    hits: Dict[str, Dict[str, Any]],
    group_field: str,
    reads_key: str = "numreads",
    sum_columns: Iterable[str] = (),
    mapq_breadth_power: float = 2.0,
    mapq_gini_power: float = 1.0,
    contig_penalty_power: float = 0.3,
    depth_concentration_power: float = 0.3,
    default_read_length: int = 150,
) -> Dict[str, Dict[str, Any]]:
    """
    Aggregate `hits` into group-level summaries keyed by `group_field`.
    """

    # sums
    SUM_FIELDS = {"numreads", "covered_bases", "length", "k2_reads"}
    SUM_FIELDS |= set(sum_columns or [])

    # weighted means (weights = numreads)
    WAVG_FIELDS = {
        "meandepth", "meanmapq", "meanbaseq",
        "minhash_score", "k2_disparity_score", "hmp_percentile",
        # REMOVED: "breadth_log_score",  # will recalculate from aggregated coverage
        "minhash_reduction", "diamond_identity", "rpkm", "rpm",
        # NOTE: gini_coefficient is NOT in WAVG_FIELDS — it is recomputed from
        # the combined coverage histogram at the group level (like breadth).
        "mapq_score", "rpm_confidence_weight", "read_fraction",
        "highmapq_fraction", "plasmid_score",
        "abundance_confidence",  # low-abundance confidence sigmoid (log-RPM based)
    }

    def _to_float(x) -> float:
        try:
            if x is None:
                return 0.0
            return float(x)
        except Exception:
            return 0.0

    def _mean_list(v) -> float:
        if isinstance(v, (list, tuple)):
            if len(v) == 0:
                return 0.0
            return sum(_to_float(xx) for xx in v) / float(len(v))
        return _to_float(v)

    def _weighted_mean(entries: List[Dict[str, Any]], field: str, wkey: str) -> float:
        num = 0.0
        den = 0.0
        for e in entries:
            w = _to_float(e.get(wkey, 0))
            if w <= 0:
                continue
            v = _mean_list(e.get(field, 0))
            num += v * w
            den += w
        return (num / den) if den > 0 else 0.0

    # 1) bucket entries by group_field
    buckets: Dict[str, List[Dict[str, Any]]] = {}
    for _id, rec in hits.items():
        g = rec.get(group_field)
        if g in (None, ""):
            g = _id
        buckets.setdefault(str(g), []).append(rec)

    # 2) aggregate each bucket
    out: Dict[str, Dict[str, Any]] = {}

    for gval, entries in buckets.items():
        if not entries:
            continue

        agg: Dict[str, Any] = {group_field: gval}

        # preserve a representative name if present
        agg["name"] = entries[0].get("name") or entries[0].get("organism") or gval
        if "subkey" not in agg and "subkey" in entries[0]:
            agg["subkey"] = entries[0].get("subkey")
        # Always propagate subkeyname from the constituent entries
        if "subkeyname" not in agg:
            agg["subkeyname"] = entries[0].get("subkeyname") or agg["name"]

        # ---------- TOPLEVEL KEY + NAME ----------
        if group_field == "toplevelkey":
            agg["toplevelkey"] = gval
            agg["toplevelname"] = (
                entries[0].get("toplevelname")
                or entries[0].get("name")
                or gval
            )
        else:
            tk = None
            tn = None
            for e in entries:
                if tk is None:
                    tk = e.get("toplevelkey")
                if tn is None:
                    tn = e.get("toplevelname")
                if tk and tn:
                    break

            agg["toplevelkey"] = tk
            agg["toplevelname"] = tn

        # sums
        sums = {f: 0.0 for f in SUM_FIELDS}
        for e in entries:
            for f in SUM_FIELDS:
                sums[f] += _to_float(e.get(f, 0))

        agg["numreads"] = sums.get("numreads", 0.0)
        agg["length"] = sums.get("length", 0.0)
        agg["covered_bases"] = sums.get("covered_bases", 0.0)
        agg["k2_reads"] = sums.get("k2_reads", 0.0)

        # counts
        agg["accessions"] = [x.get("accession") for x in entries if x.get("accession")]



        # derived coverage = covered_bases_sum / length_sum
        total_length = agg["length"]
        agg["coverage"] = (agg["covered_bases"] / total_length) if total_length > 0 else 0.0

        # weighted means
        for f in WAVG_FIELDS:
            agg[f] = _weighted_mean(entries, f, reads_key)

        # Propagate has_plasmid: true if ANY member has a plasmid
        agg['has_plasmid'] = any(bool(_e.get('has_plasmid', False)) for _e in entries)

        # ========== RECALCULATE breadth_log_score from aggregated coverage ==========
        # This ensures plasmids don't skew the score.
        # Scale by high-MAPQ fraction so organisms dominated by MAPQ=0 reads
        # get their breadth crushed (unreliable alignment ≠ real coverage).
        _hmf = float(agg.get('highmapq_fraction', 1.0))
        _mapq_scale = _hmf ** mapq_breadth_power
        agg['breadth_log_score'] = breadth_score_sigmoid(agg["coverage"]) * _mapq_scale

        # ========== RECALCULATE gini_coefficient from combined coverage ==========
        # Read-weighted averaging of per-contig Gini is misleading: a few small
        # contigs with uniform coverage can score ~1.0 even when 99.9% of the
        # total genome is uncovered.  Instead, build one combined coverage
        # histogram from ALL entries' covered_regions and compute Gini against
        # the total aggregated genome length — same philosophy as breadth.
        _combined_regions = []
        _offset = 0  # offset each contig's regions into a virtual concatenated genome
        for _e in entries:
            for _start, _end, _depth in _e.get('covered_regions', []):
                _combined_regions.append((_start + _offset, _end + _offset, _depth))
            _offset += int(_to_float(_e.get('length', 0)))
        _total_genome_len = int(agg["length"]) if agg["length"] > 0 else 1
        if _combined_regions:
            agg['gini_coefficient'] = getGiniCoeff(
                _combined_regions, _total_genome_len,
                alpha=1.8, reward_factor=2, beta=0.5)
        else:
            agg['gini_coefficient'] = 0.0
        # Preserve combined regions so the next aggregation level (strain→species)
        # can recompute Gini from the full picture instead of averaging.
        agg['covered_regions'] = _combined_regions

        # ========== MAPQ + CONTIG UTILIZATION PENALTY ON GINI ==========
        # Problem: an organism with 2000 contigs but reads on only 2 can get
        # an inflated Gini because the length-based scaling and dispersion
        # factors in getGiniCoeff over-compensate for large genomes.
        #
        # Fix 1 — Contig utilization penalty:
        #   If only 2/2000 contigs have any reads, coverage is extremely
        #   concentrated and the Gini score should reflect that.
        #   penalty = (covered_contigs / total_contigs) ^ contig_penalty_power
        #   e.g. (2/2000)^0.3 ≈ 0.063 → Gini drops to ~6% of original
        #        (50/100)^0.3  ≈ 0.81  → mild 19% penalty
        #        (100/100)     = 1.0    → no penalty
        #
        # Fix 2 — MAPQ penalty (same pattern as breadth_log_score):
        #   gini *= highmapq_fraction ^ mapq_gini_power
        #   Organisms dominated by low-MAPQ reads get their Gini crushed
        #   because the coverage is unreliable.
        _n_total_contigs = len(entries)
        _n_covered_contigs = sum(
            1 for _e in entries
            if _e.get('numreads', 0) > 0 or len(_e.get('covered_regions', [])) > 0
        )
        _contig_frac = (
            _n_covered_contigs / _n_total_contigs
            if _n_total_contigs > 0 else 1.0
        )

        # Only apply contig penalty when there are multiple contigs —
        # single-contig organisms shouldn't be penalized for having 1/1.
        if _n_total_contigs > 1 and contig_penalty_power > 0:
            _contig_penalty = _contig_frac ** contig_penalty_power
        else:
            _contig_penalty = 1.0

        # MAPQ penalty: same approach as breadth scaling
        _hmf_gini = float(agg.get('highmapq_fraction', 1.0))
        _mapq_gini_scale = _hmf_gini ** mapq_gini_power if mapq_gini_power > 0 else 1.0

        # Fix 3 — Depth-concentration penalty:
        #   Detects the "conserved human reads" pattern: many reads map to a
        #   tiny region of a large genome → absurdly high depth but negligible
        #   coverage.  Measures the ratio of actual coverage to expected coverage
        #   given the read count and genome size.
        #
        #   coverage_efficiency = actual_coverage / expected_coverage
        #   expected_coverage = (numreads * avg_read_length) / genome_length
        #
        #   e.g. Toxoplasma: 80K reads, 65Mbp genome, 0.15% actual coverage
        #        expected ~18%, efficiency = 0.15/18 = 0.008 → penalty^0.3 = 0.22
        #        Burkholderia: 575 reads, 15Mbp genome, 0.46% actual coverage
        #        expected ~0.6%, efficiency = 0.46/0.6 = 0.77 → penalty^0.3 = 0.93
        _depth_conc_penalty = 1.0
        if depth_concentration_power > 0:
            _total_numreads = float(agg.get('numreads', 0))
            _genome_len = float(agg.get('length', 1))
            _actual_cov = float(agg.get('coverage', 0))

            # Estimate avg read length from entries (weighted by numreads)
            _rl_num, _rl_den = 0.0, 0.0
            for _e in entries:
                _nr = _to_float(_e.get('numreads', 0))
                _arl = _to_float(_e.get('avg_read_length', default_read_length))
                if _nr > 0 and _arl > 0:
                    _rl_num += _arl * _nr
                    _rl_den += _nr
            _avg_rl = _rl_num / _rl_den if _rl_den > 0 else float(default_read_length)

            if _total_numreads > 0 and _genome_len > 0:
                _expected_cov = (_total_numreads * _avg_rl) / _genome_len
                _expected_cov = min(1.0, _expected_cov)  # cap at 100%

                if _expected_cov > 0:
                    _cov_efficiency = min(1.0, _actual_cov / _expected_cov)
                    _depth_conc_penalty = _cov_efficiency ** depth_concentration_power
                # else: no reads expected → no penalty

        agg['gini_coefficient'] *= _contig_penalty * _mapq_gini_scale * _depth_conc_penalty
        agg['gini_contig_frac'] = _contig_frac
        agg['gini_contig_penalty'] = _contig_penalty
        agg['gini_mapq_scale'] = _mapq_gini_scale
        agg['gini_depth_conc_penalty'] = _depth_conc_penalty
        agg['n_contigs_total'] = _n_total_contigs
        agg['n_contigs_covered'] = _n_covered_contigs

        # ========== COVERAGE-AWARE GATE ON MINHASH ==========
        # Per-contig confidence gating is defeated by weighted averaging:
        # contigs with 0 reads have minhash_reduction=0 but contribute
        # nothing to the average (weight=0), so a few high-coverage contigs
        # dominate.  E.g. Toxoplasma: 82K reads on 3/80 contigs → per-contig
        # weighted avg ≈ 0.99, but organism-level coverage is only 0.15%.
        #
        # Fix: apply the depth-concentration penalty to minhash as well.
        # This detects the "conserved human reads" pattern where many reads
        # pile onto a tiny region of a large genome, producing high minhash
        # uniqueness per-contig but negligible genome-wide coverage.
        #
        # The penalty reuses _depth_conc_penalty (already computed above for
        # gini) and the contig utilization fraction.  Together they ensure
        # that an organism like Toxoplasma (20K reads, 65Mb genome, 0.05%
        # coverage, efficiency ~0.002) gets its minhash crushed:
        #   minhash_confidence = depth_conc_penalty * contig_penalty
        #   e.g. 0.22 * 0.063 ≈ 0.014 → minhash drops from 1.0 to ~0.01
        #
        # Well-covered organisms (efficiency ≈ 1.0, full contig utilization)
        # are unaffected: penalty ≈ 1.0.
        _agg_raw_mh = float(agg.get('minhash_score', agg.get('minhash_reduction', 0)))

        # Combine depth-concentration efficiency with contig utilization
        # AND a direct breadth sigmoid so that organisms with near-zero
        # genome coverage (< ~1%) have their minhash crushed regardless
        # of the depth-concentration heuristic.
        #
        # breadth_gate uses a sigmoid centred at 1% coverage — same shape
        # as breadth_log_score but without the MAPQ scaling so we don't
        # double-penalise.  At 0.05% coverage → ~0.0;  at 2% → ~1.0.
        #
        # The three factors multiply so ALL must be healthy for the gate
        # to pass:  coverage_efficiency * contig_utilisation * breadth_gate
        _agg_cov = float(agg.get('coverage', 0))
        _breadth_gate = breadth_score_sigmoid(_agg_cov, midpoint=0.01, steepness=12_000)
        _mh_coverage_gate = _depth_conc_penalty * _contig_penalty * _breadth_gate
        agg['minhash_reduction'] = _agg_raw_mh * _mh_coverage_gate
        agg['minhash_confidence'] = _mh_coverage_gate
        agg['minhash_reduction_pre_gate'] = _agg_raw_mh
        agg['minhash_coverage_gate'] = _mh_coverage_gate
        agg['minhash_breadth_gate'] = _breadth_gate

        out[gval] = agg

    # ========== RPM-NORMALIZED MINHASH REDUCTION ==========
    # Within each parent group (toplevelkey), normalize minhash_reduction
    # by RPM share so that high-abundance organisms are boosted and
    # low-abundance siblings are dampened.
    #
    # Example: monkeypox (20 reads, RPM=X) vs 4 other orthopox (4 reads each)
    #   monkeypox rpm_share ≈ 0.56 → minhash_reduction boosted
    #   each other  rpm_share ≈ 0.11 → minhash_reduction dampened
    #
    # Formula:  minhash_reduction *= (floor + (1 - floor) * rpm_share)
    #   floor=0.3 ensures low-RPM organisms aren't zeroed out entirely

    _rpm_norm_floor = 0.3  # minimum scaling factor for lowest-RPM member

    # 1) Bucket the output entries by their parent group
    parent_buckets: Dict[str, List[str]] = {}
    for gval, agg in out.items():
        parent = agg.get("toplevelkey") or gval
        parent_buckets.setdefault(str(parent), []).append(gval)

    # 2) For each parent group, compute RPM shares and adjust minhash_reduction
    _solo_boost_strength = 0.9  # max boost toward 1.0 for a solo organism

    for parent, member_keys in parent_buckets.items():
        if len(member_keys) <= 1:
            # ----------------------------------------------------------
            # SOLO EXCLUSIVITY BOOST
            # ----------------------------------------------------------
            # Being the only organism in a group means no siblings are
            # competing for these reads — that's a positive signal.
            # Boost minhash_reduction toward 1.0, scaled by rpm_confidence_weight
            # (sigmoid on read_fraction) so the boost is proportional to how
            # confident we are the reads are real.
            #
            # Formula:
            #   gap       = 1.0 - raw_minhash        (room to grow)
            #   boost     = gap * strength * conf     (how much to fill that gap)
            #   final     = raw + boost
            #
            # At 5 reads / 200k total (read_fraction=2.5e-5, conf≈0.71):
            #   raw=0.5 → final = 0.5 + 0.5*0.4*0.71 = 0.642
            # At 5 reads / 10M total  (read_fraction=5e-7,  conf≈0.01):
            #   raw=0.5 → final = 0.5 + 0.5*0.4*0.01 = 0.502  (barely moves)
            # ----------------------------------------------------------
            mk = member_keys[0]
            raw_minhash = _to_float(out[mk].get("minhash_reduction", 0))
            conf = _to_float(out[mk].get("rpm_confidence_weight", 0))

            # ── Respect the minhash_confidence gate ──────────────────
            # The confidence gate (set earlier from coverage + gini)
            # caps how high minhash_reduction can go.  Without this,
            # the solo boost pushes low-coverage organisms (e.g.
            # Toxoplasma with 82K reads but 0.15% coverage) right
            # back to ~1.0, undoing the gate entirely.
            #
            # ceiling = minhash_confidence: the maximum the boosted
            # value should reach.  For good-coverage organisms
            # (confidence ≈ 0.94), the boost can push from gated
            # toward 0.94 — a modest helpful lift.  For low-coverage
            # noise (confidence ≈ 0.08), the boost is capped near
            # the gated value and can't inflate the score.
            _mh_conf = _to_float(out[mk].get("minhash_confidence", 1.0))
            ceiling = _mh_conf  # coverage evidence caps the max

            gap = max(0.0, ceiling - raw_minhash)
            boost = gap * _solo_boost_strength * conf
            adjusted = raw_minhash + boost

            out[mk]["minhash_reduction_raw"] = raw_minhash
            out[mk]["minhash_reduction"] = min(ceiling, adjusted)
            out[mk]["rpm_share_in_group"] = 1.0
            out[mk]["rpm_norm_scale"] = 1.0
            out[mk]["solo_exclusivity_boost"] = boost
            continue

        # ----------------------------------------------------------
        # MULTI-MEMBER: RPM-share normalization (as before)
        # ----------------------------------------------------------
        # Gather RPM values for each member in this group
        rpms = []
        for mk in member_keys:
            r = _to_float(out[mk].get("rpm", 0))
            rpms.append((mk, r))

        total_rpm = sum(r for _, r in rpms)
        if total_rpm <= 0:
            continue

        # Identify the dominant member (highest RPM in the group)
        sorted_rpms = sorted(rpms, key=lambda x: x[1], reverse=True)
        top_mk = sorted_rpms[0][0]

        for mk, r in rpms:
            rpm_share = r / total_rpm
            _mh_conf = _to_float(out[mk].get("minhash_confidence", 1.0))
            raw_minhash = _to_float(out[mk].get("minhash_reduction", 0))
            out[mk]["minhash_reduction_raw"] = raw_minhash

            if mk == top_mk:
                # ----------------------------------------------------------
                # DOMINANT MEMBER: boost minhash toward the coverage-gate
                # ceiling rather than scaling it down.
                # In high-ANI conflict groups (e.g. O104:H4 vs Shigella), the
                # minhash signal is similar for all members because they share
                # k-mers.  The organism with the most reads is the true hit —
                # reward it by filling the gap to its confidence ceiling,
                # proportional to how dominant it is.
                #
                # Example: O104:H4 rpm_share=0.60, raw=0.90, conf=0.95
                #   gap=0.05, boost=0.03 → final=0.93  (was 0.72 before)
                # ----------------------------------------------------------
                gap = max(0.0, _mh_conf - raw_minhash)
                boost = gap * rpm_share
                adjusted = min(_mh_conf, raw_minhash + boost)
                scale = (adjusted / raw_minhash) if raw_minhash > 0 else 1.0
            else:
                # ----------------------------------------------------------
                # LOSER: scale down aggressively, amplified by the coverage
                # gap vs the dominant member.
                #
                # Base formula (RPM-share):
                #   scale = floor + (1 - floor) * rpm_share
                #   floor=0.10 so Shigella with low share gets a hit.
                #
                # Coverage-ratio penalty (new):
                #   If the dominant organism has 20% coverage and this one
                #   has 8%, the coverage ratio = 8/20 = 0.4.  We blend this
                #   into the scale so that even organisms with decent RPM
                #   share get penalised when their coverage is much lower
                #   than the dominant — a strong signal that shared reads
                #   were misassigned.
                #
                #   coverage_ratio = my_coverage / top_coverage
                #   scale *= (cov_floor + (1 - cov_floor) * coverage_ratio)
                #   cov_floor=0.25 so 0% coverage → 25% of base scale
                #
                # Example: Shigella dysenteriae rpm_share=0.21, cov=8%, top_cov=20%
                #   base_scale = 0.10 + 0.90*0.21 = 0.289
                #   cov_ratio  = 8/20 = 0.4
                #   cov_scale  = 0.25 + 0.75*0.4 = 0.55
                #   final_scale = 0.289 * 0.55 = 0.159
                # ----------------------------------------------------------
                _loser_floor = 0.10
                base_scale = _loser_floor + (1.0 - _loser_floor) * rpm_share

                # Coverage-ratio penalty: compare this member's breadth to
                # the dominant member's breadth within the group.
                _top_cov = _to_float(out[top_mk].get("coverage", 0))
                _my_cov  = _to_float(out[mk].get("coverage", 0))
                if _top_cov > 0:
                    _cov_ratio = min(1.0, _my_cov / _top_cov)
                else:
                    _cov_ratio = 1.0  # can't penalise if dominant also has 0
                _cov_floor = 0.25
                _cov_scale = _cov_floor + (1.0 - _cov_floor) * _cov_ratio

                scale = base_scale * _cov_scale
                adjusted = raw_minhash * scale

            out[mk]["minhash_reduction"] = adjusted
            out[mk]["rpm_share_in_group"] = rpm_share
            out[mk]["rpm_norm_scale"] = scale

    return out


def sibling_disparity_from_group_reads(
    group_reads,
    target_key,
    reads_field="reads",
    key_field="key",
):
    """
    Compute sibling disparity metrics for `target_key` within `group_reads`.

    group_reads: list[dict], each like {"reads": <num>, "key": <id>}  (includes itself)
    target_key: the key to score (matched against dict[key_field])

    Returns dict:
      - n
      - total_reads
      - own_reads
      - proportion              = own_reads / total_reads
      - rank_index              = 0 is top
      - rank_disparity          = 1 - idx/(n-1)  (top=1, bottom=0; 0 if n<=1)
      - top_reads
      - top_ratio               = own_reads / top_reads
      - sorted                  = list of (reads, key) sorted desc
    """
    # normalize / coerce reads and keys
    items = []
    for d in group_reads:
        if key_field not in d:
            continue
        k = d[key_field]
        r = d.get(reads_field, 0) or 0
        try:
            r = float(r)
        except Exception:
            r = 0.0
        items.append((r, k))

    if not items:
        return {
            "n": 0, "total_reads": 0.0, "own_reads": 0.0,
            "proportion": 0.0, "rank_index": -1, "rank_disparity": 0.0,
            "top_reads": 0.0, "top_ratio": 0.0, "sorted": []
        }

    # sort descending by reads
    items_sorted = sorted(items, key=lambda x: x[0], reverse=True)
    n = len(items_sorted)
    total = sum(r for r, _ in items_sorted)
    top_reads = items_sorted[0][0] if n else 0.0

    # find target position + reads
    rank_index = -1
    own_reads = 0.0
    for i, (r, k) in enumerate(items_sorted):
        if k == target_key:
            rank_index = i
            own_reads = r
            break

    # disparity (rank-based)
    if n <= 1 or rank_index < 0:
        rank_disparity = 1.0
    else:
        rank_disparity = 1.0 - (rank_index / (n - 1))

    proportion = (own_reads / total) if total > 0 else 0.0
    top_ratio = (own_reads / top_reads) if top_reads > 0 else 0.0

    return {
        "n": n,
        "total_reads": total,
        "own_reads": own_reads,
        "proportion": proportion,
        "rank_index": rank_index,
        "rank_disparity": rank_disparity,
        "top_reads": top_reads,
        "top_ratio": top_ratio,
        "sorted": items_sorted,
    }

def compute_scores_per(
    data = {},
    reward_factor = 2,
    dispersion_factor = 0.9,
    alpha = 1,
    comparison_df = pd.DataFrame(),
    fallback_top = "Unknown",
    total_reads = 0,
    mapq_breadth_power = 2.0,
    breadth_midpoint = 0.01,
    breadth_steepness = 12_000,
    abundance_rpm_midpoint = 5.0,
    abundance_rpm_steepness = 2.0,
):
    data['_mapq_breadth_power'] = mapq_breadth_power
    if len(data.get('covered_regions', [])) > 0:
        gini_strain = getGiniCoeff(
            data['covered_regions'],
            data['length'],
            alpha=alpha,
            reward_factor=reward_factor,
            beta=dispersion_factor
        )
    else:
        gini_strain = 0
    # -------------------------
    # RPM / RPKM (TPKM) metrics
    # -------------------------

    reads_mapped = data.get('numreads', 0)
    ref_length_bp = data.get('length', 0)
    total_reads_millions = total_reads / 1e6 if total_reads > 0 else 0
    ref_length_kb = ref_length_bp / 1e3 if ref_length_bp > 0 else 0

    # RPM: reads per million mapped reads
    if total_reads_millions > 0:
        data['rpm'] = reads_mapped / total_reads_millions
    else:
        data['rpm'] = 0

    # RPKM / TPKM
    if total_reads_millions > 0 and ref_length_kb > 0:
        data['rpkm'] = reads_mapped / (total_reads_millions * ref_length_kb)
    else:
        data['rpkm'] = 0

    col_stat2 = 'Δ All%'
    col_stat = 'Δ^-1 Breadth'
    # -------------------------
    # Minhash block (now rpm is available)
    # -------------------------
    reads_mapped = data.get('numreads', 0)

    # Fraction of total reads hitting this reference
    read_fraction = reads_mapped / total_reads if total_reads > 0 else 0
    data['read_fraction'] = read_fraction

    rpm_weight = rpm_confidence_weight(read_fraction, k=50_000, midpoint=0.0001)
    data['rpm_confidence_weight'] = rpm_weight
    data['strainname'] = data.get('strainname', fallback_top)
    data['gini_coefficient'] = gini_strain

    # Covered bases
    data['covered_bases'] = sum(
        region[1] - region[0] + 1 for region in data.get('covered_regions', [])
    )

    # MAPQ
    mapq = data.get('meanmapq', 0)
    data['mapq_score'] = normalize_mapq(mapq)

    # Breadth — scaled by high-MAPQ read fraction
    # If most reads are MAPQ=0 (ambiguous), breadth is unreliable and gets crushed.
    # highmapq_fraction is the proportion of reads with MAPQ >= threshold (e.g. 5).
    # mapq_breadth_power controls how aggressively low-MAPQ penalizes breadth.
    coverage = data.get('coverage', 0)
    hmf = float(data.get('highmapq_fraction', 1.0))
    _mbp = float(data.get('_mapq_breadth_power', 2.0))
    mapq_scale = hmf ** _mbp  # e.g. 0.1^2 = 0.01 → breadth almost zeroed
    data['breadth_log_score'] = breadth_score_sigmoid(
        coverage, midpoint=breadth_midpoint, steepness=breadth_steepness
    ) * mapq_scale
    data['highmapq_fraction'] = hmf

    # ── Low-abundance confidence (sterile-site boost) ────────────────────
    # Uses log-RPM sigmoid to score organisms that are meaningful even at
    # low read counts.  Stored as a feature; actual weighting happens in
    # compute_tass_score / compute_tass_score_from_metrics.
    data['abundance_confidence'] = low_abundance_confidence(
        numreads=data.get('numreads', 0),
        total_reads=total_reads,
        genome_length_bp=data.get('length', 1),
        rpm_midpoint=abundance_rpm_midpoint,
        rpm_steepness=abundance_rpm_steepness,
    )
    if not comparison_df.empty:
        # Look up by subkey (species) — comparison_df is now aggregated to
        # the subkey level so all accessions in the same species share one
        # composite minhash signal.
        lookup_key = str(data.get('subkey', data.get('accession', ''))).strip()
        if lookup_key in comparison_df.index:
            c1 = float(comparison_df.loc[lookup_key, col_stat])
            d_all = float(comparison_df.loc[lookup_key, col_stat2])

            k_sig = 0.90
            x0 = -10.0
            pen = 1.0 / (1.0 + math.exp(-k_sig * (d_all - x0)))

            comparison_value = min(1.0, c1 * pen)
            data['minhash_score'] = comparison_value
        else:
            data['minhash_score'] = rpm_weight * 0.5  # default to a moderate score if no comparison available
    # ── Confidence-gated minhash_reduction ─────────────────────────────
    # Problem: minhash_score can be 1.0 even when coverage is negligible
    # (e.g. Toxoplasma at 0.001% coverage).  For sterile/blood samples
    # where minhash_weight is high, this inflates the TASS score for
    # what are likely contaminants or false positives.
    #
    # Solution: scale minhash_score by a "coverage confidence" factor
    # derived from breadth and gini.  This preserves full minhash power
    # for legitimate high-ANI conflicts (Shigella vs E. coli) that have
    # real coverage, while crushing it for noise hits.
    #
    # confidence = breadth_w * breadth_sigmoid(coverage) + gini_w * gini
    #   - breadth_sigmoid(coverage): primary gate — is this organism
    #     actually present?  At 0.001% coverage → ~0;  at 5%+ → ~1
    #   - gini: secondary signal — are reads uniformly distributed?
    #     Low gini = reads clumped in a tiny region = suspicious
    #
    # The raw coverage sigmoid is used (not breadth_log_score) to avoid
    # double-penalizing through mapq_scale, which is already captured in
    # the breadth_log_score component of TASS.
    raw_minhash = data.get('minhash_score', 1)
    # data['minhash_reduction'] = minhash_confidence * raw_minhash
    # minhash_confidence = min(1.0, max(0.0, raw_minhash))

    cov_conf = breadth_score_sigmoid(coverage)          # 0→1 based on coverage %
    gini_val = float(data.get('gini_coefficient', 0))   # 0→1 uniformity
    _mcg_breadth_w = 1   # how much coverage evidence matters
    _mcg_gini_w    = 0.0   # how much distribution uniformity matters
    minhash_confidence = (_mcg_breadth_w * cov_conf) + (_mcg_gini_w * gini_val)
    minhash_confidence = min(1.0, max(0.0, minhash_confidence))  # sanity clamp

    # ── Low-read penalty ────────────────────────────────────────────────
    # Organisms with very few reads are more likely to be false positives.
    # Scale minhash_confidence by a read-count factor so that low-read
    # organisms get a harder reduction.  rpm_weight (sigmoid on
    # read_fraction) is already computed above; we blend it in with a
    # floor so we never completely zero out an organism that has decent
    # coverage but happens to have few absolute reads.
    #
    # _read_penalty_w controls strength: 0.0 = no penalty, 1.0 = full
    # penalty (minhash_confidence *= rpm_weight).  0.35 gives a moderate
    # dampening: an organism at rpm_weight=0.1 sees confidence scaled by
    # 0.65 + 0.35*0.1 = 0.685 (≈30% reduction).
    #
    # Coverage bypass: when cov_conf is high (organism has solid breadth),
    # the low-read penalty is faded out proportionally.  This prevents
    # a dominant high-coverage organism in a conflict group from being
    # pulled down by the global read-fraction sigmoid.
    #   cov_conf=0.95 → penalty_w fades to 0.35*(1-0.95)=0.017 (nearly off)
    #   cov_conf=0.10 → penalty_w stays at 0.35*(1-0.10)=0.315 (nearly full)
    _read_penalty_w = 0.35 * (1.0 - cov_conf)
    minhash_confidence *= (1.0 - _read_penalty_w) + (_read_penalty_w * rpm_weight)

    data['minhash_reduction'] = minhash_confidence
    data['minhash_confidence'] = minhash_confidence  # store for debugging/reporting

    return data

def rpm_confidence_weight(read_fraction, k=50_000, midpoint=0.0001):
    """
    read_fraction = reads_mapped / total_reads  (value between 0 and 1)
    midpoint: fraction at which confidence = 0.5
              0.0001 = 0.01% of total reads (tune to your typical noise floor)
    k: steepness — needs to be large since fractions are tiny
    """
    return 1.0 / (1.0 + math.exp(-k * (read_fraction - midpoint)))


def breadth_score_sigmoid(coverage, midpoint=0.01, steepness=12_000):
        """Sigmoid mapping of genome coverage fraction → [0, 1].

        Parameters
        ----------
        coverage : float
            Fraction of the genome covered (0–1).
        midpoint : float
            Coverage fraction at which the sigmoid returns 0.5.
            Lower values make the curve sensitive to very low coverage
            (useful for sterile/blood sites).  Default 0.01 (1%).
        steepness : float
            Controls how sharply the sigmoid transitions.  When lowering
            *midpoint*, increase *steepness* proportionally to keep the
            curve tight (e.g. midpoint=0.001 → steepness≈120 000).
        """
        x = steepness * (coverage - midpoint)
        if x >= 50:
            return 1.0
        if x <= -50:
            return 0.0
        return 1.0 / (1.0 + math.exp(-x))

def low_abundance_confidence(numreads, total_reads, genome_length_bp,
                              rpm_midpoint=5.0, rpm_steepness=2.0):
    """Score that rewards organisms whose RPM is meaningful given sequencing
    depth, even when absolute read counts are very low.

    Operates in **log₁₀-RPM** space so the sigmoid is not crushed by the
    huge dynamic range of metagenomic read counts.

    Parameters
    ----------
    numreads : int
        Reads mapped to this organism.
    total_reads : int
        Total reads in the sample (denominator for RPM).
    genome_length_bp : int
        Reference genome length in base-pairs (used for optional RPKM
        awareness — currently kept simple with RPM only).
    rpm_midpoint : float
        RPM at which the sigmoid returns 0.5.  For sterile/blood sites
        where even 5 RPM is significant, use 1.0–5.0.  For high-biomass
        sites (gut, skin) where 50+ RPM is expected, use 50–200.
    rpm_steepness : float
        Steepness of the log₁₀-RPM sigmoid.  Default 2.0 gives a gentle
        curve that reaches ~0.95 about one order of magnitude above
        *rpm_midpoint*.

    Returns
    -------
    float in [0, 1]
    """
    if total_reads <= 0 or numreads <= 0:
        return 0.0

    rpm = (numreads / total_reads) * 1e6
    # Work in log-space so the sigmoid is not crushed by tiny fractions.
    # log10(5) ≈ 0.70,  log10(50) ≈ 1.70
    log_rpm = math.log10(max(rpm, 1e-3))
    log_mid = math.log10(max(rpm_midpoint, 1e-3))

    x = rpm_steepness * (log_rpm - log_mid)
    if x >= 50:
        return 1.0
    if x <= -50:
        return 0.0
    return 1.0 / (1.0 + math.exp(-x))
def calculate_mmbert_prob(
    mmbert_dict = {},
    taxid = None
):
    if taxid in mmbert_dict:
        avg = mmbert_dict[taxid].get('avg')
        model = mmbert_dict[taxid].get('model')

        if avg is not None:
            mmbert_avgs = avg
        if model is not None:
            mmbert_models = model
        # print(f"MicrobeRT results for {taxid}: Avg Probability: {mmbert_avgs}, Models: {(mmbert_models) if mmbert_models else 'N/A'}")
        return mmbert_avgs, mmbert_models
    return None, None

def calculate_mean(diamond_list, key):
    """
    Calculate the weighted mean of a key, weighted by 'cds' values.

    Parameters:
    diamond_list (list of dict): List of dictionaries containing the data.
    key (str): The key for which the weighted mean is to be calculated.

    Returns:
    float: The weighted mean of the values.
    """
    total_weighted_value = 0
    total_cds = 0
    # Loop through each item in the list and calculate the weighted value
    # if diamond_list is a list:
    if isinstance(diamond_list, list):
        for item in diamond_list:
            value = float(item.get(key, 0))  # Convert the key value to a float
            cds = int(item.get('cds', 0))    # Get the cds value (as the weight)

            # Add to the weighted sum
            total_weighted_value += value * cds
            total_cds += cds
    else:
        value = float(diamond_list.get(key, 0))  # Convert the key value to a float
        cds = int(diamond_list.get('cds', 0))    # Get the cds value (as the weight)

        # Add to the weighted sum
        total_weighted_value += value * cds
        total_cds += cds
    # Return the weighted mean
    return total_weighted_value / total_cds if total_cds != 0 else 0

def calculate_aggregate_scores(
    data = {},
    hmp_dists = {},
    body_sites = [],
    sampletype = "Sterile",
    k2_mapping = {},
    group_reads = defaultdict(int),
    mmbert_dict = {},
    dmnd = [],
    min_cds_found = 5,
    total_reads = 0,
):
    zscore, hmp_percentile, hmp_info = calculate_hmp_percentile(
        value = data,
        dists = hmp_dists,
        body_sites = body_sites,
        sampletype= sampletype,
        total_reads = total_reads,
    )
    data['mmbert'], data['mmbert_model'] = calculate_mmbert_prob(
        mmbert_dict = mmbert_dict,
        taxid = data.get('key', None)
    )
    data['zscore'] = zscore
    data['hmp_percentile'] = hmp_percentile
    data['hmp_norm_abundance'] = hmp_info.get('hmp_norm_abundance', 0)
    data['hmp_mean'] = hmp_info.get('hmp_mean', 0)
    data['hmp_std'] = hmp_info.get('hmp_std', 0)
    data['observed_abundance'] = hmp_info.get('observed_abundance', 0)
    data['hmp_site_count'] = hmp_info.get('hmp_site_count', 0)
    data['hmp_num_samples'] = hmp_info.get('hmp_num_samples', 0)
    data['k2_reads'] = k2_mapping.get(data.get('key', None), {}).get('clades_covered', 0)
    data['k2_disparity_score'] = calculate_k2_reads_disparity(
        data = data,
        k2_mapping = k2_mapping
    )
    # get the k2 reads using the taxid
    metrics = sibling_disparity_from_group_reads(group_reads, target_key=data.get('key', None))
    data['disparity'] = metrics.get('rank_disparity', 1.0)
    key = data.get('key', None)
    if data['key'] in dmnd:
        dmnd[key]['maxvalereached'] = dmnd[key].get('cds', 0) > min_cds_found
        dmnd_results = dmnd[key]
        mean_identity = calculate_mean(dmnd_results, 'identity')
        mean_length = calculate_mean(dmnd_results, 'lengthmedian')
        mean_mismatched = calculate_mean(dmnd_results, 'mismatchedmedian')
        mean_evalue = calculate_mean(dmnd_results, 'medianevalue')
        # mean_contigs = calculate_mean(dmnd_results, 'contigs')
        mean_cds = calculate_mean(dmnd_results, 'cds')
        max_val_reached =  mean_cds > min_cds_found
        dmnd_results = {
            'identity': mean_identity,
            'lengthmean': mean_length,
            'mismatchedmean': mean_mismatched,
            'meanevalue': mean_evalue,
            # 'contigs': mean_contigs,
            'cds': mean_cds,
            'maxvalereached': max_val_reached
        }
        data['diamond'] = dmnd_results
        data['diamond_identity'] = mean_identity
    return data



def calculate_siblings_score (
    data: dict,
    group_key: str = "toplevelkey",
    reads_key: str = "numreads",
    out_prefix: str = "siblings_",
    include_unknown: bool = False,
    unknown_label: str = "Unknown",
):
    """
    Mutates and returns `data` by adding sibling-based metrics computed from `reads_key`,
    comparing entries within the same `group_key`.

    Assumes:
      data[id] = { ..., "group": <group>, "numreads": <number>, ... }

    Adds (per entry):
      - {prefix}n
      - {prefix}total_reads
      - {prefix}proportion
      - {prefix}rank_disparity
      - {prefix}top_ratio
      - {prefix}rank_index
    """

    # 1) Build groups -> list of (id, reads)
    groups = defaultdict(list)
    for key, rec in data.items():
        grp = rec.get(group_key, unknown_label)
        if (not include_unknown) and (grp == unknown_label or grp is None or grp == ""):
            continue

        val = rec.get(reads_key, 0) or 0
        try:
            val = float(val)
        except Exception:
            val = 0.0

        groups[grp].append((key, val))

    # 2) For each group, compute sibling metrics
    for grp, items in groups.items():
        # items: list of (id, reads)
        # sort descending by reads
        items_sorted = sorted(items, key=lambda x: x[1], reverse=True)
        n = len(items_sorted)
        total = sum(v for _, v in items_sorted)
        top = items_sorted[0][1] if n else 0.0

        # map id -> rank index
        rank_index = {k: i for i, (k, _) in enumerate(items_sorted)}

        for k, own in items_sorted:
            # rank disparity: top=1, bottom=0
            if n <= 1:
                disparity = 1.0
            else:
                idx = rank_index[k]
                disparity = 1.0 - (idx / (n - 1))

            proportion = (own / total) if total > 0 else 0.0
            top_ratio = (own / top) if top > 0 else 0.0

            data[k][f"{out_prefix}n"] = n
            data[k][f"{out_prefix}total_reads"] = total
            data[k][f"{out_prefix}proportion"] = proportion
            data[k][f"{out_prefix}rank_disparity"] = disparity
            data[k][f"{out_prefix}top_ratio"] = top_ratio
            data[k][f"{out_prefix}rank_index"] = rank_index[k]

    return data

def compute_tass_score_from_metrics(metrics_df, breadth_w, minhash_w, gini_w,
                                    disparity_w=0.0, hmp_w=0.0, alpha=1.0,
                                    plasmid_bonus_w=0.0,
                                    abundance_confidence_w=0.0,
                                    abundance_gate=False,
                                    score_power=1.0,
                                    tass_mode="additive"):

    b = metrics_df["breadth_log_score"].to_numpy(float)
    m = metrics_df["minhash_reduction"].to_numpy(float)
    g = metrics_df["gini_coefficient"].to_numpy(float)
    d = metrics_df["disparity_score"].to_numpy(float) if "disparity_score" in metrics_df else 0.0
    h = metrics_df["hmp_percentile"].to_numpy(float) if "hmp_percentile" in metrics_df else 0.0

    if tass_mode == "penalized":
        # ── PENALIZED MODE ────────────────────────────────────────────────
        # Baseline 0.5; core metrics push score up (good signal) or down
        # (poor signal).  Each metric in [0,1] contributes:
        #   weight * (metric - 0.5)  →  range [-0.5*w, +0.5*w]
        # The weighted sum of deviations is scaled by alpha and added to 0.5.
        # This means an organism with ALL metrics at 0 gets penalized hard
        # (score ≈ 0.0) and one with ALL metrics at 1 gets score ≈ 1.0.
        # An organism with no signal (all zeros) scores:
        #   0.5 + alpha * sum(w_i * (0 - 0.5)) = 0.5 - 0.5*alpha ≈ 0.0

        core_signal = (breadth_w * (b - 0.5)
                       + minhash_w * (m - 0.5)
                       + gini_w * (g - 0.5)
                       + disparity_w * (d - 0.5)
                       + hmp_w * (h - 0.5))
        core_signal = alpha * core_signal

        # Raw core score (before bonuses): centered at 0.5
        score = 0.5 + core_signal

        # ── Abundance confidence: multiplicative scaler on deviation ──────
        # Instead of adding a free bonus, AC scales how far the score can
        # deviate from 0.5.  High AC (≈1) → full deviation preserved.
        # Low AC (≈0) → score collapses back toward 0.5.
        # This prevents organisms with bad core metrics from being inflated
        # by high abundance confidence.  An organism with core=0.046 gets:
        #   deviation = 0.046 - 0.5 = -0.454  (penalized)
        #   with AC=0.84: deviation stays -0.454 (already penalized, AC keeps it)
        # vs additive mode: 0.046 + 0.3*0.84 = 0.298 (inflated!)
        if abundance_confidence_w > 0 and "abundance_confidence" in metrics_df.columns:
            ac = metrics_df["abundance_confidence"].to_numpy(float)
            # Scale the deviation by a blend of 1.0 and AC, controlled by weight.
            # When abundance_confidence_w=0: deviation unchanged (blend=1.0)
            # When abundance_confidence_w=1: deviation fully scaled by AC
            # ac_scaler in [0, 1] range
            ac_scaler = (1.0 - abundance_confidence_w) + abundance_confidence_w * ac
            deviation = score - 0.5
            score = 0.5 + deviation * ac_scaler

        # ── Plasmid bonus: gated by core metric quality ───────────────────
        # Only grant plasmid bonus if the organism has meaningful core signal.
        # "Meaningful" = core weighted score (before centering) > 0.15.
        # This prevents FP organisms with gini=0.04, minhash=0.004 from
        # getting free plasmid boost.
        if plasmid_bonus_w > 0 and "plasmid_score" in metrics_df.columns:
            ps = metrics_df["plasmid_score"].to_numpy(float)
            # Core quality = raw weighted sum (same as additive mode score)
            core_quality = breadth_w * b + minhash_w * m + gini_w * g + disparity_w * d + hmp_w * h
            core_quality = alpha * core_quality
            # Gate: only apply bonus where core quality > 0.15
            plasmid_gate = np.where(core_quality > 0.15, 1.0, core_quality / 0.15)
            score = score + plasmid_bonus_w * ps * plasmid_gate

        # ── Multiplicative abundance gate (optional, same as additive) ────
        if abundance_gate and "abundance_confidence" in metrics_df.columns:
            gate = metrics_df["abundance_confidence"].to_numpy(float)
            # In penalized mode, gate dampens deviation from 0.5
            deviation = score - 0.5
            score = 0.5 + deviation * gate

        # No score_power in penalized mode — the baseline-centered approach
        # naturally produces a well-distributed range.

        return np.clip(score, 0.0, 1.0)

    # ── ADDITIVE MODE (original behavior) ─────────────────────────────────
    score = breadth_w * b + minhash_w * m + gini_w * g + disparity_w * d + hmp_w * h
    score = alpha * score

    # Plasmid bonus: additive, outside normalized weights
    if plasmid_bonus_w > 0 and "plasmid_score" in metrics_df.columns:
        ps = metrics_df["plasmid_score"].to_numpy(float)
        score = score + plasmid_bonus_w * ps

    # Abundance confidence bonus: additive, outside normalized weights.
    # Boosts organisms whose RPM is meaningful even at low absolute read
    # counts (particularly useful for sterile/blood sites).
    if abundance_confidence_w > 0 and "abundance_confidence" in metrics_df.columns:
        ac = metrics_df["abundance_confidence"].to_numpy(float)
        score = score + abundance_confidence_w * ac

    # ── Multiplicative abundance gate ─────────────────────────────────────
    # When enabled, multiply the entire score by the abundance_confidence
    # sigmoid (0→1).  Crushes noise organisms with trivially low RPM while
    # leaving real detections untouched.
    if abundance_gate and "abundance_confidence" in metrics_df.columns:
        gate = metrics_df["abundance_confidence"].to_numpy(float)
        score = score * gate

    # ── Power transform (score recalibration) ─────────────────────────────
    # When score_power < 1, compresses high scores and lifts low scores:
    #   score_power=0.5 → 0.09 becomes 0.30,  0.95 stays 0.97
    #   score_power=0.3 → 0.09 becomes 0.52,  0.95 stays 0.98
    # Preserves monotonic ordering so thresholds still separate TP/FP.
    # score_power=1.0 (default) is a no-op.
    #
    # score_power_scale (per-organism, 0→1) modulates the strength of the
    # transform.  Dominant organisms in their ANI/toplevelkey group get the
    # full boost; minor siblings sharing reads with many close relatives
    # get almost none.
    #   effective_power = 1.0 - (1.0 - score_power) * score_power_scale
    if score_power != 1.0 and score_power > 0:
        score = np.clip(score, 0.0, 1.0)
        if "score_power_scale" in metrics_df.columns:
            sp_scale = metrics_df["score_power_scale"].to_numpy(float)
            effective_power = 1.0 - (1.0 - score_power) * sp_scale
            score = np.power(score, effective_power)
        else:
            score = np.power(score, score_power)

    return np.clip(score, 0.0, 1.0)

def compute_tass_score(data = {}, weights={}):
    """
    count is a dictionary that might look like:
    {
      'normalized_disparity': <some_value>,
      'mapq_score': <some_value>,
      'meangini': <some_value>,
      ...
    }
    We apply the known formula for TASS Score using the provided weights.
    """
    _tass_mode = weights.get('tass_mode', 'additive')

    if _tass_mode == 'penalized':
        # ── PENALIZED MODE ────────────────────────────────────────────────
        # Baseline 0.5; core metrics push score up or down.
        # Metric contribution: weight * (metric - 0.5)

        _bw = float(weights.get('breadth_weight', 0))
        _mw = float(weights.get('minhash_weight', 0))
        _gw = float(weights.get('gini_weight', 0))
        _dw = float(weights.get('disparity_weight', 0))
        _hw = float(weights.get('hmp_weight', 0))

        _b = float(data.get('breadth_log_score', 0))
        _m = float(data.get('minhash_reduction', 0))
        _g = float(data.get('gini_coefficient', 0))
        _d = float(data.get('disparity', 0))
        _h = float(data.get('hmp_percentile', 0))

        core_signal = (_bw * (_b - 0.5)
                       + _mw * (_m - 0.5)
                       + _gw * (_g - 0.5)
                       + _dw * (_d - 0.5)
                       + _hw * (_h - 0.5))

        # Also include minor metrics if they have weights
        _mapq_w = float(weights.get('mapq_score', 0))
        _mapq = float(data.get('mapq_score', 0))
        if _mapq_w > 0:
            core_signal += _mapq_w * (_mapq - 0.5)

        _k2w = float(weights.get('k2_disparity_score_weight', 0))
        _k2 = float(data.get('k2_disparity_score', 0))
        if _k2w > 0:
            core_signal += _k2w * (_k2 - 0.5)

        _diw = float(weights.get('diamond_identity', 0))
        _di = float(data.get('diamond', {}).get('identity', 0))
        if _diw > 0:
            core_signal += _diw * (_di - 0.5)

        tass_score = 0.5 + core_signal

        # ── Abundance confidence: multiplicative scaler on deviation ──────
        _acw = float(weights.get('abundance_confidence_weight', 0))
        _ac = float(data.get('abundance_confidence', 0))
        if _acw > 0:
            ac_scaler = (1.0 - _acw) + _acw * _ac
            deviation = tass_score - 0.5
            tass_score = 0.5 + deviation * ac_scaler

        # ── Plasmid bonus: gated by core quality ─────────────────────────
        _plasmid_score = float(data.get('plasmid_score', 0))
        _plasmid_bonus_w = float(weights.get('plasmid_bonus_weight', 0))
        if _plasmid_score > 0 and _plasmid_bonus_w > 0:
            # Raw core quality = weighted sum without centering
            core_quality = (_bw * _b + _mw * _m + _gw * _g
                            + _dw * _d + _hw * _h)
            # Gate: ramp from 0 at core_quality=0 to 1 at core_quality=0.15
            plasmid_gate = min(1.0, core_quality / 0.15) if core_quality < 0.15 else 1.0
            tass_score += _plasmid_score * _plasmid_bonus_w * plasmid_gate

        # ── Multiplicative abundance gate (optional) ──────────────────────
        if weights.get('abundance_gate', False):
            gate = float(data.get('abundance_confidence', 0.0))
            deviation = tass_score - 0.5
            tass_score = 0.5 + deviation * gate

        # No score_power in penalized mode

        return min(1.0, max(0.0, tass_score))

    # ── ADDITIVE MODE (original behavior) ─────────────────────────────────
    tass_score = sum([
        apply_weight(data.get('disparity', 0), weights.get('disparity_weight', 0)),
        apply_weight(data.get('minhash_reduction', 0),       weights.get('minhash_weight', 0)),
        apply_weight(data.get('gini_coefficient', 0),             weights.get('gini_weight', 0)),
        apply_weight(data.get('breadth_log_score', 0),             weights.get('breadth_weight', 0)),
        apply_weight(data.get('hmp_percentile', 0),             weights.get('hmp_weight', 0)),
        apply_weight(data.get('mapq_score', 0),       weights.get('mapq_score', 0)),
        apply_weight(data.get('k2_disparity_score', 0),         weights.get('k2_disparity_score_weight', 0)),
        apply_weight(data.get('diamond', {}).get('identity', 0),
                     weights.get('diamond_identity', 0)),
        # Low-abundance confidence: boosts organisms meaningful at low read
        # counts (sterile/blood sites).  Weight = 0 by default → no effect
        # unless explicitly enabled via --abundance_confidence_weight.
        apply_weight(data.get('abundance_confidence', 0),
                     weights.get('abundance_confidence_weight', 0)),
    ])

    # ── Plasmid bonus: additive boost outside normalized weight pool ──────
    _plasmid_score = float(data.get('plasmid_score', 0))
    _plasmid_bonus_w = float(weights.get('plasmid_bonus_weight', 0))
    if _plasmid_score > 0 and _plasmid_bonus_w > 0:
        tass_score += _plasmid_score * _plasmid_bonus_w

    # ── Multiplicative abundance gate ─────────────────────────────────────
    if weights.get('abundance_gate', False):
        gate = float(data.get('abundance_confidence', 0.0))
        tass_score *= gate

    # ── Power transform (score recalibration) ─────────────────────────────
    # score_power_scale (0→1) modulates how much the power transform
    # applies to this organism.  Dominant organisms in their ANI/toplevelkey
    # group (proportion≈1) get the full boost; minor siblings sharing reads
    # with many close relatives (proportion≈0.2) get almost no boost.
    #   effective_power = 1.0 - (1.0 - score_power) * score_power_scale
    # When score_power_scale=1: effective_power = score_power (full effect)
    # When score_power_scale=0: effective_power = 1.0         (no-op)
    _score_power = float(weights.get('score_power', 1.0))
    if _score_power != 1.0 and _score_power > 0:
        _sp_scale = float(data.get('score_power_scale', 1.0))
        _effective_power = 1.0 - (1.0 - _score_power) * _sp_scale
        tass_score = max(0.0, min(1.0, tass_score))
        tass_score = tass_score ** _effective_power

    return min(1.0, max(0.0, tass_score))

def calculate_hmp_percentile(
        value = {},
        dists = {},
        sampletype = "Sterile",
        body_sites = [],
        total_reads = 0
):
    abus = []
    taxid = (
        value.get('taxid')
        or value.get('key')
        or value.get('toplevelkey')
        or value.get('ref')
    )
    for body_site in body_sites:
        if not taxid:
            continue
        try:
            # get the key where it is the (taxid, body_site)
            taxid_int = int(float(taxid))
            keys = [
                (taxid_int, body_site),
                (str(taxid_int), body_site),
                (str(taxid), body_site),
            ]
            for k in keys:
                if k in dists:
                    abus.append(
                        dict(
                            norm_abundance = dists[k].get('norm_abundance', 0),
                            norm_stdev = dists[k].get('norm_stdev', 0),
                            std = dists[k].get('std', 0),
                            mean = dists[k].get('mean', 0),
                            site_count = dists[k].get('site_count', 0),
                            num_samples = len(dists[k].get('abundances', [])),
                        )
                    )
                    break
        except Exception as e:
            print(f"Error in taxid lookup for hmp: {e}")
    # Use read fraction (0–1 scale) to match HMP distribution mean/std which are also fractions
    if total_reads > 0:
        fraction_observed = float(value.get('numreads', 0) or 0) / total_reads
    else:
        # Fallback: use pre-computed read_fraction if total_reads wasn't passed
        fraction_observed = float(value.get('read_fraction', 0) or 0)
    sum_abus_expected = sum([x.get('mean', 0) for x in abus])  # means are already fractions (0–1)
    # sum_abus_expected = sum([x.get('norm_abundance', 0) for x in abus])  # means are already fractions (0–1)

    stdsum = sum([x.get('std', 0) for x in abus])
#     stdsum = sum([x.get('norm_stdev', 0) for x in abus])
    zscore = ((fraction_observed - sum_abus_expected) / stdsum) if stdsum > 0 else 3
    # zscore = ((percent_total_reads_observed -sum_norm_abu)/ stdsum)  if stdsum > 0 else 3
    base_percentile = norm.cdf(zscore)

    # Penalize below-typical values, keep neutral near the center,
    # and boost unusually high values.
    if zscore < -1.0:
        adjusted_percentile = base_percentile ** 2
    elif zscore <= 1.0:
        adjusted_percentile = base_percentile
    elif zscore <= 2.0:
        adjusted_percentile = 1.0 - (1.0 - base_percentile) ** 1.5
    else:
        adjusted_percentile = 1.0 - (1.0 - base_percentile) ** 2

    # Aggregate site_count and num_samples across matched body sites
    total_site_count = sum(x.get('site_count', 0) for x in abus)
    total_num_samples = sum(x.get('num_samples', 0) for x in abus)

    hmp_info = {
        'hmp_norm_abundance': sum_abus_expected,
        'hmp_mean': sum([x.get('mean', 0) for x in abus]),
        'hmp_std': stdsum,
        'observed_abundance': fraction_observed,
        'hmp_site_count': total_site_count,
        'hmp_num_samples': total_num_samples,
    }
    return zscore, adjusted_percentile, hmp_info

def json_safe(x):
    """Convert numpy/pandas/scalars/containers into JSON-serializable Python types."""
    if x is None:
        return None
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating,)):
        return float(x)
    if isinstance(x, (np.bool_,)):
        return bool(x)
    if isinstance(x, (pd.Timestamp,)):
        return x.isoformat()
    if isinstance(x, dict):
        return {str(k): json_safe(v) for k, v in x.items()}
    if isinstance(x, (list, tuple, set)):
        return [json_safe(v) for v in x]
    return x


def calculate_classes(rec, ref, pathogens, sample_type="Unknown", taxdump=None):
    """
    Calculate pathogen classification for a single record (aggregate or strain).
    Now includes taxonomic lineage traversal - if a taxid isn't found in pathogens dict,
    walks up the lineage tree (species → genus → family → ...) until a match is found.

    Args:
        rec: The record dictionary (aggregate or strain data)
        ref: The reference identifier (fallback for lookups)
        pathogens: Dictionary of pathogen information keyed by taxid
        sample_type: Sample type/body site (e.g., "blood", "stool", "nasal")
                    Will be normalized automatically
        taxdump: Optional taxonomy data from load_taxdump(). If provided, enables
                lineage traversal for pathogen lookup.

    Returns:
        dict: A dictionary with the record data plus classification fields:
            - high_cons: bool
            - is_pathogen: str (category)
            - microbial_category: str (category)
            - annClass: str (Direct/Derived/Mixed)
            - is_annotated: str (Yes/No)
            - status: str
            - ref: str (the reference identifier used)
            - matched_taxid: str (the taxid that matched in pathogens, may be parent)
            - matched_rank: str (the rank of the matched taxid, if lineage used)
            - normalized_sample_site: str (the normalized body site)
    """

    CATEGORY_PRIORITY = {
        "Primary": 5,
        "Opportunistic": 4,
        "Potential": 3,
        "Commensal": 2,
        "Unknown": 1,
    }

    # Normalize the sample type/body site
    normalized_sample_type = normalize_body_site(sample_type)

    # choose identifier for pathogen lookup
    taxid = ref

    def find_pathogen_in_lineage(target_taxid):
        """
        Search for pathogen info in the lineage, starting with the taxid itself.

        This function:
        1. First checks if target_taxid is directly in pathogens dict
        2. If not found and taxdump is provided, gets the lineage
        3. Walks up the lineage (species → genus → family → ...)
        4. Returns the first match found

        Returns:
            tuple: (refpath, matched_taxid, matched_rank) or (None, None, None)
        """
        # First, try the taxid directly
        refpath = pathogens.get(str(target_taxid))
        if refpath:
            return refpath, str(target_taxid), "direct"

        # If not found and we have taxdump, walk up the lineage
        if taxdump and str(target_taxid) in taxdump:
            lineage = get_lineage(str(target_taxid), taxdump)

            # lineage is [(taxid, rank), ...] from specific to general
            # e.g., [('198214', 'species'), ('28211', 'genus'), ('267890', 'family'), ...]
            for lin_taxid, lin_rank in lineage:
                refpath = pathogens.get(str(lin_taxid))
                if refpath:
                    return refpath, str(lin_taxid), lin_rank

        # Not found in lineage either
        return None, None, None
    # if members is an attribute, and not empty then iterate through all of them
    # Members may be subkey groups (which themselves have strain-level members)
    # or flat strains.  Recursively classify nested members first.
    if "members" in rec and rec["members"]:
        # Recursively classify nested members (subkey groups with their own members)
        for member in rec["members"]:
            if "members" in member and member["members"]:
                _sub_result = calculate_classes(
                    rec=member,
                    ref=member.get("taxid") or member.get("key") or ref,
                    pathogens=pathogens,
                    sample_type=sample_type,
                    taxdump=taxdump,
                )
                # Propagate classification fields back onto the member
                for _field in ('high_cons', 'is_pathogen', 'microbial_category',
                               'annClass', 'is_annotated', 'status', 'ref',
                               'matched_taxid', 'matched_rank',
                               'normalized_sample_site', 'commensal_sites',
                               'pathogenic_sites', 'members'):
                    if _field in _sub_result:
                        member[_field] = _sub_result[_field]

        member_categories = []
        member_high_cons = []
        member_annotations = []
        member_statuses = []
        member_matched_info = []  # Store (matched_taxid, matched_rank) for each member

        for member in rec["members"]:
            member_taxid = (
                member.get("taxid")
                or member.get("key")
                or ref
            )

            # If this member was already classified recursively, use those results
            if member.get('is_annotated') is not None and member.get('microbial_category') is not None:
                st_cat = member.get('microbial_category', 'Unknown')
                st_hc = member.get('high_cons', False)
                st_ann = member.get('annClass', 'Direct')
                is_annotated = member.get('is_annotated', 'No')
                status = member.get('status', '')
                member_matched_info.append((member.get('matched_taxid'), member.get('matched_rank')))
                member_categories.append(st_cat)
                member_high_cons.append(st_hc)
                member_annotations.append((st_ann, is_annotated))
                member_statuses.append(status)
                continue

            # Try to find pathogen info in lineage FOR EACH MEMBER
            # This searches: member_taxid → genus → family → ... until match found
            refpath, matched_taxid, matched_rank = find_pathogen_in_lineage(member_taxid)

            if refpath:
                cat, direct = get_pathogen_classification(refpath, normalized_sample_type)
                st_cat = normalize_category(cat)
                st_ann = "Direct" if direct else "Derived"
                st_hc = bool(refpath.get('high_cons', False) or refpath.get("high_consequence", False))
                is_annotated = "Yes"
                status = refpath.get("status", "N/A")
                member_matched_info.append((matched_taxid, matched_rank))
                # Propagate body-site-specific flora info onto the member
                member['commensal_sites'] = refpath.get('commensal_sites', [])
                member['pathogenic_sites'] = refpath.get('pathogenic_sites', [])
            else:
                st_cat = "Unknown"
                st_ann = "Direct"
                st_hc = False
                is_annotated = "No"
                status = ""
                member_matched_info.append((None, None))
                member['commensal_sites'] = []
                member['pathogenic_sites'] = []

            member_categories.append(st_cat)
            member_high_cons.append(st_hc)
            member_annotations.append((st_ann, is_annotated))
            member_statuses.append(status)

        # Determine aggregate values based on priority
        # High cons: True if ANY member is high_cons
        st_hc = any(member_high_cons)

        # Category: highest priority among all members
        if member_categories:
            st_cat = max(member_categories, key=lambda cat: CATEGORY_PRIORITY.get(cat, 0))
        else:
            st_cat = "Unknown"

        # AnnClass: Mixed if multiple categories, else Direct/Derived from members
        unique_categories = set(member_categories)
        if len(unique_categories) > 1:
            st_ann = "Mixed"
        else:
            # Use the annotation class from the first annotated member
            st_ann = member_annotations[0][0] if member_annotations else "Direct"

        # Is annotated: Yes if ANY member is annotated
        is_annotated = "Yes" if any(ann[1] == "Yes" for ann in member_annotations) else "No"

        # Status: collect unique non-empty statuses
        unique_statuses = [s for s in set(member_statuses) if s and s != "N/A"]
        status = ", ".join(unique_statuses) if unique_statuses else "N/A"

        # For aggregate, use the first matched taxid/rank info
        matched_taxid = None
        matched_rank = None
        for m_taxid, m_rank in member_matched_info:
            if m_taxid:
                matched_taxid = m_taxid
                matched_rank = m_rank
                break

        # Union of all member commensal/pathogenic sites for group-level display
        _agg_commensal = set()
        _agg_pathogenic = set()
        for member in rec["members"]:
            for s in member.get('commensal_sites', []):
                if isinstance(s, list):
                    _agg_commensal.update(s)
                elif s:
                    _agg_commensal.add(s)
            for s in member.get('pathogenic_sites', []):
                if isinstance(s, list):
                    _agg_pathogenic.update(s)
                elif s:
                    _agg_pathogenic.add(s)
        agg_commensal_sites = sorted(_agg_commensal)
        agg_pathogenic_sites = sorted(_agg_pathogenic)

    else:
        # No members, use the aggregate taxid directly
        # Search: taxid → genus → family → ... until match found
        refpath, matched_taxid, matched_rank = find_pathogen_in_lineage(taxid)

        if refpath:
            cat, direct = get_pathogen_classification(refpath, normalized_sample_type)
            st_cat = normalize_category(cat)
            st_ann = "Direct" if direct else "Derived"
            st_hc = bool(refpath.get('high_cons', False) or refpath.get("high_consequence", False))
            is_annotated = "Yes"
            status = refpath.get("status", "N/A")
            agg_commensal_sites = refpath.get("commensal_sites", [])
            agg_pathogenic_sites = refpath.get("pathogenic_sites", [])
        else:
            st_cat = "Unknown"
            st_ann = "Direct"
            st_hc = False
            is_annotated = "No"
            status = ""
            matched_taxid = None
            matched_rank = None
            agg_commensal_sites = []
            agg_pathogenic_sites = []

    # make a shallow copy so we don't mutate the input record
    new_item = dict(rec)

    # add ONLY the requested annotation fields
    new_item.update({
        "high_cons": st_hc,
        "is_pathogen": st_cat,
        "microbial_category": st_cat,
        "annClass": st_ann,
        "is_annotated": is_annotated,
        "status": status,
        "ref": str(taxid),
        "matched_taxid": matched_taxid,  # Which taxid in lineage matched
        "matched_rank": matched_rank,    # At what rank the match occurred
        "normalized_sample_site": normalized_sample_type,
        "commensal_sites": agg_commensal_sites,
        "pathogenic_sites": agg_pathogenic_sites,
    })

    return new_item
def annotate_aggregate_dict(
    aggregate_dict,
    pathogens,
    sample_type="Unknown",
    taxdump = {}
):

    out = []

    for ref, rec in aggregate_dict.items():
        # choose identifier for pathogen lookup
        taxid = (
            rec.get("toplevelkey")
            or rec.get("key")
            or rec.get("taxid")
            or ref
        )
        result = calculate_classes(
            rec=rec,
            ref=taxid,
            pathogens=pathogens,
            sample_type=sample_type,
            taxdump = taxdump
        )
        out.append(result)

    return out



def normalize_category(label):
    """
    Normalize category labels to standard format.

    Args:
        label (str): Category label (e.g., "primary", "opportunistic")

    Returns:
        str: Capitalized standard category
    """
    if not label:
        return "Unknown"

    s = str(label).strip().lower()

    if s in ("primary", "primary (sterile)"):
        return "Primary"
    if s in ("opportunistic",):
        return "Opportunistic"
    if s in ("potential",):
        return "Potential"
    if s in ("commensal",):
        return "Commensal"

    return "Unknown"



def choose_aggregate_category_from_strains(strains):
    """
    strains: list of dicts that already contain:
      - microbial_category
      - tass_score
    Rule:
      - If only one category exists => annClass = Direct; choose it
      - If multiple:
          pick category with highest max tass_score;
          annClass = Mixed
    Returns: (category, annClass)
    """
    cats = []
    for s in strains:
        cat = normalize_category(s.get("microbial_category") or s.get("is_pathogen") or "Unknown")
        try:
            tass = float(s.get("tass_score", 0) or 0)
        except Exception:
            tass = 0.0
        cats.append((cat, tass))

    if not cats:
        return "Unknown", "Direct"

    unique_cats = sorted({c for c, _ in cats})
    if len(unique_cats) == 1:
        return unique_cats[0], "Direct"

    # Multiple categories: choose by highest max tass, break ties by severity priority
    best_by_cat = {}
    for cat, tass in cats:
        best_by_cat[cat] = max(best_by_cat.get(cat, 0.0), tass)

    # sort: highest max tass, then highest priority
    ranked = sorted(
        best_by_cat.items(),
        key=lambda kv: (kv[1], CATEGORY_PRIORITY.get(kv[0], 0)),
        reverse=True
    )
    chosen = ranked[0][0]
    return chosen, "Mixed"

def pathogen_label(rft, sample_type):
    is_pathogen = "Unknown"
    isPathi = False
    direct_match = False
    callclass = rft.get('callclass', "N/A")
    high_cons = rft.get('high_cons', False)
    pathogenic_sites = rft.get('pathogenic_sites', [])
    commensal_sites = rft.get('commensal_sites', [])
    if sample_type == "sterile":
        direct_match = True
        is_pathogen = callclass.capitalize() if callclass else "Primary (sterile)"
        isPathi = True
    elif sample_type in pathogenic_sites:
        if callclass != "commensal":
            is_pathogen = callclass.capitalize()
            isPathi = True
        else:
            is_pathogen = "Potential"
        direct_match = True
    elif sample_type in commensal_sites:
        is_pathogen = "Commensal"
        direct_match = True
    elif callclass and callclass != "":
        is_pathogen = callclass.capitalize() if callclass else "Unknown"
        isPathi = True
    return is_pathogen, isPathi, direct_match, high_cons


# ── Control comparison utilities ─────────────────────────────────────────────

def load_control_data(json_paths):
    """Load one or more match_paths.py output JSONs and index organisms by key levels.

    Each JSON is expected to have the standard ``{metadata, organisms}``
    structure *or* the legacy plain-list format.

    Returns a dict::

        {
            "by_toplevelkey": { "<tlk>": [ {tass_score, numreads, source}, ... ] },
            "by_key":         { "<key>": [ ... ] },
            "by_subkey":      { "<sk>":  [ ... ] },
        }

    where each list entry represents one control sample's value for that
    organism (allowing multiple control replicates).
    """
    import json as _json_mod
    import os

    index = {
        "by_toplevelkey": defaultdict(list),
        "by_key": defaultdict(list),
        "by_subkey": defaultdict(list),
    }

    if not json_paths:
        return index

    for fpath in json_paths:
        if not os.path.exists(fpath):
            print(f"WARNING: control file not found, skipping: {fpath}")
            continue
        with open(fpath, "r") as fh:
            raw = _json_mod.load(fh)

        # Accept both structured and legacy list format
        if isinstance(raw, dict) and "organisms" in raw:
            organisms = raw["organisms"]
        elif isinstance(raw, list):
            organisms = raw
        else:
            print(f"WARNING: unexpected control JSON format in {fpath}, skipping")
            continue

        source = os.path.basename(fpath)

        for grp in organisms:
            tlk = str(grp.get("toplevelkey", grp.get("key", "")))
            grp_name = grp.get("name") or grp.get("toplevelname") or tlk
            grp_entry = {
                "tass_score": float(grp.get("tass_score", 0) or 0),
                "numreads": float(grp.get("numreads", 0) or 0),
                "source": source,
                "name": grp_name,
            }
            if tlk:
                index["by_toplevelkey"][tlk].append(grp_entry)

            # Index individual members at key and subkey levels.
            # Supports both the 3-level hierarchy (members are subkey groups
            # with their own 'members' of strains) and the legacy flat format.
            for member in grp.get("members", []):
                m_key = str(member.get("key", ""))
                m_subkey = str(member.get("subkey", member.get("key", "")))
                m_name = member.get("name") or member.get("toplevelname") or grp_name
                m_entry = {
                    "tass_score": float(member.get("tass_score", 0) or 0),
                    "numreads": float(member.get("numreads", 0) or 0),
                    "source": source,
                    "name": m_name,
                }
                if m_subkey:
                    index["by_subkey"][m_subkey].append(m_entry)
                # If this member has nested strain-level members, index those at key level
                if "members" in member and member["members"]:
                    for strain in member["members"]:
                        s_key = str(strain.get("key", ""))
                        s_subkey = str(strain.get("subkey", strain.get("key", "")))
                        s_name = strain.get("name") or m_name
                        s_entry = {
                            "tass_score": float(strain.get("tass_score", 0) or 0),
                            "numreads": float(strain.get("numreads", 0) or 0),
                            "source": source,
                            "name": s_name,
                        }
                        if s_key:
                            index["by_key"][s_key].append(s_entry)
                        if s_subkey and s_subkey not in index["by_subkey"]:
                            index["by_subkey"][s_subkey].append(s_entry)
                else:
                    # Legacy flat structure: member IS a strain
                    if m_key:
                        index["by_key"][m_key].append(m_entry)

    return index


def compute_control_metrics(
    sample_tass,
    sample_reads,
    neg_values,
    pos_values,
    fold_threshold=2.0,
):
    """Compute fold-change metrics comparing a sample organism to controls.

    Parameters
    ----------
    sample_tass : float
        TASS score of the sample organism.
    sample_reads : float
        Read count of the sample organism.
    neg_values : list[dict]
        List of ``{tass_score, numreads, source}`` from negative controls.
    pos_values : list[dict]
        List of ``{tass_score, numreads, source}`` from positive controls.
    fold_threshold : float
        Fold-change below which the sample is considered "within_negative"
        control range.  Default 2.0.

    Returns
    -------
    dict with keys:
        neg_max_tass, neg_max_reads,
        pos_min_tass, pos_min_reads,
        tass_fold_over_neg, reads_fold_over_neg,
        control_flag ("within_negative" | "above_negative" | "no_neg_controls"),
        neg_control_values, pos_control_values
    """
    result = {
        "neg_max_tass": 0.0,
        "neg_max_reads": 0.0,
        "pos_min_tass": None,
        "pos_min_reads": None,
        "tass_fold_over_neg": None,
        "reads_fold_over_neg": None,
        "tass_fold_over_pos": None,
        "reads_fold_over_pos": None,
        "control_flag": "no_neg_controls",
        "neg_control_values": [],
        "pos_control_values": [],
    }

    # ── Negative controls ────────────────────────────────────────────────
    if neg_values:
        neg_tass_vals = [v["tass_score"] for v in neg_values]
        neg_reads_vals = [v["numreads"] for v in neg_values]
        result["neg_max_tass"] = max(neg_tass_vals)
        result["neg_max_reads"] = max(neg_reads_vals)
        result["neg_control_values"] = [
            {"tass_score": v["tass_score"], "numreads": v["numreads"],
             "source": v.get("source", "")}
            for v in neg_values
        ]

        # Fold-change: how many times higher is the sample vs the worst neg control?
        if result["neg_max_tass"] > 0:
            result["tass_fold_over_neg"] = round(
                float(sample_tass) / result["neg_max_tass"], 4)
        else:
            # Neg controls have 0 TASS → any sample signal is infinitely above
            result["tass_fold_over_neg"] = float("inf") if sample_tass > 0 else 1.0

        if result["neg_max_reads"] > 0:
            result["reads_fold_over_neg"] = round(
                float(sample_reads) / result["neg_max_reads"], 4)
        else:
            result["reads_fold_over_neg"] = float("inf") if sample_reads > 0 else 1.0

        # Flag
        tass_fold = result["tass_fold_over_neg"]
        if tass_fold is not None and tass_fold != float("inf") and tass_fold < fold_threshold:
            result["control_flag"] = "within_negative"
        else:
            result["control_flag"] = "above_negative"

    # ── Positive controls ────────────────────────────────────────────────
    if pos_values:
        pos_tass_vals = [v["tass_score"] for v in pos_values]
        pos_reads_vals = [v["numreads"] for v in pos_values]
        result["pos_min_tass"] = min(pos_tass_vals)
        result["pos_min_reads"] = min(pos_reads_vals)
        result["pos_control_values"] = [
            {"tass_score": v["tass_score"], "numreads": v["numreads"],
             "source": v.get("source", "")}
            for v in pos_values
        ]

        # Fold-change: sample vs the minimum (weakest) positive control
        pos_min_tass = result["pos_min_tass"]
        pos_min_reads = result["pos_min_reads"]
        if pos_min_tass and pos_min_tass > 0:
            result["tass_fold_over_pos"] = round(
                float(sample_tass) / pos_min_tass, 4)
        else:
            result["tass_fold_over_pos"] = float("inf") if sample_tass > 0 else 0.0

        if pos_min_reads and pos_min_reads > 0:
            result["reads_fold_over_pos"] = round(
                float(sample_reads) / pos_min_reads, 4)
        else:
            result["reads_fold_over_pos"] = float("inf") if sample_reads > 0 else 0.0

    return result


def compute_control_comparison(
    data,
    neg_index,
    pos_index,
    fold_threshold=2.0,
    level="toplevelkey",
):
    """Compute control comparison for a single organism dict at a given hierarchy level.

    Parameters
    ----------
    data : dict
        Organism or member dict containing 'toplevelkey', 'key', 'subkey',
        'tass_score', and 'numreads'.
    neg_index : dict
        Output from ``load_control_data()`` for negative controls.
    pos_index : dict
        Output from ``load_control_data()`` for positive controls.
    fold_threshold : float
        Passed through to ``compute_control_metrics()``.
    level : str
        One of "toplevelkey", "key", "subkey".

    Returns
    -------
    dict or None
        Control comparison dict, or None if no controls at all.
    """
    level_map = {
        "toplevelkey": "by_toplevelkey",
        "key": "by_key",
        "subkey": "by_subkey",
    }
    idx_key = level_map.get(level, "by_toplevelkey")

    # Determine the organism's identifier at this level
    org_id = str(data.get(level, data.get("key", "")))
    if not org_id:
        return None

    neg_vals = neg_index.get(idx_key, {}).get(org_id, [])
    pos_vals = pos_index.get(idx_key, {}).get(org_id, [])

    if not neg_vals and not pos_vals:
        return None

    sample_tass = float(data.get("tass_score", 0) or 0)
    sample_reads = float(data.get("numreads", 0) or 0)

    return compute_control_metrics(
        sample_tass=sample_tass,
        sample_reads=sample_reads,
        neg_values=neg_vals,
        pos_values=pos_vals,
        fold_threshold=fold_threshold,
    )


def find_missing_positive_controls(final_json, pos_index, levels=None):
    """Identify organisms present in the positive control(s) but absent from the sample.

    Parameters
    ----------
    final_json : list[dict]
        The sample's output organism list (each element is a group with
        ``toplevelkey`` and ``members``).
    pos_index : dict
        The positive-control index returned by :func:`load_control_data`.
    levels : list[str] or None
        Which hierarchy levels to check.  Any combination of
        ``"toplevelkey"``, ``"key"``, ``"subkey"``.  Defaults to
        ``["toplevelkey"]`` when *None*.

    Returns
    -------
    list[dict]
        One entry per missing organism, each containing:
        ``level``, ``id``, ``name``, ``pos_tass_score``, ``pos_numreads``,
        ``source``, and ``missing_control: True``.
    """
    if not pos_index:
        return []
    if levels is None:
        levels = ["toplevelkey"]

    # Collect IDs present in the sample at each level
    sample_ids = {
        "toplevelkey": set(),
        "key": set(),
        "subkey": set(),
    }
    for grp in final_json:
        tlk = str(grp.get("toplevelkey", grp.get("key", "")))
        if tlk:
            sample_ids["toplevelkey"].add(tlk)
        for sk_m in grp.get("members", []):
            # Subkey-level member
            sk = str(sk_m.get("subkey", sk_m.get("key", "")))
            if sk:
                sample_ids["subkey"].add(sk)
            # Strain-level members nested inside the subkey group
            for strain in sk_m.get("members", []):
                mk = str(strain.get("key", ""))
                msk = str(strain.get("subkey", strain.get("key", "")))
                if mk:
                    sample_ids["key"].add(mk)
                if msk:
                    sample_ids["subkey"].add(msk)

    missing = []
    _seen_ids = set()  # org IDs already emitted (dedup across levels)

    level_to_index_key = {
        "toplevelkey": "by_toplevelkey",
        "key": "by_key",
        "subkey": "by_subkey",
    }

    for lvl in levels:
        idx_key = level_to_index_key.get(lvl)
        if not idx_key:
            continue
        present = sample_ids.get(lvl, set())
        for org_id, entries in pos_index.get(idx_key, {}).items():
            if org_id in present:
                continue
            if org_id in _seen_ids:
                continue
            _seen_ids.add(org_id)
            best = max(entries, key=lambda e: e.get("tass_score", 0))
            missing.append({
                "level": lvl,
                "id": org_id,
                "name": best.get("name", org_id),
                "pos_tass_score": best.get("tass_score", 0),
                "pos_numreads": best.get("numreads", 0),
                "source": best.get("source", ""),
                "missing_control": True,
            })

    return missing

