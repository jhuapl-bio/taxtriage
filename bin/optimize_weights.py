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

        # Append the new point (x2, y2)
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
        "minhash_reduction", "gini_coefficient", "diamond_identity", "rpkm", "rpm",
        "mapq_score",
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
            agg["subkeyname"] = entries[0].get("subkeyname", agg["name"])

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

        # ========== RECALCULATE breadth_log_score from aggregated coverage ==========
        # This ensures plasmids don't skew the score
        # agg["breadth_log_score"] = agg["coverage"] ** 2  # preferred formula
        agg['breadth_log_score'] = breadth_score_sigmoid(agg["coverage"])
        # Examples:
        # agg["breadth_log_score"] = agg["coverage"] ** 3  # more aggressive
        # agg["breadth_log_score"] = breadth_score_sigmoid(agg["coverage"])  # threshold-based

        out[gval] = agg

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
    total_reads = 0
):
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

    # Breadth
    coverage = data.get('coverage', 0)
    data['breadth_log_score'] = breadth_score_sigmoid(coverage)
    if not comparison_df.empty:
        accession = str(data['accession']).strip()
        if accession in comparison_df.index:
            c1 = float(comparison_df.loc[accession, col_stat])
            d_all = float(comparison_df.loc[accession, col_stat2])

            k_sig = 0.90
            x0 = -10.0
            pen = 1.0 / (1.0 + math.exp(-k_sig * (d_all - x0)))

            comparison_value = min(1.0, c1 * pen)
            data['minhash_score'] = comparison_value
        else:
            data['minhash_score'] = rpm_weight * 0.5  # default to a moderate score if no comparison available
    data['minhash_reduction'] = data.get('minhash_score', 1)



    return data

def rpm_confidence_weight(read_fraction, k=50000, midpoint=0.0001):
    """
    read_fraction = reads_mapped / total_reads  (value between 0 and 1)
    midpoint: fraction at which confidence = 0.5
              0.0001 = 0.01% of total reads (tune to your typical noise floor)
    k: steepness — needs to be large since fractions are tiny
    """
    return 1.0 / (1.0 + math.exp(-k * (read_fraction - midpoint)))


def breadth_score_sigmoid(coverage, midpoint=0.09, steepness=20):
        return 1.0 / (1.0 + math.exp(-steepness * (coverage - midpoint)))
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

):
    zscore, hmp_percentile = calculate_hmp_percentile(
        value = data,
        dists = hmp_dists,
        body_sites = body_sites,
        sampletype= sampletype,
    )
    data['mmbert'], data['mmbert_model'] = calculate_mmbert_prob(
        mmbert_dict = mmbert_dict,
        taxid = data.get('key', None)
    )
    data['zscore'] = zscore
    data['hmp_percentile'] = hmp_percentile
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

def compute_tass_score_from_metrics(metrics_df, breadth_w, minhash_w, gini_w, alpha=1.0):

    b = metrics_df["breadth_log_score"].to_numpy(float)
    m = metrics_df["minhash_reduction"].to_numpy(float)
    g = metrics_df["gini_coefficient"].to_numpy(float)

    score = breadth_w * b + minhash_w * m + gini_w * g
    score = alpha * score
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
    # normalize score of count to 0-1, min max range is 3, if larger than 3 set to 3 first
    # convert z score to percentile

    # Summation of each sub-score * weight

    tass_score = sum([
        apply_weight(data.get('disparity', 0), weights.get('disparity_weight', 0)),
        apply_weight(data.get('minhash_reduction', 0),       weights.get('minhash_weight', 0)),
        apply_weight(data.get('gini_coefficient', 0),             weights.get('gini_weight', 0)),
        apply_weight(data.get('breadth_log_score', 0),             weights.get('breadth_weight', 0)),
        apply_weight(data.get('hmp_percentile', 0),             weights.get('hmp_weight', 0)),
        apply_weight(data.get('mapq_score', 0),       weights.get('mapq_score', 0)),
        apply_weight(data.get('k2_disparity_score', 0),         weights.get('k2_disparity_score_weight', 0)),
        apply_weight(data.get('diamond', {}).get('identity', 0),
                     weights.get('diamond_identity', 0))  ])
    return tass_score

def calculate_hmp_percentile(
        value = {},
        dists = {},
        sampletype = "Sterile",
        body_sites = [],
        total_reads = 0
):
    abus = []
    taxid = value.get('taxid', None)
    for body_site in body_sites:
        if not taxid:
            continue
        try:
            # get the key where it is the (taxid, body_site)
            if (int(taxid), body_site) in dists:
                abus.append(
                    dict(
                        norm_abundance = dists[(int(taxid), body_site)].get('norm_abundance', 0),
                        std = dists[(int(taxid), body_site)].get('std', 0),
                        mean = dists[(int(taxid), body_site)].get('mean', 0),
                    )
                )
        except Exception as e:
            print(f"Error in taxid lookup for hmp: {e}")
    percent_total_reads_observed = 100*(sum(value.get('numreads', [])) / total_reads) if total_reads > 0 else 0
    sum_abus_expected = sum([x.get('mean',0)/100 for x in abus])
    stdsum = sum([x.get('std',0) for x in abus])
    zscore = ((percent_total_reads_observed -sum_abus_expected)/ stdsum)  if stdsum > 0 else 3
    # zscore = ((percent_total_reads_observed -sum_norm_abu)/ stdsum)  if stdsum > 0 else 3
    percentile = norm.cdf(zscore)

    return zscore, percentile

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
        dict: A new dictionary with the record data plus classification fields:
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
    if "members" in rec and rec["members"]:
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

            # Try to find pathogen info in lineage FOR EACH MEMBER
            # This searches: member_taxid → genus → family → ... until match found
            refpath, matched_taxid, matched_rank = find_pathogen_in_lineage(member_taxid)

            if refpath:
                # Use the new context-aware classification

                cat, direct = get_pathogen_classification(refpath, normalized_sample_type)
                st_cat = normalize_category(cat)
                st_ann = "Direct" if direct else "Derived"
                st_hc = bool(refpath.get('high_cons', False) or refpath.get("high_consequence", False))
                is_annotated = "Yes"
                status = refpath.get("status", "N/A")
                member_matched_info.append((matched_taxid, matched_rank))
            else:
                st_cat = "Unknown"
                st_ann = "Direct"
                st_hc = False
                is_annotated = "No"
                status = ""
                member_matched_info.append((None, None))

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

    else:
        # No members, use the aggregate taxid directly
        # Search: taxid → genus → family → ... until match found
        refpath, matched_taxid, matched_rank = find_pathogen_in_lineage(taxid)

        if refpath:
            # Use the new context-aware classification
            cat, direct = get_pathogen_classification(refpath, normalized_sample_type)
            st_cat = normalize_category(cat)
            st_ann = "Direct" if direct else "Derived"
            st_hc = bool(refpath.get('high_cons', False) or refpath.get("high_consequence", False))
            is_annotated = "Yes"
            status = refpath.get("status", "N/A")
        else:
            st_cat = "Unknown"
            st_ann = "Direct"
            st_hc = False
            is_annotated = "No"
            status = ""
            matched_taxid = None
            matched_rank = None

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

