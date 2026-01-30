from collections import defaultdict
from utils import format_non_zero_decimals, logarithmic_weight, apply_weight

def compute_tass_score(count, weights):
    """
    count is a dictionary that might look like:
    {
      'normalized_disparity': <some_value>,
      'alignment_score': <some_value>,
      'meangini': <some_value>,
      ...
    }
    We apply the known formula for TASS Score using the provided weights.
    """
    # normalize score of count to 0-1, min max range is 3, if larger than 3 set to 3 first
    # convert z score to percentile

    # Summation of each sub-score * weight



    tass_score = sum([
        apply_weight(count.get('normalized_disparity', 0), weights.get('disparity_score', 0)),
        apply_weight(count.get('alignment_score', 0),       weights.get('mapq_score', 0)),
        apply_weight(count.get('meanminhash_reduction', 0),       weights.get('minhash_weight', 0)),
        apply_weight(count.get('meangini', 0),             weights.get('gini_coefficient', 0)),
        apply_weight(count.get('hmp_percentile', 0),             weights.get('hmp_weight', 0)),
        apply_weight(count.get('log_weight_breadth', 0),             weights.get('breadth_weight', 0)),
        apply_weight(count.get('k2_disparity', 0),         weights.get('k2_disparity_weight', 0)),
        apply_weight(count.get('diamond', {}).get('identity', 0),
                     weights.get('diamond_identity', 0))  ])

    return tass_score

def calculate_scores(
        aggregated_stats,
        pathogens,
        sample_name="No_Name",
        sample_type="Unknown",
        total_reads=0,
        aligned_total = 0,
        weights={}
    ):
    """
    Write reference hits and pathogen information to a TSV file.

    Args:
    reference_hits (dict): Dictionary with reference names as keys and counts as values.
    pathogens (dict): Dictionary with reference names as keys and dictionaries of additional attributes as values.
    output_file_path (str): Path to the output TSV file.
    """
    final_scores = []

    # Write the header row
    for ref, count in aggregated_stats.items():
        strainlist = count.get('strainslist', [])
        is_pathogen = "Unknown"
        callfamclass = ""
        annClass = "None"
        refpath = pathogens.get(ref)
        # check if the sample type is in the pathogenic sites
        direct_match = False
        high_cons = False
        def pathogen_label(rft):
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

        formatname = count.get('name', "N/A")

        if refpath:
            is_pathogen, isPathi, direct_match, high_cons = pathogen_label(refpath)
            is_annotated = "Yes"
            status = refpath.get('status', "N/A")
            formatname = refpath.get('name', formatname)
            if is_pathogen == "Commensal":
                callfamclass = "Commensal Listing"
        else:
            is_annotated = "No"
            taxid = count.get(ref, "")
            status = ""
        if direct_match:
            annClass = "Direct"
        else:
            annClass = "Derived"
            # take the species_taxid and see if it is a pathogen
            if ref != count.get('species_taxid'):
                ref_spec = pathogens.get(count.get('species_taxid'), None)
                if ref_spec:
                    is_pathogen_spec, _, dir_spec, high_cons_spec = pathogen_label(ref_spec)
                    is_pathogen = is_pathogen_spec
                    high_cons = high_cons_spec
                    if dir_spec:
                        annClass = "Direct"
            # iterate through strainslist if it exists, and get pathogen label until dir_spec is true or all strains are checked
            elif not direct_match and strainlist:
                for strain in strainlist:
                    strain_taxid = strain.get('taxid', None)

                    if strain_taxid:
                        ref_strain = pathogens.get(strain_taxid, None)

                        if ref_strain:
                            is_pathogen_strain, _, dir_strain, high_cons_strain = pathogen_label(ref_strain)
                            if dir_strain:
                                is_pathogen = is_pathogen_strain
                                high_cons = high_cons_strain
                                annClass = "Direct"
        listpathogensstrains = []
        fullstrains = []
        callclasses = set()
        pathogenic_sites = refpath.get('pathogenic_sites', []) if refpath else []
        print(formatname, ref, annClass, sample_type)
        if strainlist:
            pathogenic_reads = 0
            merged_strains = defaultdict(dict)

            for x in strainlist:
                strainname = x.get('strainname', None)
                taxid = x.get('taxid', None)
                if taxid:
                    keyx = taxid
                elif not taxid and strainname:
                    keyx = strainname
                else:
                    keyx = None
                if keyx in merged_strains:
                    merged_strains[keyx]['numreads'] += x.get('numreads', 0)
                    merged_strains[keyx]['subkeys'].append(x.get('subkey', ""))
                else:
                    merged_strains[keyx] = x
                    merged_strains[keyx]['numreads'] = x.get('numreads', 0)
                    merged_strains[keyx]['subkeys'] = [x.get('subkey', "")]

            for xref, x in merged_strains.items():
                pathstrain = None
                if x.get('taxid'):
                    fullstrains.append(f"{x.get('strainname', 'N/A')} ({x.get('taxid', '')}: {x.get('numreads', 0)} reads)")
                else:
                    fullstrains.append(f"{x.get('strainname', 'N/A')} ({x.get('numreads', 0)} reads)")

                if x.get('taxid') in pathogens:
                    pathstrain = pathogens.get(x.get('taxid'))
                elif x.get('fullname') in pathogens:
                    pathstrain = pathogens.get(x.get('fullname'))
                if pathstrain:
                    taxx = x.get('taxid', "")
                    if pathstrain.get('callclass') not in ["commensal", "Unknown", 'unknown', '', None]:
                        callclasses.add(pathstrain.get('callclass').capitalize())

                    if pathstrain.get('high_cons', False):
                        high_cons = True
                    if sample_type in pathstrain.get('pathogenic_sites', []) or sample_type == pathstrain.get('general_classification', ''):
                        pathogenic_reads += x.get('numreads', 0)
                        percentreads = f"{x.get('numreads', 0)*100/aligned_total:.1f}" if aligned_total > 0 and x.get('numreads', 0) > 0 else "0"
                        listpathogensstrains.append(f"{x.get('strainname', 'N/A')} ({percentreads}%)")

            if callfamclass == "" or len(listpathogensstrains) > 0:
                callfamclass = f"{', '.join(listpathogensstrains)}" if listpathogensstrains else ""
        if len(callclasses) > 0:
            # if Primary set is_pathogen to primary, if opposite set to opportunistic if potential set to potential
            if "Primary" in callclasses:
                is_pathogen = "Primary"
            elif "Opportunistic" in callclasses:
                is_pathogen = "Opportunistic"
            elif "Potential" in callclasses:
                is_pathogen = "Potential"
        breadth_total = count.get('breadth_total', 0)
        countreads = sum(count['numreads'])
        if aligned_total == 0:
            percent_aligned = 0
        else:
            percent_aligned = format_non_zero_decimals(100*countreads / aligned_total)
        if total_reads == 0:
            percent_total = 0
        else:
            percent_total = format_non_zero_decimals(100*countreads / total_reads)
        if len(pathogenic_sites) == 0:
            pathogenic_sites = ""

        # Apply weights to the relevant scores
        # Example usage:
        breadth = count.get('breadth_total')  # your raw value between 0 and 1
        # print the breadth and the name
        log_weight_breadth = 1-logarithmic_weight(breadth)
        # get log weight of the breadth_total
        count['log_weight_breadth'] = log_weight_breadth
        tass_score = compute_tass_score(
            count,
            weights,
        )
        final_scores.append(
            dict(
                tass_score=tass_score,
                formatname=formatname,
                sample_name=sample_name,
                sample_type=sample_type,
                fullstrains=fullstrains,
                listpathogensstrains=listpathogensstrains,
                percent_total=percent_total,
                percent_aligned=percent_aligned,
                callfamclass=callfamclass,
                log_breadth_weight = count.get('log_weight_breadth', 0),
                total_reads=sum(count['numreads']),
                reads_aligned=countreads,
                hmp_percentile = count.get('hmp_percentile', 0),
                zscore= count.get('zscore', 0),
                meancoverage=count.get('meancoverage', 0),
                breadth_total = breadth_total,
                total_length = count.get('total_length', 0),
                is_annotated=is_annotated,
                is_pathogen=is_pathogen,
                ref=ref,
                status=status,
                annClass=annClass,
                high_cons = high_cons,
                mmbert = count.get('mmbert', None),
                mmbert_model = count.get('mmbert_model', None),
                pathogenic_reads = pathogenic_reads,
                gini_coefficient=count.get('meangini', 0),
                meanbaseq=count.get('meanbaseq', 0),
                meanmapq=count.get('meanmapq', 0),
                alignment_score=count.get('alignment_score', 0),
                disparity_score=count.get('normalized_disparity', 0),
                covered_regions=count.get('covered_regions', 0),
                siblings_score=count.get('raw_disparity', 0),
                k2_reads=count.get('k2_numreads', 0),
                k2_disparity=count.get('k2_disparity', 0),
                gtcov=count.get('gtcov', 0),
                minhash_score=count.get('meanminhash_reduction', 0),

            )
        )
    return final_scores
def compute_cost(
    aggregated_stats,
    pathogens,
    sample_name,
    sample_type,
    total_reads,
    aligned_total,
    weights
):
    """
    Runs calculate_scores(...), then sums up penalties based on meancoverage & tass_score.
    """
    total_weight = sum(weights.values())
    if total_weight != 1:
        for key in weights:
            if total_weight != 0:
                weights[key] = weights[key] / total_weight
            else:
                weights[key] = 0
    # First, get final_scores from your existing pipeline:
    final_scores = calculate_scores(
        aggregated_stats=aggregated_stats,
        pathogens=pathogens,
        sample_name=sample_name,
        sample_type=sample_type,
        total_reads=total_reads,
        aligned_total=aligned_total,
        weights=weights,
    )

    cost = 0.0
    # print("__________")
    perrefcost = {}
    for rec in final_scores:
        ref = rec.get('ref')
        tass_score = rec['tass_score']
        actual_coverage = rec.get('meancoverage', 0.0)

        # The "ground-truth" coverage for this ref (already placed in aggregated_stats)
        # e.g. aggregated_stats[ref]['gtcov'] may be 0 or >0
        coverage = aggregated_stats.get(ref, {}).get('gtcov', 0.0)
        # CASE 1: coverage == 0 => we want TASS ~ 0,
        #         also penalize if actual_coverage is large (contradiction).
        if coverage <= 0:
            # The bigger the TASS (wrongly claiming positivity),
            # or the bigger the actual coverage (contradiction),
            # the larger this penalty.
            cost += (1.0 + actual_coverage) * (tass_score ** 2)
            perrefcost[ref] = (1.0 + actual_coverage) * (tass_score ** 2)
        else:
            # CASE 2: coverage > 0 => we want TASS near 1 for big coverage.
            #  (A) penalize if TASS is low => cost ~ coverage * (1 - TASS)^2
            #  (B) penalize if TASS is high but actual coverage is too small
            #      => e.g. coverage_short = (coverage - actual_coverage)
            #             cost2 = (tass_score^2) * coverage_short
            #      meaning if ground-truth coverage is large but the pipeline
            #      found less, it’s contradictory to have TASS=high.
            cost1 = coverage * ((1.0 - tass_score) ** 2)

            coverage_short = max(0.0, coverage - actual_coverage)
            cost2 = (tass_score ** 2) * coverage_short

            cost += (cost1 + cost2)
            perrefcost[ref] =  (cost1 + cost2)
        # print("\t",ref, actual_coverage, coverage, tass_score, rec['gini_coefficient'])
    return cost, perrefcost


def optimize_three_weights(
    aggregated_stats,
    pathogens,
    sample_name="No_Name",
    sample_type="Unknown",
    aligned_total=0,
    total_reads=0,
    initial_weights=None,
    max_iterations=500
):
    """
    Optimize only the three weights:
      - 'breadth_weight'
      - 'gini_coefficient'
      - 'minhash_weight'

    The function returns a dict mapping the three weight names to optimized values
    that sum to 1 (compute_cost normalizes weights if necessary).
    """

    # default initial weights (equal)
    if initial_weights is None:
        initial_weights = [1.0, 1.0, 1.0]

    # Ensure we have three elements
    if len(initial_weights) != 3:
        raise ValueError("initial_weights must be a list/tuple of length 3")

    # bounds: each weight in [0,1]
    bounds = [(0.0, 1.0), (0.0, 1.0), (0.0, 1.0)]

    # constraint: sum(weights) == 1
    constraints = {'type': 'eq', 'fun': lambda w: w[0] + w[1] + w[2] - 1.0}

    # objective: use compute_cost; we negate because minimize() minimizes the objective.
    # The previous optimize_weights used -cost to try to maximize the metric used internally.
    def objective(w):
        # Build a weights dictionary compatible with compute_cost / calculate_scores.
        # Put the three weights in the canonical names used by the pipeline.
        weight_dict = {
            'breadth_weight': float(w[0]),
            'gini_coefficient': float(w[1]),
            'minhash_weight': float(w[2]),
            # leave other weights absent (calculate_scores / compute_cost use .get(..., 0))
        }
        cost, _ = compute_cost(
            aggregated_stats=aggregated_stats,
            pathogens=pathogens,
            sample_name=sample_name,
            sample_type=sample_type,
            total_reads=total_reads,
            aligned_total=aligned_total,
            weights=weight_dict
        )
        # Negate to match earlier pattern (if you want to maximize cost); change sign if you want to minimize.
        return -cost

    # use SLSQP like your other optimizer
    from scipy.optimize import minimize
    result = minimize(
        objective,
        x0=initial_weights,
        method='SLSQP',
        bounds=bounds,
        constraints=constraints,
        options={'maxiter': max_iterations, 'disp': True}
    )

    if not result.success:
        raise RuntimeError(f"3-weight optimization failed: {result.message}")

    optimized = {
        'breadth_weight': float(result.x[0]),
        'gini_coefficient': float(result.x[1]),
        'minhash_weight': float(result.x[2]),
    }

    # compute final cost (optional) and print summary
    final_cost, perref = compute_cost(
        aggregated_stats=aggregated_stats,
        pathogens=pathogens,
        sample_name=sample_name,
        sample_type=sample_type,
        total_reads=total_reads,
        aligned_total=aligned_total,
        weights=optimized
    )

    print("Optimized 3 weights (breadth, gini, minhash):")
    print(f"\t breadth_weight  = {optimized['breadth_weight']:.4f}")
    print(f"\t gini_coefficient = {optimized['gini_coefficient']:.4f}")
    print(f"\t minhash_weight  = {optimized['minhash_weight']:.4f}")
    print(f"\t final cost = {final_cost:.6f}")

    return optimized
