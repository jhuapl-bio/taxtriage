import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from io import BytesIO
import math
def body_site_map(x):
    body_sites = {
        "stool": "gut",
        "nose": "nasal",
        "vagina": "vaginal",
        "teeth": "oral",
        "resp": "nasal",
        "abscess": "skin",
        "absscess": "skin",
        "sputum": "oral",
        "mouth": "oral",
        "urinary tract": "urine",
        "ear": "skin",
        "lung": ["oral", "nasal"],
        "eye": "eye",
        "sinus": "nasal",
        "urogenital": "skin",
        "cornea": "eye",
    }
    if body_sites.get(x):
        return body_sites[x]
    else:
        return x




def import_distributions(
        distribution_data,
        column_id,
        body_sites = [],
    ):
    site_counts = dict()
    # convert any empty or NaN body site to "unknown"
    distribution_data  = str(distribution_data)
    # Load the full dataset for comparison
    #if distribution_data is a .gz file, decompress otherwise read as is
    stats_dict_new = {}
    saved_amts = []
    if distribution_data.endswith('.gz'):
        stats_dict = pd.read_csv(distribution_data, sep='\t', compression='gzip').to_dict(orient='records')
    else:
        stats_dict = pd.read_csv(distribution_data, sep='\t').to_dict(orient='records')
    for value in stats_dict:
        value['body_site'] = body_site_map(value['body_site'])
    uniquesites = set([x['body_site'] for x in stats_dict])
    if len(body_sites) > 0:
        for value in stats_dict:
            if value['body_site'] not in site_counts:
                site_counts[value['body_site']] = value['site_count']

            if value['rank'] == "species" or value['rank'] == "subspecies":
                vals ={
                    'site_counts': value['site_count'],
                    'body_site':value['body_site'],
                    'number_samples': (len(value['abundances'].split(","))),
                    'name': value['name'],
                    "rank": value['rank'],
                    "proportion total": (len(value['abundances'].split(","))) / value['site_count']
                }
                saved_amts.append(vals)

            if 'body_site' in value and value['body_site'] in body_sites:
                stats_dict_new[(value[column_id], value['body_site'])] = value
    else:
        stats_dict_new = {(x[column_id], body_site_map(x['body_site'])): x for x in stats_dict}
    i = 0
    # Sort the saved_amts list of dicts on the "number_samples" key
    saved_amts = sorted(saved_amts, key=lambda x: x['proportion total'], reverse=False)
    # for val in saved_amts:
    #     print(f"{val['name']}, {val['site_counts']} - {val['number_samples']} - {val['rank']} - {val['body_site']}, Found in {(100*val['proportion total']):.2f}% of samples")



    # convert all abundances to a list of floats
    def custom_stddev(row):
        total_count = row['site_count']
        abundances = row['abundances']
        norm_abundance = row['norm_abundance']
        summation = 0
        for x in abundances:
            summation += (x - norm_abundance)**2
        stdev =  (summation / total_count)**0.5
        return stdev
    for key in stats_dict_new:
        stats_dict_new[key]['abundances'] = [float(x) for x in stats_dict_new[key]['abundances'].split(",")]
        stats_dict_new[key]['norm_abundance'] = sum(stats_dict_new[key]['abundances'])/stats_dict_new[key]['site_count']
        # print(f"Original std dev: {stats_dict_new[key]['std']}, Original Mean: {stats_dict_new[key]['mean']}")
        stats_dict_new[key]['std'] = custom_stddev(stats_dict_new[key])
        # print(f"New std dev: {stats_dict_new[key]['std']} ")
    return stats_dict_new, site_counts


def make_vplot(taxid, stats, column, result_df, percentile_column=None, abu_col='abundances'):
    import statistics
    import numpy as np
    import math
    import matplotlib.pyplot as plt
    import seaborn as sns
    from io import BytesIO
    from scipy.stats import norm


    # Pre-calculations
    std_dev = stats['std']
    mean = stats.get('norm_abundance', 0)
    abundances = stats['abundances']
    mean_recalc = statistics.mean(abundances) if len(abundances) > 1 else 0
    std_dev_recalc = statistics.stdev(abundances) if len(abundances) > 1 else 0
    std_dev_recalc = std_dev
    # Create figure in memory
    fig, ax = plt.subplots(figsize=(4, 2))

    # if stats.get('name') == "Alistipes shahii":
    #     print(stats['name'], stats['mean'], stats['std'], stats['norm_abundance'])
    #     print(mean_recalc, max(abundances))
    #     print(f"75th percentile: {np.percentile(abundances, 75)}")
    #     print(f"25th percentile: {np.percentile(abundances, 25)}")
    #     print(f"50th percentile: {np.percentile(abundances, 50)}")
    #     print(f"93rd percentile: {np.percentile(abundances, 93)}")
    #     print(f"10th percentile: {np.percentile(abundances, 10)}")
    #     print(f"40th percentile: {np.percentile(abundances, 40)}")
    #     print(f"4th percentile: {np.percentile(abundances, 4)}")
    #     print(f"23rd percentile: {np.percentile(abundances, 23)}")
    #     print(f"3rd percentile: {np.percentile(abundances, 3)}")
    #     print(f"22nd percentile: {np.percentile(abundances, 22)}")

    filtered_data = abundances
    if taxid in result_df[column].values:
        specific_row = result_df[result_df[column] == taxid].iloc[0]
        specific_abundance = specific_row[abu_col]

        # Calculate the z-score based on raw abundances
        if std_dev_recalc:
            zscore = (specific_abundance - mean_recalc) / std_dev_recalc
            specific_percentile_new = 100 * (0.5 * (1 + math.erf(zscore / math.sqrt(2))))
        else:
            zscore = 0
            specific_percentile_new = 0
        # Compute empirical percentile as before
        if not percentile_column:
            if len(filtered_data) == 0:
                specific_percentile = 100
            else:
                specific_percentile = int(np.sum(np.array(filtered_data) <= specific_abundance) / len(filtered_data) * 100)
        else:
            specific_percentile = math.floor(specific_row[percentile_column] * 100)
            # if specific_row['name'] == "Alistipes shahii":
            #     print("\tspecific percentile (empirical):", specific_percentile,
            #         "\n\tspecific_abundance:", specific_abundance,
            #         '\n\tzscore:', zscore,
            #         '\n\tstd-dev:', std_dev_recalc,
            #         '\n\tstd-dev original:', std_dev,
            #         '\n\tmean-recalc:', mean_recalc,
            #         "\n\tspecific percentile (parametric):", specific_percentile_new)

        # Boxplot with extended whiskers (3 * IQR)
        sns.boxplot(x=filtered_data, showfliers=False, whis=3, ax=ax)

        # Dot color logic remains unchanged
        color = 'r'
        if std_dev_recalc:
            try:
                z = abs((specific_abundance - mean_recalc) / std_dev_recalc)
                if z <= 1:
                    color = 'g'
                elif z <= 2:
                    color = 'y'
                elif z <= 3:
                    color = 'orange'
            except Exception as ex:
                print(ex)
        plt.plot(specific_abundance, 0, color, marker="o", markersize=15,
                markeredgecolor="black", linewidth=3)

        # Position text above the dot
        text_y_position = 0.27 * ax.get_ylim()[1]
        if len(abundances) == 0:
            plt.text(specific_abundance, -2.4 * text_y_position, 'Not Prev. in HHS',
                    ha='center', va='bottom', color='black',
                    bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2'),
                    fontsize=22)
        plt.text(specific_abundance, text_y_position, f'{specific_percentile}%', ha='center', va='bottom', color='black',
                bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2'), fontsize=22)

    # Save the figure to a BytesIO object and return it
    fig_buffer = BytesIO()
    plt.savefig(fig_buffer, format='png')
    plt.close(fig)
    fig_buffer.seek(0)
    return fig_buffer
