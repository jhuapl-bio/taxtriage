import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from io import BytesIO
import math
def body_site_map(x):
    body_sites = {
        'gut': "stool",
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
        # stats_dict_new[key]['std'] = custom_stddev(stats_dict_new[key])
        # print(f"New std dev: {stats_dict_new[key]['std']} ")
    return stats_dict_new, site_counts

def make_vplot(taxid, stats, column, result_df, percentile_column =None):
    # Calculate the standard deviation and mean
    std_dev = stats['std']
    mean = stats['mean']
    abundances = stats['abundances']
    # Create figure in memory, not on disk
    fig, ax = plt.subplots(figsize=(4, 2))

    # Continue with your data filtering and plotting logic
    filtered_data = [x for x in abundances if mean - 3 * std_dev < x < mean + 3 * std_dev]

    if taxid in result_df[column].values:
        specific_row = result_df[result_df[column] == taxid].iloc[0]
        specific_abundance = specific_row['abundance']
        if not percentile_column:
            if len(filtered_data) == 0:
                specific_percentile = 100
            else:
                specific_percentile = int(np.sum(np.array(filtered_data) <= specific_abundance) / len(filtered_data) * 100)
        else:
            specific_percentile = math.floor(specific_row[percentile_column]*100)
        # Boxplot
        sns.boxplot(x=filtered_data, showfliers=False, ax=ax,)

        # Dot color logic remains the same
        color = 'r'
        if std_dev != 0:
            try:
                if abs((specific_abundance - mean) / std_dev) <= 1:
                    color = 'g'
                elif abs((specific_abundance - mean) / std_dev) <= 2:
                    color = 'y'
                elif abs((specific_abundance - mean) / std_dev) <= 3:
                    color = 'orange'
            except Exception as ex :
                print(ex)

        plt.plot(specific_abundance, 0, color, marker="o", markersize=15, markeredgecolor="black", linewidth=3)

        # Set text position y coordinate slightly above or below the dot based on your plot's range or data
        text_y_position = 0.27 * ax.get_ylim()[1]  # Example: 10% of the way up the y-axis range
        if len(abundances) == 0:
            plt.text(specific_abundance, -2.4*text_y_position, f'Not Prev. in HHS', ha='center', va='bottom', color='black',
                bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2'), fontsize=22)
        plt.text(specific_abundance, text_y_position, f'{specific_percentile}%', ha='center', va='bottom', color='black',
                bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2'), fontsize=22)
    # plt.title("Percentile of " + str(name) + " to microbiome of " + body_site)
    # plt.xlabel("Abundance")

    # Instead of saving, return the figure as a BytesIO object
    fig_buffer = BytesIO()
    plt.savefig(fig_buffer, format='png')
    plt.close(fig)  # Make sure to close the plot
    fig_buffer.seek(0)  # Move to the beginning of the BytesIO buffer for reading

    return fig_buffer
