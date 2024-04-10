import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sys
from io import BytesIO
import random

def import_distributions(
        distribution_data,
        column_id,
        body_sites = [],
    ):
    # convert any empty or NaN body site to "unknown"
    distribution_data  = str(distribution_data)
    # Load the full dataset for comparison
    #if distribution_data is a .gz file, decompress otherwise read as is
    stats_dict_new = {}
    if distribution_data.endswith('.gz'):
        stats_dict = pd.read_csv(distribution_data, sep='\t', compression='gzip').to_dict(orient='records')
    else:
        stats_dict = pd.read_csv(distribution_data, sep='\t').to_dict(orient='records')
    if len(body_sites) > 0:

        for value in stats_dict:
            if 'body_site' in value and value['body_site'] in body_sites:
                stats_dict_new[(value[column_id], value['body_site'])] = value
    else:
        stats_dict_new = {(x[column_id], x['body_site']): x for x in stats_dict}
    i = 0
    # convert all abundances to a list of floats
    for key in stats_dict_new:
        stats_dict_new[key]['abundances'] = [float(x) for x in stats_dict_new[key]['abundances'].split(",")]
    return stats_dict_new
def make_vplot(taxid, stats, column, result_df):
    # Calculate the standard deviation and mean
    std_dev = stats['std']
    mean = stats['mean']
    filtered_data = stats['abundances']

    # Create figure in memory, not on disk
    fig, ax = plt.subplots(figsize=(4, 2))

    # Continue with your data filtering and plotting logic
    filtered_data = [x for x in filtered_data if mean - 3 * std_dev < x < mean + 3 * std_dev]

    if taxid in result_df[column].values:
        specific_row = result_df[result_df[column] == taxid].iloc[0]
        specific_abundance = specific_row['abundance']

        if len(filtered_data) == 0:
            specific_percentile = 100
        else:
            specific_percentile = int(np.sum(np.array(filtered_data) <= specific_abundance) / len(filtered_data) * 100)

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
        if len(filtered_data) == 0:
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
