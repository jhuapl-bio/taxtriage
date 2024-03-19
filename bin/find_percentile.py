#!/usr/bin/env python
##############################################################################################
# Copyright 2022 The Johns Hopkins University Applied Physics Laboratory LLC
# All rights reserved.
# Permission is hereby granted, free of charge, to any person obtaining a copy of this
# software and associated documentation files (the "Software"), to deal in the Software
# without restriction, including without limitation the rights to use, copy, modify,
# merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
# PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
# OR OTHER DEALINGS IN THE SOFTWARE.

import pandas as pd
import argparse
import matplotlib.pyplot as plt
import sys
import seaborn as sns
import numpy as np
from scipy.stats import entropy


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Calculate the percentile rank for given abundance values based on type.",
        epilog="Example: script.py -i parameters.tsv -t name -d data.pkl",
    )
    parser.add_argument("-i", "--input", metavar="INPUT", required=True,
                        help="Path to input TSV file containing ID, body site, and abundance for which to calculate percentiles.")
    parser.add_argument("-d", "--data", metavar="DATA", required=True,
                        help="Path to input .pkl file containing all abundance data.")
    parser.add_argument("-t", "--type", metavar="TYPE", required=False, default='tax_id',
                        help="What type of data is being processed. Options: 'tax_id' or 'name'.",
                        choices=['tax_id', 'name'])
    parser.add_argument("-a", "--abundance_col", metavar="ABU", required=False, default='abundance',
                        help="Name of abundance column, default is abundance")
    parser.add_argument("-x", "--id_col", metavar="IDCOL", required=False, default='id',
                        help="Name of id column, default is id")
    parser.add_argument("-s", "--sitecol", metavar="SCOL", required=False, default='body_site',
                        help="Name of site column, default is body_site")
    return parser.parse_args(argv)

mapRank = {
    "S1": "subspecies",
    "S2": "strain",
}


def make_vplot(taxid, body_site, stats, args, result_df):

    # Calculate the standard deviation and mean
    std_dev = stats['std']
    mean = stats['mean']
    name = taxid
    filtered_data = stats['abundances']

    # filter the filtered_data to remove tobe above or below -3 and 3 std dev
    filtered_data = [x for x in filtered_data if (x > (mean - 3 * std_dev)) and (x < (mean + 3 * std_dev))]
    # Check if tax_id 626931 is present in result_df and plot a red dot and text if it is
    if taxid in result_df[args.type].values:
        specific_row = result_df[result_df[args.type] == taxid].iloc[0]
        specific_abundance = specific_row['abundance']
        # specific_percentile = round(specific_row['percentile'])
        # calcualte the percentile from stddev and mean
        # if the filtered_Data is empty then set the percentile to 100
        if len(filtered_data) == 0:
            specific_percentile = 100
        else:
            specific_percentile = int(np.sum(filtered_data <= specific_abundance) / len(filtered_data) * 100)
        # Calculate the difference between the specific abundance and the mean
        difference = specific_abundance - mean
        # Determine the number of standard deviations away from the mean
        num_std_devs = abs(difference) / std_dev

        # Calculate alpha based on the distance from the 50th percentile, ensuring it's within the range [0, 1]
        alpha = abs(specific_percentile - 50) / 50

        # Ensure that alpha is at least slightly visible and not greater than 1
        alpha = max(0.1, min(alpha, 1))
        plt.figure(figsize=(5, 2))

        # Plot the boxplot without showing outliers
        boxplot=sns.boxplot(x=filtered_data, showfliers=False,  )
        # Calculate alpha based on the distance from the 50th percentile, ensuring it's within the range [0, 1]
        alpha = abs(specific_percentile - 50) / 50
        alpha = max(0.1, min(alpha, 1))  # Ensure alpha is in a sensible range

        # Get the boxes (patches) from the boxplot and set their colors
        # Iterate over the boxes in the plot and set their facecolor and alpha
        # Get the patch objects (boxes) and modify their properties
        for patch in boxplot.patches:
            patch.set_facecolor("#666792")
            patch.set_alpha(alpha)


        # Plotting the red dot
        # if equal to or less than 1st std deviation color gree, if equal to or less than 2nd std deviation color orange, if equal to or less than 3rd std deviation color red
        color = "r"  # Default color green for within 1 standard deviation

        # Determine the color based on the number of standard deviations
        if num_std_devs <= 1:
            color = "g"  # Green for within ±1 standard deviation
        elif num_std_devs <= 2:
            color = "y"  # Yellow for within ±2 standard deviations
        elif num_std_devs <= 3:
            color = "orange"  # Orange for within ±3 standard deviations


        plt.plot(specific_abundance, 0, color, marker="o", markersize=6, markeredgecolor="black", linewidth=3 )  # 'ro' plots a red dot

        # Adding text annotation above the red dot
        plt.text(specific_abundance, 0.19, f'{specific_percentile}%',
                ha='center', va='bottom', color='black',
                bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'))

    plt.title("Percentile of "+str(name) + " to microbiome of "+body_site)
    plt.xlabel("Abundance")
    plt.show()


def main():
    args = parse_args()


    # Filter full_df to contain only the rows for the relevant body site(s) and calculate mean abundance
    abundance_data = pd.read_csv(args.input, sep="\t")

    # change column "id" in avbundance_data to "tax_id" if args.type is "name"
    abundance_data = abundance_data.rename(columns={args.id_col: args.type})
    abundance_data = abundance_data.rename(columns={args.sitecol: 'body_site'})
    abundance_data = abundance_data.rename(columns={args.abundance_col: 'abundance'})
    abundance_data = abundance_data[[args.type, 'body_site', 'abundance']]

    abundance_data = abundance_data.dropna(subset=[args.type])
    # Map several names for common groups for body_site
    body_site_map = {
        "gut": "stool",
        "nose": "nasal",
        "vagina": "vaginal",
        "teeth": "oral"
    }
    # convert all body_site with map
    abundance_data['body_site'] = abundance_data['body_site'].map(lambda x: body_site_map[x] if x in body_site_map else x)
    body_sites = abundance_data['body_site'].unique()
    # convert any empty or NaN body site to "unknown"


    # Load the full dataset for comparison
    full_df = pd.read_csv(args.data, sep='\t', nrows=1000000 )
    full_df['tax_id'] = full_df['tax_id'].astype(int)
    filtered_data = full_df[full_df['body_site'].isin(body_sites)]
    # print unique values of name
    stats_dict = {}
    for (taxid, body_site), group_data in filtered_data.groupby([args.type, 'body_site']):
        # get the first value of rank
        rank = group_data['rank'].iloc[0]
        stats_dict[(taxid, body_site)] = {
            'mean': group_data['abundance'].mean(),
            'std': group_data['abundance'].std(),
            "rank": rank,
            "abundances": list(group_data['abundance'])
        }
    for index, row in abundance_data.iterrows():
        # taxid, body_site, stats, args, result_df
        # if taxid and body site not in stats dict then make it empty or 0
        if (row[args.type], row['body_site']) not in stats_dict:
            stats = {
                'mean': 0,
                'std': 0,
                "rank": "NA",
                "abundances": []
            }
        else:
            stats = stats_dict[(row[args.type], row['body_site'])]
        rank = stats['rank']
        if rank in ['species', 'subspecies', 'strain']:
            make_vplot(
                row[args.type],
                row['body_site'],
                stats,
                args,
                abundance_data

            )
            break




if __name__ == "__main__":
    sys.exit(main())
