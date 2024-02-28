#!/usr/bin/env python3

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
#

import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os
import numpy as np

def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Plot coverage depth across chromosomes.",
        epilog="Example: python coverage_plot.py -i coverage.bedgraph -o plot.png",
    )
    parser.add_argument(
        "-i",
        "--input",
        metavar="INPUT",
        help="Path to coverage BEDGRAPH file",
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="OUTPUT",
        help="Output plot file (e.g., plot.png)",
    )
    return parser.parse_args(argv)


args = parse_args()

# Load the coverage data
coverage_data = pd.read_csv(args.input, sep='\t', header=None, names=['chrom', 'start', 'end', 'depth'])

# Normalize the coverage depth
max_depth = coverage_data['depth'].max()
coverage_data['normalized_depth'] = coverage_data['depth'] / max_depth

# Determine the number of chromosomes and create a subplot grid
chromosomes = coverage_data['chrom'].unique()
n_chromosomes = len(chromosomes)
n_cols = int(np.ceil(np.sqrt(n_chromosomes)))
n_rows = int(np.ceil(n_chromosomes / n_cols))

# Adjust the figure size to make plots less tall
fig_width = 15  # You can adjust this as needed
fig_height = fig_width / n_cols * n_rows  # Adjust the height based on the number of rows and columns

# Create a grid of sparkline plots
plt.figure(figsize=(fig_width, fig_height))


for i, chrom in enumerate(chromosomes, 1):
    ax = plt.subplot(n_rows, n_cols, i)
    data = coverage_data[coverage_data['chrom'] == chrom]
    format_chrom = chrom.split('|')[0]
    ax.plot(data['start'], data['normalized_depth'], color='blue', linewidth=1)
    ax.set_title(format_chrom, fontsize=8)
    ax.axis('off')  # Sparklines typically do not have axes

plt.tight_layout()
plt.savefig(args.output)
plt.close()
