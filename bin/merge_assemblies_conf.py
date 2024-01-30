#!/usr/bin/env python3
import sys
import csv
import argparse
from collections import defaultdict
import math

def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Merge a paf to confidence output mqc.tsv with the kraken report it corresponds to include tax info",
        epilog="Example: ",
    )
    parser.add_argument(
        "-i",
        "--conf_file",
        metavar="CONF",
        type=str,
        required=True,
        help="Confidences file",
    )
    parser.add_argument(
        "-o",
        "--file_out",
        metavar="FILE_OUT",
        type=str,
        required=True,
        help="Name of the output tsv file containing a mapping of top n organisms at individual taxa levels",
    )

    return parser.parse_args(argv)

def weighted_average(values, weights):
    weighted_sum = sum(value * weight for value, weight in zip(values, weights))
    total_weight = sum(weights)
    if total_weight == 0:
        return 0
    return weighted_sum / total_weight

def calculate_stdev(depths, ref_size, weighted_mean_depth):
    sum_squares = sum((depth - weighted_mean_depth)**2 * size for depth, size in zip(depths, ref_size))
    count = len(depths)
    variance = sum_squares / count if count > 0 else 0
    return math.sqrt(variance)

def coverage_thresholds(depths, thresholds, ref_size):
    coverage_counts = {threshold: 0 for threshold in thresholds}
    for depth, size in zip(depths, ref_size):
        for threshold in thresholds:
            if depth >= threshold:
                coverage_counts[threshold] += size
    return coverage_counts

def main(argv=None):
    args = parse_args(argv)

    # Initialize data storage
    data = defaultdict(lambda: {'Mean Depth': [], 'Average Coverage': [], 'Ref Size': [], 'Total Reads Aligned': 0, 'Individual Reads Aligned': []})

    # Read the confidences file
    with open(args.conf_file, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        total_reads_aligned = sum(int(row['Total Reads Aligned']) for row in reader)
        file.seek(0)  # Reset file pointer to the beginning
        reader = csv.DictReader(file, delimiter='\t')  # Reinitialize reader
        for row in reader:
            name = row['Name']
            data[name]['Mean Depth'].append(float(row['Mean Depth']))
            data[name]['Average Coverage'].append(float(row['Average Coverage']))
            data[name]['Ref Size'].append(float(row['Ref Size']))
            data[name]['Total Reads Aligned'] += int(row['Total Reads Aligned'])
            data[name]['Individual Reads Aligned'].append(int(row['Individual Reads Aligned']))
    file.close()
    # Process data
    aggregated = []
    thresholds = [1, 10, 50, 100, 300]
    for name, values in data.items():
        mean_depth = weighted_average(values['Mean Depth'], values['Ref Size'])
        avg_coverage = weighted_average(values['Average Coverage'], values['Ref Size'])
        stdev = calculate_stdev(values['Mean Depth'], values['Ref Size'], mean_depth)
        abu_total_aligned = sum(values['Individual Reads Aligned']) / total_reads_aligned
        coverage_counts = coverage_thresholds(values['Mean Depth'], thresholds, values['Ref Size'])
        coverage_percentages = [coverage_counts[t] / sum(values['Ref Size']) * 100 for t in thresholds]

        aggregated.append({
            'Name': name,
            'Weighted Mean Depth': mean_depth,
            'Weighted Avg Coverage': avg_coverage,
            'Stdev': stdev,
            '% Reads Aligned': abu_total_aligned * 100,
            '1X Cov.': coverage_percentages[0],
            '10X Cov.': coverage_percentages[1],
            '50X Cov.': coverage_percentages[2],
            '100X Cov.': coverage_percentages[3],
            '300X Cov.': coverage_percentages[4],
            'Total Ref Size': sum(values['Ref Size']),
            'Total Reads Aligned': values['Total Reads Aligned']
        })

    # Output the aggregated data
    print(aggregated)
    exit()
    output_file = args.file_out
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=aggregated[0].keys(), delimiter='\t')
        writer.write
