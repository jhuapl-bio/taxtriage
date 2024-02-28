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
    coverage_counts = [0] * len(thresholds)
    print(depths, thresholds, ref_size)
    for depth, size in zip(depths, ref_size):
        for i, threshold in enumerate(thresholds):
            if depth >= threshold:
                coverage_counts[i] += size
    return coverage_counts

def main(argv=None):
    args = parse_args(argv)

    # Initialize data storage
    data = defaultdict(lambda: {'Mean Depth': [], 'Average Coverage': [], 'Ref Size': [], 'Total Reads Aligned': 0, '% Reads Aligned': []})

    # Read the confidences file
    with open(args.conf_file, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        total_reads_aligned = sum(int(row['Total Reads Aligned']) for row in reader)
        file.seek(0)  # Reset file pointer to the beginning
        reader = csv.DictReader(file, delimiter='\t')  # Reinitialize reader
        for row in reader:
            name = row['Organism']
            data[name]['Mean Depth'].append(float(row['Mean Depth']))
            data[name]['Average Coverage'].append(float(row['Average Coverage']))
            data[name]['Ref Size'].append(float(row['Ref Size']))
            data[name]['Total Reads Aligned'] += float(row['Total Reads Aligned'])
            # data[name]['Individual Reads Aligned'].append(float(row['Individual Reads Aligned']))
            # split 1:10:50:100:300X Cov. column into list on colons, setting thresholds and making it as a dict
            data[name]['Coverage Thresholds'] = {int(t): float(p) for t, p in zip([1, 10, 50, 100, 300], row['1:10:50:100:300X Cov.'].split(':'))}

    file.close()
    # Process data
    aggregated = []
    thresholds = [1, 10, 50, 100, 300]
    for name, values in data.items():
        mean_depth = weighted_average(values['Mean Depth'], values['Ref Size'])
        avg_coverage = weighted_average(values['Average Coverage'], values['Ref Size'])
        stdev = calculate_stdev(values['Mean Depth'], values['Ref Size'], mean_depth)
        if values['Total Reads Aligned'] == 0:
            abu_total_aligned = 0
        else:
            abu_total_aligned = values['Total Reads Aligned'] / total_reads_aligned
        # use the Coverage Thresholds column to find weighted threhsolds for merged row
        cov_thresholds = values['Coverage Thresholds']

        # coverage_counts =  coverage_thresholds(cov_thresholds, thresholds, values['Ref Size'])
        # print(coverage_counts)
        # exit()
        aggregated.append( {
            'Organism': name,
            'Weighted Mean Depth': mean_depth,
            'Weighted Avg Coverage': avg_coverage,
            'Total Ref Size': sum(values['Ref Size']),
            'Total Reads Aligned': values['Total Reads Aligned'],
            '% Reads Aligned': abu_total_aligned * 100,
            'Stdev': stdev,
            "1X:10x:50x:100x:300x Cov.": "Disabled"
            # '1X:10x:50x:100x:300x Cov.': f"{coverage_percentages[0]:.2f}:{coverage_percentages[1]:.2f}:{coverage_percentages[2]:.2f}:{coverage_percentages[3]:.2f}:{coverage_percentages[4]:.2f}",
        })

    # Output the aggregated data
    output_file = args.file_out
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=aggregated[0].keys(), delimiter='\t')
        writer.writeheader()
        writer.writerows(aggregated)
    f.close()

if __name__ == "__main__":
    main()
