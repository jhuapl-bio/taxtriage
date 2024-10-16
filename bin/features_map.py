#!/usr/bin/env python
import argparse
import csv
from collections import defaultdict

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Count unique feature occurrences in BAM file regions and report counts both per accession and aggregated per organism.",
        epilog="Example usage: python count_features.py --cov my_features.bed --map my_mapping.txt --output feature_counts.txt --group_output organism_feature_counts.txt",
    )
    parser.add_argument("-c", "--cov", required=True, help="BED file coverage input file")
    parser.add_argument("-m", "--map", required=True, help="Mapping file correlating accessions to organisms")
    parser.add_argument("-o", "--output", required=True, help="Output file to write feature counts per accession")
    parser.add_argument("-g", "--group_output", help="Output file to write aggregated feature counts per organism")
    return parser.parse_args(argv)

def read_mapping(map_file):
    """
    Read the accession to organism mapping file and return a dictionary of mappings.
    """
    accession_to_organism = {}
    with open(map_file, "r") as mapf:
        reader = csv.reader(mapf, delimiter='\t')
        for row in reader:
            accession = row[0]
            organism = row[2]
            taxid = row[4]  # Assuming taxid is in the 5th column (index 4)
            accession_to_organism[accession] = (organism, taxid)
    return accession_to_organism

def process_coverage_file(cov_file, accession_to_organism):
    """
    Process the BED coverage file and return counts per accession and per organism.
    """
    counts_per_accession = defaultdict(lambda: defaultdict(lambda: [0, 0.0]))  # {accession: {feature_name: [count, sum_fraction]}}
    counts_per_organism = defaultdict(lambda: defaultdict(lambda: [0, 0.0]))   # {organism: {feature_name: [count, sum_fraction]}}

    with open(cov_file, "r") as covf:
        reader = csv.reader(covf, delimiter='\t')
        for row in reader:
            accession = row[0]
            start = row[1]
            end = row[2]
            feature_name = row[3]
            count = int(row[4])
            fraction = float(row[7])

            if accession in accession_to_organism:
                organism, taxid = accession_to_organism[accession]

                # Update counts per accession
                counts_per_accession[accession][feature_name][0] += count
                counts_per_accession[accession][feature_name][1] += fraction

                # Update counts per organism
                counts_per_organism[organism][feature_name][0] += count
                counts_per_organism[organism][feature_name][1] += fraction

    return counts_per_accession, counts_per_organism

def write_output(output_file, counts_per_accession, accession_to_organism):
    """
    Write the per accession output to a file, skipping features with count == 0.
    """
    with open(output_file, "w") as outf:
        writer = csv.writer(outf, delimiter='\t')
        writer.writerow(["Organism", "TaxID", "Accession", "Feature", "Count", "Fraction"])

        for accession, features in counts_per_accession.items():
            organism, taxid = accession_to_organism[accession]
            for feature, (count, fraction) in features.items():
                if count > 0:  # Skip if count is 0
                    writer.writerow([organism, taxid, accession, feature, count, f"{fraction:.2f}"])

def write_grouped_output(group_output_file, counts_per_organism):
    """
    Write the aggregated per organism output to a file, skipping features with count == 0.
    """
    with open(group_output_file, "w") as outf:
        writer = csv.writer(outf, delimiter='\t')
        writer.writerow(["Organism", "Feature", "Count", "Fraction"])

        for organism, features in counts_per_organism.items():
            for feature, (count, fraction) in features.items():
                if count > 0:  # Skip if count is 0
                    writer.writerow([organism, feature, count, f"{fraction:.2f}"])

def main(argv=None):
    args = parse_args(argv)

    # Read the accession to organism mapping file
    accession_to_organism = read_mapping(args.map)

    # Process the coverage file
    counts_per_accession, counts_per_organism = process_coverage_file(args.cov, accession_to_organism)

    # Write output per accession
    write_output(args.output, counts_per_accession, accession_to_organism)

    # Write aggregated output per organism if -g is provided
    if args.group_output:
        write_grouped_output(args.group_output, counts_per_organism)

if __name__ == "__main__":
    main()
