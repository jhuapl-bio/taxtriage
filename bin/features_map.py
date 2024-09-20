#!/usr/bin/env python3
import argparse
import pysam
import csv
import os

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Count unique feature occurrences in BAM file regions and report counts both per accession and aggregated per organism.",
        epilog="Example usage: python count_features.py --bed my_features.bed --bam my_alignments.bam --map my_mapping.txt --output feature_counts.txt --group_output organism_feature_counts.txt",
    )
    parser.add_argument("-b", "--bed", required=True, help="BED file with feature annotations")
    parser.add_argument("-a", "--bam", required=True, help="Indexed BAM file with alignment data")
    parser.add_argument("-m", "--map", required=True, help="Mapping file correlating accessions to organisms")
    parser.add_argument("-o", "--output", required=True, help="Output file to write feature counts per accession")
    parser.add_argument("-i", "--index", required=False, help="BAM Bai file")
    parser.add_argument("-g", "--group_output", required=True, help="Output file to write aggregated feature counts per organism")
    return parser.parse_args(argv)

def determine_paired_single(bam_file):
    max_reads_to_check = 5
    read_count = 0

    for read in bam_file.head(max_reads_to_check):
        if read.is_paired:
            return True
        read_count += 1

    return False

def count_features(args):
    if not args.index and not os.path.exists(args.bam+".bai"):
        pysam.index(args.bam)
    else:
        print(f"Index file found  ")
    # Load mapping file
    map_dict = {}
    mapname = dict()
    with open(args.map, 'r') as mapfile:
        reader = csv.reader(mapfile, delimiter='\t')
        for row in reader:
            map_dict[row[0]] = row[2]  # Assuming the third column is organism name
            if len(row) > 4:
                mapname[row[0]] = row[4]
    # Read BED file and initialize count structures
    feature_counts = {}
    with open(args.bed, 'r') as bedfile:
        reader = csv.reader(bedfile, delimiter='\t')
        for row in reader:
            feature_counts[row[0]] = { 'features': {}, 'total_reads': 0 }

    bam_file = pysam.AlignmentFile(args.bam, "rb")
    is_paired = determine_paired_single(bam_file)
    bam_file.close()
    bam_file = pysam.AlignmentFile(args.bam, "rb")

    total_fragments = 0
    read_names = set()

    if is_paired:
        for read in bam_file.fetch(until_eof=True):
            if (read.is_read1 or read.is_read2) and not read.is_secondary and not read.is_supplementary and read.query_name not in read_names:
                total_fragments += 1
                read_names.add(read.query_name)
    else:
        for read in bam_file.fetch(until_eof=True):
            if not read.is_secondary and not read.is_supplementary:
                total_fragments += 1

    print(f"Total fragments: {total_fragments}")

    # Count reads for each feature
    for acc in feature_counts.keys():
        with open(args.bed, 'r') as bedfile:
            reader = csv.reader(bedfile, delimiter='\t')
            for row in reader:
                if row[0] == acc:
                    try:
                        count = bam_file.count(contig=row[0], start=int(row[1]), end=int(row[2]))
                        feature_counts[acc]['total_reads'] += count
                        feature_counts[acc]['features'][row[3]] = count
                    except ValueError:
                        continue

    # Write the output files
    with open(args.output, 'w', newline='') as outfile:
        writer = csv.writer(outfile, delimiter='\t')
        writer.writerow([ 'Assembly Name', 'Mapped Value', 'Accession', 'Feature', '# reads align', 'Percent Read Align'])
        for acc, data in feature_counts.items():
            # sort feature items by count
            mapv = mapname.get(acc, 'Unknown')
            data['features'] = dict(sorted(data['features'].items(), key=lambda item: item[1], reverse=True))
            for feature, count in data['features'].items():
                if count > 0:
                    percent_read_align = (100 * count / total_fragments) if total_fragments > 0 else 0
                    writer.writerow([  map_dict.get(acc, 'Unknown'), mapv, acc, feature, count, f"{percent_read_align:.2f}"])
    # Aggregate per organism and feature
    organism_feature_counts = {}

    for acc, data in feature_counts.items():
        organism = map_dict.get(acc, 'Unknown')
        for feature, count in data['features'].items():
            if count > 0:
                # Check if we have this organism already
                if organism not in organism_feature_counts:
                    organism_feature_counts[organism] = {}

                # Check if we have this feature for the organism already
                if feature in organism_feature_counts[organism]:
                    organism_feature_counts[organism][feature] += count
                else:
                    organism_feature_counts[organism][feature] = count

    # Now, write the organism and feature-specific counts to the group output file
    with open(args.group_output, 'w', newline='') as groupfile:
        writer = csv.writer(groupfile, delimiter='\t')
        writer.writerow(['Organism Name', 'Mapped Value', 'Feature', '# reads align', 'Percent Read Align'])
        mapv = mapname.get(acc, 'Unknown')

        for organism, features in organism_feature_counts.items():
            # sort feature items9 by count
            features = dict(sorted(features.items(), key=lambda item: item[1], reverse=True))
            for feature, count in features.items():
                percent_read_align = (100 * count / total_fragments) if total_fragments > 0 else 0
                writer.writerow([organism,  mapv, feature, count, "{:.2f}".format(percent_read_align)])

if __name__ == "__main__":
    args = parse_args()
    count_features(args)
