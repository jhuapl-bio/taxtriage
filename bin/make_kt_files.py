#!/usr/bin/env python3

import argparse


def count_leading_spaces(s):
    return int((len(s) - len(s.lstrip()) )  / 2)

def parse_kraken_report(file_path, root="root"):
    taxonomy_tree = []
    # lineage_by_indent = {0: ['root']}
    current_lineage = []

    with open(file_path, 'r') as file:
        i=0
        for line in file:

            i+=1
            # if (i > 18):
            #     break;
            parts = line.strip().split('\t')
            if len(parts) >= 6:
                count = int(parts[2])  # Read count
                tax_name_with_indent = parts[5]  # Taxonomic name with leading spaces
                # tax_name = tax_name_with_indent.strip()
                indent = count_leading_spaces(tax_name_with_indent)

                # Check if the length of the lineage is greater than index, if so, truncate up to that index and replace index at indent with current taxname
                if (len(current_lineage) > indent):
                    current_lineage = current_lineage[:indent]
                current_lineage.append(tax_name_with_indent.strip())
                if root != current_lineage[0] and current_lineage[0] != "unclassified":
                    current_lineage[0] = root
                if (count > 0):
                    taxonomy_tree.append({
                        'count': count,
                        'lineage': current_lineage.copy()
                    })
    # for key, val in enumerate(taxonomy_tree):
    #     print(val)
    return taxonomy_tree



def write_krona_input(taxonomy_tree, output_file_path):
    with open(output_file_path, 'w') as out_file:
        for entry in taxonomy_tree:
            lineage_str = '\t'.join(entry['lineage'])  # Exclude 'root' from the lineage
            out_file.write(f"{entry['count']}\t{lineage_str}\n")


def main():
    parser = argparse.ArgumentParser(description="Convert Kraken report to Krona input format")
    parser.add_argument("-i", "--input", help="Path to the Kraken report file")
    parser.add_argument("-o", "--output", help="Path for the output file")
    parser.add_argument("-r", "--root", default="root", help="Root name: default is: root")

    args = parser.parse_args()

    taxonomy_tree = parse_kraken_report(args.input, args.root)
    write_krona_input(taxonomy_tree, args.output)
    print(f"Krona input file created at {args.output}")

if __name__ == "__main__":
    main()
