#!/usr/bin/env python3

from Bio import Entrez
import argparse
import time

def parse_args():
    parser = argparse.ArgumentParser(description="Fetch assembly accessions for a list of nuccore accessions.")
    parser.add_argument("-i", "--input", required=True, help="Path to a file containing nuccore accessions, one per line.")
    parser.add_argument("-r", "--restart", action='store_true', help="Email for entrez")
    parser.add_argument("-o", "--file_out", required=True, help="Output file to append the assembly accessions.")
    parser.add_argument("-e", "--email", required=False, help="Email for entrez")
    parser.add_argument("-b", "--batch_count", required=False, default=100, help="Batch count for fetching assembly accessions. Default is 100.")
    parser.add_argument("-f", "--failed", default="failed_accessions.txt", help="File to save failed accessions.")
    return parser.parse_args()

def read_processed_accessions(output_file):
    processed = set()
    try:
        with open(output_file, 'r') as f:
            for line in f:
                parts = line.strip().split("\t")
                if parts:
                    processed.add(parts[0])
    except FileNotFoundError:
        pass  # File doesn't exist, which is fine for the first run
    return processed

def fetch_assembly_accession(nuccore_id, retries=2):
    try:
        handle = Entrez.elink(dbfrom="nuccore", db="assembly", id=nuccore_id, retmax=1)
        link_record = Entrez.read(handle)
        handle.close()
        if link_record[0]["LinkSetDb"]:
            assembly_id = link_record[0]["LinkSetDb"][0]["Link"][0]["Id"]
            handle = Entrez.esummary(db="assembly", id=assembly_id)
            summary_record = Entrez.read(handle)
            handle.close()
            return summary_record["DocumentSummarySet"]["DocumentSummary"][0]["AssemblyAccession"]
    except Exception as e:
        if retries > 0:
            print(f"Failed to fetch assembly accession for {nuccore_id}. Retrying {retries} more times.")
            time.sleep(1)  # Wait for 5 seconds before retrying
            return fetch_assembly_accession(nuccore_id, retries-1)
    return None

def main():
    args = parse_args()
    if args.email:
        Entrez.email = args.email
    processed_accessions = read_processed_accessions(args.file_out)
    failed_accessions = []

    with open(args.input, 'r') as infile:
        for line in infile:
            nuccore_id = line.strip()
            print(nuccore_id)
            if nuccore_id and nuccore_id not in processed_accessions:
                assembly_accession = fetch_assembly_accession(nuccore_id)
                print(assembly_accession)
                if assembly_accession:
                    mode="w" if args.restart else "a"
                    with open(args.file_out, mode) as outfile:
                        print(f"{nuccore_id}\t{assembly_accession}\n")
                        outfile.write(f"{nuccore_id}\t{assembly_accession}\n")
                else:
                    failed_accessions.append(nuccore_id)
    if failed_accessions:
        with open(args.failed, 'w') as failed_file:
            for acc in failed_accessions:
                failed_file.write(acc + "\n")

if __name__ == "__main__":
    main()
