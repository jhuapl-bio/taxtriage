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
##############################################################################################
import gzip
import argparse
import time
import re
import csv
from Bio import SeqIO
import os
# import function determine_priority_assembly from fil in same dir
from determine_priority_assembly import determine_priority_assembly, format_description, parse_fasta_header

def parse_args():
    parser = argparse.ArgumentParser(description="Fetch assembly accessions for a list of nuccore accessions.")
    parser.add_argument("-i", "--input", required=True, help="Path to a file containing nuccore accessions, one per line.")
    parser.add_argument("-a", "--assembly", required=True, help="Path to a file of assembly refseq information")
    parser.add_argument("-p", "--pathogens", required=False, help="OPTIONAL - if assembly mapping fails for an organism, attempt to query using the  pathogen sheet's name and or taxid columns")
    parser.add_argument("-o", "--output", required=True, help="Output File for GCFs")
    return parser.parse_args()

def search_in_dicts(organism_name, strain, assemblies, gcfs):
    if organism_name in assemblies:
        if strain and strain in assemblies[organism_name]:
            # print(f"Match found for {organism_name} and strain {strain}")
            # You can return or process the matched data here
            return organism_name, gcfs[organism_name][strain]
        else:
            # print(f"Organism {organism_name} found, but strain {strain} is not present.")
            return None, None
    else:
        # print(f"Organism {organism_name} not found.")
        return None, None

def long_match_names(organism_name, names):
    split_name = organism_name.split(" ")
    match = False
    # iterate through split_name combining index each time until full combine, for each combined iteration check if it is names. If so , print and return
    # move backwards from the full length to only first word until first match
    for i in range(len(split_name), 0, -1):
        name = " ".join(split_name[:i])
        if  name in names:
            match = name
            break
    return match

def extract_after_colon(text):
    # Regular expression to match everything after ': '
    match = re.search(r':\s*(.*)', text)
    if match:
        return match.group(1)  # Return the matched group, which is everything after ': '
    else:
        return text  # Return the original text if no ': ' is found

from collections import defaultdict

def main():
    args = parse_args()
    assemblies = dict()
    gcfs = dict()
    backup_mapping = defaultdict(list )

    if args.pathogens:
        # import the csv file and read the columns using csv writer/reader
        with open(args.pathogens, 'r') as csvfile:
            # read file as dictionary
            csvreader = csv.DictReader(csvfile)
            for row in csvreader:
                backup_mapping[row['name']].append(
                    dict(taxid=row['taxid'], basename=row['name'])
                )
                alternative_names = row.get('alternative_names', []).split(";")
                if len(alternative_names) > 0:
                    for altname in alternative_names:
                        backup_mapping[altname].append(
                            dict(taxid=row['taxid'], basename=row['name'])
                        )
    parent_taxids = dict()
    taxids = dict()
    with open(args.assembly, 'r') as assembly_file:
        lines = assembly_file.readlines()
        for line in lines:
            if line.startswith("#"):
                pass
            cols = line.split("\t")
            if len(cols) > 7:
                name = cols[7]
                strain = cols[8]
                strain = strain.split("=")[1] if "=" in strain else strain
                if strain == "na":
                    strain = "None"
                ## set 2d list for [name] in the assemblies dict
                if name not in gcfs:
                    gcfs[name] = dict()
                if  name in assemblies and assemblies[name] == 0:
                    continue
                priority = determine_priority_assembly(line)
                if name not in assemblies or (name in assemblies and priority < assemblies[name]):
                    assemblies[name] = priority
                    taxid = cols[5]
                    # extract value of strain=value. If strain=value is not present, set to just the same value for strain
                    gcfs[name] = cols[0]
                    taxids[taxid] = cols[0]
                    parent_taxids[taxid] = cols[6]
    assembly_file.close()
    accessions = dict()
    total = 0
    count = 0
    final_list = []
    missing_elements = []
    if os.path.exists(args.input):
        # Determine the file type based on the extension
        file_mode = 'rt'  # text mode for reading

        # Use gzip.open if it's a gzipped file, otherwise open normally
        if args.input.endswith('.gz'):
            file_handler = gzip.open(args.input, file_mode)
        else:
            file_handler = open(args.input, file_mode)
        # if the file is compressed decompress and read otherwise just read
        with file_handler as file:
            for seq_record in SeqIO.parse(file, "fasta"):
                id, desc = format_description(seq_record.id, seq_record.description)
                desc = desc.replace("TPA_inf: ", "")
                desc = desc.replace("MAG: ", "")
                desc = desc.replace("MAGL: ", "")
                # if there is something like MAG: or TPA_inf: in the description, remove it but use fuzzy match so that if it is up : and a space then remove it, dont hardcode the values.
                # desc = extract_after_colon(desc)
                total +=1
                name = long_match_names(desc, assemblies.keys())
                # organism, gcf = search_in_dicts(organism_name, explicit_strain, assemblies, gcfs)
                if name:
                    count+=1
                    final_list.append(f"{seq_record.id}\t{gcfs[name]}\t{name}\t{desc}")
                else:
                    # print(f"No match found in teh assembly refseq file for {seq_record.id}\t{desc}")
                    # print(f"Attempting to match through backup sheet of name->taxid")
                    name = long_match_names(desc, backup_mapping.keys())
                    if name:
                        count+=1
                        # try to get GCFS from the backup sheet
                        gcf = None
                        for taxid_names in backup_mapping.get(name, []):
                            taxid = taxid_names.get('taxid', None)
                            gcf = taxids.get(taxid, None)
                            basename = taxid_names.get('basename', name)
                            if gcf:
                                break
                        if gcf:
                            final_list.append(f"{seq_record.id}\t{gcf}\t{basename}\t{desc}")
                        else:
                            missing_elements.append(f"{seq_record.id}\t{desc}")
                    else:
                        # print(f"No match found in the backup sheet for {seq_record.id}\t{desc}")
                        missing_elements.append(f"{seq_record.id}\t{desc}")
    if total > 0:
        print(f"Found {count} out of {total} FASTA accessions ({count/total*100:.2f}%) in the assembly file.")


    if len(final_list) == 0:
        print("No assemblies found")
    else:
        with open(args.output, 'w') as infile:
            for line in final_list:
                infile.write(line + "\n")
        infile.close()
    if len(missing_elements) > 0:
        dirname = os.path.dirname(args.output)
        with open(os.path.join(dirname, "missing_assemblies.txt"), 'w') as missing_file:
            for line in missing_elements:
                missing_file.write(line + "\n")
        missing_file.close()




if __name__ == "__main__":
    main()
