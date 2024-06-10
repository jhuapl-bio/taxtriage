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

"""Provide a command line tool to fetch a list of refseq genome ids to a single file, useful for kraken2 database building or alignment purposes"""
from collections import defaultdict
import sys
import os
import gzip
import argparse
import csv
import pysam
from math import log2



def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and transform a tabular samplesheet.",
        epilog="Example: python check_samplesheet.py samplesheet.csv samplesheet.valid.csv",
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        metavar="INPUT",
        help="BAM File to process",
    )
    parser.add_argument(
        "-d",
        "--depth",
        metavar="DEPTH",

        help="Depth File (from samtools) corresponding to your sample",
    )
    parser.add_argument(
        "-m",
        "--match",
        metavar="MATCH",
        help="Reference accessions match, 2 cols min: accession (NZ/NC), name",
    )
    parser.add_argument(
        "-r",
        "--min_reads_align",
        metavar="MINREADSALIGN",
        default=3,
        help="Filter for minimum reads aligned to reference per organism. Default is 1",
    )
    parser.add_argument(
        "-v",
        "--mincoverage",
        metavar="Minimum coverage value to consider acceptable cutoff for confidence. Anything above == confidence",
        default=1,
        type=int,
    )
    parser.add_argument(
        "-c",
        "--capval",
        metavar="Cap val for depth for entropy calculation. Default: 15",
        default=15,
        type=int,
        help="At what threshold to cutoff for determining shannon entropy stats for organism depth of cov. Default is 15, -1 for no cap value",
    )
    parser.add_argument(
        "-a",
        "--accessioncol",
        metavar="ACCCOL",
        default=0,
        help="Index of the column in mapfile (if specified) to match to the reference accession. 0 index start",
    )
    parser.add_argument(
        "-q",
        "--taxcol",
        metavar="TAXCOL",
        default=4, required =False,
        help="Index of the column in mapfile (if specified) to add taxid col",
    )
    parser.add_argument(
        "-n",
        "--namecol",
        default=2,
        metavar="NAMECOL",
        help="Index of the column in mapfile (if specified) to match to the name. 0 index start",
    )
    parser.add_argument(
        "-x",
        "--coverage",
        metavar="COVERAGEFILE",
        default=None,
        help="Samtools coverage file",
    )
    parser.add_argument(
        "-j",
        "--assembly",
        metavar="ASSEMBLY",
        default=None, required=False,
        help="Assembly refseq file",
    )
    parser.add_argument(
        "-k",  "--compress_species", default=False,  help="Compress species to species level",  action='store_true'
    )
    parser.add_argument(
        "-t",
        "--sampletype",
        metavar="SAMPLETYPE",
        default='Unknown',
        help="Sample Type to process. If Empty, defaults to null",
    )
    parser.add_argument(
        "-s",
        "--samplename",
        metavar="SAMPLENAME",
        default="No_Name",
        help="Name of the sample to process. If Empty, defaults to 'No_Name'",
    )
    parser.add_argument(
        "-p",
        "--pathogens",
        metavar="PATHOGENS",
        help="TXT File to process. Must be in header format: ",
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="OUTPUT",
        type=str,
        help="Name of output directory for feature table files",
    )


    return parser.parse_args(argv)

def lorenz_curve(depths):
    """Compute the Lorenz curve for a list of depths."""
    sorted_depths = sorted(depths)
    cumulative_depths = [0]
    total = sum(sorted_depths)
    for depth in sorted_depths:
        cumulative_depths.append(cumulative_depths[-1] + depth)
    lorenz_curve = [x / total if x >0 else 0 for x in cumulative_depths]
    return lorenz_curve


def gini_coefficient(depths):
    """Calculate the Gini coefficient for a list of depths."""
    lorenz = lorenz_curve(depths)
    n = len(lorenz)
    area_under_lorenz = sum((lorenz[i] + lorenz[i + 1]) / 2 for i in range(n - 1)) / (n - 1)
    gini = 1 - 2 * area_under_lorenz
    return gini

def breadth_of_coverage(depths, genome_length):
    """Calculate the breadth of coverage for a list of depths."""
    non_zero_positions = len([depth for depth in depths if depth > 0])
    return non_zero_positions / genome_length if genome_length > 0 else 0

def adjusted_fair_distribution_score(depths, genome_length):
    """Calculate the adjusted score for fair distribution of depths."""
    gini = gini_coefficient(depths)
    breadth = breadth_of_coverage(depths, genome_length)

    total_depth = sum(depths)
    if total_depth == 0:
        return 0  # If there's no coverage, we consider it the lowest score

    # Compute penalty based on the distribution of depths
    max_depth = max(depths)
    avg_depth = sum(depths) / len(depths) if len(depths) > 0 else 0

    # Define a penalty factor that reduces the impact of high depths
    penalty_factor = 0.5  # You can adjust this factor based on desired sensitivity

    if avg_depth > 0:
        penalty = penalty_factor * (max_depth / avg_depth - 1) / (max_depth / avg_depth + 1)
    else:
        penalty = penalty_factor * max_depth

    # Adjust the score: 1 - Gini coefficient, penalized by the penalty factor
    gini_score = (1 - gini) * (1 - penalty)

    # Ensure both scores are between 0 and 1
    gini_score = max(0, min(1, gini_score))
    breadth = max(0, min(1, breadth))
    # log transform breadth to 0 and 1, more weight closer to 1
    if breadth > 0:
        breadth = log2(breadth + 1) / log2(2)

    # Combine the Gini score and breadth of coverage with equal weight
    final_score = 0.1 * gini_score + 0.9 * breadth
    return final_score


def get_fair_distribution_score(data):
    # Assuming 'data' is a dictionary with 'depths' as a list and 'total_length' as an int
    depths = data.get('depths', [])
    genome_length = data.get('total_length', len(depths))  # Use provided total length or length of depths

    # Calculate the adjusted fair distribution score
    fair_distribution_score = adjusted_fair_distribution_score(depths, genome_length)

    return fair_distribution_score


def import_pathogens(pathogens_file):
    """Import the pathogens from the input CSV file, correctly handling commas in quoted fields."""
    pathogens_dict = {}
    # Open the file using the `with` statement
    with open(pathogens_file, 'r', newline='', encoding='ISO-8859-1') as file:
        # Create a CSV reader object that handles commas inside quotes automatically
        reader = csv.reader(file, delimiter=',', quotechar='"')

        # Iterate over each row in the CSV file
        for row in reader:
            # Assign each part of the row to variables if available
            pathogen_name = row[0] if len(row) > 0 else None
            taxid = row[1] if len(row) > 1 else None
            call_class = row[2] if len(row) > 2 else None

            pathogenic_sites = row[3] if len(row) > 3 else None
            commensal_sites = row[4] if len(row) > 4 else None
            status = row[5] if len(row) > 5 else None
            pathology = row[6] if len(row) > 6 else None

            # Store the data in the dictionary, keyed by pathogen name
            pathogens_dict[pathogen_name] = {
                'taxid': taxid,
                'callclass': call_class,
                'pathogenic_sites': pathogenic_sites,
                'name': pathogen_name,
                'commensal_sites': commensal_sites,
                'status': status,
                'pathology': pathology
            }
            pathogens_dict[taxid] = pathogens_dict[pathogen_name]


    # No need to explicitly close the file, `with` statement handles it.
    return pathogens_dict

def identify_pathogens(inputfile, pathogens):
    """Identify the pathogens in the input file"""
    with open(inputfile, 'r') as f:
        for line in f:
            if line.startswith('pathogen'):
                pathogens = line.split('\t')
                return pathogens
    f.close()
    return None

def calculate_entropy(values):
    """Calculate the Shannon entropy of the given values without NumPy."""
    total = sum(values)
    probabilities = [value / total for value in values]
    return -sum(p * log2(p) for p in probabilities if p > 0)



# Function to calculate weighted mean
def calculate_weighted_mean(data, numreads):
    if len(numreads) > 0:
        total_weight = sum(numreads)  # Sum of all weights (numreads)
        weighted_sum = []
        for x in range(0, len(numreads)):
            if x >= len(data):
                break
            if numreads[x] > 0:
                weighted_sum.append(data[x] * numreads[x]) # Sum of weight*value
            else:
                weighted_sum.append(0)
        weighted_sum = sum(weighted_sum)
        if total_weight > 0:
            weighted_mean = weighted_sum / total_weight
        else:
            return 0
    else:
        weighted_mean = 0
    return weighted_mean


def count_reference_hits(bam_file_path, depthfile, covfile, matchdct):
    """
    Count the number of reads aligned to each reference in a BAM file.

    Args:
    bam_file_path (str): Path to the BAM file.

    Returns:
    dict: A dictionary with reference names as keys and counts of aligned reads as values.
    """
    # Initialize a dictionary to hold the count of reads per reference
    reference_counts = {}
    reference_coverage = defaultdict(lambda: defaultdict(dict))
    reference_lengths = {}
    notseen = set()
    unaligned = 0
    aligned_reads = 0
    total_reads = 0
    amount_pre_read = dict()
    # Open the BAM file for reading
    testacc = None


    with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
        # get total reads
        total_reads = sum(1 for _ in bam_file)
        for ref in bam_file.header.references:
            reference_lengths[ref] = bam_file.get_reference_length(ref)
            reference_coverage[ref] = dict(
                length = reference_lengths[ref],
                depths =  defaultdict(int),
                gini_coefficient = 0,
                mapqs = [],
                baseqs = [],
                meanmapq = 0,
                numreads = 0,
                meanbaseq = 0,
                coverage = 0,
                meandepth = 0,
                name = ref,
                accession = ref,
                isSpecies = False,
            )
        if covfile:
            print("Reading coverage information from coverage file")
            with open(covfile, 'r') as f:
                for line in f:
                    line = line.strip()
                    splitline = line.split('\t')
                    # if not header line
                    if "#rname"  in splitline:
                        continue
                    reference_name = splitline[0]
                    numreads = int(splitline[3])
                    covbases = float(splitline[4])
                    coverage = float(splitline[5])
                    meandepth = float(splitline[6])
                    meanbaseq = float(splitline[7])
                    meanmapq = float(splitline[8])
                    if reference_name in reference_coverage:
                        reference_coverage[reference_name]['numreads'] = numreads
                        reference_coverage[reference_name]['coverage'] = coverage
                        reference_coverage[reference_name]['meandepth'] = meandepth
                        reference_coverage[reference_name]['meanbaseq'] = meanbaseq
                        reference_coverage[reference_name]['meanmapq'] = meanmapq
            f.close()
        if not depthfile or not covfile:
            print("No depthfile or covfile supplied, reading input from bam file")
            for read in bam_file.fetch(until_eof=True):
                if not read.is_unmapped:  # Check if the read is aligned
                    reference_name = bam_file.get_reference_name(read.reference_id)
                    aligned_reads +=1
                    # Accumulate coverage data
                    if not covfile:
                        # get the baseq of the read
                        # Get the base qualities
                        base_qualities = read.query_qualities
                        # Convert base qualities to ASCII
                        base_qualities_ascii = [chr(qual + 33) for qual in base_qualities]
                        # Convert ASCII to integers
                        base_qualities_int = [ord(qual) - 33 for qual in base_qualities_ascii]
                        baseq = sum(base_qualities) / len(base_qualities)
                        # get the mapq of the read
                        mapq = float(read.mapping_quality)
                        reference_coverage[reference_name]['mapqs'].append(mapq)
                        reference_coverage[reference_name]['baseqs'].append(baseq)
                        reference_coverage[reference_name]['numreads']+=1
                    if not depthfile:
                        for pos in range(read.reference_start, read.reference_end):
                            reference_coverage[reference_name]['depths'][pos] += 1
                else:
                    unaligned += 1

        if not covfile:
            # make meanbaseq and meanmapq
            for ref, data in reference_coverage.items():
                if data['numreads'] > 0:
                    data['meanbaseq'] = sum(data['baseqs']) / len(data['baseqs'])
                    data['meanmapq'] = sum(data['mapqs']) / len(data['mapqs'])

    bam_file.close()
    if depthfile:
        print("Reading depth information from depth file")
        # read in depthfile and reference_coverage[reference_name][pos]  as value
        with open(depthfile, 'r') as f:
            for line in f:
                splitline = line.split('\t')
                reference_name = splitline[0]
                pos = int(splitline[1])
                depth = int(splitline[2])
                reference_coverage[reference_name]['depths'][pos] = depth
        f.close()
     # make a new dict which is name then count, percentage across total reads and percentage across aligned reads
    i=0
    return reference_coverage, total_reads

def main():
    args = parse_args()
    inputfile = args.input
    pathogenfile = args.pathogens
    covfile = args.coverage
    output = args.output
    matcher = args.match
    matchdct = dict()
    header = True
    i =0
    capval  = args.capval
    mincov = args.mincoverage
    reference_hits, total_reads = count_reference_hits(
        inputfile,
        args.depth,
        covfile,
        matchdct
    )
    testacc = "CP000253.1"
    for key, value in reference_hits.items():
        if not testacc:
            testacc = key
        break
    assembly_to_accession = defaultdict(set)
    if args.match:
        # open the match file and import the match file
        header = True
        i=0
        with open (matcher, 'r') as f:
            accindex = args.accessioncol
            nameindex = args.namecol
            taxcol = args.taxcol
            for line in f:
                line = line.strip()
                if header and i==0:
                    header = False
                    continue
                taxid=None
                assembly = None
                splitline = line.split('\t')
                if len(splitline) > 0:
                    accession = splitline[accindex]
                else:
                    accession = None
                if len(splitline) > 0:
                    assembly = splitline[1]
                if len(splitline) > taxcol:
                    taxid = splitline[taxcol]
                if len(splitline) > nameindex:
                    name = splitline[nameindex]
                else:
                    name = None
                if accession in reference_hits:
                    reference_hits[accession]['taxid'] = taxid
                    reference_hits[accession]['name'] = name
                    reference_hits[accession]['assembly'] = assembly
                if accession not in assembly_to_accession[assembly]:
                    assembly_to_accession[assembly].add(accession)
                i+=1

        f.close()
    # print(f"Step2: {testacc}, {reference_hits[testacc]}","\n\n")
    spectaxidmatch = dict()
    collect_subspecies = defaultdict(list)
    if args.assembly:
        i=0
        with open(args.assembly, 'r') as f:
            for line in f:
                line = line.strip()
                i+=1
                if i<2:
                    continue
                splitline = line.split('\t')
                if len(splitline) > 0:
                    accession = splitline[0]
                    taxid = splitline[5]
                    species_taxid = splitline[6]
                    name = splitline[7]
                    strain = splitline[8].replace("strain=", "")
                    isSpecies = False if species_taxid != taxid else True
                    # fine value where assembly == accession from reference_hits
                    if accession in assembly_to_accession:
                        for acc in assembly_to_accession[accession]:
                            if acc in reference_hits:
                                reference_hits[acc]['isSpecies'] = isSpecies
                                reference_hits[acc]['toplevelkey'] = species_taxid
                                reference_hits[acc]['strain'] = strain
                                reference_hits[acc]['assemblyname'] = name
                                reference_hits[acc]['name'] = name

        f.close()
    # print(f"Step3: {testacc}, {reference_hits[testacc]}","\n\n")

    final_format = defaultdict(dict)

    if args.compress_species:
        # Need to convert reference_hits to only species species level in a new dictionary
        for key, value in reference_hits.items():
            value['accession'] = key
            if 'toplevelkey' in value and value['toplevelkey']:
                valtoplevel = value['toplevelkey']
            else:
                valtoplevel = key
            valkey = value['accession']

            final_format[valtoplevel][valkey]= value

    else:
        # We don't aggregate, so do final format on the organism name only
        for key, value in reference_hits.items():
            if value['toplevelkey']:
                valtoplevel = value['toplevelkey']
            else:
                valtoplevel = key
            valkey = key
            final_format[valtoplevel][valkey] = value
    # for key, value in final_format.items():
    #     for key2, value2 in value.items():
    #         print(value2.keys())
            # if value2['accession'] == testacc:
            #     print(f"Step4: {key}, {key2}","\n\n")
            # break
    # Dictionary to store aggregated species-level data
    species_aggregated = {}
    def getGiniCoeff(data, acc_length=0):
        # Assuming 'data' is a dictionary with 'depths' as a list and 'total_length' as an int.
        depths = list(data.values()) if data.values() else [0]
        # Get get proportion of depths in increments 10 or higher like 0-10 x 11-20 x etc


        # Gini Coefficient
        gini_coefficient = get_fair_distribution_score({"depths": depths, 'total_length': acc_length})
        return gini_coefficient
    # Aggregate data at the species level
    for top_level_key, entries in final_format.items():
        for val_key, data in entries.items():

            if top_level_key not in species_aggregated:
                species_aggregated[top_level_key] = {
                    'key': top_level_key,
                    'numreads': [],
                    'mapqs': [],
                    'depths': [],
                    "accs": [],
                    "coverages": [],
                    "coeffs": [],
                    'baseqs': [],
                    "isSpecies": True if  args.compress_species  else data['isSpecies'],
                    'strainslist': [],
                    'name': data['name'],  # Assuming the species name is the same for all strains
                }
            gini_strain = getGiniCoeff(data['depths'], data['length'])
            species_aggregated[top_level_key]['coeffs'].append(gini_strain)
            species_aggregated[top_level_key]['numreads'].append(data['numreads'])
            species_aggregated[top_level_key]['coverages'].append(data['coverage'])
            species_aggregated[top_level_key]['baseqs'].append(data['meanbaseq'])
            species_aggregated[top_level_key]['accs'].append(data['accession'])
            species_aggregated[top_level_key]['mapqs'].append(data['meanmapq'])
            species_aggregated[top_level_key]['depths'].append(data['meandepth'])
            if 'strain' in data:
                species_aggregated[top_level_key]['strainslist'].append({
                    "strainname":data['strain'],
                    "fullname":data['name'],
                    "subkey": val_key,
                    "numreads": data['numreads'],
                    "taxid": data['taxid'] if "taxid" in data else None,
                })
            else:
                species_aggregated[top_level_key]['strainslist'].append({
                    "strainname":val_key,
                    "fullname":val_key,
                    "subkey": val_key,
                    "numreads": data['numreads'],
                    "taxid": data['taxid'] if "taxid" in data else None,
                })


            # species_aggregated[top_level_key]['mapqs'].append(20)
            # species_aggregated[top_level_key]['numreads'].append(data['numreads'])
    # Calculate weighted means for aggregated data
    for top_level_key, aggregated_data in species_aggregated.items():
        numreads = aggregated_data['numreads']
        aggregated_data['meanmapq'] = calculate_weighted_mean(aggregated_data['mapqs'], numreads)
        aggregated_data['meanbaseq'] = calculate_weighted_mean(aggregated_data['baseqs'],numreads)
        aggregated_data['meandepth'] = calculate_weighted_mean(aggregated_data['depths'],numreads)
        aggregated_data['meancoverage'] = calculate_weighted_mean(aggregated_data['coverages'],numreads)
        aggregated_data['meangini'] = calculate_weighted_mean(aggregated_data['coeffs'],numreads)

    # Print the final aggregated data
    for top_level_key, aggregated_data in species_aggregated.items():
        print(f"Entry Top Key: {top_level_key}")
        print(f"\tNum Reads: {aggregated_data['numreads']}")
        print(f"\tMean MapQ: {aggregated_data['meanmapq']}")
        print(f"\tMean BaseQ: {aggregated_data['meanbaseq']}")
        print(f"\tStrains List: {aggregated_data['strainslist']}")
        print()

    for key, value in species_aggregated.items():
        if testacc in value['accs']:
            print(f"Step5: {key}, {value['name']}","\n\n")
            break

    pathogens = import_pathogens(pathogenfile)
    # Next go through the BAM file (inputfile) and see what pathogens match to the reference, use biopython
    # to do this

    if args.min_reads_align:
        # filter the reference_hits based on the minimum number of reads aligned
        print(f"Filtering for minimum reads aligned: {args.min_reads_align}")
        species_aggregated = {k: v for k, v in species_aggregated.items() if sum(v['numreads']) >= int(args.min_reads_align)}
    # print("Step 6:\n")
    # for key, value in species_aggregated.items():
    #     print(f"\t{key}, {value}","\n\n")
    # Create a new dictionary to store aggregated species-level data
    species_references = defaultdict(dict)
    for top_level_key, data in species_aggregated.items():
        print(f"Reference: {data['name']} (Taxid: {data['key']})")
        strainnames = [f"{strain['strainname']} ({strain['taxid']})" for strain in data['strainslist']]
        # print(f"\tStrains seen: {', '.join(strainnames)}")
        print(f"\tNumber of Reads: {data['numreads']}")
        print(f"\tMean Coverage: {data['meancoverage']}")
        print(f"\Alignment Conf: {data['meangini']}")
        print(f"\tDepth of Coverage: {data['meandepth']}")
        print(f"\tMean BaseQ: {data['meanbaseq']}")
        print(f"\tMean MapQ: {data['meanmapq']}")
        print(f"\tisSpecies: {data['isSpecies']}")
        print()



    write_to_tsv(
        aggregated_stats=species_aggregated,
        pathogens=pathogens,
        output_file_path=output,
        sample_name=args.samplename,
        sample_type = args.sampletype,
        total_reads = total_reads
    )
def format_non_zero_decimals(number):
    # Convert the number to a string
    num_str = str(number)
    if '.' not in num_str:
        # If there's no decimal point, return the number as is
        return number
    else:
        # Split into integer and decimal parts
        integer_part, decimal_part = num_str.split('.')
        non_zero_decimals = ''.join([d for d in decimal_part if d != '0'][:2])
        # Count leading zeros in the decimal part
        leading_zeros = len(decimal_part) - len(decimal_part.lstrip('0'))
        # Construct the new number with two significant decimal digits
        formatted_number = f"{integer_part}.{('0' * leading_zeros) + non_zero_decimals}"
        # Convert back to float, then to string to remove trailing zeros
        return str(float(formatted_number))



def write_to_tsv(aggregated_stats, pathogens, output_file_path, sample_name="No_Name", sample_type="Unknown", total_reads=0):
    """
    Write reference hits and pathogen information to a TSV file.

    Args:
    reference_hits (dict): Dictionary with reference names as keys and counts as values.
    pathogens (dict): Dictionary with reference names as keys and dictionaries of additional attributes as values.
    output_file_path (str): Path to the output TSV file.
    """
    with open(output_file_path, 'w') as file:
        # Write the header row

        header =  "Name\tSample\tSample Type\t% Reads\t% Aligned Reads\t# Aligned\tIsAnnotated\tPathogenic Sites\tType\tTaxid\tStatus\tGini Coefficient\tMean BaseQ\tMean MapQ\tMean Coverage\tMean Depth\tAnnClass\tisSpecies\tPathogenic Subsp/Strains\n"
        file.write(f"{header}")
        print("________________________________________")
        total_reads_aligned = 0
        for ref, count in aggregated_stats.items():
            strainlist = count['strainslist']
            isSpecies = count['isSpecies']
            is_pathogen = "Unknown"
            callfamclass = ""
            derived_pathogen = False
            isPath = False
            annClass = "None"
            pathogenic_sites = []
            total_reads_aligned += sum(count['numreads'])
            if ref in pathogens:
                refpath = pathogens[ref]
                pathogenic_sites = refpath['pathogenic_sites'] if 'pathogenic_sites' in refpath else []
            elif ref in pathogens:
                refpath = pathogens[ref]
                pathogenic_sites = refpath['pathogenic_sites'] if 'pathogenic_sites' in refpath else []
            else:
                refpath = None
            def pathogen_label(ref):
                is_pathogen = "Unknown"
                isPathi = False
                callclass= ref['callclass']
                pathogenic_sites = ref['pathogenic_sites']
                commensal_sites = ref['commensal_sites']
                if sample_type in pathogenic_sites:
                    if callclass != "commensal":
                        is_pathogen = callclass.capitalize()
                        isPathi = True
                    else:
                        is_pathogen = "Potential"

                elif sample_type in commensal_sites:
                    is_pathogen = "Commensal"
                elif callclass and callclass != "":
                    is_pathogen = callclass.capitalize() if callclass and callclass !="" else "Unknown"
                    isPathi = True
                return is_pathogen, isPathi
            formatname = count['name']

            if refpath:
                is_pathogen, isPathi = pathogen_label(refpath)
                isPath = isPathi
                if isPathi:
                    annClass = "Direct"
                if ref in refpath and refpath[ref]:
                    taxid = refpath[ref]
                else:
                    taxid = count[ref] if ref in count and count[ref] else ""
                is_annotated = "Yes"
                commsites = refpath['commensal_sites']
                status = refpath['status']
                formatname = refpath['name']
                if is_pathogen == "Commensal":
                    callfamclass = "Commensal Listing"
                # elif is_pathogen != "N/A":
                #     callfamclass = "Listed Pathogen"
                # else:
                #     callfamclass = "Unknown Listing"
            else:
                is_annotated = "No"
                taxid = count[ref] if ref in count and count[ref]  else ""
                sites = ""
                status = ""
            listpathogensstrains = []
            fullstrains = []
            if count['strainslist']:
                pathogenic_reads = 0
                lenstrains = len(count['strainslist'])
                # fullstrains = [x['fullname'] for x in count['strainslist']]
                # check, for each strainslist, if it is listed as a pathogen in pathogens[taxid], get count
                # merge all strainslist on the strainname attribute, sum the numreads
                # if taxid is present, use that, otherwise use the strainname
                merged_strains = defaultdict(dict)
                for x in count['strainslist']:
                    keyx = x['strainname'] if not x['taxid'] else x['taxid']
                    if formatname != x['strainname']:
                        formatname= formatname.replace(x['strainname'], "")
                    if keyx in merged_strains:
                        merged_strains[keyx]['numreads'] += x['numreads']
                        merged_strains[keyx]['subkeys'].append(x['subkey'])
                    else:
                        merged_strains[keyx] = x
                        merged_strains[keyx]['subkeys'] = [x['subkey']]

                for xref, x in merged_strains.items():
                    pathstrain = None
                    if x['taxid']:
                        fullstrains.append(f"{x['strainname']} ({x['taxid']}: {x['numreads']} reads)")
                    else:
                        fullstrains.append(f"{x['strainname']} ({x['numreads']} reads)")
                    if x['taxid'] in pathogens:
                        pathstrain = pathogens[x['taxid']]
                    elif x['fullname'] in pathogens:
                        pathstrain = pathogens[x['fullname']]
                    if pathstrain:
                        taxx  = f"{x['taxid']}" if x['taxid'] else ""
                        if sample_type in pathstrain['pathogenic_sites'] or pathstrain['callclass'] != "commensal":
                            pathogenic_reads+= x['numreads']
                            percentreads = f"{x['numreads']/total_reads:.1f}" if total_reads > 0 and x['numreads'] > 0 else 0
                            listpathogensstrains.append(f"{x['strainname']} ({percentreads}%)")
                            # listpathogensstrains.append(f"{x['strainname']} ({taxx}: {x['numreads']} reads - {percentreads}%)")
                if callfamclass == "":
                    callfamclass = f"{len(listpathogensstrains) if len(listpathogensstrains) > 0 else ''}"
                    if len(listpathogensstrains) > 0:
                        callfamclass = f"{', '.join(listpathogensstrains)}"
                if (is_pathogen == "N/A" or is_pathogen == "Unknown") and len(listpathogensstrains) > 0:
                    is_pathogen = "Potential"
                    annClass = "Derived"

            print(f"Reference: {ref} - {formatname}")
            print(f"\tIsPathogen: {is_pathogen}")
            print(f"\tCallClass: {callfamclass}")
            print(f"\tPathogenic Strains: {listpathogensstrains}")
            percentreads = f"{pathogenic_reads/total_reads:.1f}" if total_reads > 0 and pathogenic_reads>0 else 0
            print(f"\tPathogenic Reads: {pathogenic_reads} - {percentreads}%")
            print(f"\tAligned Strains: {fullstrains}")
            print(f"\tTotal reads: {sum(count['numreads'])}")
            print(f"\tAlignment Conf: {count['meangini']}")
            print()

            meanbaseq = format_non_zero_decimals(count['meanbaseq'])
            gini_coefficient = format_non_zero_decimals(count['meangini'])
            meanmapq = format_non_zero_decimals(count['meanmapq'])
            meancoverage = format_non_zero_decimals(count['meancoverage'])
            meandepth = format_non_zero_decimals(count['meandepth'])
            # if total_reads_aligned is 0 then set percent_aligned to 0
            countreads = sum(count['numreads'])
            if total_reads_aligned == 0:
                percent_aligned = 0
            else:
                percent_aligned = format_non_zero_decimals(100*countreads / total_reads_aligned)
            if total_reads == 0:
                percent_total = 0
            else:
                percent_total = format_non_zero_decimals(100*countreads / total_reads)
            if len(pathogenic_sites) == 0:
                pathogenic_sites = ""
            # Assuming 'count' is a simple value; if it's a dictionary or complex structure, adjust accordingly.
            file.write(f"{formatname}\t{sample_name}\t{sample_type}\t{percent_aligned}\t{percent_total}\t{countreads}\t{is_annotated}\t{pathogenic_sites}\t{is_pathogen}\t{ref}\t{status}\t{gini_coefficient}\t{meanbaseq}\t{meanmapq}\t{meancoverage}\t{meandepth}\t{annClass}\t{isSpecies}\t{callfamclass}\n")

if __name__ == "__main__":
    sys.exit(main())


