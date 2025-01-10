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

"""Provide a command line tool to validate and transform tabular samplesheets."""


import argparse
import csv
import logging
import sys
from collections import Counter
from pathlib import Path
import re
import os
import numpy as np
from distributions import import_distributions, body_site_map
import pandas as pd
logger = logging.getLogger()


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and transform a tabular samplesheet.",
        epilog="Example: python check_samplesheet.py samplesheet.csv samplesheet.valid.csv",
    )
    parser.add_argument(
        "-i",
        "--file_in",
        metavar="FILE_IN",
        type=Path,
        help="Tabular input Kraken2 report",
    )
    parser.add_argument(
        "-f",
        "--filter_ranks",
        metavar="TOP_HITS_STRING",
        type=str, nargs="+", default=['S5', 'S4', 'S3', 'S2', 'S1', 'S', 'C', 'O', 'F', 'G', 'P', 'K', 'D', 'U'],
        help="Filter on only showing specific ranks ",
    )
    parser.add_argument(
        "--remove_commensals", default=False,  help="Remove any and all commensals that are listed in the pathogen sheet",  action='store_true'
    )
    parser.add_argument(
        "-s",
        "--top_hits_string",
        metavar="TOP_HITS_STRINGT",
        type=str, nargs="+", default=[],
        help="Top Hits ",
    )
    parser.add_argument(
        "-d",
        "--distributions",
        metavar="DISTRIBUTIONS",
        type=Path,
        default = None,
        help="Path to OPTIONAL file that contains distributions of zscores for each taxid in the input file.",
    )
    parser.add_argument(
        "-p",
        "--pathogens",
        metavar="PATHOGENS",
        type=Path,
        default = None,
        help="Path to OPTIONAL file that contains pathogens that have been annotated. Match taxids only",
    )
    parser.add_argument(
        "-z",
        "--zscore",
        metavar="ZSCORE",
        type=float,
        default = 1.5,
        help="OPTIONAL: If providing a distributions file, max zscore to consider for top hits",
    )
    parser.add_argument(
        "-o",
        "--file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Name of the output tsv file containing a mapping of top n organisms at individual taxa levels",
    )
    parser.add_argument(
        "-t",
        "--top_per_rank",
        default=5,
        metavar="TOP_PER_RANK",
        type=int,
        help="Max Top per rank code",
    )
    parser.add_argument(
        "-b",
        "--body_site",
        default="Unknown",
        metavar="BODYSITE",
        type=str,
        help="What body site to loook at distributions for. If empty and listed distributions enabled, this is assumed as Unknown",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def import_file(input, filter_ranks):
    tsv_file = open(input, newline='')
    read_tsv = csv.reader(tsv_file, delimiter="\t")
    mapping = []

    taxids = dict()
    header = ['abundance', 'clade_fragments_covered',
        'number_fragments_assigned', 'rank', 'taxid', 'name', 'parents']
    for row in read_tsv:
        entry = dict()
        for x in range(0, len(header)):

            if (header[x] != 'name' and header[x] != 'rank' and header[x] != 'taxid' and x < len(row)):
                entry[header[x]] = float(row[x])
            elif x < len(row):
                if header[x] == 'taxid':
                    entry[header[x]] = int(row[x])
                else:
                    entry[header[x]] = row[x]
            else:
                entry[header[x]] = ""
        mapping.append(entry)

        taxids[entry['taxid']] = entry['name']
    # k2_regex = re.compile(r"^\s{0,2}(\d{1,3}\.\d{1,2})\t(\d+)\t(\d+)\t([\dUDKRPCOFGS-]{1,3})\t(\d+)(\s+)(.+)")
    k2_regex = re.compile(r"^(\s+)(.+)")
    depth = dict()
    lastparents = dict()
    for i, l in enumerate(mapping):
        match = k2_regex.search(l['name'])

        if match:
            depth[l['taxid']] = int(len(match.group(1))/2)

            l['name'] = match.group(2)
        else:
            depth[l['taxid']] = 0
        parents = []

        for i in range(depth[l['taxid']]-1, 0, -1):
            parents.append(int(lastparents[i]))
        lastdepth = depth[l['taxid']]
        lastparents[depth[l['taxid']]] = l['taxid']

        l['depth'] = depth[l['taxid']]
        l['parents'] = parents
    if len(filter_ranks) > 0:
        mapping = [m for m in mapping if m['rank'] in filter_ranks]
        return mapping
    else:
        return mapping


def top_hit(mapping, specific_limits, top_per_rank, dist_orgs = []):
    uniq_ranks = list([x['rank'] for x in mapping])
    sorted_mapping = dict()
    for x in mapping:
        if not x['rank'] in sorted_mapping:
            sorted_mapping[x['rank']] = []
        sorted_mapping[x['rank']].append(x)
    i = 0
    for rank in uniq_ranks:
        sorted_specific_rank = sorted(
            sorted_mapping[rank], key=lambda d: d['abundance'], reverse=True)
        sorted_mapping[rank] = sorted_specific_rank
    newdata_seen = dict()
    countranks = dict()
    newdata = dict()
    for org in dist_orgs:
        if not org in newdata:
            # find the index and get value of the org in mapping
            row = None
            for row in mapping:
                rank = row['rank']
                if row['taxid'] == org:
                    newdata[row['taxid']] = row
                    break
                if not rank in countranks:
                    countranks[rank] = 1
                else:
                    countranks[rank] += 1
    for specifics, value in specific_limits.items():
        rank = value['rank']
        limit = value['limit']
        if rank in sorted_mapping:
            i = 0
            for row in sorted_mapping[rank]:
                if specifics in row['parents']:
                    if i >= limit:
                        break
                    else:
                        newdata_seen[row['taxid']] = True
                        if not rank in countranks:
                            countranks[rank] = 1
                        else:
                            countranks[rank] += 1
                        i += 1
                        if not row['taxid'] in newdata:
                            newdata[row['taxid']] = row
    for key, value in sorted_mapping.items():
        for row in value:
            if not key in countranks:
                countranks[key] = 0
            if countranks[key] >= top_per_rank:
                break
            else:
                if not row['taxid'] in newdata:
                    newdata[row['taxid']] = row
                    countranks[key] += 1
    return newdata


def make_files(mapping, outpath):
    header = ['abundance', 'clade_fragments_covered',
        'number_fragments_assigned', 'rank', 'taxid', 'name']
    path = str(outpath)
    path = open(path, "w")
    writer = csv.writer(path, delimiter='\t')
    writer.writerow(header)
    for taxid, row in mapping.items():
        out = []
        for head_item in header:
            if head_item == 'name':
                out.append(row[head_item].strip())
            elif head_item == 'parents':
                parentt = ";".join([str(x) for x in row[head_item]])
                out.append(parentt)
            else:
                out.append(row[head_item])
        writer.writerow(out)
    path.close()
    print("Done exporting the top hits report ")

def get_stats(row, stats_dict):
    taxid = row['tax_id']
    if (taxid, row['body_site']) in stats_dict:
        return stats_dict[(taxid, row['body_site'])]
    else:
        return None
def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level,
                        format="[%(levelname)s] %(message)s")
    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    specific_limits = dict()
    if args.top_hits_string:
        for x in args.top_hits_string:
            fulllist = x.split(":")
            if len(fulllist) >= 3:
                specific_limits[int(fulllist[0])] = dict(
                    limit=int(fulllist[1]), rank=fulllist[2])


    mapping = import_file(args.file_in, args.filter_ranks)
    # find mapping where taxid == 335341
    seentaxids = dict()
    for val in mapping:
        seentaxids[val['taxid']] = True
    # rename df_full['taxid'] to 'tax_id'
    dist_orgs = None
    pathogen_orgs = []
    extra_orgs = []
    if args.body_site == "Unknown" or not args.body_site:
        body_sites = []
    else:
        body_sites = [body_site_map(args.body_site.lower())]
    remove_taxids = []
    if args.pathogens:
        try:
            with open(args.pathogens, 'r', encoding='utf-8', errors='replace') as f:
                pathogen_sheet = pd.read_csv(f, sep=',')
            f.close()
            # A function to apply the mapping and remove duplicates
            def translate_and_deduplicate_sites(sites):
                # Split the string by comma and space, and remove empty strings if any
                sites_list = [site.strip() for site in sites.split(',') if site.strip()]
                # Initialize an empty set for the translated sites
                translated_sites = set()
                for site in sites_list:
                    # Get the mapped value from the dictionary
                    mapped_value = body_site_map(site.lower())
                    # If the mapped value is a list, add all its items to the set
                    if isinstance(mapped_value, list):
                        translated_sites.update(mapped_value)
                    else:
                        translated_sites.add(mapped_value)
                # Join the unique sites back into a string
                return ', '.join(sorted(translated_sites))
            # convert all pathogen_sites nan to "Unknown"
            pathogen_sheet['pathogenic_sites'].fillna("Unknown", inplace=True)
            pathogen_sheet['pathogenic_sites'] = pathogen_sheet['pathogenic_sites'].apply(translate_and_deduplicate_sites)
            # check if (lowercase) args.body_site is anywhere in pathogen_orgs body_site column, if not filter
            if len(body_sites) > 0:
                pathogen_sheet = pathogen_sheet[pathogen_sheet['pathogenic_sites'].str.lower().isin(body_sites)]

            # filter out where general_classification is pathogen or opportunistic pathogen
            pathogen_orgs = pathogen_sheet[pathogen_sheet['general_classification'].isin(["primary", "opportunistic", "potential", "oportunistic"])]['taxid']

            # remove all Nan values
            pathogen_orgs = pathogen_orgs.dropna()
            pathogen_orgs = pathogen_orgs.astype(int).tolist()

            # do the same for commensals and sites
            pathogen_sheet['commensal_sites'].fillna("", inplace=True)
            pathogen_sheet['commensal_sites'] = pathogen_sheet['commensal_sites'].apply(translate_and_deduplicate_sites)
            pathogen_sheet['commensal_sites'] = pathogen_sheet['commensal_sites'].str.lower()


            # Check if `remove_commensals` is True, and if so, remove taxids from the mapping
            if args.remove_commensals:
                # get all taxids where general classification is commensal
                commensal_orgs = pathogen_sheet[pathogen_sheet['general_classification'].isin(["commensal"])]['taxid']
                remove_taxids.extend(commensal_orgs)
            # even if general classification is commensal and it is a pathogen, add it back in since it is a pathogen for that body site likely
            for orgn in pathogen_orgs:
                if orgn in seentaxids:
                    extra_orgs.append(orgn)
        except Exception as e:
            print(e)
            print("Error reading pathogens file")



    if args.distributions:
        dists, site_counts = import_distributions(
            args.distributions,
            "tax_id",
            body_sites
        )
        df_full = pd.DataFrame(mapping)
        # only get dists where args.body_site is in body_site column
        # only keep rank has S in it
        df_full = df_full[df_full['rank'].str.contains('S')]
        df_full['body_site'] = body_sites[0] if len(body_sites) >= 1 else None
        df_full.rename(columns={'taxid': 'tax_id'}, inplace=True)
        # convert body_site to lowercase
        df_full['body_site'] = df_full['body_site'].str.lower()
        # only get rows where args.body_site is in body_site column
        # filter out 2d dict of taxid, body_site of df_full where taxid AND body_site is in dists
        datanew = []
        for index, row in df_full.iterrows():
            taxidsonly = [key[0] for key in dists.keys()]
            if (row['tax_id'], row['body_site']) in dists or row['tax_id'] not in taxidsonly or len(body_sites)== 0:
                datanew.append(row)
        df_full = pd.DataFrame(datanew)
        if len(df_full) == 0:
            print("No data to process")
        else:
            df_full['stats'] = df_full.apply(lambda x: get_stats(x, dists), axis=1)
            # get the mean, std, and zscore for each row
            df_full['mean'] = df_full['stats'].apply(lambda x: x['mean'] if x else None)
            df_full['norm_abundance'] = df_full['stats'].apply(lambda x: x['norm_abundance'] if x else None)
            df_full['std'] = df_full['stats'].apply(lambda x: x['std'] if x else None)
            ## calculate the stddev by considering both the "abundances" list and including empty values or 0 for missing length between site_counts and abundances
            df_full['zscore'] = (df_full['abundance'] - df_full['norm_abundance']) / df_full['std']

            # if zscore is NaN convert to -1
            df_full['zscore'] = df_full['zscore'].fillna(-1)
            # filter any row with zscore outside of absolute value of args.zscore
            df_full = df_full[  ( df_full['zscore'] > args.zscore)   ]
            # get all where zscore is less than or equal to args.zscore and is not -1
            if args.remove_commensals:
                taxa_non = df_full[  ( df_full['zscore'] <= args.zscore)   ]['tax_id'].tolist()
                remove_taxids.extend(taxa_non)

            # Removed -1 from distributions - need to figure out a better way as too many things are downloaded currently

            # get taxid is 2 stool body_site
            # Get percentile of abundance relative to all abundances
            # df_full['percentile'] = df_full['abundance'].rank(pct=True) * 100
            dist_orgs = df_full['tax_id'].tolist()
            extra_orgs = extra_orgs + dist_orgs
    extra_orgs = list(set(extra_orgs))
    # remove all taxids from mapping that are in remove_taxids
    mapping_taxids = [x['taxid'] for x in mapping]
    mark_for_removal = []
    mark_for_removal_extra = []
    if len(remove_taxids) > 0:
        for taxid in remove_taxids:
            # find the index and get value of the org in mapping
            if taxid in mapping_taxids:
                idx = mapping_taxids.index(taxid)
                mark_for_removal.append(idx)
                print(f"removing {taxid} commensal organisms from the mapping")
            if taxid in extra_orgs:
                idx_extra = extra_orgs.index(taxid)
                mark_for_removal_extra.append(idx_extra)
    # Remove from mapping_taxids
    for idx in sorted(mark_for_removal, reverse=True):
        del mapping_taxids[idx]
    # Remove from extra_orgs
    for idx_extra in sorted(mark_for_removal_extra, reverse=True):
        del extra_orgs[idx_extra]
    mapping = top_hit(mapping, specific_limits, args.top_per_rank, extra_orgs)
    make_files(mapping, args.file_out)


if __name__ == "__main__":
    sys.exit(main())
