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
from xmlrpc.client import Boolean
import sys
import os
import gzip
import argparse
import logging
from pathlib import Path
import urllib.request as request
import ssl

ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE


# import requests
logger = logging.getLogger()


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and transform a tabular samplesheet.",
        epilog="Example: python check_samplesheet.py samplesheet.csv samplesheet.valid.csv",
    )
    parser.add_argument(
        "-i",
        "--input",
        metavar="INPUT",
        help="List of taxIDs or names",
    )
    parser.add_argument(
        "-a",
        "--assembly_file",
        type=Path,
        nargs="+",
        help="Assembly refseq file to download from, with taxid and accession in the header",
    )
    parser.add_argument(
        "-m",
        "--mapfile",
        type=Path,
        help="Optional mapping file to map accessions to taxids",
    )
    parser.add_argument(
        "-d",
        "--decompress",
        action='store_true',
        default=False,
        help="Assembly refseq file to download from, with taxid and accession in the header",
    )
    parser.add_argument(
        "-f",
        "--force",
        action='store_true',
        default=False,
        help="Force/overwrite download even if it exists",
    )
    parser.add_argument(
        "-p",
        "--get_aas",
        action='store_true',
        default=False,
        help="Download the corresponding amino acid sequences for the CDS features",
    )
    parser.add_argument(
        "-b",
        "--bed_out",
        type=Path,
        default=None,
        help="--mask_products: also write the product regions of every assembly as one BED "
             "(0-based, half-open). Taken from the assembly's GFF3 (one CDS row per exon, so "
             "introns are left out); if the GFF cannot be fetched, from its feature table (CDS only)",
    )
    parser.add_argument(
        "-t",
        "--bed_types",
        default="CDS",
        help="Feature types written to --bed_out (GFF column 3), comma-separated. Default CDS",
    )
    parser.add_argument(
        "-o",
        "--dir_out",
        metavar="DIR_OUT",
        type=Path,
        help="Name of output directory for feature table files",
    )

    return parser.parse_args(argv)

import shutil
from contextlib import closing
from mimetypes import guess_type
from functools import partial

def download_feature_file(dir_out: str, url: str,  outfile: str, featurefile: bool = True):
    """Download the assembly file from NCBI"""
    try:
        accessions = set()
        response = request.urlopen(url, context=ctx)
        encoding = guess_type(url)[1]   # uses file extension
        _open = partial(
            gzip.open, mode='rt') if encoding == 'gzip' else open
        outpath = os.path.join(dir_out, outfile)
        with closing(request.urlopen(url, context=ctx)) as r:
            with open(outpath, 'wb') as f:
                shutil.copyfileobj(r, f)
            f.close()
        r.close()

        if encoding == 'gzip':
            print("Encoding is gzip")
            # decompress the outpath to dir_out, write to file
            outpathdecompress = outpath.replace(".gz", "")
            with gzip.open(outpath, 'rt') as f:
                #write contents to outpath
                contents = f.readlines()
                # remove the gz extension
                with open(outpathdecompress, 'w') as r:
                    # filter lines with CDS in beginning of line only
                    for line in contents:
                        if line.startswith("CDS") and featurefile:
                            r.write(line)
                            parseline = line.split("\t")
                            accession = parseline[10]
                            if accession:
                                #add to set
                                accessions.add(accession)
                        elif not featurefile:
                            r.write(line)
                            # get the header if starts with ">"
                            if line.startswith(">"):
                                # get the accession
                                accession = line.split(" ")[0].replace(">", "")
                                if accession:
                                    accessions.add(accession)
                r.close()
                # remove the gz file
                os.remove(outpath)
            f.close()
        print(f"DONE: Downloaded feature file")
        return accessions
    except Exception as e:
        logger.error(f"Failed to download assembly file from {url} with error {e}")
        return None
def _gff_products(path, types):
    """(seqid, start0, end, label) for the GFF3 rows whose type is in `types`."""
    rows = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in types:
                continue
            attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            label = attrs.get("gene") or attrs.get("Name") or attrs.get("locus_tag") or attrs.get("ID") or f[2]
            rows.append((f[0], int(f[3]) - 1, int(f[4]), f"product:{f[2]}:{label}"[:120].replace(" ", "_")))
    return rows


def _feature_table_products(path, types):
    """Fallback: CDS rows of a downloaded feature table (no header; spliced CDS = one span)."""
    rows = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[0] not in types:
                continue
            try:
                a, b = int(f[7]) - 1, int(f[8])
            except ValueError:
                continue
            label = (f[14] if len(f) > 14 and f[14] else (f[13] if len(f) > 13 else f[0])) or f[0]
            rows.append((f[6], a, b, f"product:{f[0]}:{label}"[:120].replace(" ", "_")))
    return rows


def write_products_bed(gcfs, assembly_urls, dir_out, bed_out, types):
    """One BED of product regions for all assemblies; a log line per assembly."""
    rows = []
    for gcf in gcfs:
        if gcf not in assembly_urls:
            print(f"[products] {gcf}: not in the assembly summary - no products")
            continue
        gcf_url = assembly_urls[gcf].rstrip("/")
        base = os.path.basename(gcf_url)
        gff = os.path.join(dir_out, f"{gcf}_genomic.gff.gz")
        got = []
        try:
            if not os.path.exists(gff) or os.path.getsize(gff) == 0:
                with closing(request.urlopen(f"{gcf_url}/{base}_genomic.gff.gz", context=ctx)) as r, open(gff, "wb") as w:
                    shutil.copyfileobj(r, w)
            got = _gff_products(gff, types)
            src = "GFF3"
        except Exception as e:
            if os.path.exists(gff):
                os.remove(gff)  # never reuse a partial download
            print(f"[products] {gcf}: GFF3 unavailable ({e}); using the feature table")
            ft = os.path.join(dir_out, f"{gcf}_feature_table.txt")
            got = _feature_table_products(ft, types) if os.path.exists(ft) else []
            src = "feature table"
        print(f"[products] {gcf}: {len(got)} {','.join(sorted(types))} regions from the {src}")
        rows.extend(got)
    rows.sort()
    with open(bed_out, "w") as fh:
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    print(f"[products] {len(rows)} regions -> {bed_out}")


def main():
    args = parse_args()
    logger.setLevel(logging.INFO)
    inputfile = args.input
    assembly_file = args.assembly_file
    dir_out = args.dir_out
    if not dir_out.exists():
        dir_out.mkdir()
    with open(inputfile, 'r') as f:
        # read in lines, remove newline and get rid of empty lines
        gcfs = [line.strip() for line in f if line.strip()]
    f.close()
    assembly_urls = dict()
    maptaixd = dict()
    for assembly_path in assembly_file:
        with open(assembly_path, 'r') as f:
            for line in f:
                if line.startswith("#"):
                    continue
                # Strip only the line ending, then split strictly on tab so that
                # fields containing multiple spaces (e.g. asm_submitter) are not
                # broken across columns.
                line = line.rstrip("\r\n").split("\t")
                gcf = line[0]
                url = line[19]
                taxid = line[5]
                maptaixd[gcf] = taxid
                assembly_urls[gcf] = url

    if args.mapfile:
        mfa = open(args.mapfile, 'w')
        mfa.write(f"accession.version\ttaxid\tassembly\n")
    else:
        mfa = None
    accessions_map = dict()
    for gcf in gcfs:
        if gcf in assembly_urls:
            accessions_map[gcf] = set()
            try:
                gcf_url = assembly_urls[gcf]
                # strip "/" if at end of the gcf_url
                if gcf_url.endswith("/"):
                    gcf_url = gcf_url.rstrip("/")
                baseurl = os.path.basename(gcf_url)
                url = f"{gcf_url}/{baseurl}_feature_table.txt.gz"

                if args.force or not os.path.exists(f"{dir_out}/{gcf}_feature_table.txt") or os.path.getsize(f"{dir_out}/{gcf}_feature_table.txt") == 0:
                    print(f"Downloading feature file for {gcf} at {url}")
                    download_feature_file(dir_out, url, f"{gcf}_feature_table.txt.gz", featurefile=True)
                else:
                    print(f"Feature file for {gcf} at {url} already exists, skipping")
                    # read in the file and get the accessions
                if args.get_aas:
                    url_aa = f"{gcf_url}/{baseurl}_protein.faa.gz"
                    print(f"Downloading amino acid sequence file for {gcf} at {url_aa}")
                    # Check if file exists, if so then skip and file isnt empty
                    if args.force or not os.path.exists(f"{dir_out}/{gcf}_protein.faa") or os.path.getsize(f"{dir_out}/{gcf}_protein.faa") == 0:
                        accessions = download_feature_file(dir_out, url_aa, f"{gcf}_protein.faa.gz", featurefile=False)
                        accessions_map[gcf] = accessions
                    else:
                        print(f"Amino acid sequence file for {gcf} at {dir_out}/{gcf}_protein.faa already exists, skipping")
                        # read in only header lines, parse and add to accessions set
                        with open(f"{dir_out}/{gcf}_protein.faa", 'r') as f:
                            for line in f:
                                if line.startswith(">"):
                                    accession = line.split(" ")[0].replace(">", "")
                                    if accession:
                                        accessions_map[gcf].add(accession)
                            f.close()
            except Exception as e:
                print(f"Failed to write feature file for {gcf} at {url} with error {e}")
    if args.bed_out:
        types = {t.strip() for t in args.bed_types.split(",") if t.strip()}
        write_products_bed(gcfs, assembly_urls, dir_out, args.bed_out, types)
    if mfa:
        # get only accessions in the amino acid file that are not in the feature file accessions
        seen = dict()
        for gcf, accessions in accessions_map.items():
            taxid = maptaixd.get(gcf, None)
            if taxid and accessions:
                try:
                    for acc_n in accessions:
                        if acc_n not in seen:
                            mfa.write(f"{acc_n}\t{taxid}\t{gcf}\n")
                            seen[acc_n] = True
                except Exception as ex:
                    print(ex)
if __name__ == "__main__":
    sys.exit(main())


