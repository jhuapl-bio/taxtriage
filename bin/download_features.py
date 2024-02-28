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
        help="Assembly refseq file to download from, with taxid and accession in the header",
    )
    parser.add_argument(
        "-d",
        "--decompress",
        action='store_true',
        default=False,
        help="Assembly refseq file to download from, with taxid and accession in the header",
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

def download_feature_file(dir_out: str, url: str,  outfile: str):
    """Download the assembly file from NCBI"""
    try:
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
            # decompress the outpath to dir_out, write to file
            outpathdecompress = outpath.replace(".gz", "")
            with gzip.open(outpath, 'rt') as f:
                #write contents to outpath
                contents = f.readlines()
                # remove the gz extension
                with open(outpathdecompress, 'w') as r:
                    # filter lines with CDS in beginning of line only
                    for line in contents:
                        if line.startswith("CDS"):
                            r.write(line)
                r.close()
                # remove the gz file
                os.remove(outpath)
            f.close()
    except Exception as e:
        logger.error(f"Failed to download assembly file from {url} with error {e}")
        return None
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
    with open(assembly_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            line = line.strip().split("\t")
            gcf = line[0]
            url = line[19]
            assembly_urls[gcf] = url
    f.close()
    for gcf in gcfs:
        if gcf in assembly_urls:
            try:
                url = assembly_urls[gcf]
                baseurl = os.path.basename(url)
                url = f"{url}/{baseurl}_feature_table.txt.gz"
                download_feature_file(dir_out, url, f"{gcf}_feature_table.txt.gz")
                print(f"Downloaded feature file for {gcf} at {url} to {dir_out}/{gcf}_feature_table.txt")
            except Exception as e:
                print(f"Failed to write feature file for {gcf} at {url} with error {e}")

if __name__ == "__main__":
    sys.exit(main())


