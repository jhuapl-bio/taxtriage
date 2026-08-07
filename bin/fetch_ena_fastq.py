#!/usr/bin/env python

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

"""
Download ENA-hosted FASTQ files with retries and md5 verification.

Written for the SRA_FETCH_ENA module so the download does not depend on curl
or wget being present in the container.  Files are written to a temporary
name and only renamed into place once the md5 matches, so a killed task can
never leave a truncated FASTQ behind in a storeDir cache.
"""

import argparse
import hashlib
import logging
import os
import shutil
import sys
import time
import urllib.error
import urllib.request

logger = logging.getLogger()

CHUNK = 1024 * 1024

# ENA's file host serves /vol1 over plain HTTP and FTP but NOT HTTPS — an https://
# request to ftp.sra.ebi.ac.uk returns 403 Forbidden. Rather than hard-coding one
# scheme and breaking if that ever changes, every candidate scheme is tried in
# turn, starting with whatever the manifest supplied.
SCHEME_FALLBACKS = ("http", "https", "ftp")


def md5sum(path):
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(CHUNK), b""):
            digest.update(chunk)
    return digest.hexdigest()


def candidate_urls(url):
    """
    Return *url* followed by the same location under the other schemes.

    Order matters: the supplied URL is tried first, so a manifest that already
    names a working scheme costs nothing extra.
    """
    parsed = urllib.parse.urlsplit(url)
    if not parsed.scheme or not parsed.netloc:
        return [url]

    rest = urllib.parse.urlunsplit(("", "") + parsed[2:])
    candidates = [url]
    for scheme in SCHEME_FALLBACKS:
        if scheme == parsed.scheme:
            continue
        candidates.append("{}://{}{}".format(scheme, parsed.netloc, rest))
    return candidates


def _fetch_once(url, dest, expected_md5, timeout):
    """
    Download *url* to *dest*, verifying md5 before the file is put in place.

    Writes to a `.partial` file and only renames on success, so an interrupted
    or corrupt download can never leave a truncated FASTQ in the storeDir cache
    where a later run would mistake it for a complete one.
    """
    tmp = dest + ".partial"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "taxtriage-sra-fetch"})
        with urllib.request.urlopen(req, timeout=timeout) as response, open(tmp, "wb") as out:
            shutil.copyfileobj(response, out, CHUNK)

        if os.path.getsize(tmp) == 0:
            raise IOError("downloaded file is empty")

        if expected_md5:
            actual = md5sum(tmp)
            if actual != expected_md5:
                raise IOError(
                    "md5 mismatch for {}: expected {}, got {}".format(dest, expected_md5, actual)
                )
            logger.info("md5 verified for %s", dest)

        os.replace(tmp, dest)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def download(url, dest, expected_md5=None, retries=4, backoff=5.0, timeout=300):
    candidates = candidate_urls(url)
    last = None
    for attempt in range(1, retries + 1):
        for candidate in candidates:
            try:
                logger.info("Downloading %s -> %s (attempt %s/%s)", candidate, dest, attempt, retries)
                _fetch_once(candidate, dest, expected_md5, timeout)
                logger.info("Downloaded %s", dest)
                return
            except (urllib.error.URLError, urllib.error.HTTPError, IOError, OSError) as err:
                last = err
                logger.warning("Failed %s: %s", candidate, err)
        if attempt < retries:
            time.sleep(backoff * attempt)

    raise SystemExit(
        "ERROR: failed to download {} after {} attempts across {}: {}".format(
            dest, retries, ", ".join(candidates), last
        )
    )


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Download ENA FASTQ files with md5 verification.")
    parser.add_argument("--url", required=True, help="URL of read 1 (or the only read file).")
    parser.add_argument("--out", required=True, help="Output filename for read 1.")
    parser.add_argument("--md5", default=None, help="Expected md5 of read 1.")
    parser.add_argument("--url2", default=None, help="URL of read 2, for paired runs.")
    parser.add_argument("--out2", default=None, help="Output filename for read 2.")
    parser.add_argument("--md52", default=None, help="Expected md5 of read 2.")
    parser.add_argument("--retries", type=int, default=4)
    parser.add_argument("-l", "--log-level", default="INFO",
                        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"))
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")

    download(args.url, args.out, expected_md5=args.md5 or None, retries=args.retries)
    if args.url2 and args.out2:
        download(args.url2, args.out2, expected_md5=args.md52 or None, retries=args.retries)
    return 0


if __name__ == "__main__":
    sys.exit(main())
