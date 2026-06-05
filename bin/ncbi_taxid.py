#!/usr/bin/env python3
##############################################################################################
# Copyright 2024 The Johns Hopkins University Applied Physics Laboratory LLC
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
"""Backup lookup of a NCBI taxid for a nuccore accession (e.g. OR833055.1).

Used when an accession supplied as a local reference FASTA (or pulled by a
classifier) cannot be matched to an entry in the RefSeq / GenBank assembly
summary files. In that case there is no GCF/GCA -> taxid mapping available, so
we query NCBI E-utilities directly by accession to recover the taxid and pass
it downstream in the map files.

The helper only uses the Python standard library so it works in any of the
pipeline containers regardless of whether Biopython is installed.
"""

import json
import ssl
import sys
import time
import urllib.parse
import urllib.request

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

# NCBI hosts a valid certificate, but several pipeline containers ship outdated
# CA bundles; fall back to an unverified context so the backup lookup still
# succeeds rather than hard-failing the whole run.
_CTX = ssl.create_default_context()
_CTX.check_hostname = False
_CTX.verify_mode = ssl.CERT_NONE


def _get(url, timeout=30):
    req = urllib.request.Request(url, headers={"User-Agent": "taxtriage-ncbi-backup"})
    with urllib.request.urlopen(req, context=_CTX, timeout=timeout) as resp:
        return resp.read().decode("utf-8", errors="replace")


def fetch_taxid(accession, email=None, api_key=None, retries=3, timeout=30, pause=0.4):
    """Return the taxid (as a string) for a nuccore accession, or None.

    Resolves the accession to an internal UID with esearch, then reads the
    taxid from esummary. Network/parse failures are swallowed and return None
    so callers can degrade gracefully (the taxid is simply left blank).
    """
    accession = (accession or "").strip()
    if not accession:
        return None

    common = {"db": "nuccore"}
    if email:
        common["email"] = email
    if api_key:
        common["api_key"] = api_key

    last_err = None
    for attempt in range(retries):
        try:
            # 1) accession -> internal UID
            search_params = dict(common)
            search_params.update({"term": accession, "retmode": "json"})
            search_url = EUTILS_BASE + "esearch.fcgi?" + urllib.parse.urlencode(search_params)
            search = json.loads(_get(search_url, timeout=timeout))
            idlist = search.get("esearchresult", {}).get("idlist", [])
            if not idlist:
                return None
            uid = idlist[0]

            time.sleep(pause)

            # 2) UID -> taxid via document summary
            summ_params = dict(common)
            summ_params.update({"id": uid, "retmode": "json"})
            summ_url = EUTILS_BASE + "esummary.fcgi?" + urllib.parse.urlencode(summ_params)
            summary = json.loads(_get(summ_url, timeout=timeout))
            record = summary.get("result", {}).get(uid, {})
            taxid = record.get("taxid")
            if taxid in (None, "", 0, "0"):
                return None
            return str(taxid)
        except Exception as err:  # noqa: BLE001 - degrade gracefully
            last_err = err
            time.sleep(pause * (attempt + 1))
    if last_err is not None:
        sys.stderr.write(
            "NCBI taxid backup lookup failed for {}: {}\n".format(accession, last_err)
        )
    return None


if __name__ == "__main__":
    # Simple CLI: ncbi_taxid.py ACCESSION [EMAIL]
    if len(sys.argv) < 2:
        sys.stderr.write("Usage: ncbi_taxid.py ACCESSION [EMAIL]\n")
        sys.exit(2)
    acc = sys.argv[1]
    mail = sys.argv[2] if len(sys.argv) > 2 else None
    result = fetch_taxid(acc, email=mail)
    if result:
        print(result)
    else:
        sys.stderr.write("No taxid found for {}\n".format(acc))
        sys.exit(1)
