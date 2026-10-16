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
import os
import random
import ssl
import sys
import threading
import time
import urllib.error
import urllib.parse
import urllib.request

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

# NCBI E-utilities allows 3 requests/second without an API key and 10/second
# with one. Stay comfortably under those ceilings - the published limit is
# enforced per-IP and bursts get answered with HTTP 429.
RATE_NO_KEY = float(os.environ.get("NCBI_RATE_NO_KEY", "2.0"))
RATE_WITH_KEY = float(os.environ.get("NCBI_RATE_WITH_KEY", "7.0"))

# Accessions resolved per esearch/esummary pair. Batching is the single biggest
# win: 200 accessions cost 2 requests instead of 400.
BATCH_SIZE = int(os.environ.get("NCBI_BATCH_SIZE", "100"))

# NCBI hosts a valid certificate, but several pipeline containers ship outdated
# CA bundles; fall back to an unverified context so the backup lookup still
# succeeds rather than hard-failing the whole run.
_CTX = ssl.create_default_context()
_CTX.check_hostname = False
_CTX.verify_mode = ssl.CERT_NONE


class _RateLimiter:
    """Process-wide minimum-interval throttle shared by every request."""

    def __init__(self, rate_per_sec):
        self._lock = threading.Lock()
        self._interval = 1.0 / rate_per_sec if rate_per_sec > 0 else 0.0
        self._next_allowed = 0.0
        # Extra delay imposed after a 429, decayed on success.
        self._penalty = 0.0

    def set_rate(self, rate_per_sec):
        with self._lock:
            self._interval = 1.0 / rate_per_sec if rate_per_sec > 0 else 0.0

    def wait(self):
        while True:
            with self._lock:
                now = time.monotonic()
                if now >= self._next_allowed:
                    self._next_allowed = now + self._interval + self._penalty
                    return
                sleep_for = self._next_allowed - now
            time.sleep(sleep_for)

    def penalize(self):
        """Back off globally after a 429 so sibling calls slow down too."""
        with self._lock:
            self._penalty = min(5.0, (self._penalty * 2) if self._penalty else 0.5)

    def relax(self):
        with self._lock:
            if self._penalty:
                self._penalty = max(0.0, self._penalty / 2.0)


_LIMITER = _RateLimiter(RATE_NO_KEY)


def configure_rate(api_key=None, rate_per_sec=None):
    """Set the shared request rate. Called once before a batch of lookups."""
    if rate_per_sec is None:
        rate_per_sec = RATE_WITH_KEY if api_key else RATE_NO_KEY
    _LIMITER.set_rate(rate_per_sec)


def _get(url, timeout=30, retries=6):
    """Rate-limited GET with exponential backoff that honours HTTP 429."""
    last_err = None
    for attempt in range(retries):
        _LIMITER.wait()
        try:
            req = urllib.request.Request(
                url, headers={"User-Agent": "taxtriage-ncbi-backup"}
            )
            with urllib.request.urlopen(req, context=_CTX, timeout=timeout) as resp:
                body = resp.read().decode("utf-8", errors="replace")
            _LIMITER.relax()
            return body
        except urllib.error.HTTPError as err:
            last_err = err
            if err.code in (429, 500, 502, 503, 504):
                if err.code == 429:
                    _LIMITER.penalize()
                retry_after = err.headers.get("Retry-After") if err.headers else None
                try:
                    delay = float(retry_after)
                except (TypeError, ValueError):
                    delay = (2.0 ** attempt) + random.uniform(0, 0.5)
                time.sleep(min(delay, 60.0))
                continue
            raise
        except Exception as err:  # noqa: BLE001 - transient network issues
            last_err = err
            time.sleep(min((2.0 ** attempt) + random.uniform(0, 0.5), 60.0))
    raise last_err


def _common_params(email, api_key):
    common = {"db": "nuccore"}
    if email:
        common["email"] = email
        common["tool"] = "taxtriage"
    if api_key:
        common["api_key"] = api_key
    return common


def _esearch_uids(accessions, common, timeout):
    """Resolve a list of accessions to UIDs in one esearch call."""
    term = " OR ".join("{}[ACCN]".format(a) for a in accessions)
    params = dict(common)
    params.update({"term": term, "retmode": "json", "retmax": str(len(accessions) * 2)})
    url = EUTILS_BASE + "esearch.fcgi"
    # Terms get long fast; POST avoids URL length limits.
    data = urllib.parse.urlencode(params).encode("utf-8")
    for attempt in range(6):
        _LIMITER.wait()
        try:
            req = urllib.request.Request(
                url, data=data, headers={"User-Agent": "taxtriage-ncbi-backup"}
            )
            with urllib.request.urlopen(req, context=_CTX, timeout=timeout) as resp:
                payload = json.loads(resp.read().decode("utf-8", errors="replace"))
            _LIMITER.relax()
            return payload.get("esearchresult", {}).get("idlist", [])
        except urllib.error.HTTPError as err:
            if err.code in (429, 500, 502, 503, 504):
                if err.code == 429:
                    _LIMITER.penalize()
                retry_after = err.headers.get("Retry-After") if err.headers else None
                try:
                    delay = float(retry_after)
                except (TypeError, ValueError):
                    delay = (2.0 ** attempt) + random.uniform(0, 0.5)
                time.sleep(min(delay, 60.0))
                continue
            raise
        except Exception:  # noqa: BLE001
            time.sleep(min((2.0 ** attempt) + random.uniform(0, 0.5), 60.0))
    return []


def _esummary_taxids(uids, common, timeout):
    """Return {accession(-version and bare): taxid} for a list of UIDs."""
    params = dict(common)
    params.update({"id": ",".join(uids), "retmode": "json"})
    url = EUTILS_BASE + "esummary.fcgi?" + urllib.parse.urlencode(params)
    summary = json.loads(_get(url, timeout=timeout))
    result = summary.get("result", {})
    out = {}
    for uid in result.get("uids", uids):
        record = result.get(uid, {})
        taxid = record.get("taxid")
        if taxid in (None, "", 0, "0"):
            continue
        taxid = str(taxid)
        for key in (record.get("accessionversion"), record.get("caption")):
            if key:
                out[str(key)] = taxid
                out[str(key).split(".")[0]] = taxid
    return out


def fetch_taxids(accessions, email=None, api_key=None, timeout=30, batch_size=None):
    """Return {accession: taxid or None} for many accessions, batched.

    Requests are throttled process-wide and retried with exponential backoff on
    HTTP 429, so large unmapped-accession lists no longer trip NCBI's rate
    limit. Failures degrade to None rather than aborting the run.
    """
    accessions = [str(a).strip() for a in accessions if a and str(a).strip()]
    accessions = list(dict.fromkeys(accessions))
    results = {a: None for a in accessions}
    if not accessions:
        return results

    api_key = api_key or os.environ.get("NCBI_API_KEY") or None
    configure_rate(api_key)
    common = _common_params(email, api_key)
    size = batch_size or BATCH_SIZE

    for start in range(0, len(accessions), size):
        chunk = accessions[start:start + size]
        try:
            uids = _esearch_uids(chunk, common, timeout)
            if not uids:
                continue
            lookup = _esummary_taxids(uids, common, timeout)
        except Exception as err:  # noqa: BLE001 - degrade gracefully
            sys.stderr.write(
                "NCBI taxid backup lookup failed for batch of {} starting at {}: {}\n".format(
                    len(chunk), chunk[0], err
                )
            )
            continue
        for acc in chunk:
            taxid = lookup.get(acc) or lookup.get(acc.split(".")[0])
            if taxid:
                results[acc] = taxid
    return results


def fetch_taxid(accession, email=None, api_key=None, retries=6, timeout=30, pause=None):
    """Return the taxid (as a string) for a single nuccore accession, or None.

    Thin wrapper over :func:`fetch_taxids` kept for backwards compatibility;
    `retries` and `pause` are accepted but handled by the shared limiter.
    """
    accession = (accession or "").strip()
    if not accession:
        return None
    return fetch_taxids([accession], email=email, api_key=api_key, timeout=timeout).get(
        accession
    )


if __name__ == "__main__":
    # Simple CLI: ncbi_taxid.py ACCESSION [ACCESSION ...] [--email EMAIL]
    argv = sys.argv[1:]
    mail = None
    if "--email" in argv:
        i = argv.index("--email")
        mail = argv[i + 1] if len(argv) > i + 1 else None
        del argv[i:i + 2]
    if not argv:
        sys.stderr.write("Usage: ncbi_taxid.py ACCESSION [ACCESSION ...] [--email EMAIL]\n")
        sys.exit(2)
    found = fetch_taxids(argv, email=mail)
    missing = [a for a, t in found.items() if not t]
    for acc in argv:
        if found.get(acc):
            print("{}\t{}".format(acc, found[acc]) if len(argv) > 1 else found[acc])
    for acc in missing:
        sys.stderr.write("No taxid found for {}\n".format(acc))
    sys.exit(1 if len(missing) == len(argv) else 0)
