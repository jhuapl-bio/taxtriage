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
Resolve NCBI SRA / ENA accessions into a concrete download manifest.

Accepts run (SRR/ERR/DRR), experiment (SRX/ERX/DRX), sample (SRS/ERS/DRS,
SAMN/SAMEA/SAMD) and study/project (SRP/ERP/DRP, PRJNA/PRJEB/PRJDB) accessions.
Non-run accessions are expanded into every child run, producing one pipeline
sample per run.

Resolution order
----------------
1. ENA Portal API (``filereport?result=read_run``).  Gives ready-made
   ``fastq.gz`` URLs plus ``library_layout`` and ``instrument_platform``, so
   paired-vs-single is known before anything is downloaded.
2. NCBI eutils ``runinfo`` (esearch + efetch on ``db=sra``).  Used when ENA
   returns nothing for the accession.  Yields the run list and layout but no
   FASTQ URLs, so those runs are marked ``source=sratools`` and are fetched
   with ``prefetch`` + ``fasterq-dump`` downstream.

Output: a CSV manifest, one row per run, consumed by ``input_check.nf``.
"""

import argparse
import csv
import io
import logging
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

logger = logging.getLogger()

ENA_PORTAL = "https://www.ebi.ac.uk/ena/portal/api/filereport"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

ENA_FIELDS = [
    "run_accession",
    "experiment_accession",
    "sample_accession",
    "study_accession",
    "instrument_platform",
    "instrument_model",
    "library_layout",
    "library_strategy",
    "library_source",
    "read_count",
    "fastq_ftp",
    "fastq_md5",
    "fastq_bytes",
]

MANIFEST_FIELDS = [
    "sample",           # pipeline sample name (unique per run)
    "request",          # the accession the user actually asked for
    "run_accession",
    "source",           # 'ena' (direct fastq URLs) or 'sratools' (prefetch/fasterq-dump)
    "platform",         # ILLUMINA / OXFORD / PACBIO / ...  (TaxTriage vocabulary)
    "single_end",       # 'true' / 'false'
    "fastq_1_url",
    "fastq_2_url",
    "fastq_1_md5",
    "fastq_2_md5",
    "instrument_platform",
    "instrument_model",
    "library_strategy",
    "read_count",
]

# ── Accession patterns ──────────────────────────────────────────────────────
RUN_RE = re.compile(r"^(SRR|ERR|DRR)\d{5,}$", re.IGNORECASE)
EXPERIMENT_RE = re.compile(r"^(SRX|ERX|DRX)\d{5,}$", re.IGNORECASE)
SAMPLE_RE = re.compile(r"^(SRS|ERS|DRS)\d{5,}$|^SAM(N|EA|EG|D)\d+$", re.IGNORECASE)
STUDY_RE = re.compile(r"^(SRP|ERP|DRP)\d{5,}$|^PRJ(NA|EB|EA|DA|DB)\d+$", re.IGNORECASE)

ACCESSION_RES = (RUN_RE, EXPERIMENT_RE, SAMPLE_RE, STUDY_RE)

# ENA instrument_platform / NCBI runinfo Platform → TaxTriage platform vocabulary.
# Anything unmapped is passed through uppercased; downstream minimap2 preset
# selection falls back to map-ont for unrecognised platforms.
PLATFORM_MAP = {
    "ILLUMINA": "ILLUMINA",
    "OXFORD_NANOPORE": "OXFORD",
    "OXFORD NANOPORE": "OXFORD",
    "NANOPORE": "OXFORD",
    "PACBIO_SMRT": "PACBIO",
    "PACBIO SMRT": "PACBIO",
    "PACBIO": "PACBIO",
    "ION_TORRENT": "ILLUMINA",  # short reads; use the short-read presets
    "ION TORRENT": "ILLUMINA",
    "BGISEQ": "ILLUMINA",
    "DNBSEQ": "ILLUMINA",
    "CAPILLARY": "ILLUMINA",
}


def is_accession(value):
    """Return True when *value* looks like an SRA/ENA accession rather than a path."""
    if not value:
        return False
    value = value.strip()
    # A path-ish string is never an accession
    if any(c in value for c in "/\\") or "." in value:
        return False
    return any(rx.match(value) for rx in ACCESSION_RES)


def normalise_platform(raw):
    if not raw:
        return ""
    key = raw.strip().upper()
    return PLATFORM_MAP.get(key, key.replace(" ", "_"))


def http_get(url, retries=4, backoff=3.0, timeout=120):
    """GET *url*, retrying transient failures with linear backoff."""
    last = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "taxtriage-sra-resolver"})
            with urllib.request.urlopen(req, timeout=timeout) as handle:
                return handle.read().decode("utf-8", errors="replace")
        except (urllib.error.URLError, urllib.error.HTTPError, OSError) as err:
            last = err
            if isinstance(err, urllib.error.HTTPError) and err.code in (400, 404):
                # Definitive "not here" — no point retrying
                logger.debug("%s returned HTTP %s", url, err.code)
                return ""
            logger.warning("Attempt %s/%s failed for %s: %s", attempt, retries, url, err)
            if attempt < retries:
                time.sleep(backoff * attempt)
    logger.warning("Giving up on %s: %s", url, last)
    return ""


# ── ENA ─────────────────────────────────────────────────────────────────────
def query_ena(accession):
    """Return a list of ENA read_run records for *accession* (possibly empty)."""
    params = urllib.parse.urlencode(
        {
            "accession": accession,
            "result": "read_run",
            "fields": ",".join(ENA_FIELDS),
            "format": "tsv",
            "limit": "0",
        }
    )
    body = http_get("{}?{}".format(ENA_PORTAL, params))
    if not body.strip():
        return []
    reader = csv.DictReader(io.StringIO(body), delimiter="\t")
    return [row for row in reader if row.get("run_accession")]


def pick_ena_fastqs(record):
    """
    Split an ENA ``fastq_ftp`` field into (read1, read2) URLs plus their md5s.

    ENA lists 1, 2 or 3 files.  Three files means the submitter provided a
    paired run plus an orphan/unpaired file: ``RUN.fastq.gz;RUN_1.fastq.gz;
    RUN_2.fastq.gz``.  In that case the ``_1``/``_2`` pair is used and the
    orphan file is dropped, which is what every downstream paired-end tool
    expects.
    """
    urls = [u.strip() for u in (record.get("fastq_ftp") or "").split(";") if u.strip()]
    md5s = [m.strip() for m in (record.get("fastq_md5") or "").split(";") if m.strip()]
    if len(md5s) != len(urls):
        md5s = [""] * len(urls)

    if not urls:
        return "", "", "", ""

    paired = [(u, m) for u, m in zip(urls, md5s) if re.search(r"_[12]\.f(ast)?q\.gz$", u)]
    if len(paired) >= 2:
        paired.sort(key=lambda pair: pair[0])
        (u1, m1), (u2, m2) = paired[0], paired[1]
        return _url(u1), _url(u2), m1, m2

    return _url(urls[0]), "", md5s[0], ""


def _url(url):
    """
    Add a scheme to the bare host paths ENA returns in ``fastq_ftp``.

    HTTP, not HTTPS: ftp.sra.ebi.ac.uk answers 403 Forbidden over TLS. The
    downloader falls back through https:// and ftp:// if this ever stops
    working, so this is a preference rather than a hard requirement.
    """
    if not url:
        return ""
    if url.startswith(("http://", "https://", "ftp://")):
        return url
    return "http://" + url.lstrip("/")


def ena_records_to_rows(accession, records):
    rows = []
    for rec in records:
        url1, url2, md5_1, md5_2 = pick_ena_fastqs(rec)
        layout = (rec.get("library_layout") or "").strip().upper()
        # Trust the file listing over the declared layout: a run declared PAIRED
        # that only ships one FASTQ is single-end in practice.
        single_end = not bool(url2) if url1 else layout != "PAIRED"
        rows.append(
            {
                "request": accession,
                "run_accession": rec["run_accession"],
                "source": "ena" if url1 else "sratools",
                "platform": normalise_platform(rec.get("instrument_platform")),
                "single_end": "true" if single_end else "false",
                "fastq_1_url": url1,
                "fastq_2_url": url2,
                "fastq_1_md5": md5_1,
                "fastq_2_md5": md5_2,
                "instrument_platform": (rec.get("instrument_platform") or "").strip(),
                "instrument_model": (rec.get("instrument_model") or "").strip(),
                "library_strategy": (rec.get("library_strategy") or "").strip(),
                "read_count": (rec.get("read_count") or "").strip(),
            }
        )
    return rows


# ── NCBI eutils fallback ────────────────────────────────────────────────────
def query_ncbi_runinfo(accession, api_key=None):
    """
    Expand *accession* to runs using NCBI eutils, returning runinfo CSV records.

    Used only when ENA has no record of the accession.  No FASTQ URLs are
    available on this path, so every run is marked ``source=sratools``.
    """
    search_params = {
        "db": "sra",
        "term": accession,
        "retmax": "10000",
        "usehistory": "y",
        "retmode": "json",
    }
    if api_key:
        search_params["api_key"] = api_key
    body = http_get("{}/esearch.fcgi?{}".format(EUTILS, urllib.parse.urlencode(search_params)))
    if not body.strip():
        return []

    ids = re.findall(r'"(\d+)"', body.split('"idlist"')[-1]) if '"idlist"' in body else []
    if not ids:
        return []

    fetch_params = {
        "db": "sra",
        "id": ",".join(ids[:500]),
        "rettype": "runinfo",
        "retmode": "csv",
    }
    if api_key:
        fetch_params["api_key"] = api_key
    csv_body = http_get("{}/efetch.fcgi?{}".format(EUTILS, urllib.parse.urlencode(fetch_params)))
    if not csv_body.strip():
        return []
    reader = csv.DictReader(io.StringIO(csv_body))
    return [row for row in reader if row.get("Run")]


def ncbi_records_to_rows(accession, records):
    rows = []
    for rec in records:
        layout = (rec.get("LibraryLayout") or "").strip().upper()
        rows.append(
            {
                "request": accession,
                "run_accession": rec["Run"].strip(),
                "source": "sratools",
                "platform": normalise_platform(rec.get("Platform")),
                "single_end": "false" if layout == "PAIRED" else "true",
                "fastq_1_url": "",
                "fastq_2_url": "",
                "fastq_1_md5": "",
                "fastq_2_md5": "",
                "instrument_platform": (rec.get("Platform") or "").strip(),
                "instrument_model": (rec.get("Model") or "").strip(),
                "library_strategy": (rec.get("LibraryStrategy") or "").strip(),
                "read_count": (rec.get("spots") or "").strip(),
            }
        )
    return rows


# ── Driver ──────────────────────────────────────────────────────────────────
def resolve(accession, force_sratools=False, api_key=None):
    """Resolve one user-supplied accession into a list of manifest rows."""
    if not is_accession(accession):
        raise ValueError(
            "'{}' is not a recognised SRA/ENA accession. Expected a run "
            "(SRR/ERR/DRR), experiment (SRX/ERX/DRX), sample (SRS/ERS/DRS, "
            "SAMN/SAMEA) or study (SRP/ERP/DRP, PRJNA/PRJEB) accession.".format(accession)
        )

    rows = ena_records_to_rows(accession, query_ena(accession))
    if not rows:
        logger.warning("ENA has no read_run records for %s; falling back to NCBI eutils.", accession)
        rows = ncbi_records_to_rows(accession, query_ncbi_runinfo(accession, api_key=api_key))

    if not rows:
        raise ValueError(
            "Could not resolve accession '{}' via ENA or NCBI. Check the accession "
            "is public and spelled correctly.".format(accession)
        )

    if force_sratools:
        for row in rows:
            row["source"] = "sratools"

    return rows


def name_samples(base_name, rows):
    """
    Assign a unique pipeline sample name to each run.

    One run  → the requested sample name is used as-is.
    Many runs → names are suffixed with the run accession so a project
    accession fans out into ``<sample>_SRR000001``, ``<sample>_SRR000002``, …
    """
    single = len(rows) == 1
    for row in rows:
        row["sample"] = base_name if single else "{}_{}".format(base_name, row["run_accession"])
    return rows


def read_requests(path):
    """
    Read the ``sample,accession`` request file written by input_check.nf.

    The file is normally headerless (Nextflow's collectFile just concatenates
    one ``sample,accession`` line per row), but a ``sample,accession`` header
    line is tolerated and skipped so the script is also usable by hand.
    """
    requests = []
    with open(path, newline="", encoding="utf-8-sig") as handle:
        for fields in csv.reader(handle):
            if not fields:
                continue
            sample = fields[0].strip()
            accession = fields[1].strip() if len(fields) > 1 else ""
            if not accession:
                continue
            if sample.lower() == "sample" and accession.lower() == "accession":
                continue  # header line
            requests.append((sample or accession, accession))
    if not requests:
        raise ValueError("No accessions found in {}".format(path))
    return requests


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Resolve SRA/ENA accessions into a downloadable FASTQ manifest.",
        epilog="Example: resolve_sra.py -i sra_requests.csv -o sra_manifest.csv",
    )
    parser.add_argument("-i", "--input", required=True,
                        help="CSV of requests with columns 'sample,accession'.")
    parser.add_argument("-o", "--output", required=True, help="Output manifest CSV.")
    parser.add_argument("--force-sratools", action="store_true",
                        help="Ignore ENA FASTQ URLs and route every run through sra-tools.")
    parser.add_argument("--api-key", default=None, help="Optional NCBI eutils API key.")
    parser.add_argument("-l", "--log-level", default="INFO",
                        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"))
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")

    try:
        requests = read_requests(args.input)
    except (OSError, ValueError) as err:
        logger.critical(str(err))
        return 2

    all_rows = []
    seen_samples = set()
    seen_runs = set()
    for base_name, accession in requests:
        try:
            rows = name_samples(base_name, resolve(accession,
                                                   force_sratools=args.force_sratools,
                                                   api_key=args.api_key))
        except ValueError as err:
            logger.critical(str(err))
            return 1

        # A run may be reachable from more than one request (e.g. the user lists
        # both PRJNA000000 and one of its runs). Keep the first occurrence only:
        # downstream the run accession is the join key and the storeDir cache key,
        # both of which must be unique.
        deduped = []
        for row in rows:
            if row["run_accession"] in seen_runs:
                logger.warning("Run %s already resolved from an earlier accession; skipping duplicate from %s.",
                               row["run_accession"], accession)
                continue
            seen_runs.add(row["run_accession"])
            deduped.append(row)
        rows = deduped

        for row in rows:
            # Guard against two requests colliding on one sample name
            name = row["sample"]
            suffix = 2
            while name in seen_samples:
                name = "{}_{}".format(row["sample"], suffix)
                suffix += 1
            row["sample"] = name
            seen_samples.add(name)

        logger.info("%s -> %s run(s): %s", accession, len(rows),
                    ", ".join(r["run_accession"] for r in rows))
        all_rows.extend(rows)

    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, MANIFEST_FIELDS, extrasaction="ignore")
        writer.writeheader()
        for row in all_rows:
            writer.writerow(row)

    logger.info("Wrote %s run(s) to %s", len(all_rows), args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
