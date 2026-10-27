#!/usr/bin/env python3
##############################################################################################
# Copyright 2026 The Johns Hopkins University Applied Physics Laboratory LLC
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
"""Build a reference FASTA + taxid mapping from a plain list of contig accessions.

TaxTriage pulls whole assemblies for taxids/names (--organisms / --organisms_file /
--get_pathogens). When you only care about a fixed set of CONTIGS (nuccore
accessions such as NC_045512.2, NZ_CP020549.1, OR833055.1), this script builds the
two files the pipeline needs to use exactly those sequences:

    <prefix>.fasta               ->  --reference_fasta
    <prefix>.accession_map.tsv   ->  --custom_accession_map   (accession, taxid, name)

plus the same mapping in the pipeline's own format and a full metadata table.

Where sequences come from (in order):
  1. --source-fasta  local FASTA(s), e.g. the prebuilt pathogens.fasta(.gz). Plain
                     FASTAs with a samtools .fai are read by random access; anything
                     else (incl. gzip/bgzip) is streamed once.
  2. NCBI nuccore    efetch, batched by total sequence length (skipped with --offline).

Where metadata comes from (taxid, organism, assembly):
  1. --source-map    existing mapping file(s): a pipeline mapping
                     (Acc/Assembly/Organism_Name/Description/Mapped_Value, e.g.
                     pathogens.mapping.taxid.tsv), a 4-column gcfmapping/localmap, or a
                     custom accession map (accession/taxid/name). Wins when it has a taxid.
  2. NCBI esummary   nuccore docsum (taxid, organism, title, length, AssemblyAcc).
  3. NCBI elink      nuccore -> assembly, for contigs whose docsum has no assembly
                     (skipped with --no-elink / --offline).
  4. --assembly      assembly_summary_*.txt: fills a missing taxid from the assembly
                     and adds organism / species_taxid / refseq_category columns.
  5. --taxdump       the taxonomy the pipeline runs with (names.dmp + nodes.dmp, and
                     merged.dmp / delnodes.dmp when present). Merged taxids are
                     remapped, deleted / unknown ones flagged, the scientific name and
                     species-level taxid are taken from it.

Resume: every batch of sequences is fsync'd to <prefix>.fasta and checkpointed in
<outdir>/.<prefix>.state/; metadata and assembly lookups are cached there too.
Re-running the same command picks up where it stopped (a half-written batch is
truncated away). --restart discards the state and starts over. Accessions added to
the list later are simply fetched on the next run; ones removed are dropped from
the final FASTA.

Outputs (in --outdir, prefix defaults to the input file's stem):
  <prefix>.fasta / .fasta.fai     sequences, 80-col wrapped, headers ">ACC.V description"
  <prefix>.accession_map.tsv      accession, taxid, name  (for --custom_accession_map)
  <prefix>.mapping.taxid.tsv      Acc, Assembly, Organism_Name, Description, Mapped_Value
  <prefix>.metadata.tsv           everything: requested id, record id, assembly, taxid,
                                  original taxid, taxid status, species, rank, length,
                                  sources, notes
  <prefix>.assemblies.txt         unique assembly accessions
  <prefix>.missing.tsv            requested id + reason, for anything not delivered
  <prefix>.summary.json           counts for the run
  <prefix>.log                    timestamped log (appended across resumes)

Examples
  # online, normalise taxids against the taxdump the pipeline uses
  pull_contigs.py -i contigs.txt -o refs/ --taxdump taxdump/ -e you@org.org

  # offline, from the prebuilt pathogen FASTA + its mapping
  pull_contigs.py -i contigs.txt -o refs/ --offline \\
      --source-fasta pathogens.fasta --source-map pathogens.mapping.taxid.tsv \\
      --assembly assembly_summary_refseq.txt --taxdump taxdump/

  # then
  nextflow run ... --reference_fasta refs/contigs.fasta \\
      --custom_accession_map refs/contigs.accession_map.tsv

Only the Python standard library is used, so it runs in any pipeline container.
"""

import argparse
import array
import csv
import gzip
import json
import logging
import os
import random
import re
import shutil
import ssl
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import OrderedDict, defaultdict

EUTILS_BASE = os.environ.get("NCBI_EUTILS_BASE", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/")
WRAP = 80
LOG = logging.getLogger("pull_contigs")

HEADER_WORDS = {"accession", "accessions", "contig", "contigs", "acc", "id", "seqid", "seq_id",
                "accession_version", "nuccore", "name", "sequence", "sequence_id"}
DB_TAGS = {"ref", "gb", "emb", "dbj", "tpg", "tpe", "tpd", "pdb", "sp", "tr", "lcl", "gnl"}


# ─────────────────────────────────────────────────────────────────────────────
# Accession helpers
# ─────────────────────────────────────────────────────────────────────────────
def base_acc(acc):
    return re.sub(r"\.\d+$", "", acc)


def has_version(acc):
    return re.search(r"\.\d+$", acc) is not None


def clean_id(token):
    """Normalise one id token: strip '>' and pick the accession out of pipe ids
    (gi|123|ref|NC_045512.2|, kraken:taxid|9606|NC_000001.11)."""
    tok = token.strip().lstrip(">").strip()
    if "|" in tok:
        parts = [p for p in tok.split("|") if p]
        for i, p in enumerate(parts):
            if p.lower() in DB_TAGS and i + 1 < len(parts):
                return parts[i + 1]
        cands = [p for p in parts if re.search(r"[A-Za-z]", p) and re.search(r"\d", p)
                 and ":" not in p]
        if cands:
            return cands[-1]
    return tok


def read_contig_list(path):
    """First column of each line; blank / '#' lines skipped; a header row is ignored."""
    out, seen, dups = [], set(), 0
    opener = gzip.open if _is_gzip(path) else open
    with opener(path, "rt") as fh:
        first = True
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            tok = re.split(r"[\t,;\s]+", line)[0]
            if first:
                first = False
                if tok.lower().lstrip("#>") in HEADER_WORDS:
                    continue
            acc = clean_id(tok)
            if not acc:
                continue
            if acc in seen:
                dups += 1
                continue
            seen.add(acc)
            out.append(acc)
    return out, dups


def _is_gzip(path):
    with open(path, "rb") as fh:
        return fh.read(2) == b"\x1f\x8b"


# ─────────────────────────────────────────────────────────────────────────────
# NCBI E-utilities (stdlib, rate limited, retry with backoff)
# ─────────────────────────────────────────────────────────────────────────────
class Eutils:
    def __init__(self, email=None, api_key=None, rate=None, timeout=120, retries=6):
        self.email = email
        self.api_key = api_key or os.environ.get("NCBI_API_KEY") or None
        rate = rate or (8.0 if self.api_key else 2.5)
        self.interval = 1.0 / rate
        self.next_ok = 0.0
        self.timeout = timeout
        self.retries = retries
        self.ctx = ssl.create_default_context()
        self.insecure = False
        self.requests = 0

    def _common(self, db):
        p = {"db": db, "tool": "taxtriage_pull_contigs"}
        if self.email:
            p["email"] = self.email
        if self.api_key:
            p["api_key"] = self.api_key
        return p

    def _wait(self):
        now = time.monotonic()
        if now < self.next_ok:
            time.sleep(self.next_ok - now)
        self.next_ok = time.monotonic() + self.interval

    def request(self, endpoint, params, retries=None, timeout=None):
        data = urllib.parse.urlencode(params, doseq=True).encode()
        url = EUTILS_BASE.rstrip("/") + "/" + endpoint
        retries = retries or self.retries
        timeout = timeout or self.timeout
        last = None
        for attempt in range(retries):
            self._wait()
            try:
                req = urllib.request.Request(url, data=data,
                                             headers={"User-Agent": "taxtriage-pull-contigs"})
                with urllib.request.urlopen(req, context=self.ctx, timeout=timeout) as r:
                    body = r.read().decode("utf-8", errors="replace")
                self.requests += 1
                return body
            except urllib.error.HTTPError as e:
                last = e
                if e.code in (429, 500, 502, 503, 504):
                    ra = e.headers.get("Retry-After") if e.headers else None
                    try:
                        delay = float(ra)
                    except (TypeError, ValueError):
                        delay = 2.0 ** attempt + random.uniform(0, 0.5)
                    if e.code == 429:
                        self.interval = min(self.interval * 1.5, 2.0)
                    LOG.warning("NCBI %s HTTP %s, retry %d/%d in %.1fs",
                                endpoint, e.code, attempt + 1, retries, delay)
                    time.sleep(min(delay, 60))
                    continue
                raise
            except urllib.error.URLError as e:
                last = e
                if isinstance(getattr(e, "reason", None), ssl.SSLError) and not self.insecure:
                    # Several pipeline containers ship stale CA bundles (same
                    # fallback as bin/ncbi_taxid.py).
                    LOG.warning("SSL verification failed (%s); retrying without verification",
                                e.reason)
                    self.ctx = ssl.create_default_context()
                    self.ctx.check_hostname = False
                    self.ctx.verify_mode = ssl.CERT_NONE
                    self.insecure = True
                    continue
                delay = 2.0 ** attempt + random.uniform(0, 0.5)
                LOG.warning("NCBI %s network error (%s), retry %d/%d in %.1fs",
                            endpoint, e.reason, attempt + 1, retries, delay)
                time.sleep(min(delay, 60))
            except Exception as e:  # noqa: BLE001 - socket timeouts, resets
                last = e
                delay = 2.0 ** attempt + random.uniform(0, 0.5)
                LOG.warning("NCBI %s error (%s), retry %d/%d in %.1fs",
                            endpoint, e, attempt + 1, retries, delay)
                time.sleep(min(delay, 60))
        raise RuntimeError(f"NCBI {endpoint} failed after {retries} attempts: {last}")

    def esearch_accn(self, accs):
        p = self._common("nuccore")
        p.update(term=" OR ".join(f"{a}[ACCN]" for a in accs), retmode="json",
                 retmax=str(len(accs) * 3 + 10))
        js = json.loads(self.request("esearch.fcgi", p))
        return js.get("esearchresult", {}).get("idlist", [])

    def esummary(self, db, ids):
        p = self._common(db)
        p.update(id=",".join(ids), retmode="json")
        js = json.loads(self.request("esummary.fcgi", p))
        res = js.get("result", {}) or {}
        return [res[u] for u in res.get("uids", []) if isinstance(res.get(u), dict)
                and "error" not in res[u]]

    def efetch_fasta(self, uids):
        p = self._common("nuccore")
        p.update(id=",".join(uids), rettype="fasta", retmode="text")
        return self.request("efetch.fcgi", p)

    def elink_assembly(self, uids, retries=2, timeout=60):
        """{nuccore uid: [assembly uid, ...]} (one linkset per id).

        Each repeated id= is a separate link computation on NCBI's side, so a large
        batch can outlast NCBI's backend limit (~15 s) and the connection is dropped
        with an empty body (IncompleteRead(0 bytes)). Retrying the same batch does
        not help, so fail fast here and let ncbi_assemblies() split it."""
        p = self._common("assembly")
        p.update(dbfrom="nuccore", id=list(uids), retmode="json", linkname="nuccore_assembly")
        js = json.loads(self.request("elink.fcgi", p, retries=retries, timeout=timeout))
        out = {}
        for ls in js.get("linksets", []):
            ids = [str(i) for i in ls.get("ids", [])]
            links = []
            for db in ls.get("linksetdbs", []) or []:
                if db.get("dbto") == "assembly":
                    links.extend(str(x) for x in db.get("links", []))
            for i in ids:
                out[i] = links
        return out


# ─────────────────────────────────────────────────────────────────────────────
# State: checkpointed FASTA + jsonl caches
# ─────────────────────────────────────────────────────────────────────────────
class State:
    """<prefix>.fasta is appended batch by batch; <state>/fasta.ckpt holds one line per
    record that is fully on disk (acc, requested, length, header offset, seq offset,
    end offset, source). meta.jsonl caches NCBI/source-map metadata per requested id;
    assembly.jsonl caches nuccore->assembly links (incl. 'none found')."""

    def __init__(self, outdir, prefix, restart=False):
        self.fasta = os.path.join(outdir, f"{prefix}.fasta")
        self.dir = os.path.join(outdir, f".{prefix}.state")
        if restart:
            for p in (self.fasta, self.fasta + ".fai"):
                if os.path.exists(p):
                    os.remove(p)
            shutil.rmtree(self.dir, ignore_errors=True)
        os.makedirs(self.dir, exist_ok=True)
        self.ckpt = os.path.join(self.dir, "fasta.ckpt")
        self.meta_path = os.path.join(self.dir, "meta.jsonl")
        self.asm_path = os.path.join(self.dir, "assembly.jsonl")
        self.records = OrderedDict()   # record acc -> dict
        self.by_requested = {}         # requested id -> record acc
        self.meta = {}                 # requested id -> dict
        self.assembly = {}             # record acc -> {"assembly", "source"}
        self._recover()

    def _recover(self):
        entries = []
        if os.path.exists(self.ckpt):
            with open(self.ckpt) as fh:
                for line in fh:
                    if not line.endswith("\n"):
                        break  # torn last line
                    f = line.rstrip("\n").split("\t")
                    if len(f) < 7:
                        continue
                    entries.append(dict(acc=f[0], requested=f[1].split(",") if f[1] else [],
                                        length=int(f[2]), hdr_off=int(f[3]), seq_off=int(f[4]),
                                        end=int(f[5]), source=f[6]))
        size = os.path.getsize(self.fasta) if os.path.exists(self.fasta) else 0
        good = [e for e in entries if e["end"] <= size]
        good_end = max((e["end"] for e in good), default=0)
        dropped = len(entries) - len(good)
        if size > good_end:
            LOG.info("Resume: truncating %s from %d to %d bytes (partial batch from an "
                     "interrupted run)", os.path.basename(self.fasta), size, good_end)
            with open(self.fasta, "r+b") as fh:
                fh.truncate(good_end)
        if dropped:
            LOG.warning("Resume: %d checkpoint entries point past the end of the FASTA; "
                        "they will be fetched again", dropped)
        if dropped or size > good_end or len(good) != len(entries):
            self._rewrite_ckpt(good)
        for e in good:
            self.records[e["acc"]] = e
            for r in e["requested"]:
                self.by_requested[r] = e["acc"]
        for path, target, key in ((self.meta_path, self.meta, "requested"),
                                  (self.asm_path, self.assembly, "acc")):
            if os.path.exists(path):
                with open(path) as fh:
                    for line in fh:
                        try:
                            d = json.loads(line)
                        except ValueError:
                            continue  # torn line
                        target[d[key]] = d
        if self.records or self.meta:
            LOG.info("Resume: %d sequences already on disk, %d metadata entries and %d "
                     "assembly lookups cached", len(self.records), len(self.meta),
                     len(self.assembly))

    def _rewrite_ckpt(self, entries):
        tmp = self.ckpt + ".tmp"
        with open(tmp, "w") as fh:
            for e in entries:
                fh.write(self._ckpt_line(e))
        os.replace(tmp, self.ckpt)

    @staticmethod
    def _ckpt_line(e):
        return "\t".join([e["acc"], ",".join(e["requested"]), str(e["length"]), str(e["hdr_off"]),
                          str(e["seq_off"]), str(e["end"]), e["source"]]) + "\n"

    def write_batch(self, recs, source):
        """recs: list of (acc, description, sequence, [requested ids]). Durable on return."""
        if not recs:
            return 0
        new = []
        with open(self.fasta, "ab") as fh:
            for acc, desc, seq, requested in recs:
                if acc in self.records:
                    # same record reached via another requested id (X vs X.2)
                    e = self.records[acc]
                    for r in requested:
                        if r not in e["requested"]:
                            e["requested"].append(r)
                            self.by_requested[r] = acc
                    new.append(e)
                    continue
                hdr_off = fh.tell()
                fh.write((f">{acc} {desc}".rstrip() + "\n").encode())
                seq_off = fh.tell()
                for i in range(0, len(seq), WRAP):
                    fh.write(seq[i:i + WRAP].encode() + b"\n")
                e = dict(acc=acc, requested=list(requested), length=len(seq), hdr_off=hdr_off,
                         seq_off=seq_off, end=fh.tell(), source=source)
                self.records[acc] = e
                for r in requested:
                    self.by_requested[r] = acc
                new.append(e)
            fh.flush()
            os.fsync(fh.fileno())
        with open(self.ckpt, "a") as ck:
            for e in new:
                ck.write(self._ckpt_line(e))
            ck.flush()
            os.fsync(ck.fileno())
        return len(new)

    def put_meta(self, items):
        with open(self.meta_path, "a") as fh:
            for d in items:
                self.meta[d["requested"]] = d
                fh.write(json.dumps(d) + "\n")
            fh.flush()
            os.fsync(fh.fileno())

    def put_assembly(self, items):
        with open(self.asm_path, "a") as fh:
            for d in items:
                self.assembly[d["acc"]] = d
                fh.write(json.dumps(d) + "\n")
            fh.flush()
            os.fsync(fh.fileno())


# ─────────────────────────────────────────────────────────────────────────────
# Matching requested ids to record ids
# ─────────────────────────────────────────────────────────────────────────────
class Wanted:
    """Requested ids indexed exactly and by version-less base. An unversioned request
    matches any version; a versioned one prefers the exact version and falls back to
    another version of the same accession (flagged as version_changed)."""

    def __init__(self, ids):
        self.exact = set(ids)
        self.by_base = defaultdict(list)
        for a in ids:
            self.by_base[base_acc(a)].append(a)

    def match(self, record_id):
        """requested ids satisfied by record_id: (exact matches, loose matches)."""
        exact = [record_id] if record_id in self.exact else []
        loose = [r for r in self.by_base.get(base_acc(record_id), []) if r != record_id]
        return exact, loose

    def remove(self, ids):
        for a in ids:
            self.exact.discard(a)
            b = base_acc(a)
            if a in self.by_base.get(b, []):
                self.by_base[b].remove(a)
                if not self.by_base[b]:
                    del self.by_base[b]

    def __len__(self):
        return len(self.exact)


# ─────────────────────────────────────────────────────────────────────────────
# Local FASTA source
# ─────────────────────────────────────────────────────────────────────────────
def _read_fai(path):
    out = OrderedDict()
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 5:
                out[f[0]] = (int(f[1]), int(f[2]), int(f[3]), int(f[4]))
    return out


def _fai_fetch(fh, length, offset, linebases, linewidth):
    nlines = (length + linebases - 1) // linebases if linebases else 0
    span = length + nlines * (linewidth - linebases)
    fh.seek(offset)
    raw = fh.read(span).decode("ascii", errors="replace")
    return re.sub(r"\s+", "", raw)[:length]


def _pick(wanted, record_id, allow_other_version):
    """Requested ids this record satisfies. An unversioned request takes any version;
    a versioned request takes a different version only when allow_other_version
    (offline: the local copy is the best there is; online: NCBI has the exact one)."""
    exact, loose = wanted.match(record_id)
    loose = [r for r in loose if allow_other_version or not has_version(r)]
    return exact + loose


def scan_local_fastas(paths, wanted, state, batch_bp, allow_other_version):
    """Copy wanted records from local FASTA(s) into the output. Returns count."""
    found = 0
    for path in paths:
        if not len(wanted):
            break
        if not os.path.exists(path):
            LOG.warning("Source FASTA %s not found, skipping", path)
            continue
        gz = _is_gzip(path)
        fai = path + ".fai"
        t0 = time.time()
        n_before = found
        batch, bp = [], 0

        def flush():
            nonlocal batch, bp, found
            if batch:
                found += state.write_batch(batch, f"local:{os.path.basename(path)}")
                batch, bp = [], 0

        use_fai = not gz and os.path.exists(fai)
        if use_fai and os.path.getmtime(fai) < os.path.getmtime(path):
            LOG.warning("Source %s: %s is older than the FASTA (stale index; re-run "
                        "'samtools faidx %s'); streaming instead", path,
                        os.path.basename(fai), path)
            use_fai = False
        stale = False
        if use_fai:
            LOG.info("Source %s: random access through %s", path, os.path.basename(fai))
            idx = _read_fai(fai)
            # map wanted -> fai names (by exact id or base)
            fai_by_base = defaultdict(list)
            for name in idx:
                fai_by_base[base_acc(clean_id(name))].append(name)
            hits = []
            for req in sorted(wanted.exact):
                names = fai_by_base.get(base_acc(req), [])
                exact = [n for n in names if clean_id(n) == req]
                if exact:
                    pick = exact[0]
                elif names and (not has_version(req) or allow_other_version):
                    pick = names[-1]  # other version: noted as version_changed
                else:
                    continue
                hits.append((idx[pick][1], pick, req))
            # descriptions live on header lines; read them while seeking in file order
            with open(path, "rb") as fh:
                for off, name, req in sorted(hits):
                    if req not in wanted.exact:
                        continue
                    L, o, lb, lw = idx[name]
                    # header line ends right before the sequence offset
                    back = min(o, 4096)
                    fh.seek(o - back)
                    chunk = fh.read(back).decode("ascii", errors="replace")
                    hdr = chunk.rstrip("\n").rsplit("\n", 1)[-1]
                    hdr = hdr[hdr.rfind(">") + 1:] if ">" in hdr else name
                    parts = hdr.split(None, 1)
                    rid = clean_id(parts[0]) if parts else ""
                    desc = parts[1] if len(parts) > 1 else ""
                    seq = _fai_fetch(fh, L, o, lb, lw)
                    # Trust the index only if it agrees with the file: header name at
                    # the offset, no '>' inside the span, next record starts right after.
                    nlines = (L + lb - 1) // lb if lb else 0
                    fh.seek(o + L + nlines * (lw - lb))
                    nxt = fh.read(1)
                    if rid != clean_id(name) or ">" in seq or nxt not in (b">", b""):
                        stale = True
                        LOG.warning("Source %s: %s does not match the FASTA at %s (stale "
                                    "index; re-run 'samtools faidx %s'). Falling back to "
                                    "streaming.", path, os.path.basename(fai), name, path)
                        break
                    reqs = _pick(wanted, rid, allow_other_version)
                    if not reqs:
                        continue
                    wanted.remove(reqs)
                    batch.append((rid, desc, seq, reqs))
                    bp += L
                    if bp >= batch_bp:
                        flush()
            flush()
        if (not use_fai or stale) and len(wanted):
            LOG.info("Source %s: streaming (%s)", path,
                     "gzip" if gz else ("stale .fai" if (stale or os.path.exists(fai))
                                        else "no .fai index"))
            opener = gzip.open if gz else open
            with opener(path, "rt") as fh:
                cur, desc, seqparts, reqs = None, "", [], []
                for line in fh:
                    if line.startswith(">"):
                        if reqs:
                            seq = "".join(seqparts)
                            batch.append((cur, desc, seq, reqs))
                            bp += len(seq)
                            if bp >= batch_bp:
                                flush()
                        if not len(wanted):
                            reqs = []
                            break
                        parts = line[1:].strip().split(None, 1)
                        cur = clean_id(parts[0]) if parts else ""
                        desc = parts[1] if len(parts) > 1 else ""
                        reqs = _pick(wanted, cur, allow_other_version)
                        if reqs:
                            wanted.remove(reqs)
                        seqparts = []
                    elif reqs:
                        seqparts.append(line.strip())
                if reqs:
                    batch.append((cur, desc, "".join(seqparts), reqs))
            flush()
        LOG.info("Source %s: %d record(s) taken in %.1fs", os.path.basename(path),
                 found - n_before, time.time() - t0)
    return found


# ─────────────────────────────────────────────────────────────────────────────
# Source maps (offline metadata)
# ─────────────────────────────────────────────────────────────────────────────
_ACC_H = {"acc", "accession", "accession_version", "seqid", "seq_id", "id", "nuccore", "contig"}
_TAX_H = {"taxid", "tax_id", "taxonomy_id", "ncbi_taxid", "mapped_value"}
_ASM_H = {"assembly", "assembly_accession", "gcf", "gca", "assembly_acc"}
_NAME_H = {"organism_name", "organism", "name", "species", "custom_name"}
_DESC_H = {"description", "title", "desc"}


def _norm_h(h):
    return re.sub(r"[^a-z0-9]+", "_", str(h or "").strip().lower()).strip("_")


def _clean_taxid(v):
    t = re.sub(r"\.0+$", "", str(v or "").strip())
    return t if re.fullmatch(r"\d+", t) and t != "0" else ""


def load_source_maps(paths):
    """{acc and base acc: {taxid, assembly, name, description}} from mapping files."""
    out = {}
    for path in paths or []:
        opener = gzip.open if _is_gzip(path) else open
        with opener(path, "rt", newline="", encoding="utf-8-sig") as fh:
            lines = [l for l in fh if l.strip() and not l.startswith("# ")]
        if not lines:
            continue
        if "\t" in lines[0]:
            delim = "\t"   # descriptions contain commas; never sniff a TSV as CSV
        else:
            try:
                delim = csv.Sniffer().sniff(lines[0], delimiters=",;").delimiter
            except csv.Error:
                delim = "\t"
        rows = list(csv.reader(lines, delimiter=delim))
        header = [_norm_h(h) for h in rows[0]]

        def col(names, hdr=header):
            for i, h in enumerate(hdr):
                if h in names:
                    return i
            return None

        ai = col(_ACC_H)
        if ai is not None:
            ti, si, ni, di = col(_TAX_H), col(_ASM_H), col(_NAME_H), col(_DESC_H)
            body = rows[1:]
        else:
            # headerless gcfmapping / localmap: Acc, Assembly, Organism_Name, Description[, taxid]
            ai, si, ni, di = 0, 1, 2, 3
            ti = 4 if len(rows[0]) > 4 else None
            body = rows
        n = 0
        for r in body:
            if ai >= len(r) or not r[ai].strip():
                continue
            acc = clean_id(r[ai])
            get = lambda i: r[i].strip() if i is not None and i < len(r) else ""  # noqa: E731
            d = dict(taxid=_clean_taxid(get(ti)), assembly=get(si), name=get(ni),
                     description=get(di), source=os.path.basename(path))
            if d["assembly"] and not re.match(r"^GC[AF]_\d+", d["assembly"]):
                d["assembly"] = ""
            for k in (acc, base_acc(acc)):
                prev = out.get(k)
                # keep the most complete entry
                if prev is None or (not prev["taxid"] and d["taxid"]):
                    out[k] = d
            n += 1
        LOG.info("Source map %s: %d accession rows", path, n)
    return out


# ─────────────────────────────────────────────────────────────────────────────
# Assembly summaries and taxdump
# ─────────────────────────────────────────────────────────────────────────────
def load_assembly_summaries(paths, wanted):
    """{GCF/GCA: {taxid, species_taxid, organism_name, refseq_category, assembly_level}}"""
    out = {}
    if not wanted:
        return out
    wanted_base = {base_acc(a): a for a in wanted}
    for path in paths or []:
        with open(path, "r", errors="replace") as fh:
            cols = None
            for line in fh:
                if line.startswith("#"):
                    h = line.lstrip("#").strip().split("\t")
                    if "assembly_accession" in [x.strip() for x in h]:
                        cols = [x.strip() for x in h]
                    continue
                acc = line.split("\t", 1)[0]
                if acc not in wanted and base_acc(acc) not in wanted_base:
                    continue
                f = line.rstrip("\r\n").split("\t")
                c = cols or ["assembly_accession", "bioproject", "biosample", "wgs_master",
                             "refseq_category", "taxid", "species_taxid", "organism_name",
                             "infraspecific_name", "isolate", "version_status", "assembly_level"]
                rec = {k: (f[i] if i < len(f) else "") for i, k in enumerate(c)}
                for key in (acc, base_acc(acc)):
                    if key not in out:
                        out[key] = rec
    return out


class Taxdump:
    def __init__(self, path):
        self.path = path
        nodes = os.path.join(path, "nodes.dmp")
        names = os.path.join(path, "names.dmp")
        for p in (nodes, names):
            if not os.path.exists(p):
                raise SystemExit(f"ERROR: --taxdump {path} has no {os.path.basename(p)}")
        t0 = time.time()
        self.merged, self.deleted = {}, set()
        mp = os.path.join(path, "merged.dmp")
        if os.path.exists(mp):
            with open(mp) as fh:
                for line in fh:
                    f = line.split("\t|\t")
                    self.merged[int(f[0])] = int(f[1].split("\t")[0])
        dp = os.path.join(path, "delnodes.dmp")
        if os.path.exists(dp):
            with open(dp) as fh:
                for line in fh:
                    self.deleted.add(int(line.split("\t")[0]))
        # compact arrays indexed by taxid: parent + rank code
        self.parent = array.array("i")
        self.rank = array.array("b")
        self.rank_names, rank_code = [], {}
        with open(nodes) as fh:
            for line in fh:
                f = line.split("\t|\t", 3)
                tid, par, rk = int(f[0]), int(f[1]), f[2].split("\t|")[0].strip()
                if tid >= len(self.parent):
                    grow = tid + 1 - len(self.parent) + 100000
                    self.parent.extend([-1] * grow)
                    self.rank.extend([-1] * grow)
                if rk not in rank_code:
                    rank_code[rk] = len(self.rank_names)
                    self.rank_names.append(rk)
                self.parent[tid] = par
                self.rank[tid] = rank_code[rk]
        self.names = {}
        LOG.info("Taxdump %s: nodes loaded (%d merged, %d deleted) in %.1fs", path,
                 len(self.merged), len(self.deleted), time.time() - t0)

    def exists(self, t):
        return 0 < t < len(self.parent) and self.parent[t] != -1

    def resolve(self, taxid):
        """(current taxid or None, status)"""
        try:
            t = int(taxid)
        except (TypeError, ValueError):
            return None, "missing"
        if self.exists(t):
            return t, "ok"
        if t in self.merged and self.exists(self.merged[t]):
            return self.merged[t], "merged"
        if t in self.deleted:
            return None, "deleted"
        return None, "not_in_taxdump"

    def rank_of(self, t):
        return self.rank_names[self.rank[t]] if self.exists(t) else ""

    def ancestor_at(self, t, rank="species"):
        seen = 0
        while self.exists(t) and seen < 100:
            if self.rank_of(t) == rank:
                return t
            p = self.parent[t]
            if p == t:
                break
            t, seen = p, seen + 1
        return None

    def load_names(self, taxids):
        want = {int(t) for t in taxids if t}
        if not want:
            return
        with open(os.path.join(self.path, "names.dmp")) as fh:
            for line in fh:
                if "scientific name" not in line:
                    continue
                f = line.split("\t|\t")
                t = int(f[0])
                if t in want:
                    self.names[t] = f[1]

    def name(self, t):
        return self.names.get(int(t), "") if t else ""


# ─────────────────────────────────────────────────────────────────────────────
# NCBI steps
# ─────────────────────────────────────────────────────────────────────────────
def _docsum_to_meta(ds):
    return dict(acc=ds.get("accessionversion") or ds.get("caption") or "",
                uid=str(ds.get("uid") or ""), taxid=_clean_taxid(ds.get("taxid")),
                organism=ds.get("organism") or "", title=ds.get("title") or "",
                length=int(ds.get("slen") or 0), strain=ds.get("strain") or "",
                assembly=ds.get("assemblyacc") or "", moltype=ds.get("moltype") or "",
                completeness=ds.get("completeness") or "", source="ncbi_esummary")


def ncbi_metadata(eu, queries, state, batch_size):
    """queries: {requested id: accession to look up}. Caches results in state."""
    items = list(queries.items())
    nb = (len(items) + batch_size - 1) // batch_size
    found = 0
    for b in range(nb):
        chunk = items[b * batch_size:(b + 1) * batch_size]
        accs = list(dict.fromkeys(q for _, q in chunk))
        docs = []
        try:
            uids = eu.esearch_accn(accs)
            if uids:
                docs = eu.esummary("nuccore", uids)
        except Exception as e:  # noqa: BLE001
            LOG.error("Metadata batch %d/%d failed: %s (will retry on the next run)", b + 1, nb, e)
            continue
        by_ver, by_base = {}, defaultdict(list)
        for ds in docs:
            m = _docsum_to_meta(ds)
            if m["acc"]:
                by_ver[m["acc"]] = m
                by_base[base_acc(m["acc"])].append(m)
        # versioned requests whose exact version esearch did not return (older
        # versions): ask esummary for them directly.
        retry = [q for _, q in chunk if has_version(q) and q not in by_ver]
        if retry:
            try:
                for ds in eu.esummary("nuccore", list(dict.fromkeys(retry))):
                    m = _docsum_to_meta(ds)
                    if m["acc"]:
                        by_ver[m["acc"]] = m
                        by_base[base_acc(m["acc"])].append(m)
            except Exception as e:  # noqa: BLE001
                LOG.debug("direct esummary for %d versioned ids failed: %s", len(retry), e)
        out = []
        for req, q in chunk:
            m = by_ver.get(q)
            note = ""
            if m is None and by_base.get(base_acc(q)):
                cands = sorted(by_base[base_acc(q)],
                               key=lambda x: int(x["acc"].rsplit(".", 1)[-1]) if has_version(x["acc"]) else 0)
                m = cands[-1]
                if has_version(q):
                    note = f"version_changed:{q}->{m['acc']}"
            if m:
                d = dict(m, requested=req, note=note)
                out.append(d)
        state.put_meta(out)
        found += len(out)
        LOG.info("Metadata batch %d/%d: %d/%d resolved at NCBI", b + 1, nb, len(out), len(chunk))
    return found


def parse_fasta_text(text):
    recs, cur, desc, parts = [], None, "", []
    for line in text.splitlines():
        if line.startswith(">"):
            if cur is not None:
                recs.append((cur, desc, "".join(parts)))
            p = line[1:].strip().split(None, 1)
            cur = clean_id(p[0]) if p else ""
            desc = p[1] if len(p) > 1 else ""
            parts = []
        elif cur is not None:
            parts.append(line.strip())
    if cur is not None:
        recs.append((cur, desc, "".join(parts)))
    return recs


def ncbi_sequences(eu, todo, state, max_bp, max_ids):
    """todo: list of (requested, meta). Fetch by uid, grouped by total length."""
    groups, cur, bp = [], [], 0
    for req, m in todo:
        L = m.get("length") or 0
        if cur and (bp + L > max_bp or len(cur) >= max_ids):
            groups.append(cur)
            cur, bp = [], 0
        cur.append((req, m))
        bp += L
    if cur:
        groups.append(cur)
    got = 0
    for gi, grp in enumerate(groups, 1):
        uids = list(dict.fromkeys(m["uid"] for _, m in grp if m.get("uid")))
        t0 = time.time()
        try:
            recs = parse_fasta_text(eu.efetch_fasta(uids))
        except Exception as e:  # noqa: BLE001
            LOG.error("Sequence batch %d/%d failed: %s (will retry on the next run)",
                      gi, len(groups), e)
            continue
        by_id = {r[0]: r for r in recs}
        by_base = {base_acc(r[0]): r for r in recs}
        batch, short = [], 0
        for req, m in grp:
            r = by_id.get(m["acc"]) or by_base.get(base_acc(m["acc"]))
            if r is None:
                continue
            if m.get("length") and len(r[2]) != m["length"]:
                short += 1
                LOG.warning("%s: got %d bp, docsum says %d bp; not written (will retry)",
                            r[0], len(r[2]), m["length"])
                continue
            batch.append((r[0], r[1], r[2], [req]))
        n = state.write_batch(batch, "ncbi_efetch")
        got += n
        LOG.info("Sequence batch %d/%d: %d/%d records, %.1f Mbp in %.1fs", gi, len(groups), n,
                 len(grp), sum(len(x[2]) for x in batch) / 1e6, time.time() - t0)
    return got


def _elink_split(eu, chunk, depth=0):
    """elink a chunk; on failure split it in half until single ids.
    Returns ({nuccore uid: [assembly uids]}, [failed (acc, uid)], n_splits)."""
    try:
        return eu.elink_assembly([u for _, u in chunk]), [], 0
    except Exception as e:  # noqa: BLE001
        if len(chunk) == 1:
            LOG.warning("elink failed for %s (uid %s): %s", chunk[0][0], chunk[0][1], e)
            return {}, list(chunk), 0
        LOG.info("elink batch of %d failed (%s); splitting", len(chunk),
                 str(e).split(":")[-1].strip())
        mid = len(chunk) // 2
        l1, f1, s1 = _elink_split(eu, chunk[:mid], depth + 1)
        l2, f2, s2 = _elink_split(eu, chunk[mid:], depth + 1)
        l1.update(l2)
        return l1, f1 + f2, s1 + s2 + 1


def ncbi_assemblies(eu, records, state, batch_size):
    """records: list of (record acc, uid). nuccore->assembly via elink.

    Batches are small (--elink-batch-size) and a batch NCBI drops is split in half
    recursively, so one slow id no longer costs the whole batch. Ids that still
    fail are not cached and are retried on the next run."""
    items = [(a, u) for a, u in records if u]
    nb = (len(items) + batch_size - 1) // batch_size
    all_failed = []
    for b in range(nb):
        chunk = items[b * batch_size:(b + 1) * batch_size]
        links, failed, splits = _elink_split(eu, chunk)
        asm_uids = sorted({x for v in links.values() for x in v})
        asm = {}
        try:
            for i in range(0, len(asm_uids), 200):
                for ds in eu.esummary("assembly", asm_uids[i:i + 200]):
                    asm[str(ds.get("uid"))] = ds.get("assemblyaccession", "")
        except Exception as e:  # noqa: BLE001
            LOG.error("Assembly batch %d/%d: assembly esummary failed: %s (will retry on the "
                      "next run)", b + 1, nb, e)
            all_failed.extend(chunk)
            continue
        failed_uids = {u for _, u in failed}
        out = []
        for acc, uid in chunk:
            if uid in failed_uids:
                continue
            accs = [asm[x] for x in links.get(uid, []) if asm.get(x)]
            pref = "GCF_" if re.match(r"^(N[CZTGW]|AC|NW)_", acc) else "GCA_"
            accs.sort(key=lambda a: (not a.startswith(pref), a))
            out.append(dict(acc=acc, assembly=accs[0] if accs else "", source="ncbi_elink",
                            n_linked=len(accs)))
        state.put_assembly(out)
        all_failed.extend(failed)
        LOG.info("Assembly batch %d/%d: %d/%d linked to an assembly%s%s", b + 1, nb,
                 sum(1 for d in out if d["assembly"]), len(chunk),
                 f", {splits} split(s)" if splits else "",
                 f", {len(failed)} failed" if failed else "")
    if all_failed:
        LOG.warning("Assembly lookup failed for %d record(s) (%s%s); the FASTA and taxids are "
                    "unaffected, only their Assembly column is empty. Re-run to retry, or "
                    "pass --no-elink to skip.", len(all_failed),
                    ", ".join(a for a, _ in all_failed[:5]), ", ..." if len(all_failed) > 5 else "")


# ─────────────────────────────────────────────────────────────────────────────
# Finalise
# ─────────────────────────────────────────────────────────────────────────────
def finalize_fasta(state, keep_accs):
    """Drop records no longer requested (list edited between runs) and write .fai."""
    extras = [a for a in state.records if a not in keep_accs]
    if extras:
        LOG.info("Removing %d record(s) no longer in the input list from the FASTA", len(extras))
        tmp = state.fasta + ".tmp"
        kept = []
        with open(state.fasta, "rb") as src, open(tmp, "wb") as dst:
            for acc, e in state.records.items():
                if acc not in keep_accs:
                    continue
                src.seek(e["hdr_off"])
                blob = src.read(e["end"] - e["hdr_off"])
                hdr_off = dst.tell()
                dst.write(blob)
                ne = dict(e, hdr_off=hdr_off, seq_off=hdr_off + (e["seq_off"] - e["hdr_off"]),
                          end=dst.tell())
                kept.append(ne)
            dst.flush()
            os.fsync(dst.fileno())
        os.replace(tmp, state.fasta)
        state._rewrite_ckpt(kept)
        state.records = OrderedDict((e["acc"], e) for e in kept)
    with open(state.fasta + ".fai", "w") as fh:
        for acc, e in state.records.items():
            fh.write(f"{acc}\t{e['length']}\t{e['seq_off']}\t{WRAP}\t{WRAP + 1}\n")


def write_tsv(path, header, rows):
    tmp = path + ".tmp"
    with open(tmp, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        if header:
            w.writerow(header)
        w.writerows(rows)
    os.replace(tmp, path)


# ─────────────────────────────────────────────────────────────────────────────
def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Pull a list of contig accessions into a FASTA plus the accession->taxid "
                    "mapping TaxTriage needs (--reference_fasta / --custom_accession_map).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split("Examples", 1)[1] if "Examples" in __doc__ else None)
    p.add_argument("-i", "--input", required=True,
                   help="File with one contig accession per line (first column used; "
                        "'#' lines and a header row are ignored; .gz ok)")
    p.add_argument("-o", "--outdir", default=".", help="Output directory (default: .)")
    p.add_argument("-p", "--prefix", default=None,
                   help="Output prefix (default: input file name without extension)")
    p.add_argument("--source-fasta", nargs="+", default=[],
                   help="Local FASTA(s) to take sequences from before NCBI (plain + .fai "
                        "for random access, otherwise streamed; gzip ok)")
    p.add_argument("--source-map", nargs="+", default=[],
                   help="Existing mapping file(s) with accession/taxid/assembly/name for "
                        "offline metadata (pipeline mapping.taxid.tsv, gcfmapping, or "
                        "custom accession map)")
    p.add_argument("-a", "--assembly", nargs="+", default=[],
                   help="assembly_summary_refseq.txt / _genbank.txt for assembly-level "
                        "taxid / organism / species_taxid")
    p.add_argument("-t", "--taxdump", default=None,
                   help="Directory with names.dmp + nodes.dmp (+ merged.dmp, delnodes.dmp): "
                        "the taxonomy the pipeline runs with")
    p.add_argument("--map-rank", choices=["taxid", "species"], default="taxid",
                   help="Taxid written to the accession map: the record's own taxid "
                        "(default) or its species-level ancestor (needs --taxdump)")
    p.add_argument("--no-map-names", action="store_true",
                   help="Leave the 'name' column out of the accession map (the pipeline "
                        "then names organisms from its own taxonomy)")
    p.add_argument("-e", "--email", default=None, help="Email for NCBI E-utilities")
    p.add_argument("-k", "--api-key", default=None,
                   help="NCBI API key (or $NCBI_API_KEY): raises the rate limit")
    p.add_argument("--offline", action="store_true", help="Never contact NCBI")
    p.add_argument("--no-elink", action="store_true",
                   help="Skip the nuccore->assembly elink lookup for contigs whose docsum "
                        "has no assembly")
    p.add_argument("--batch-size", type=int, default=200,
                   help="Accessions per esearch/esummary metadata request (default 200)")
    p.add_argument("--elink-batch-size", type=int, default=20,
                   help="Records per nuccore->assembly elink request (default 20). Each id "
                        "is a separate link computation at NCBI; large batches time out "
                        "and are split automatically")
    p.add_argument("--max-bp", type=float, default=50e6,
                   help="Max total bp per efetch request / local write batch (default 50e6)")
    p.add_argument("--max-ids", type=int, default=500,
                   help="Max ids per efetch request (default 500)")
    p.add_argument("--restart", action="store_true",
                   help="Discard previous output/state for this prefix and start over")
    p.add_argument("--strict", action="store_true",
                   help="Exit 1 if any requested contig is missing or has no taxid")
    p.add_argument("-v", "--verbose", action="store_true")
    return p.parse_args(argv)


def setup_logging(path, verbose):
    LOG.setLevel(logging.DEBUG if verbose else logging.INFO)
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%Y-%m-%d %H:%M:%S")
    fh = logging.FileHandler(path)
    fh.setFormatter(fmt)
    sh = logging.StreamHandler(sys.stderr)
    sh.setFormatter(fmt)
    LOG.handlers[:] = [fh, sh]


def main(argv=None):
    args = parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)
    prefix = args.prefix or re.sub(r"(\.(txt|tsv|csv|list|lst))?(\.gz)?$", "",
                                   os.path.basename(args.input))
    out = lambda ext: os.path.join(args.outdir, f"{prefix}{ext}")  # noqa: E731
    setup_logging(out(".log"), args.verbose)
    LOG.info("=" * 70)
    LOG.info("pull_contigs.py %s", " ".join(sys.argv[1:] if argv is None else argv))
    if args.map_rank == "species" and not args.taxdump:
        raise SystemExit("ERROR: --map-rank species needs --taxdump")
    if args.offline and not (args.source_fasta):
        LOG.warning("--offline without --source-fasta: no sequences can be obtained")
    t_start = time.time()

    requested, dups = read_contig_list(args.input)
    if not requested:
        raise SystemExit(f"ERROR: no accessions in {args.input}")
    LOG.info("Input: %d unique accession(s)%s", len(requested),
             f" ({dups} duplicate line(s) ignored)" if dups else "")

    state = State(args.outdir, prefix, restart=args.restart)
    smap = load_source_maps(args.source_map)
    eu = None if args.offline else Eutils(args.email, args.api_key)
    if eu and not args.email and not eu.api_key:
        LOG.info("Tip: pass --email / --api-key for faster, NCBI-compliant requests")

    # 1) sequences from local FASTA(s)
    pending = [r for r in requested if r not in state.by_requested]
    if pending and args.source_fasta:
        n = scan_local_fastas(args.source_fasta, Wanted(pending), state, args.max_bp,
                              allow_other_version=args.offline)
        LOG.info("Local sources: %d new record(s)", n)

    # 2) metadata: source map first, NCBI for the rest
    def lookup_id(req):
        acc = state.by_requested.get(req)
        return acc or req

    def smap_get(acc):
        return smap.get(acc) or smap.get(base_acc(acc))

    need_meta = [r for r in requested if r not in state.meta]
    from_map, need_ncbi = [], {}
    for r in need_meta:
        q = lookup_id(r)
        sm = smap_get(q)
        rec = state.records.get(q)
        # A source-map taxid saves the NCBI round trip only when the sequence is
        # already on disk (or we are offline); otherwise NCBI is needed for the uid
        # to fetch the sequence, and the source-map taxid still wins at finalise.
        if sm and sm["taxid"] and (rec or not eu):
            from_map.append(dict(requested=r, acc=q if rec else r, uid="", taxid=sm["taxid"],
                                 organism=sm["name"], title=sm["description"],
                                 length=rec["length"] if rec else 0, strain="",
                                 assembly=sm["assembly"], moltype="", completeness="",
                                 source=f"source_map:{sm['source']}", note=""))
        else:
            need_ncbi[r] = q
    if from_map:
        state.put_meta(from_map)
        LOG.info("Metadata: %d from source map(s)", len(from_map))
    if need_ncbi and eu:
        LOG.info("Metadata: querying NCBI for %d accession(s)", len(need_ncbi))
        ncbi_metadata(eu, need_ncbi, state, args.batch_size)
    elif need_ncbi:
        LOG.info("Metadata: %d accession(s) have no source-map taxid (offline)", len(need_ncbi))

    # 3) sequences from NCBI
    todo = [(r, state.meta[r]) for r in requested
            if r not in state.by_requested and r in state.meta and state.meta[r].get("uid")]
    # requested ids that resolve to a record already on disk (X and X.2 both listed)
    for r, m in list(todo):
        if m["acc"] in state.records:
            state.write_batch([(m["acc"], "", "", [r])], "dup")
    todo = [(r, m) for r, m in todo if r not in state.by_requested]
    if todo and eu:
        LOG.info("Sequences: fetching %d record(s) (%.1f Mbp) from NCBI", len(todo),
                 sum(m.get("length", 0) for _, m in todo) / 1e6)
        ncbi_sequences(eu, todo, state, args.max_bp, args.max_ids)

    # 4) assemblies
    def assembly_for(req):
        m = state.meta.get(req, {})
        acc = state.by_requested.get(req, m.get("acc", ""))
        if m.get("assembly"):
            return m["assembly"], m["source"].split(":")[0]
        sm = smap_get(acc) if acc else None
        if sm and sm["assembly"]:
            return sm["assembly"], "source_map"
        a = state.assembly.get(acc)
        if a and a.get("assembly"):
            return a["assembly"], a["source"]
        return "", ""

    if eu and not args.no_elink:
        need = []
        for r in requested:
            acc = state.by_requested.get(r)
            if not acc or acc in state.assembly or assembly_for(r)[0]:
                continue
            uid = state.meta.get(r, {}).get("uid")
            if uid:
                need.append((acc, uid))
        need = list(dict.fromkeys(need))
        if need:
            LOG.info("Assemblies: elink for %d record(s) without one in their docsum", len(need))
            ncbi_assemblies(eu, need, state, args.elink_batch_size)

    # 5) finalise FASTA + taxonomy
    keep = {state.by_requested[r] for r in requested if r in state.by_requested}
    finalize_fasta(state, keep)

    asm_ids = {assembly_for(r)[0] for r in requested} - {""}
    asm_info = load_assembly_summaries(args.assembly, asm_ids) if args.assembly else {}
    if args.assembly:
        LOG.info("Assembly summaries: %d/%d assemblies found", len({k for k in asm_info if k in asm_ids}),
                 len(asm_ids))

    tax = Taxdump(args.taxdump) if args.taxdump else None

    rows, missing, by_record = [], [], OrderedDict()
    for r in requested:
        acc = state.by_requested.get(r)
        m = state.meta.get(r, {})
        if not acc:
            if not m:
                reason = ("not_in_source_fasta(offline)" if args.offline
                          else "not_found_at_ncbi")
            else:
                reason = "sequence_not_retrieved(rerun to retry)" if eu else \
                         "metadata_only(offline; not in source FASTA)"
            missing.append((r, reason))
            continue
        if acc in by_record:
            by_record[acc]["requested"].append(r)
            continue
        rec = state.records[acc]
        sm = smap_get(acc) or {}
        asm, asm_src = assembly_for(r)
        ai = asm_info.get(asm) or asm_info.get(base_acc(asm)) or {}
        # curated source map > NCBI docsum > assembly summary
        if sm.get("taxid"):
            orig_taxid, tax_src = sm["taxid"], "source_map"
        elif m.get("taxid"):
            orig_taxid, tax_src = m["taxid"], m.get("source", "").split(":")[0]
        else:
            orig_taxid, tax_src = "", ""
        if not orig_taxid and ai.get("taxid"):
            orig_taxid, tax_src = _clean_taxid(ai["taxid"]), "assembly_summary"
        desc = m.get("title") or sm.get("description") or ""
        if not desc:
            with open(state.fasta, "rb") as fh:
                fh.seek(rec["hdr_off"])
                hdr = fh.readline().decode(errors="replace").strip()[1:].split(None, 1)
                desc = hdr[1] if len(hdr) > 1 else ""
        notes = [m["note"]] if m.get("note") else []
        if has_version(r) and acc != r and not any(n.startswith("version_changed") for n in notes):
            notes.append(f"version_changed:{r}->{acc}")
        by_record[acc] = dict(requested=[r], acc=acc, assembly=asm, assembly_source=asm_src,
                              orig_taxid=orig_taxid, taxid_source=tax_src,
                              organism=m.get("organism") or sm.get("name") or ai.get("organism_name", ""),
                              description=desc, length=rec["length"], seq_source=rec["source"],
                              asm_organism=ai.get("organism_name", ""),
                              asm_species_taxid=_clean_taxid(ai.get("species_taxid")),
                              refseq_category=ai.get("refseq_category", ""),
                              assembly_level=ai.get("assembly_level", ""), notes=notes)

    if tax:
        cur = {}
        for d in by_record.values():
            t, st = tax.resolve(d["orig_taxid"]) if d["orig_taxid"] else (None, "missing")
            sp = tax.ancestor_at(t, "species") if t else None
            cur[d["acc"]] = (t, st, sp)
        tax.load_names({x for v in cur.values() for x in (v[0], v[2]) if x})
    for d in by_record.values():
        if tax:
            t, st, sp = cur[d["acc"]]
            d.update(taxid=str(t or ""), taxid_status=st, rank=tax.rank_of(t) if t else "",
                     species_taxid=str(sp or ""), species=tax.name(sp) if sp else "",
                     sci_name=tax.name(t) if t else "")
            if st == "merged":
                d["notes"].append(f"taxid_merged:{d['orig_taxid']}->{t}")
            elif st in ("deleted", "not_in_taxdump"):
                d["notes"].append(f"taxid_{st}:{d['orig_taxid']}")
        else:
            d.update(taxid=d["orig_taxid"], taxid_status="ok" if d["orig_taxid"] else "missing",
                     rank="", species_taxid=d["asm_species_taxid"], species="",
                     sci_name="")
        d["name"] = d["sci_name"] or d["organism"] or d["asm_organism"] or d["description"]
        if args.map_rank == "species":
            d["map_taxid"] = d["species_taxid"] or d["taxid"]
            d["map_name"] = d["species"] or d["name"]
        else:
            d["map_taxid"], d["map_name"] = d["taxid"], d["name"]

    # outputs
    recs = list(by_record.values())
    if args.no_map_names:
        write_tsv(out(".accession_map.tsv"), ["accession", "taxid"],
                  [[d["acc"], d["map_taxid"]] for d in recs if d["map_taxid"]])
    else:
        write_tsv(out(".accession_map.tsv"), ["accession", "taxid", "name"],
                  [[d["acc"], d["map_taxid"], d["map_name"]] for d in recs if d["map_taxid"]])
    write_tsv(out(".mapping.taxid.tsv"),
              ["Acc", "Assembly", "Organism_Name", "Description", "Mapped_Value"],
              [[d["acc"], d["assembly"], d["map_name"], d["description"], d["map_taxid"]]
               for d in recs])
    write_tsv(out(".metadata.tsv"),
              ["requested", "accession", "assembly", "assembly_source", "taxid", "taxid_status",
               "original_taxid", "taxid_source", "rank", "name", "species_taxid", "species",
               "ncbi_organism", "assembly_organism", "refseq_category", "assembly_level",
               "length", "sequence_source", "description", "notes"],
              [[",".join(d["requested"]), d["acc"], d["assembly"], d["assembly_source"], d["taxid"],
                d["taxid_status"], d["orig_taxid"], d["taxid_source"], d["rank"], d["name"],
                d["species_taxid"], d["species"], d["organism"], d["asm_organism"],
                d["refseq_category"], d["assembly_level"], d["length"], d["seq_source"],
                d["description"], ";".join(d["notes"])] for d in recs])
    with open(out(".assemblies.txt"), "w") as fh:
        for a in sorted({d["assembly"] for d in recs if d["assembly"]}):
            fh.write(a + "\n")
    no_taxid = [d for d in recs if not d["map_taxid"]]
    write_tsv(out(".missing.tsv"), ["requested", "reason"],
              missing + [(",".join(d["requested"]), f"no_taxid({d['taxid_status']})")
                         for d in no_taxid])

    summ = OrderedDict(
        input=os.path.abspath(args.input), requested=len(requested),
        duplicate_lines=dups, sequences_written=len(recs),
        total_bp=sum(d["length"] for d in recs), missing_sequence=len(missing),
        with_taxid=len(recs) - len(no_taxid), without_taxid=len(no_taxid),
        with_assembly=sum(1 for d in recs if d["assembly"]),
        taxid_merged=sum(1 for d in recs if d["taxid_status"] == "merged"),
        taxid_deleted_or_unknown=sum(1 for d in recs if d["taxid_status"] in ("deleted", "not_in_taxdump")),
        version_changed=sum(1 for d in recs if any(n.startswith("version_changed") for n in d["notes"])),
        distinct_taxids=len({d["map_taxid"] for d in recs if d["map_taxid"]}),
        sequence_sources=dict(sorted(
            ((s, sum(1 for d in recs if d["seq_source"] == s)) for s in {d["seq_source"] for d in recs}))),
        ncbi_requests=eu.requests if eu else 0, taxdump=args.taxdump, map_rank=args.map_rank,
        seconds=round(time.time() - t_start, 1),
        complete=not missing and not no_taxid)
    with open(out(".summary.json"), "w") as fh:
        json.dump(summ, fh, indent=2)

    LOG.info("-" * 70)
    LOG.info("Requested %d | written %d (%.1f Mbp) | missing %d | without taxid %d | "
             "with assembly %d | distinct taxids %d", summ["requested"], summ["sequences_written"],
             summ["total_bp"] / 1e6, summ["missing_sequence"], summ["without_taxid"],
             summ["with_assembly"], summ["distinct_taxids"])
    if summ["taxid_merged"] or summ["taxid_deleted_or_unknown"]:
        LOG.info("Taxdump: %d merged taxid(s) remapped, %d deleted/unknown",
                 summ["taxid_merged"], summ["taxid_deleted_or_unknown"])
    if summ["version_changed"]:
        LOG.warning("%d contig(s) delivered at a different version than requested "
                    "(see notes in %s)", summ["version_changed"], os.path.basename(out(".metadata.tsv")))
    if missing or no_taxid:
        LOG.warning("See %s for the %d problem accession(s)%s",
                    os.path.basename(out(".missing.tsv")), len(missing) + len(no_taxid),
                    "; re-run the same command to retry" if eu else "")
    LOG.info("FASTA:          %s", out(".fasta"))
    LOG.info("Accession map:  %s", out(".accession_map.tsv"))
    LOG.info("Use with: --reference_fasta %s --custom_accession_map %s",
             out(".fasta"), out(".accession_map.tsv"))
    if args.strict and not summ["complete"]:
        sys.exit(1)


if __name__ == "__main__":
    main()
