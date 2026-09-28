#!/usr/bin/env python3
"""
Resolve the taxid of a local spike-in FASTA from the accessions in its headers.

A spike-in row is scored in the report against ONE taxid, so the FASTA must hold
one organism (any number of contigs / chromosomes / plasmids). Each record's
leading id is looked up with NCBI esummary (nuccore):

Output: <accession> <taxid> <organism> <match keys>

A taxid is NOT required. Each record gets a MATCH KEY: its taxid when one is
known, else its own accession — the same fallback match_paths.py uses to key a
reference with no taxid, so the report pairs the spike with that detection.
Match keys are written as  key|n_records|bp;key|n_records|bp  so the mixing
step can split a row's reads across them (ISS: per record; NanoSim: per bp).

  * all records -> one taxid  : taxid column = that taxid, one match key
  * several taxids / some or  : taxid column empty, one match key per organism
    all records unresolved      (or per unresolved accession), with a warning

Precedence: sheet taxid (--taxid) > --custom-map (--custom_accession_map; keyed
by the sheet accession or by each record id) > NCBI esummary. With
--only-if-mapped the script is a pure custom-map override: it rewrites --output
only when the map resolves the organism, and leaves it untouched otherwise
(used for accession rows whose taxid NCBI already supplied).
"""
import argparse
import gzip
import json
import os
import sys
import time
import urllib.parse
import urllib.request
from collections import OrderedDict

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from custom_accession_map import CustomAccessionMap

ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"


def headers(path):
    """[(id, description, length_bp)] for every record."""
    op = gzip.open if path.endswith(".gz") else open
    out = []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                parts = line[1:].strip().split(None, 1)
                out.append([parts[0] if parts else "", (parts[1] if len(parts) > 1 else ""), 0])
            elif out:
                out[-1][2] += len(line.strip())
    return [tuple(r) for r in out]


def keys_str(pairs):
    """[(key, id, bp)] -> 'key|n|bp;...' aggregated per key, first-seen order."""
    agg = OrderedDict()
    for key, _rid, bp in pairs:
        n, b = agg.get(key, (0, 0))
        agg[key] = (n + 1, b + bp)
    return ";".join(f"{k}|{n}|{b}" for k, (n, b) in agg.items())


def esummary(ids, api_key=None, batch=100):
    """{id: (taxid, organism)} for the ids NCBI recognises."""
    out = {}
    for i in range(0, len(ids), batch):
        chunk = ids[i:i + batch]
        q = {"db": "nuccore", "id": ",".join(chunk), "retmode": "json"}
        if api_key:
            q["api_key"] = api_key
        url = ESUMMARY + "?" + urllib.parse.urlencode(q)
        data = None
        for attempt in range(3):
            try:
                with urllib.request.urlopen(url, timeout=60) as r:
                    data = json.load(r)
                break
            except Exception as e:  # network / rate limit
                print(f"[resolve_taxids] esummary attempt {attempt + 1} failed: {e}", file=sys.stderr)
                time.sleep(2 * (attempt + 1))
        if not data:
            continue
        res = data.get("result", {})
        for uid in res.get("uids", []):
            d = res.get(uid, {})
            tid = str(d.get("taxid") or "").strip()
            acc = d.get("accessionversion") or d.get("caption") or uid
            if tid and tid != "0":
                out[acc] = (tid, d.get("organism") or "")
                out.setdefault(d.get("caption") or acc, out[acc])
        time.sleep(0.12 if api_key else 0.4)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--accession", required=True, help="Spike-in id this FASTA is keyed by")
    ap.add_argument("--taxid", default="", help="Taxid from the sheet; skips the lookup")
    ap.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY", ""))
    ap.add_argument("--custom-map", default=None, help="--custom_accession_map file")
    ap.add_argument("--no-ncbi", action="store_true", help="Never query NCBI")
    ap.add_argument("--only-if-mapped", action="store_true",
                    help="Write --output only when the custom map resolves the organism")
    ap.add_argument("--output", required=True)
    a = ap.parse_args()

    recs = headers(a.fasta)
    # Label from the header only when the file is one record; a multi-record file
    # (several organisms, or chromosomes named per-chromosome) is labelled by its id.
    first_desc = recs[0][1].split(",")[0] if len(recs) == 1 else a.accession

    def whole_keys(tid):
        return keys_str([(tid, rid, bp) for rid, _d, bp in recs])

    if a.taxid:
        with open(a.output, "w") as fh:
            fh.write(f"{a.accession}\t{a.taxid}\t{first_desc}\t{whole_keys(a.taxid)}\n")
        print(f"[resolve_taxids] {a.accession}: taxid {a.taxid} from the sheet", file=sys.stderr)
        return

    cmap = CustomAccessionMap.load(a.custom_map)

    # The sheet accession itself may be in the map (e.g. a GCF_ id, or the id of a
    # custom genome): that names the organism outright.
    whole = cmap.get(a.accession)
    if whole:
        with open(a.output, "w") as fh:
            fh.write(f"{a.accession}\t{whole[0]}\t{whole[1] or first_desc}\t{whole_keys(whole[0])}\n")
        print(f"[resolve_taxids] {a.accession}: taxid {whole[0]} from the custom accession map",
              file=sys.stderr)
        return

    ids = [rid for rid, _d, _b in recs]
    found = {}
    for rid in ids:
        hit = cmap.get(rid)
        if hit:
            found[rid] = hit
    n_map = len(found)
    prior = None   # --only-if-mapped: the taxid NCBI already gave the whole row
    if a.only_if_mapped:
        if not n_map:
            return
        try:
            with open(a.output) as fh:
                f = fh.readline().rstrip("\n").split("\t")
            if len(f) > 2 and f[1].strip():
                prior = (f[1].strip(), f[2].strip())
        except OSError:
            pass
        if n_map < len(ids):
            print(f"[resolve_taxids] {a.accession}: custom map covers {n_map}/{len(ids)} records; "
                  f"the rest keep {'taxid ' + prior[0] if prior else 'their accession'}", file=sys.stderr)
        for rid in ids:
            if rid not in found and prior:
                found[rid] = prior
    rest = [rid for rid in ids if rid not in found]
    if rest and not a.no_ncbi:
        found.update(esummary(rest, a.api_key or None))
    if n_map:
        print(f"[resolve_taxids] {a.accession}: {n_map}/{len(ids)} record(s) resolved by the custom map",
              file=sys.stderr)

    # esummary may key on the version-less caption: match either form
    def look(rid):
        return found.get(rid) or found.get(rid.split(".")[0])

    groups = OrderedDict()          # taxid -> {"org": str, "ids": [..]}
    unresolved = []
    pairs = []                      # (match key, record id, bp)
    for rid, desc, bp in recs:
        hit = look(rid)
        if not hit:
            unresolved.append(rid)
            pairs.append((rid, rid, bp))       # no taxid -> keyed by its accession
            continue
        g = groups.setdefault(hit[0], {"org": hit[1], "ids": []})
        g["ids"].append(rid)
        pairs.append((hit[0], rid, bp))

    if unresolved:
        print(f"[resolve_taxids] {a.accession}: {len(unresolved)}/{len(recs)} record(s) have no taxid "
              "and are reported by accession: "
              + ", ".join(unresolved[:10]) + (" ..." if len(unresolved) > 10 else ""), file=sys.stderr)

    if len(groups) == 1 and not unresolved:
        tid, g = next(iter(groups.items()))
        org = g["org"] or first_desc
    else:
        tid, org = "", first_desc
        if len(groups) + len(unresolved) > 1:
            print(f"[resolve_taxids] NOTE: '{a.accession}' holds {len(groups) + len(unresolved)} organisms/"
                  "accessions; its spiked reads are split across them in the report. To spike one,\n"
                  "                 add a `record` column (e.g. record=" + (recs[0][0] if recs else "ID") +
                  ") or give each organism its own row.", file=sys.stderr)
    with open(a.output, "w") as fh:
        fh.write(f"{a.accession}\t{tid}\t{org}\t{keys_str(pairs)}\n")
    print(f"[resolve_taxids] {a.accession}: taxid {tid or '-'} ({org}); match keys {keys_str(pairs)}",
          file=sys.stderr)

if __name__ == "__main__":
    main()
