#!/usr/bin/env python3
"""
Shared, dependency-free reader for --custom_accession_map.

One file, used everywhere an accession is turned into a taxid:
  * MAP_TAXID_ASSEMBLY (append_taxid.py)   -> the general pipeline + insilico
  * FETCH_SPIKEIN_REFS (resolve_fasta_taxids.py) -> spike-in series

Format: delimited text (tab, comma or semicolon — sniffed), header required.
  accession   accession | acc | accession_version | seqid | id | nuccore
  taxid       taxid | tax_id | taxonomy_id | ncbi_taxid
  name        name | organism | organism_name        (optional)
Lines starting with '#' are ignored.

Keys are matched exactly first, then without the version suffix, so a map
written as OR833055 also covers OR833055.1 (and vice versa).
"""
import csv
import re
import sys

ACC_COLS = ["accession", "acc", "accession_version", "seqid", "seq_id", "id", "nuccore", "contig"]
TAX_COLS = ["taxid", "tax_id", "taxonomy_id", "ncbi_taxid", "mapped_value"]
NAME_COLS = ["name", "organism", "organism_name", "species"]


def _norm(h):
    return re.sub(r"[^a-z0-9]+", "_", str(h or "").strip().lower()).strip("_")


def _clean_taxid(v):
    t = str(v or "").strip()
    t = re.sub(r"\.0+$", "", t)
    return t if re.fullmatch(r"\d+", t) and t != "0" else ""


def _base(acc):
    return re.sub(r"\.\d+$", "", acc)


class CustomAccessionMap:
    def __init__(self, entries=None):
        self.exact = {}   # acc -> (taxid, name)
        self.base = {}    # acc without version -> (taxid, name)
        for acc, tid, name in entries or []:
            self.exact[acc] = (tid, name)
            self.base.setdefault(_base(acc), (tid, name))

    def __len__(self):
        return len(self.exact)

    def get(self, acc):
        """(taxid, name) for an accession, or None."""
        if not acc:
            return None
        acc = str(acc).strip().lstrip(">")
        return self.exact.get(acc) or self.base.get(_base(acc))

    @classmethod
    def load(cls, path, log=sys.stderr):
        import os
        if not path or os.path.basename(str(path)).startswith("NO_FILE"):
            return cls()   # the pipeline's "not provided" placeholder
        try:
            with open(path, newline="", encoding="utf-8-sig") as fh:
                lines = [l for l in fh if l.strip() and not l.lstrip().startswith("#")]
        except OSError as e:
            print(f"[custom_map] WARNING: cannot read {path}: {e}; ignoring it.", file=log)
            return cls()
        if not lines:
            return cls()   # NO_FILE placeholder / empty map
        try:
            delim = csv.Sniffer().sniff(lines[0], delimiters="\t,;").delimiter
        except Exception:
            delim = "\t"
        rows = list(csv.reader(lines, delimiter=delim))
        header = [_norm(h) for h in rows[0]]

        def col(cands):
            for c in cands:
                if c in header:
                    return header.index(c)
            return None

        ai, ti, ni = col(ACC_COLS), col(TAX_COLS), col(NAME_COLS)
        if ai is None or ti is None:
            raise SystemExit(
                f"ERROR: --custom_accession_map {path} needs a header with an accession column "
                f"({' | '.join(ACC_COLS[:4])}) and a taxid column ({' | '.join(TAX_COLS[:4])}). "
                f"Found: {', '.join(rows[0])}"
            )
        entries, bad = [], 0
        for r in rows[1:]:
            acc = r[ai].strip().lstrip(">") if ai < len(r) else ""
            tid = _clean_taxid(r[ti]) if ti < len(r) else ""
            name = r[ni].strip() if ni is not None and ni < len(r) else ""
            if acc and tid:
                entries.append((acc, tid, name))
            elif acc:
                bad += 1
        if bad:
            print(f"[custom_map] WARNING: {bad} row(s) in {path} have no valid taxid; skipped.", file=log)
        print(f"[custom_map] loaded {len(entries)} accession -> taxid entries from {path}", file=log)
        return cls(entries)
