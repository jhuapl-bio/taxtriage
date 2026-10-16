#!/usr/bin/env python3
# ##############################################################################################
# # Copyright 2022 The Johns Hopkins University Applied Physics Laboratory LLC
# # All rights reserved.
# # Permission is hereby granted, free of charge, to any person obtaining a copy of this
# # software and associated documentation files (the "Software"), to deal in the Software
# # without restriction, including without limitation the rights to use, copy, modify,
# # merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
# # permit persons to whom the Software is furnished to do so.
# #
# # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# # INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
# # PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# # LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# # TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
# # OR OTHER DEALINGS IN THE SOFTWARE.
# #
"""
Reference-FASTA helpers shared by the TaxTriage Python scripts.

Why this exists: ``pysam.FastaFile`` (htslib faidx) can open a plain FASTA or a
**bgzip**-compressed one, but not an ordinary ``gzip`` file -- which is what
NCBI ships and what most people produce. Opening one fails with::

    OSError: error when opening file `x.fasta.gz`

``open_pysam_fasta()`` returns an object with the pysam.FastaFile API the
scripts use (``references``, ``lengths``, ``fetch(ref[, start, end])``,
``get_reference_length``, ``close``) for any accepted reference:

  plain FASTA           pysam.FastaFile (a .fai is built next to it, or next to
                        a symlink in a cache dir if that directory is read-only)
  bgzip FASTA (.gz)     pysam.FastaFile when htslib can index it
  gzip FASTA (.gz)      GzFastaReader: READ DIRECTLY FROM THE .gz, nothing is
                        decompressed to disk. One streaming pass records each
                        contig's name, length and uncompressed offset (cached
                        as a few-KB JSON sidecar so worker processes skip it);
                        fetch() then streams to that contig and keeps the last
                        contigs in memory, so whole-contig and in-order region
                        access -- what TaxTriage does -- costs ~1 pass.

Set TAXTRIAGE_GZ_FASTA=decompress to fall back to the old behaviour (one
plain copy ``<name>.taxtriage.fasta`` + .fai) for heavily random access to
very large gzip references.

Accepted extensions (case-insensitive), each optionally followed by .gz or
.bgz: .fa .fasta .fna .fas .fsa .ffn .frn .mfa (see FASTA_EXTENSIONS).
"""
from __future__ import annotations

import gzip
import os
import re
import shutil
import sys
import tempfile

FASTA_EXTENSIONS = (".fa", ".fasta", ".fna", ".fas", ".fsa", ".ffn", ".frn", ".mfa")
COMPRESSED_SUFFIXES = (".gz", ".bgz")
_FASTA_RE = re.compile(r"\.(fa|fasta|fna|fas|fsa|ffn|frn|mfa)(\.(gz|bgz))?$", re.IGNORECASE)

_resolved: dict = {}


def is_fasta_name(path: str) -> bool:
    """True when the file name carries an accepted FASTA extension."""
    return bool(_FASTA_RE.search(os.path.basename(str(path))))


def check_fasta_path(path: str, what: str = "reference FASTA") -> str:
    """Raise a clear error unless ``path`` exists and looks like a FASTA."""
    p = str(path)
    if not is_fasta_name(p):
        raise ValueError(
            f"{what} '{p}' does not have a FASTA extension. Accepted: "
            f"{', '.join(FASTA_EXTENSIONS)} (each optionally followed by .gz or .bgz).")
    if not os.path.exists(p):
        raise FileNotFoundError(f"{what} '{p}' not found")
    return p


def is_gzip(path: str) -> bool:
    with open(path, "rb") as fh:
        return fh.read(2) == b"\x1f\x8b"


def is_bgzf(path: str) -> bool:
    """BGZF = gzip member with the 'BC' extra subfield (what bgzip writes)."""
    with open(path, "rb") as fh:
        head = fh.read(18)
    return (len(head) >= 18 and head[:2] == b"\x1f\x8b" and head[3] & 4
            and head[12:14] == b"BC")


def open_fasta_text(path: str):
    """Open a FASTA for text reading, gzip/bgzip or plain."""
    return gzip.open(path, "rt") if is_gzip(path) else open(path, "r")


def _cache_dir() -> str:
    d = os.environ.get("TAXTRIAGE_FASTA_CACHE") or os.path.join(tempfile.gettempdir(), "taxtriage_fasta_cache")
    os.makedirs(d, exist_ok=True)
    return d


def _writable_dir(path: str) -> str:
    # Use the directory of the path AS GIVEN (not resolved): in a Nextflow task
    # that is the work dir, where we may write, even when the input is a
    # symlink into a read-only store.
    d = os.path.dirname(os.path.abspath(path)) or "."
    return d if os.access(d, os.W_OK) else _cache_dir()


class _Lock:
    """Exclusive inter-process lock (fcntl) around decompress + .fai creation.

    pysam builds a missing .fai itself on open; when several worker processes
    open the same new FASTA at once, one can read another's half-written index
    and see zero-length contigs. Doing the first open under a lock avoids that.
    """

    def __init__(self, target: str):
        d = os.path.dirname(target) or "."
        if not os.access(d, os.W_OK):
            d = _cache_dir()
        self.path = os.path.join(d, os.path.basename(target) + ".taxtriage.lock")
        self.fh = None

    def __enter__(self):
        self.fh = open(self.path, "a")
        try:
            import fcntl
            fcntl.flock(self.fh, fcntl.LOCK_EX)
        except (ImportError, OSError):
            pass
        return self

    def __exit__(self, *exc):
        try:
            import fcntl
            fcntl.flock(self.fh, fcntl.LOCK_UN)
        except (ImportError, OSError):
            pass
        self.fh.close()
        return False


def _fresh(out: str, src: str) -> bool:
    return os.path.exists(out) and os.path.getsize(out) > 0 and os.path.getmtime(out) >= os.path.getmtime(src)


def _decompress(path: str) -> str:
    base = os.path.basename(path)
    for suf in COMPRESSED_SUFFIXES:
        if base.lower().endswith(suf):
            base = base[: -len(suf)]
            break
    out = os.path.join(_writable_dir(path), base + ".taxtriage.fasta")
    if _fresh(out, path) and _fresh(out + ".fai", out):
        return out
    with _Lock(out):
        if not _fresh(out, path):              # re-check: another process may have done it
            tmp = f"{out}.tmp.{os.getpid()}"
            with gzip.open(path, "rb") as src, open(tmp, "wb") as dst:
                shutil.copyfileobj(src, dst, 1 << 20)
            os.replace(tmp, out)               # atomic: nobody sees a partial file
            print(f"[fasta] decompressed gzip reference {os.path.basename(path)} -> {out} "
                  f"(pysam needs plain or bgzip FASTA)", file=sys.stderr)
        if not _fresh(out + ".fai", out):
            _try_faidx(out)                    # build the .fai while still holding the lock
    return out


def _try_faidx(path: str) -> bool:
    import pysam
    try:
        fa = pysam.FastaFile(path)
        fa.close()
        return True
    except Exception:
        return False


def _indexable(path: str) -> bool:
    """Open once under the lock so a missing index is built by one process only."""
    if os.path.exists(path + ".fai") and os.path.getmtime(path + ".fai") >= os.path.getmtime(path):
        return _try_faidx(path)
    with _Lock(path):
        return _try_faidx(path)


def pysam_fasta_path(path: str) -> str:
    """Return a path ``pysam.FastaFile`` can open for this reference (cached)."""
    key = os.path.abspath(str(path))
    hit = _resolved.get(key)
    if hit and os.path.exists(hit):
        return hit
    p = str(path)
    if not os.path.exists(p):
        raise FileNotFoundError(f"reference FASTA '{p}' not found")
    if is_gzip(p):
        out = p if (is_bgzf(p) and _indexable(p)) else _decompress(p)
    else:
        out = p
        if not _indexable(p):
            # .fai cannot be written next to the file: index a symlink instead
            link = os.path.join(_cache_dir(), f"{abs(hash(key))}_{os.path.basename(p)}")
            if not os.path.exists(link):
                try:
                    os.symlink(os.path.abspath(p), link)
                except FileExistsError:
                    pass
            _indexable(link)
            out = link
    _resolved[key] = out
    return out


class GzFastaReader:
    """Read-only, pysam.FastaFile-compatible access to an ordinary gzip FASTA.

    No decompressed copy is written. Index = [(name, length, seq_offset)] in
    uncompressed coordinates, built by one streaming pass (or loaded from the
    sidecar cache). ``fetch`` seeks the gzip stream forward to the contig
    (re-opening only to go backwards) and keeps the ``cache_contigs`` most
    recently used sequences in memory.
    """

    def __init__(self, path: str, cache_contigs: int = 2):
        self.filename = str(path)
        self._index = _gz_index(self.filename)
        self.references = [n for n, _, _ in self._index]      # list, as pysam returns
        self.lengths = [l for _, l, _ in self._index]
        self.nreferences = len(self.references)
        self._pos = {n: i for i, (n, _, _) in enumerate(self._index)}
        self._fh = None
        self._cache = {}
        self._order = []
        self._max = max(1, int(cache_contigs))

    # -- pysam-compatible API -------------------------------------------------
    def get_reference_length(self, ref: str) -> int:
        return self._index[self._pos[ref]][1]

    def __contains__(self, ref):
        return ref in self._pos

    def fetch(self, reference=None, start=None, end=None, region=None):
        if region is not None and reference is None:
            m = re.match(r"^(.+?)(?::(\d[\d,]*)-(\d[\d,]*))?$", region)
            reference = m.group(1)
            if m.group(2):
                start, end = int(m.group(2).replace(",", "")) - 1, int(m.group(3).replace(",", ""))
        if reference not in self._pos:
            raise KeyError(f"sequence '{reference}' not present")
        seq = self._contig(reference)
        if start is None and end is None:
            return seq
        return seq[(start or 0):(end if end is not None else len(seq))]

    def close(self):
        if self._fh is not None:
            self._fh.close()
            self._fh = None
        self._cache.clear()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False

    # -- internals ------------------------------------------------------------
    def _contig(self, ref: str) -> str:
        hit = self._cache.get(ref)
        if hit is not None:
            self._order.remove(ref)
            self._order.append(ref)
            return hit
        _, length, off = self._index[self._pos[ref]]
        if self._fh is None or self._fh.tell() > off:
            if self._fh is not None:
                self._fh.close()
            self._fh = gzip.open(self.filename, "rb")
        self._fh.seek(off)                       # forward seek: streams, no disk
        parts, got = [], 0
        for line in self._fh:
            if line.startswith(b">"):
                break
            line = line.strip()
            parts.append(line)
            got += len(line)
            if got >= length:
                break
        seq = b"".join(parts).decode("ascii")
        self._cache[ref] = seq
        self._order.append(ref)
        while len(self._order) > self._max:
            self._cache.pop(self._order.pop(0), None)
        return seq


def _gz_index_scan(path: str):
    """One streaming pass: [(name, length, uncompressed offset of sequence)]."""
    out = []
    name, length, off, pos = None, 0, 0, 0
    with gzip.open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if name is not None:
                    out.append((name, length, off))
                name = line[1:].split(None, 1)[0].decode() if len(line) > 1 else ""
                length = 0
                pos += len(line)
                off = pos
                continue
            length += len(line.strip())
            pos += len(line)
    if name is not None:
        out.append((name, length, off))
    return out


def _gz_index(path: str):
    """Index from the sidecar cache when it matches the file, else scan + save."""
    import json
    st = os.stat(path)
    sig = [os.path.abspath(path), st.st_size, int(st.st_mtime)]
    side = os.path.join(_cache_dir(), f"{abs(hash(sig[0]))}_{os.path.basename(path)}.gzindex.json")
    try:
        with open(side) as fh:
            data = json.load(fh)
        if data.get("sig") == sig:
            return [tuple(x) for x in data["index"]]
    except (OSError, ValueError, KeyError):
        pass
    idx = _gz_index_scan(path)
    tmp = f"{side}.tmp.{os.getpid()}"
    try:
        with open(tmp, "w") as fh:
            json.dump({"sig": sig, "index": idx}, fh)
        os.replace(tmp, side)
    except OSError:
        pass                                   # read-only cache: just rescan next time
    return idx


def open_pysam_fasta(path: str):
    """Open any accepted reference with the pysam.FastaFile API.

    plain / bgzip -> pysam.FastaFile; ordinary gzip -> GzFastaReader (read
    directly, no decompressed copy) unless TAXTRIAGE_GZ_FASTA=decompress.
    """
    import pysam
    p = str(path)
    if not os.path.exists(p):
        raise FileNotFoundError(f"reference FASTA '{p}' not found")
    if is_gzip(p) and not (is_bgzf(p) and _indexable(p)):
        if os.environ.get("TAXTRIAGE_GZ_FASTA", "stream").lower() == "decompress":
            return pysam.FastaFile(_decompress(p))
        return GzFastaReader(p)
    return pysam.FastaFile(pysam_fasta_path(p))
