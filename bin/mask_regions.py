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
Region masking for TASS scoring.

TaxTriage normally scores every reference over its WHOLE length. A region mask
restricts scoring to part of each reference:

  mode = "exclude"   the BED lists regions to IGNORE (e.g. NCBI FCS contaminant
                     spans, adaptor/vector, low-complexity or shared mobile
                     elements). Everything else is scored.
  mode = "target"    the BED lists the ONLY regions to score (e.g. marker
                     genes, an amplicon panel, a toxin locus). References with
                     no target interval are dropped from scoring unless
                     untargeted="keep".

Masking works by COORDINATE COMPRESSION. For every reference the "kept"
intervals (complement of the exclude set, or the targets themselves) are laid
end to end into a new coordinate space whose length is the kept bp. Each read
is clipped to the kept intervals and its surviving pieces are translated into
that space. Every downstream metric (breadth, mean depth, Gini evenness,
positional histograms, RPKM, the contig penalty...) then operates on a genome
that simply does not contain the masked bases, so no metric needs to know about
masking and a masked base can neither help nor hurt the score.

A read is dropped entirely when at least ``read_overlap`` of its aligned span
falls in masked space (default 0.5). Kept reads contribute only their
unmasked bases to coverage.

The BED loader accepts 3+ column BED (0-based, half-open). Column 4, when
present, is kept as a label (e.g. "FCS:EXCLUDE:prok:proteobacteria").
Accession matching is tolerant: versioned/unversioned IDs and the RefSeq
"NZ_" prefix are reconciled automatically (NZ_CP012345.1 <-> CP012345.1).

Library only -- no third-party imports so it can be used from any container.
"""
from __future__ import annotations

import bisect
import gzip
import os
import sys
from collections import defaultdict
from typing import Dict, Iterable, List, Optional, Tuple

Interval = Tuple[int, int]

VALID_MODES = ("exclude", "target")


# ─────────────────────────────────────────────────────────────────────────────
# Accession helpers
# ─────────────────────────────────────────────────────────────────────────────
def accession_aliases(acc: str) -> List[str]:
    """Return lookup aliases for an accession, most specific first.

    NZ_CP012345.1 -> [NZ_CP012345.1, CP012345.1, NZ_CP012345, CP012345]
    CP012345.1    -> [CP012345.1, NZ_CP012345.1, CP012345, NZ_CP012345]
    """
    acc = acc.strip()
    out = [acc]
    if acc.startswith("NZ_"):
        out.append(acc[3:])
    elif "_" not in acc:
        # GenBank / ENA / DDBJ style -- RefSeq copies of WGS & complete genomes
        # are the same accession with an NZ_ prefix.
        out.append("NZ_" + acc)
    unversioned = []
    for a in out:
        if "." in a:
            unversioned.append(a.rsplit(".", 1)[0])
    out.extend(unversioned)
    seen = set()
    return [a for a in out if not (a in seen or seen.add(a))]


def open_text(path: str):
    """Open plain or gzip text transparently."""
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt")
    return open(path, "r")


# ─────────────────────────────────────────────────────────────────────────────
# Interval helpers
# ─────────────────────────────────────────────────────────────────────────────
def merge_intervals(ivals: Iterable[Interval]) -> List[Interval]:
    """Sort + merge overlapping/abutting half-open intervals."""
    s = sorted((int(a), int(b)) for a, b in ivals if int(b) > int(a))
    out: List[List[int]] = []
    for a, b in s:
        if out and a <= out[-1][1]:
            if b > out[-1][1]:
                out[-1][1] = b
        else:
            out.append([a, b])
    return [(a, b) for a, b in out]


def complement(ivals: List[Interval], length: int) -> List[Interval]:
    """Complement of merged intervals over [0, length)."""
    out = []
    pos = 0
    for a, b in ivals:
        a = max(0, min(a, length))
        b = max(0, min(b, length))
        if a > pos:
            out.append((pos, a))
        pos = max(pos, b)
    if pos < length:
        out.append((pos, length))
    return out


def subtract(a: List[Interval], b: List[Interval]) -> List[Interval]:
    """a minus b; both merged & sorted."""
    out = []
    j = 0
    for s, e in a:
        cur = s
        while j < len(b) and b[j][1] <= cur:
            j += 1
        k = j
        while k < len(b) and b[k][0] < e:
            bs, be = b[k]
            if bs > cur:
                out.append((cur, min(bs, e)))
            cur = max(cur, be)
            if cur >= e:
                break
            k += 1
        if cur < e:
            out.append((cur, e))
    return out


def clip_to_length(ivals: List[Interval], length: int) -> List[Interval]:
    out = []
    for a, b in ivals:
        a = max(0, min(a, length))
        b = max(0, min(b, length))
        if b > a:
            out.append((a, b))
    return out


# ─────────────────────────────────────────────────────────────────────────────
# BED loading
# ─────────────────────────────────────────────────────────────────────────────
def load_bed(paths: Iterable[str]) -> Tuple[Dict[str, List[Interval]], Dict[str, List[Tuple[int, int, str]]]]:
    """Load one or more BED files.

    Returns (merged, labelled):
      merged   {chrom: [(start, end), ...]} merged per chrom
      labelled {chrom: [(start, end, label), ...]} raw rows for reporting
    """
    raw: Dict[str, List[Interval]] = defaultdict(list)
    labelled: Dict[str, List[Tuple[int, int, str]]] = defaultdict(list)
    for p in paths:
        if not p:
            continue
        if not os.path.exists(p):
            raise FileNotFoundError(f"mask BED not found: {p}")
        with open_text(p) as fh:
            for ln, line in enumerate(fh, 1):
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 3:
                    f = line.split()
                if len(f) < 3:
                    print(f"[mask] WARNING {p}:{ln}: fewer than 3 columns, skipped", file=sys.stderr)
                    continue
                try:
                    s, e = int(f[1]), int(f[2])
                except ValueError:
                    print(f"[mask] WARNING {p}:{ln}: non-integer coordinates, skipped", file=sys.stderr)
                    continue
                if e <= s:
                    continue
                chrom = f[0]
                raw[chrom].append((s, e))
                labelled[chrom].append((s, e, f[3] if len(f) > 3 else "mask"))
    merged = {c: merge_intervals(v) for c, v in raw.items()}
    return merged, dict(labelled)


# ─────────────────────────────────────────────────────────────────────────────
# The mask
# ─────────────────────────────────────────────────────────────────────────────
class RefMask:
    """Kept-interval index for ONE reference with coordinate compression."""

    __slots__ = ("length", "kept", "_starts", "_cum", "kept_bp", "masked")

    def __init__(self, length: int, kept: List[Interval], masked: List[Interval]):
        self.length = int(length)
        self.kept = kept
        self.masked = masked
        self._starts = [a for a, _ in kept]
        cum = []
        tot = 0
        for a, b in kept:
            cum.append(tot)          # compressed coord of interval start
            tot += b - a
        self._cum = cum
        self.kept_bp = tot

    @property
    def masked_bp(self) -> int:
        return self.length - self.kept_bp

    def clip(self, start: int, end: int) -> Tuple[List[Interval], List[Interval], int]:
        """Clip [start, end) to kept space.

        Returns (orig_pieces, compressed_pieces, kept_bases).
        """
        if end <= start or not self.kept:
            return [], [], 0
        i = bisect.bisect_right(self._starts, start) - 1
        if i < 0:
            i = 0
        orig, comp = [], []
        kb = 0
        n = len(self.kept)
        while i < n:
            a, b = self.kept[i]
            if a >= end:
                break
            s = max(a, start)
            e = min(b, end)
            if e > s:
                orig.append((s, e))
                cs = self._cum[i] + (s - a)
                comp.append((cs, cs + (e - s)))
                kb += e - s
            i += 1
        return orig, comp, kb

    def to_compressed(self, pos: int) -> int:
        """Map an ORIGINAL coordinate inside kept space to compressed space.
        Positions inside a masked gap snap to the next kept base."""
        i = bisect.bisect_right(self._starts, pos) - 1
        if i < 0:
            return 0
        a, b = self.kept[i]
        if pos >= b:
            return self._cum[i] + (b - a)
        return self._cum[i] + (pos - a)


class RegionMask:
    """Collection of per-reference masks built from a BED file."""

    def __init__(
        self,
        bed_paths: Iterable[str],
        mode: str = "exclude",
        read_overlap: float = 0.5,
        untargeted: str = "drop",
        min_region: int = 0,
        pad: int = 0,
        exclude_paths: Iterable[str] = (),
    ):
        mode = (mode or "exclude").lower()
        if mode not in VALID_MODES:
            raise ValueError(f"mask mode must be one of {VALID_MODES}, got {mode!r}")
        if not (0.0 < float(read_overlap) <= 1.0):
            raise ValueError("read_overlap must be in (0, 1]")
        self.mode = mode
        self.read_overlap = float(read_overlap)
        self.untargeted = (untargeted or "drop").lower()
        self.bed_paths = [p for p in bed_paths if p]
        self.pad = int(pad or 0)
        self.min_region = int(min_region or 0)
        merged, labelled = load_bed(self.bed_paths)
        if self.pad or self.min_region:
            adj = {}
            for c, iv in merged.items():
                iv2 = [(max(0, a - self.pad), b + self.pad) for a, b in iv if (b - a) >= self.min_region]
                adj[c] = merge_intervals(iv2)
            merged = adj
        self.exclude_paths = [p for p in (exclude_paths or []) if p]
        excl, excl_labels = load_bed(self.exclude_paths) if self.exclude_paths else ({}, {})
        if self.mode == "exclude" and excl:
            # exclude mode: extra exclusions simply join the mask
            for c, iv in excl.items():
                merged[c] = merge_intervals(list(merged.get(c, [])) + list(iv))
                labelled.setdefault(c, []).extend(excl_labels.get(c, []))
            excl = {}
        self._excl = excl          # target mode only: subtracted from targets
        self._excl_alias: Dict[str, str] = {}
        for c in excl:
            for a in accession_aliases(c):
                self._excl_alias.setdefault(a, c)
        self._bed = merged
        self._labels = labelled
        # alias -> bed chrom name
        self._alias: Dict[str, str] = {}
        for c in merged:
            for a in accession_aliases(c):
                self._alias.setdefault(a, c)
        self._refs: Dict[str, Optional[RefMask]] = {}
        self.matched_refs: set = set()

    # -- lookup ---------------------------------------------------------------
    def bed_chrom_for(self, ref: str) -> Optional[str]:
        for a in accession_aliases(ref):
            c = self._alias.get(a)
            if c is not None:
                return c
        return None

    def labels_for(self, ref: str) -> List[Tuple[int, int, str]]:
        c = self.bed_chrom_for(ref)
        return list(self._labels.get(c, [])) if c else []

    def build(self, ref: str, length: int) -> Optional[RefMask]:
        """Build (and cache) the RefMask for a reference.

        Returns None when the reference is NOT affected by the mask (exclude
        mode with no intervals, or target mode with untargeted='keep'), so the
        caller can take the unmasked fast path.
        """
        if ref in self._refs:
            return self._refs[ref]
        c = self.bed_chrom_for(ref)
        ivals = clip_to_length(self._bed.get(c, []), length) if c else []
        if self.mode == "exclude":
            if not ivals:
                rm = None
            else:
                rm = RefMask(length, complement(ivals, length), ivals)
                self.matched_refs.add(ref)
        else:  # target
            ex = []
            if self._excl:
                for a in accession_aliases(ref):
                    ec = self._excl_alias.get(a)
                    if ec is not None:
                        ex = clip_to_length(self._excl[ec], length)
                        break
            if not ivals:
                if self.untargeted == "keep":
                    # untargeted refs are scored whole -- minus any exclusions
                    rm = RefMask(length, complement(ex, length), ex) if ex else None
                else:
                    rm = RefMask(length, [], [(0, length)] if length > 0 else [])
            else:
                kept = subtract(ivals, ex) if ex else ivals
                rm = RefMask(length, kept, complement(kept, length))
                self.matched_refs.add(ref)
        self._refs[ref] = rm
        return rm

    def read_is_masked(self, aln_len: int, kept_bases: int) -> bool:
        if aln_len <= 0:
            return kept_bases <= 0
        if kept_bases <= 0:
            return True
        return (1.0 - kept_bases / aln_len) >= self.read_overlap

    def describe(self) -> dict:
        return {
            "enabled": True,
            "mode": self.mode,
            "bed": [os.path.basename(p) for p in self.bed_paths],
            "exclude_bed": [os.path.basename(p) for p in self.exclude_paths] or None,
            "read_overlap": self.read_overlap,
            "untargeted": self.untargeted if self.mode == "target" else None,
            "pad": self.pad,
            "min_region": self.min_region,
            "bed_sequences": len(self._bed),
            "bed_intervals": sum(len(v) for v in self._bed.values()),
            "bed_bp": sum(b - a for v in self._bed.values() for a, b in v),
        }


def mask_fraction_bins(rm: Optional[RefMask], length: int, offset: int, bin_size: int, n_bins: int, acc_bins: List[float]):
    """Accumulate masked-bp per bin (concatenated coords) into acc_bins."""
    if rm is None or not rm.masked or bin_size <= 0:
        return
    for a, b in rm.masked:
        a += offset
        b += offset
        bi = a // bin_size
        while a < b and bi < n_bins:
            be = min(b, (bi + 1) * bin_size)
            acc_bins[bi] += be - a
            a = be
            bi += 1
