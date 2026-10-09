#!/usr/bin/env python3
# ##############################################################################################
# # Copyright 2025 The Johns Hopkins University Applied Physics Laboratory LLC
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
Compact, cross-sample-comparable depth profiles for the per-sample JSON.

A strain's ``depth_profile`` is written by match_paths.py and read by
alignment_trends.py and the report's Alignment Trends tab:

    {
      "v": 1,
      "window": 1600,                # bp per window (100 * 2**k ladder)
      "enc": "log2p1x16",            # depth byte = round(log2(1 + depth) * 16)
      "read_len": 150,               # mean aligned read length (optional)
      "reads": 1234,                 # reads aligned to the strain (optional)
      "contigs": [
        ["NC_000001.1", 185000, "<b64 depth bytes>", "<b64 breadth bytes>"],
        ["NC_000002.1", 4200]        # no reads in this sample -> all-zero windows
      ],
      "total_len": 189200,           # every accession of the reference
      "n_contigs": 2,
      "n_uncovered": 0,              # only present when zero-read accessions were
      "uncovered_bp": 0              #   too many (> MAX_LISTED_UNCOVERED) to list
    }

Design choices that make samples line up:
  * Every accession of the strain is listed, covered or not, in sorted order,
    each in its OWN coordinates. (The breadth_histogram concatenates only the
    covered accessions, so its coordinates shift from sample to sample.)
  * The window size is picked from a fixed 100 * 2**k ladder, so two samples
    with different reference sets still have integer window ratios and can be
    re-binned onto a common grid.
  * Depth is a log-scaled byte (about 4% resolution, saturating ~61,000x);
    breadth is the percent of the window covered by >=1 read (0-100).
"""
import base64
import math

PROFILE_VERSION = 1
LADDER_BASE = 100
DEPTH_ENC = "log2p1x16"
# Zero-read accessions are listed one by one only up to this many. A draft
# assembly with tens of thousands of scaffolds would otherwise add megabytes to
# every (indented) per-sample JSON; past the cap they are summarised by count and
# bp ("n_uncovered" / "uncovered_bp"), and "total_len" keeps the full size.
MAX_LISTED_UNCOVERED = 500


def choose_window(total_len, target_windows):
    """Smallest 100*2**k window giving <= target_windows windows over total_len."""
    w = LADDER_BASE
    if target_windows <= 0 or total_len <= 0:
        return w
    while total_len / w > target_windows:
        w *= 2
    return w


def encode_depth(d):
    if d <= 0:
        return 0
    return min(255, int(round(math.log2(1.0 + d) * 16)))


def decode_depth(b):
    return 0.0 if b <= 0 else (2.0 ** (b / 16.0)) - 1.0


def _window_profile_np(regions, length, window):
    """Vectorised window_profile (numpy). Runs inside one window — almost all of
    them on a deeply covered genome — are summed with bincount; the few that
    straddle a window edge are split in a short Python loop."""
    import numpy as np
    length = int(length)
    n = max(1, (length + window - 1) // window)
    arr = np.asarray(regions, dtype=np.float64).reshape(-1, 3)
    s = np.clip(arr[:, 0], 0, length).astype(np.int64)
    e = np.clip(arr[:, 1], 0, length).astype(np.int64)
    d = arr[:, 2]
    keep = e > s
    s, e, d = s[keep], e[keep], d[keep]
    b0 = s // window
    b1 = (e - 1) // window
    one = b0 == b1
    ln = (e - s)[one].astype(np.float64)
    cov = np.bincount(b0[one], weights=ln, minlength=n)[:n]
    dsum = np.bincount(b0[one], weights=ln * d[one], minlength=n)[:n]
    for si, ei, di in zip(s[~one].tolist(), e[~one].tolist(), d[~one].tolist()):
        for b in range(si // window, min(n - 1, (ei - 1) // window) + 1):
            lo = b * window
            ov = min(ei, lo + window, length) - max(si, lo)
            if ov > 0:
                cov[b] += ov
                dsum[b] += ov * di
    lo = np.arange(n, dtype=np.int64) * window
    blen = np.maximum(1, np.minimum(lo + window, length) - lo).astype(np.float64)
    return (dsum / blen).tolist(), np.minimum(100.0, cov * 100.0 / blen).tolist()


def window_profile(regions, length, window):
    """Mean depth and percent-breadth per window from non-overlapping
    (start, end, depth) runs on one accession."""
    if len(regions) > 2000:
        try:
            return _window_profile_np(regions, length, window)
        except ImportError:
            pass
    n = max(1, (int(length) + window - 1) // window)
    dsum = [0.0] * n
    cov = [0] * n
    for s, e, d in regions or ():
        s = max(0, int(s)); e = min(int(length), int(e))
        if e <= s:
            continue
        b0 = s // window
        b1 = min(n - 1, (e - 1) // window)
        for b in range(b0, b1 + 1):
            lo = b * window
            hi = min(lo + window, int(length))
            ov = min(e, hi) - max(s, lo)
            if ov > 0:
                cov[b] += ov
                dsum[b] += ov * float(d)
    depth, breadth = [], []
    for b in range(n):
        lo = b * window
        blen = max(1, min(lo + window, int(length)) - lo)
        depth.append(dsum[b] / blen)
        breadth.append(min(100.0, cov[b] * 100.0 / blen))
    return depth, breadth


def build_profile(accessions, target_windows=400, read_len=None, reads=None):
    """accessions: iterable of (accession, length, covered_regions).
    read_len (optional) is stored so readers can tell how many reads a window
    should hold at a given depth (i.e. whether an empty window is meaningful)."""
    accs = sorted((str(a), int(l), r) for a, l, r in accessions if int(l or 0) > 0)
    if not accs:
        return None
    total = sum(l for _, l, _ in accs)
    w = choose_window(total, target_windows)
    uncovered = [(a, l) for a, l, r in accs if not r]
    list_uncovered = len(uncovered) <= MAX_LISTED_UNCOVERED
    contigs = []
    for acc, length, regs in accs:
        if not regs:
            if list_uncovered:
                contigs.append([acc, length])
            continue
        dep, br = window_profile(regs, length, w)
        db = bytes(encode_depth(x) for x in dep)
        bb = bytes(min(100, int(round(x))) for x in br)
        contigs.append([acc, length,
                        base64.b64encode(db).decode(),
                        base64.b64encode(bb).decode()])
    prof = {"v": PROFILE_VERSION, "window": w, "enc": DEPTH_ENC, "contigs": contigs,
            "total_len": total, "n_contigs": len(accs)}
    if not list_uncovered:
        prof["n_uncovered"] = len(uncovered)
        prof["uncovered_bp"] = sum(l for _, l in uncovered)
    if read_len:
        prof["read_len"] = int(round(read_len))
    if reads:
        prof["reads"] = int(reads)
    return prof


def decode_profile(profile):
    """-> (window, {accession: (length, [depth...], [breadth...])})"""
    if not profile or not profile.get("contigs"):
        return None, {}
    w = int(profile.get("window") or LADDER_BASE)
    out = {}
    for entry in profile["contigs"]:
        acc, length = str(entry[0]), int(entry[1])
        n = max(1, (length + w - 1) // w)
        if len(entry) >= 4 and entry[2]:
            dep = [decode_depth(b) for b in base64.b64decode(entry[2])]
            br = [float(b) for b in base64.b64decode(entry[3])]
        else:
            dep, br = [0.0] * n, [0.0] * n
        out[acc] = (length, dep[:n], br[:n])
    return w, out


def rebin(values, factor, length, window, mode="mean"):
    """Merge `factor` consecutive windows (window -> window*factor), weighting
    by each window's true bp length (the last window of a contig is short)."""
    if factor <= 1:
        return list(values)
    out = []
    n = len(values)
    for i in range(0, n, factor):
        tot, wsum = 0.0, 0.0
        for j in range(i, min(n, i + factor)):
            lo = j * window
            bl = max(1, min(lo + window, length) - lo)
            tot += values[j] * bl
            wsum += bl
        out.append(tot / wsum if wsum else 0.0)
    return out
