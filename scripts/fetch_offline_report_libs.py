#!/usr/bin/env python3
"""
Download the CDN libraries that assets/heatmap.html references into a local
folder, so the report can be built fully offline with --offline_report_files.

The CDN URLs are read straight out of assets/heatmap.html, so this stays in
sync if a library version is bumped in the template. For each stylesheet, the
fonts / marker images it references via url(...) are downloaded too, so Font
Awesome icons and Leaflet markers work with no network access.

Usage:
    python scripts/fetch_offline_report_libs.py
    python scripts/fetch_offline_report_libs.py -o /path/to/libs
    python scripts/fetch_offline_report_libs.py -t assets/heatmap.html

Then build an offline report with the downloaded folder:
    nextflow run . ... --offline_report_files assets/offline_report_libs
    # or, to test the embedding directly:
    python bin/report_template.py -t assets/heatmap.html \\
        -o report_offline.html --offline_report_files assets/offline_report_libs
"""
import argparse
import os
import re
import sys
from pathlib import Path
from urllib.parse import urljoin, urlparse
from urllib.request import Request, urlopen

_REPO_ROOT = Path(__file__).resolve().parent.parent
_DEFAULT_TEMPLATE = _REPO_ROOT / "assets" / "heatmap.html"
_DEFAULT_OUTDIR = _REPO_ROOT / "assets" / "offline_report_libs"

# CDN <script src> and <link href> tags in the template head.
_CDN_SCRIPT_RE = re.compile(
    r'<script\b[^>]*\bsrc="(?P<url>https?://[^"]+)"', re.IGNORECASE)
_CDN_LINK_RE = re.compile(
    r'<link\b[^>]*\bhref="(?P<url>https?://[^"]+\.css)"', re.IGNORECASE)
# url(...) references inside CSS (fonts, marker images); skips data: URIs.
_CSS_URL_RE = re.compile(r'url\(\s*["\']?(?P<u>[^"\')]+)["\']?\s*\)')


def _download(url):
    req = Request(url, headers={"User-Agent": "taxtriage-offline-report"})
    with urlopen(req, timeout=60) as resp:  # nosec B310 - fixed https CDN URLs
        return resp.read()


def _save(data, outdir, name):
    dest = outdir / name
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open(dest, "wb") as fh:
        fh.write(data)
    return dest


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-t", "--template", default=str(_DEFAULT_TEMPLATE),
                    help="report template to read CDN URLs from "
                         "(default: assets/heatmap.html)")
    ap.add_argument("-o", "--outdir", default=str(_DEFAULT_OUTDIR),
                    help="directory to download into "
                         "(default: assets/offline_report_libs)")
    args = ap.parse_args(argv)

    template = Path(args.template)
    outdir = Path(args.outdir)
    if not template.is_file():
        sys.exit(f"ERROR: template not found: {template}")

    html = template.read_text(encoding="utf-8", newline="")
    scripts = [m.group("url") for m in _CDN_SCRIPT_RE.finditer(html)]
    stylesheets = [m.group("url") for m in _CDN_LINK_RE.finditer(html)]

    if not scripts and not stylesheets:
        sys.exit(f"ERROR: no CDN <script>/<link> references found in {template}")

    outdir.mkdir(parents=True, exist_ok=True)
    print(f"Reading CDN references from : {template}")
    print(f"Downloading into           : {outdir.resolve()}\n")

    n_files = 0
    n_assets = 0

    for url in scripts:
        name = os.path.basename(urlparse(url).path)
        print(f"  JS  {name:<30} <- {url}")
        _save(_download(url), outdir, name)
        n_files += 1

    for url in stylesheets:
        name = os.path.basename(urlparse(url).path)
        print(f"  CSS {name:<30} <- {url}")
        data = _download(url)
        _save(data, outdir, name)
        n_files += 1
        # Pull the fonts / images the stylesheet references via url(...).
        css_text = data.decode("utf-8", "replace")
        seen = set()
        for m in _CSS_URL_RE.finditer(css_text):
            raw = m.group("u").strip()
            if not raw or raw.startswith("data:"):
                continue
            clean = raw.split("?", 1)[0].split("#", 1)[0]
            asset_url = urljoin(url, clean)
            asset_name = os.path.basename(urlparse(asset_url).path)
            if not asset_name or asset_name in seen:
                continue
            seen.add(asset_name)
            try:
                _save(_download(asset_url), outdir, asset_name)
            except Exception as exc:  # noqa: BLE001 - report and continue
                print(f"      ! skipped {asset_name}: {exc}")
                continue
            print(f"      asset {asset_name}")
            n_assets += 1

    rel = os.path.relpath(outdir.resolve(), Path.cwd())
    print(f"\nDone. {n_files} libraries + {n_assets} referenced assets saved to:")
    print(f"  {outdir.resolve()}")
    print("\nBuild a fully offline report by pointing --offline_report_files at it:")
    print(f"  nextflow run . ... --offline_report_files {rel}")
    print("\nOr test the embedding directly, without running the pipeline:")
    print(f"  python bin/report_template.py -t {args.template} \\")
    print(f"      -o report_offline.html --offline_report_files {rel}")
    print("\n(Note: --offline_report_files embeds these local copies with no "
          "network.\n --offline_report instead downloads the same libraries at "
          "build time.)")


if __name__ == "__main__":
    main()
