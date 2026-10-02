#!/usr/bin/env python3
"""
Reassemble the self-contained interactive-report HTML from the thin
assets/heatmap.html shell plus its external parts under assets/src/.

assets/heatmap.html references its CSS and JS as external files:

    <link rel="stylesheet" href="src/css/main.css" />
    <script src="src/js/05_tab_heatmap.js"></script>

That keeps the committed template small and editable per tab/plot. The two
downstream builders each need ONE portable file, so they call inline_template()
to fold every local part back inline. CDN links and the heatmap_boot.js data
anchor are left untouched.

    bin/make_report.py          -> all.comparison.report.html  (nextflow)
    scripts/inline_boot_json.py -> _site/index.html            (GitHub Pages)

Part resolution: references are looked up relative to the template's directory
and, as a fallback, relative to the template's real path (Nextflow stages the
template as a symlink into the task workdir, so the real path points back at
assets/ where src/ lives even if src/ wasn't staged).

Library use:
    from report_template import inline_template
    html = inline_template("assets/heatmap.html")

CLI (inspect / diff the assembled output):
    python bin/report_template.py -o out.html
"""
import argparse
import base64
import json
import os
import re
import sys
from pathlib import Path
from urllib.parse import urljoin, urlparse
from urllib.request import Request, urlopen

_LINK_RE = re.compile(
    r'[ \t]*<link\b[^>]*\bhref="(?P<href>src/[^"]+\.css)"[^>]*>[ \t]*\n?',
    re.IGNORECASE,
)
_SCRIPT_RE = re.compile(
    r'[ \t]*<script\b[^>]*\bsrc="(?P<src>src/[^"]+\.js)"[^>]*>\s*</script>[ \t]*\n?',
    re.IGNORECASE,
)
_STYLE_ID_RE = re.compile(r'\bdata-style-id="([^"]+)"')

# ── CDN parts (http/https) — inlined only in offline builds ───────────────────
# Left untouched by default so the report still loads libraries from the CDN at
# view time (the historical behaviour). When --offline_report /
# --offline_report_files is requested, these are folded inline instead.
_CDN_LINK_RE = re.compile(
    r'[ \t]*<link\b[^>]*\bhref="(?P<href>https?://[^"]+\.css)"[^>]*>[ \t]*\n?',
    re.IGNORECASE,
)
_CDN_SCRIPT_RE = re.compile(
    r'[ \t]*<script\b[^>]*\bsrc="(?P<src>https?://[^"]+)"[^>]*>\s*</script>[ \t]*\n?',
    re.IGNORECASE,
)
_CDN_ANY_RE = re.compile(
    _CDN_SCRIPT_RE.pattern + "|" + _CDN_LINK_RE.pattern, re.IGNORECASE)
# url(...) references inside CSS (fonts, marker images). Skips data: URIs.
_CSS_URL_RE = re.compile(r'url\(\s*(?P<q>["\']?)(?P<u>[^"\')]+)(?P=q)\s*\)')

# Extensions -> MIME for data-URI embedding of CSS-referenced assets.
_ASSET_MIME = {
    ".woff2": "font/woff2", ".woff": "font/woff", ".ttf": "font/ttf",
    ".otf": "font/otf", ".eot": "application/vnd.ms-fontobject",
    ".svg": "image/svg+xml", ".png": "image/png", ".gif": "image/gif",
    ".jpg": "image/jpeg", ".jpeg": "image/jpeg", ".webp": "image/webp",
}


# ── Data the report fetches at VIEW time (not <script>/<link> tags) ──────────
# The Metadata choropleth / map-group views and the map's offline basemap load
# Natural Earth boundaries with d3.json(). In an offline build they are embedded
# as window.TT_OFFLINE.geo so the report never reaches for the network. Keep
# these lists in step with _geoCountries()/_geoAdmin1() in
# assets/src/js/33_categorical_metadata.js (same URLs, same order).
OFFLINE_GEO_SOURCES = {
    "countries": [
        "https://cdn.jsdelivr.net/gh/nvkelso/natural-earth-vector@master/geojson/ne_110m_admin_0_countries.geojson",
        "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_110m_admin_0_countries.geojson",
        "https://cdn.jsdelivr.net/gh/johan/world.geo.json@master/countries.geo.json",
    ],
    "admin1": [
        "https://cdn.jsdelivr.net/gh/nvkelso/natural-earth-vector@master/geojson/ne_50m_admin_1_states_provinces.geojson",
        "https://cdn.jsdelivr.net/gh/nvkelso/natural-earth-vector@v5.1.2/geojson/ne_50m_admin_1_states_provinces.geojson",
        "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_50m_admin_1_states_provinces.geojson",
    ],
}
# Feature properties the report reads (name matching + country drill-down);
# everything else is dropped to keep the embedded copy small.
_GEO_KEEP_PROPS = {
    "name", "NAME", "name_en", "name_long", "NAME_LONG", "name_alt", "gn_name", "woe_name",
    "abbrev", "postal", "admin", "ADMIN", "geounit", "geonunit", "GEOUNIT",
    "sovereignt", "SOVEREIGNT", "iso_a2", "ISO_A2", "iso_a3", "ISO_A3", "iso_3166_2",
}

# Inline-script safety: a library body (or JSON) containing the literal text
# "</script" would end the <script> element early and spill the rest of the
# file into the page. "<\/script" means the same thing to JavaScript.
_SCRIPT_CLOSE_RE = re.compile(r"</(script)", re.IGNORECASE)


def _script_safe(text):
    text = _SCRIPT_CLOSE_RE.sub(r"<\\/\1", text)
    # "<!--" followed later by "<script" puts the HTML parser in its
    # double-escaped script state, where the closing tag is not recognised.
    # Only touch a body that actually has both ("<\\!--" reads the same in JS
    # strings and non-unicode regexes).
    if "<!--" in text and re.search(r"<script", text, re.IGNORECASE):
        text = text.replace("<!--", "<\\!--")
    return text


# Library bundle shipped with the pipeline (scripts/fetch_offline_report_libs.py
# writes it; manifest.json maps each CDN URL to its file). --offline_report
# uses it first and only downloads what it lacks, so an offline build works even
# when the build host has no internet.
BUNDLED_LIBS_DIR = Path(__file__).resolve().parent.parent / "assets" / "offline_report_libs"
_MANIFEST_CACHE = {}


def _manifest(offline_dir):
    key = str(offline_dir)
    if key not in _MANIFEST_CACHE:
        mf = Path(offline_dir) / "manifest.json"
        try:
            _MANIFEST_CACHE[key] = json.loads(mf.read_text(encoding="utf-8")).get("files", {})
        except Exception:  # noqa: BLE001 - no / unreadable manifest: basename matching only
            _MANIFEST_CACHE[key] = None
    return _MANIFEST_CACHE[key]


def _local_lookup(url, offline_dir):
    """(bytes, exact) for a URL from a local library folder, or (None, False).

    exact=False means the file was matched by basename only although the folder
    has a manifest that does not list this URL -- i.e. probably a different
    library version than the template asks for.
    """
    mf = _manifest(offline_dir)
    if mf and url in mf:
        path = Path(offline_dir) / mf[url]
        if path.is_file():
            return path.read_bytes(), True
    name = os.path.basename(urlparse(url).path)
    for root, _dirs, files in os.walk(offline_dir):
        if name in files:
            with open(os.path.join(root, name), "rb") as fh:
                return fh.read(), not mf
    return None, False


def _download(url):
    req = Request(url, headers={"User-Agent": "taxtriage-offline-report"})
    with urlopen(req, timeout=60) as resp:  # nosec B310 - fixed https CDN URLs
        return resp.read()


def _fetch_bytes(url, offline_dir, allow_download=False):
    """Return the raw bytes for a CDN asset.

    offline_dir set, allow_download False -> local folder only (fully offline):
        exact manifest entry, else any file with the same basename (searched
        recursively, so the folder layout does not matter).
    offline_dir set, allow_download True  -> local folder first; a file that is
        missing, or present only as a different version, is downloaded instead.
    offline_dir None                      -> download the URL at build time.
    """
    if not offline_dir:
        return _download(url)
    data, exact = _local_lookup(url, offline_dir)
    if data is not None and exact:
        return data
    if allow_download:
        try:
            return _download(url)
        except Exception:  # noqa: BLE001
            if data is None:
                raise
    if data is not None:
        print(f"[report_template] WARNING: {url} is not in {offline_dir}/manifest.json; "
              f"embedding the local file with the same name (possibly another version). "
              f"Re-run scripts/fetch_offline_report_libs.py to refresh the folder.",
              file=sys.stderr)
        return data
    raise FileNotFoundError(
        f"offline asset '{os.path.basename(urlparse(url).path)}' not found under "
        f"{offline_dir} (for {url})")


def _inline_css_urls(css_text, base_url, offline_dir, allow_download=False):
    """Replace url(...) refs in a stylesheet with base64 data URIs.

    Fonts (Font Awesome) and marker images (Leaflet) are pulled the same way as
    the stylesheet itself. Unresolvable refs are left as-is so a partial offline
    bundle still produces a working (if icon-less) report.
    """
    # Font Awesome lists every font twice (woff2 + a truetype fallback). Every
    # browser that can run this report reads woff2, so embedding the .ttf copies
    # would only double the size of an offline build.
    css_text = re.sub(
        r',\s*url\(\s*["\']?[^"\')]+\.ttf["\']?\s*\)\s*format\(\s*["\']truetype["\']\s*\)',
        "", css_text)

    def repl(m):
        raw = m.group("u").strip()
        if not raw or raw.startswith(("data:", "#")):
            return m.group(0)  # data URI, or an IE behavior hook like url(#default#VML)
        clean = raw.split("?", 1)[0].split("#", 1)[0]
        if not clean:
            return m.group(0)
        try:
            data = _fetch_bytes(urljoin(base_url, clean), offline_dir, allow_download)
        except Exception:  # noqa: BLE001 - best effort, keep original ref
            return m.group(0)
        ext = os.path.splitext(urlparse(clean).path)[1].lower()
        mime = _ASSET_MIME.get(ext, "application/octet-stream")
        b64 = base64.b64encode(data).decode("ascii")
        return f"url(data:{mime};base64,{b64})"

    return _CSS_URL_RE.sub(repl, css_text)


def inline_cdn_assets(html, offline_dir=None, allow_download=False) -> str:
    """Fold CDN <script src>/<link href> tags inline for an offline report.

    offline_dir None -> download each library at build time.
    offline_dir set  -> read local copies (by basename) from that directory.
    """
    def _get(url):
        try:
            return _fetch_bytes(url, offline_dir, allow_download)
        except Exception as exc:  # noqa: BLE001
            how = (f"no usable local copy in {offline_dir}"
                   + (" and the download failed" if allow_download else "")
                   if offline_dir else "download failed (no network at build time?)")
            raise RuntimeError(
                f"offline report: could not embed {url} -- {how}: {exc}. "
                "Prepare a library folder with scripts/fetch_offline_report_libs.py on a "
                "machine with internet access and pass it via --offline_report_files."
            ) from exc

    def _script(m):
        body = _get(m.group("src")).decode("utf-8", "replace")
        return f"    <script>\n{_script_safe(body)}\n    </script>\n"

    def _css(m):
        text = _get(m.group("href")).decode("utf-8", "replace")
        text = _inline_css_urls(text, m.group("href"), offline_dir, allow_download)
        text = re.sub(r"</(style)", r"<\\/\1", text, flags=re.IGNORECASE)
        return f"    <style>\n{text}\n    </style>\n"

    # View-time data (boundaries) + the offline flag the report reads. Placed
    # just before the first CDN tag of the template head -- located BEFORE the
    # libraries are inlined, since library bodies contain strings like
    # "</head>" that a text search would otherwise land inside.
    geo = offline_geo_payload(offline_dir, allow_download)
    boot = {"build": True, "geo": geo}
    tag = ("    <script>\nwindow.TT_OFFLINE = "
           + _script_safe(json.dumps(boot, separators=(",", ":"))) + ";\n    </script>\n")
    first = _CDN_SCRIPT_RE.search(html) or _CDN_LINK_RE.search(html)
    head = re.search(r"<head\b[^>]*>\s*", html, re.IGNORECASE)
    at = first.start() if first else (head.end() if head else 0)
    html = html[:at] + tag + html[at:]

    # ONE pass over the template: a library body that has already been inlined
    # (jsPDF / SheetJS carry HTML snippets) must never be rescanned and have
    # tags inside its own strings "inlined".
    def _either(m):
        return _script(m) if m.group("src") else _css(m)

    html = _CDN_ANY_RE.sub(_either, html)

    verify_offline(html)
    return html


def _round_coords(c, nd):
    if isinstance(c, (int, float)):
        return round(c, nd)
    out = [_round_coords(x, nd) for x in c]
    # drop consecutive duplicate points created by rounding (rings stay closed)
    if out and isinstance(out[0], list) and out[0] and isinstance(out[0][0], (int, float)):
        dedup = [out[0]]
        for pt in out[1:]:
            if pt != dedup[-1]:
                dedup.append(pt)
        if len(dedup) >= 4 or len(dedup) == len(out):
            out = dedup
    return out


def _slim_geojson(fc, nd=3):
    feats = []
    for f in (fc or {}).get("features", []):
        g = f.get("geometry")
        if not g:
            continue
        props = {k: v for k, v in (f.get("properties") or {}).items()
                 if k in _GEO_KEEP_PROPS and v not in (None, "", -99, "-99")}
        feats.append({"type": "Feature", "properties": props,
                      "geometry": {"type": g["type"],
                                   "coordinates": _round_coords(g["coordinates"], nd)}})
    return {"type": "FeatureCollection", "features": feats}


def offline_geo_payload(offline_dir=None, allow_download=False):
    """Boundary GeoJSON for an offline build: {"countries": fc, "admin1": fc}.

    Best effort -- a missing boundary set only disables the views that need it
    (the report already falls back to ranked lists), so it warns, never fails.
    """
    out = {}
    for key, urls in OFFLINE_GEO_SOURCES.items():
        last = None
        for url in urls:
            try:
                fc = json.loads(_fetch_bytes(url, offline_dir, allow_download).decode("utf-8"))
                out[key] = _slim_geojson(fc)
                break
            except Exception as exc:  # noqa: BLE001
                last = exc
        if key not in out:
            print(f"[report_template] WARNING: offline report: '{key}' boundaries not embedded "
                  f"({last}); the metadata choropleth / outline basemap will be unavailable.",
                  file=sys.stderr)
    return out


_REMOTE_TAG_RE = re.compile(
    r'<(?:script\b[^>]*\bsrc|link\b[^>]*\bhref)\s*=\s*["\'](?P<u>(?:https?:)?//[^"\']+)["\']',
    re.IGNORECASE)


def verify_offline(html):
    """Fail loudly if an offline build still loads anything over the network."""
    left = [m.group("u") for m in _REMOTE_TAG_RE.finditer(html)]
    if left:
        raise RuntimeError("offline report still references remote resources: " + ", ".join(left))
    css_remote = re.findall(r'url\(\s*["\']?(https?://[^"\')]+)', html)
    if css_remote:
        print("[report_template] WARNING: offline report: CSS still points at "
              + ", ".join(sorted(set(css_remote))[:5]), file=sys.stderr)


def _default_template():
    return Path(__file__).resolve().parent.parent / "assets" / "heatmap.html"


def inline_template(template_path=None, offline=False, offline_dir=None) -> str:
    """Return the fully self-contained report HTML as a string.

    Local parts under src/ are always folded inline. CDN parts (d3, xlsx, jspdf,
    Leaflet, Font Awesome) are left as external links by default so the report
    loads them at view time. Set offline=True to download and embed them, or
    offline_dir=<path> to embed local copies without any network access
    (offline_dir takes precedence over offline).
    """
    template_path = Path(template_path) if template_path else _default_template()
    with open(template_path, "r", encoding="utf-8", newline="") as fh:
        html = fh.read()

    # Bases to search for referenced parts, in priority order.
    bases = []
    for cand in (template_path.parent,
                 Path(os.path.realpath(template_path)).parent):
        if cand not in bases:
            bases.append(cand)

    def _read_part(rel):
        for base in bases:
            part = base / rel
            if part.is_file():
                with open(part, "r", encoding="utf-8", newline="") as fh:
                    return fh.read()
        sys.exit(f"ERROR: referenced part not found: {rel} "
                 f"(looked in {', '.join(str(b) for b in bases)})")

    def _css(m):
        body = _read_part(m.group("href"))
        sid = _STYLE_ID_RE.search(m.group(0))
        open_tag = f'    <style id="{sid.group(1)}">' if sid else "    <style>"
        return f"{open_tag}\n{body}    </style>\n"

    def _js(m):
        return f'    <script>\n{_read_part(m.group("src"))}    </script>\n'

    html = _LINK_RE.sub(_css, html)
    html = _SCRIPT_RE.sub(_js, html)

    leftover = re.search(r'(?:href|src)="src/[^"]+"', html)
    if leftover:
        sys.exit(f"ERROR: un-inlined local reference remains: {leftover.group(0)}")

    # Offline builds: fold the CDN libraries inline too. Default leaves them as
    # external links (loaded from the CDN when the report is opened).
    #   --offline_report_files DIR                 -> DIR only, no network
    #   --offline_report_files DIR --offline_report -> DIR first, download the rest
    #   --offline_report                            -> bundled assets/offline_report_libs
    #                                                  first (if present), download the rest
    if offline_dir is not None:
        html = inline_cdn_assets(html, offline_dir=offline_dir, allow_download=bool(offline))
    elif offline:
        bundled = BUNDLED_LIBS_DIR if BUNDLED_LIBS_DIR.is_dir() else None
        html = inline_cdn_assets(html, offline_dir=bundled, allow_download=True)
    return html


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-t", "--template", default=None,
                    help="thin template to assemble (default: assets/heatmap.html)")
    ap.add_argument("-o", "--output", help="write here (default: stdout)")
    ap.add_argument("--offline_report", action="store_true",
                    help="download the CDN libraries and embed them inline")
    ap.add_argument("--offline_report_files", default=None, metavar="DIR",
                    help="directory of local CDN library copies to embed inline "
                         "(no network; takes precedence over --offline_report)")
    args = ap.parse_args()
    html = inline_template(args.template, offline=args.offline_report,
                           offline_dir=args.offline_report_files)
    if args.output:
        with open(args.output, "w", encoding="utf-8", newline="") as fh:
            fh.write(html)
        print(f"Wrote {args.output} ({len(html)} bytes)")
    else:
        sys.stdout.write(html)


if __name__ == "__main__":
    main()
