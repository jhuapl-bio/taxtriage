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
# url(...) references inside CSS (fonts, marker images). Skips data: URIs.
_CSS_URL_RE = re.compile(r'url\(\s*(?P<q>["\']?)(?P<u>[^"\')]+)(?P=q)\s*\)')

# Extensions -> MIME for data-URI embedding of CSS-referenced assets.
_ASSET_MIME = {
    ".woff2": "font/woff2", ".woff": "font/woff", ".ttf": "font/ttf",
    ".otf": "font/otf", ".eot": "application/vnd.ms-fontobject",
    ".svg": "image/svg+xml", ".png": "image/png", ".gif": "image/gif",
    ".jpg": "image/jpeg", ".jpeg": "image/jpeg", ".webp": "image/webp",
}


def _fetch_bytes(url, offline_dir):
    """Return the raw bytes for a CDN asset.

    offline_dir set  -> resolve by basename under that directory (fully offline,
                        no network). Searched recursively so fonts/images placed
                        anywhere inside the directory are found.
    offline_dir None -> download the URL at build time.
    """
    if offline_dir:
        name = os.path.basename(urlparse(url).path)
        for root, _dirs, files in os.walk(offline_dir):
            if name in files:
                with open(os.path.join(root, name), "rb") as fh:
                    return fh.read()
        raise FileNotFoundError(
            f"offline asset '{name}' not found under {offline_dir} (for {url})"
        )
    req = Request(url, headers={"User-Agent": "taxtriage-offline-report"})
    with urlopen(req, timeout=60) as resp:  # nosec B310 - fixed https CDN URLs
        return resp.read()


def _inline_css_urls(css_text, base_url, offline_dir):
    """Replace url(...) refs in a stylesheet with base64 data URIs.

    Fonts (Font Awesome) and marker images (Leaflet) are pulled the same way as
    the stylesheet itself. Unresolvable refs are left as-is so a partial offline
    bundle still produces a working (if icon-less) report.
    """
    def repl(m):
        raw = m.group("u").strip()
        if not raw or raw.startswith("data:"):
            return m.group(0)
        clean = raw.split("?", 1)[0].split("#", 1)[0]
        try:
            data = _fetch_bytes(
                clean if offline_dir else urljoin(base_url, clean), offline_dir
            )
        except Exception:  # noqa: BLE001 - best effort, keep original ref
            return m.group(0)
        ext = os.path.splitext(urlparse(clean).path)[1].lower()
        mime = _ASSET_MIME.get(ext, "application/octet-stream")
        b64 = base64.b64encode(data).decode("ascii")
        return f"url(data:{mime};base64,{b64})"

    return _CSS_URL_RE.sub(repl, css_text)


def inline_cdn_assets(html, offline_dir=None) -> str:
    """Fold CDN <script src>/<link href> tags inline for an offline report.

    offline_dir None -> download each library at build time.
    offline_dir set  -> read local copies (by basename) from that directory.
    """
    def _script(m):
        body = _fetch_bytes(m.group("src"), offline_dir).decode("utf-8", "replace")
        return f"    <script>\n{body}\n    </script>\n"

    def _css(m):
        text = _fetch_bytes(m.group("href"), offline_dir).decode("utf-8", "replace")
        text = _inline_css_urls(text, m.group("href"), offline_dir)
        return f"    <style>\n{text}\n    </style>\n"

    html = _CDN_SCRIPT_RE.sub(_script, html)
    html = _CDN_LINK_RE.sub(_css, html)
    return html


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
    if offline_dir is not None:
        html = inline_cdn_assets(html, offline_dir=offline_dir)
    elif offline:
        html = inline_cdn_assets(html, offline_dir=None)
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
