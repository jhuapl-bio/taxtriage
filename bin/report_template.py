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
import os
import re
import sys
from pathlib import Path

_LINK_RE = re.compile(
    r'[ \t]*<link\b[^>]*\bhref="(?P<href>src/[^"]+\.css)"[^>]*>[ \t]*\n?',
    re.IGNORECASE,
)
_SCRIPT_RE = re.compile(
    r'[ \t]*<script\b[^>]*\bsrc="(?P<src>src/[^"]+\.js)"[^>]*>\s*</script>[ \t]*\n?',
    re.IGNORECASE,
)
_STYLE_ID_RE = re.compile(r'\bdata-style-id="([^"]+)"')


def _default_template():
    return Path(__file__).resolve().parent.parent / "assets" / "heatmap.html"


def inline_template(template_path=None) -> str:
    """Return the fully self-contained report HTML as a string."""
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
    return html


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-t", "--template", default=None,
                    help="thin template to assemble (default: assets/heatmap.html)")
    ap.add_argument("-o", "--output", help="write here (default: stdout)")
    args = ap.parse_args()
    html = inline_template(args.template)
    if args.output:
        with open(args.output, "w", encoding="utf-8", newline="") as fh:
            fh.write(html)
        print(f"Wrote {args.output} ({len(html)} bytes)")
    else:
        sys.stdout.write(html)


if __name__ == "__main__":
    main()
