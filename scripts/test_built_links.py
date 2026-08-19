#!/usr/bin/env python3
"""
Link-integrity check for the *built* site.

`mkdocs build --strict` validates Markdown links, but it does not see raw HTML
(iframe/img src) or URLs constructed in JavaScript. Those resolve against the
page's directory URL, which is exactly where relative paths go wrong.

This walks every generated HTML page, resolves each local reference relative to
that page's real URL, and asserts the target exists on disk.

Run:  python scripts/test_built_links.py [--site site]
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from urllib.parse import unquote, urljoin, urlsplit

# href/src on any tag.
REF = re.compile(r'(?:href|src)\s*=\s*"([^"]+)"', re.IGNORECASE)

# Same-origin URLs the front-end fetches at runtime; each entry is
# (page the code runs on, path as the browser will resolve it). The pathogen
# sheet is deliberately absent here: it is fetched cross-origin from
# raw.githubusercontent.com at page load and is not part of the built site.
RUNTIME_FETCHES: list[tuple[str, str]] = []

SKIP_PREFIXES = ("http://", "https://", "//", "mailto:", "data:", "javascript:", "#")


def site_base(config: Path) -> str:
    """Extract the URL path from site_url without needing a YAML dependency."""
    if config.is_file():
        for line in config.read_text(encoding="utf-8").splitlines():
            m = re.match(r"\s*site_url:\s*(\S+)", line)
            if m:
                return urlsplit(m.group(1)).path or "/"
    return "/"


def page_url(site: Path, page: Path, base: str) -> str:
    """The URL path a page is served at, e.g. '/taxtriage/demo-report/index.html'."""
    return base + page.relative_to(site).as_posix()


def resolve(site: Path, page: Path, ref: str, base: str) -> Path | None:
    """Resolve a reference the way a browser would; None if it escapes the site."""
    ref = ref.split("#", 1)[0].split("?", 1)[0]
    if not ref:
        return None
    target = urljoin(page_url(site, page, base), ref)
    # Absolute refs include the deployment base path (e.g. /taxtriage/css/...);
    # strip it to get back to a path inside the built tree.
    if not target.startswith(base):
        return None  # escapes the deployment root
    rel = unquote(target[len(base) :])
    path = site / rel
    if rel.endswith("/") or path.is_dir():
        path = path / "index.html"
    return path


def main() -> int:
    here = Path(__file__).resolve().parent.parent
    ap = argparse.ArgumentParser()
    ap.add_argument("--site", type=Path, default=here / "site")
    ap.add_argument(
        "--base",
        default=None,
        help="URL path the site is served from (default: taken from site_url in mkdocs.yml)",
    )
    args = ap.parse_args()

    base = args.base if args.base is not None else site_base(here / "mkdocs.yml")
    if not base.startswith("/"):
        base = "/" + base
    if not base.endswith("/"):
        base += "/"

    site: Path = args.site
    if not site.is_dir():
        print(f"error: built site not found at {site}. Run `mkdocs build` first.", file=sys.stderr)
        return 1

    pages = sorted(site.rglob("*.html"))
    if not pages:
        print(f"error: no HTML pages under {site}", file=sys.stderr)
        return 1

    broken: list[tuple[str, str, str]] = []
    checked = 0

    for page in pages:
        # The embedded report is a 4MB third-party artifact with its own
        # internal asset graph; it is validated by its own build, not here.
        if page.relative_to(site).as_posix().startswith("demo/"):
            continue

        html = page.read_text(encoding="utf-8", errors="ignore")
        for ref in REF.findall(html):
            if ref.startswith(SKIP_PREFIXES):
                continue
            target = resolve(site, page, ref, base)
            if target is None:
                continue
            checked += 1
            if not target.exists():
                broken.append((page_url(site, page, base), ref, str(target)))

    # Runtime fetches, which never appear as href/src in the markup.
    for page_rel, ref in RUNTIME_FETCHES:
        page = site / page_rel
        if not page.exists():
            broken.append((page_rel, "(page missing)", page_rel))
            continue
        target = resolve(site, page, ref, base)
        checked += 1
        if not target.exists():
            broken.append((page_url(site, page, base), ref + "  [runtime fetch]", str(target)))

    print(f"checked {checked} local references across {len(pages)} pages (base {base})")
    if broken:
        print(f"\n{len(broken)} BROKEN:")
        for src, ref, resolved in broken:
            print(f"  {src}\n    ref:      {ref}\n    resolves: {resolved}  (missing)")
        return 1

    print("all local references resolve")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
