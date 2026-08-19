#!/usr/bin/env python3
"""
Migrate the TaxTriage GitHub wiki into the MkDocs `docs/` tree.

Usage:
    python scripts/migrate_wiki.py [--wiki PATH] [--docs PATH]

Idempotent: safe to re-run. Rewrites GitHub-wiki-style links
(`[Text](Page-Name)`, `[Text](Page-Name#anchor)`) into MkDocs relative
Markdown links (`[Text](page-name.md)`), copies images, and drops
wiki-only files (_Sidebar.md, _Footer.md).
"""

from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path

# Wiki page name -> destination filename in docs/
PAGE_MAP = {
    "Home": "index.md",
    "Installation": "installation.md",
    "Quick-Start": "quick-start.md",
    "Samplesheet": "samplesheet.md",
    "Pre-Aligned-Input": "pre-aligned-input.md",
    "Running-the-Pipeline": "running-the-pipeline.md",
    "CLI-Parameters": "cli-parameters.md",
    "Pipeline-Modules": "pipeline-modules.md",
    "Output": "output.md",
    "Interactive-Report": "interactive-report.md",
    "TASS-Scoring": "tass-scoring.md",
    "Microbial-Categories": "microbial-categories.md",
    "Novelty-Detection": "novelty-detection.md",
    "Detection-Rescue": "detection-rescue.md",
    "In-Silico": "in-silico.md",
    "Cloud-and-Seqera": "cloud-and-seqera.md",
    "Geneious-Plugin": "geneious-plugin.md",
    "Troubleshooting": "troubleshooting.md",
    "Contributing": "contributing.md",
    "Citations": "citations.md",
}

# Stale / renamed targets found in the wiki that don't match PAGE_MAP.
ALIASES = {
    "insilico_simulation.md": "in-silico.md",
    "Insilico": "In-Silico",
    "In-Silico-Simulation": "In-Silico",
    "Cloud-&-Seqera": "Cloud-and-Seqera",
}

SKIP_FILES = {"_Sidebar.md", "_Footer.md", "_Header.md"}

# Links the wiki wrote as if the pipeline repo were a sibling directory.
# On the docs site these must become absolute GitHub URLs.
REPO_BLOB = "https://github.com/jhuapl-bio/taxtriage/blob/main/"
REPO_RELATIVE_PREFIXES = ("../assets/", "../bin/", "../conf/", "../docs/", "../modules/")

# Anything with a scheme, a leading slash, or a fragment-only target is left alone.
EXTERNAL = re.compile(r"^(?:[a-z][a-z0-9+.-]*:|//|/|#|mailto:)", re.IGNORECASE)

LINK = re.compile(r"(?<!\!)\[([^\]]*)\]\(([^)\s]+)(\s+\"[^\"]*\")?\)")

# Raw HTML in Markdown is passed through untouched by MkDocs — unlike Markdown
# links, its src is NOT rewritten. Pages are served at directory URLs
# (/output/), so "images/x.png" would resolve to /output/images/x.png. Docs
# pages are flat, so stepping up one level is always correct.
HTML_ASSET = re.compile(r'((?:src|href)\s*=\s*")(images/)', re.IGNORECASE)


def rewrite_target(target: str) -> str:
    """Map a single link target to its MkDocs equivalent."""
    target = ALIASES.get(target, target)

    if EXTERNAL.match(target):
        return target
    for prefix in REPO_RELATIVE_PREFIXES:
        if target.startswith(prefix):
            return REPO_BLOB + target[len("../") :]
    # Already a real file reference (image, asset, .md with a path) -> leave it.
    if target.startswith(("images/", "assets/", "./", "../")):
        return target

    page, sep, anchor = target.partition("#")
    page = ALIASES.get(page, page)

    if page in PAGE_MAP:
        return PAGE_MAP[page] + (sep + anchor if sep else "")
    if page.endswith(".md"):
        return target
    return target


def rewrite_links(text: str, unresolved: set[str]) -> str:
    def repl(m: re.Match) -> str:
        label, target, title = m.group(1), m.group(2), m.group(3) or ""
        new = rewrite_target(target)
        if (
            new == target
            and not EXTERNAL.match(target)
            and not target.startswith(("images/", "assets/", "./", "../"))
            and not target.endswith((".md", ".png", ".svg", ".jpg", ".pdf", ".csv"))
        ):
            unresolved.add(target)
        return f"[{label}]({new}{title})"

    return LINK.sub(repl, text)


def strip_leading_h1(text: str, keep: bool) -> str:
    """MkDocs Material shows the nav title; duplicate H1s are fine, so keep by default."""
    if keep:
        return text
    lines = text.splitlines()
    for i, line in enumerate(lines):
        if line.startswith("# "):
            return "\n".join(lines[:i] + lines[i + 1 :]).lstrip("\n")
        if line.strip():
            break
    return text


def main() -> int:
    here = Path(__file__).resolve().parent.parent
    ap = argparse.ArgumentParser()
    ap.add_argument("--wiki", type=Path, default=here.parent / "taxtriage.wiki")
    ap.add_argument("--docs", type=Path, default=here / "docs")
    args = ap.parse_args()

    wiki: Path = args.wiki
    docs: Path = args.docs

    if not wiki.is_dir():
        print(f"error: wiki not found at {wiki}", file=sys.stderr)
        return 1

    docs.mkdir(parents=True, exist_ok=True)

    unresolved: set[str] = set()
    written = 0

    for md in sorted(wiki.glob("*.md")):
        if md.name in SKIP_FILES:
            continue
        stem = md.stem
        dest_name = PAGE_MAP.get(stem)
        if dest_name is None:
            dest_name = stem.lower().replace("_", "-") + ".md"
            print(f"  note: {md.name} not in PAGE_MAP, writing as {dest_name}")

        text = md.read_text(encoding="utf-8")
        text = rewrite_links(text, unresolved)
        text = HTML_ASSET.sub(r"\1../\2", text)
        text = text.replace("\r\n", "\n").rstrip() + "\n"

        (docs / dest_name).write_text(text, encoding="utf-8")
        written += 1

    # Images
    src_images = wiki / "images"
    if src_images.is_dir():
        dst_images = docs / "images"
        # dirs_exist_ok rather than rmtree+copy: the docs tree may live on a
        # mount that disallows unlink, and overwriting in place works there.
        shutil.copytree(
            src_images,
            dst_images,
            dirs_exist_ok=True,
            ignore=shutil.ignore_patterns(".DS_Store", "Thumbs.db"),
        )
        n = sum(1 for _ in dst_images.rglob("*") if _.is_file())
        print(f"copied {n} image files")

    print(f"wrote {written} pages to {docs}")
    if unresolved:
        print("\nunresolved link targets (check these):")
        for t in sorted(unresolved):
            print(f"  - {t}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
