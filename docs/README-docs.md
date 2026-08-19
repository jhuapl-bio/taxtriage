# Documentation site

Source for <https://jhuapl-bio.github.io/taxtriage/>, built with
[MkDocs](https://www.mkdocs.org/) + [Material](https://squidfunk.github.io/mkdocs-material/)
and published with [mike](https://github.com/jimporter/mike), which keeps one
rendered copy per version on the `gh-pages` branch.

This replaces the old GitHub wiki. Docs live in the repo and change through
pull requests, so docs and code review together.

## Versions

| Version         | Built from      | Trigger                     |
| --------------- | --------------- | --------------------------- |
| `latest`        | `main`          | every push to `main`        |
| `3.3`, `3.2`, … | the release tag | publishing a GitHub Release |

`latest` tracks `main` and is the default, so the site root redirects there.
Releasing `v3.3.9` publishes as **`3.3`**, overwriting `v3.3.8`'s docs — patch
releases update their minor series in place, so the version selector stays
short. Releases never move the default.

`latest` is a real version, not an alias. mike keeps versions and aliases in a
single namespace, so a release cannot also claim `latest` as an alias — it
fails with `alias 'latest' already specified as a version`. The workflow
guards against a tag that would resolve to `latest` and fails loudly.

`gh-pages` contains **only rendered HTML/CSS/JS** — no source, no `mkdocs.yml`,
no scripts. It is entirely machine-generated; never commit to it by hand.

## Branch model

Work happens on `dev` via pull requests, then `dev` merges into `main`.
`.github/workflows/docs-check.yml` builds and tests the site on PRs into `dev`
and `main` without deploying, so breakage is caught in review.
`.github/workflows/docs.yml` does the actual deploying.

## The three tabs

| Tab            | Source                                                  | Refreshes              |
| -------------- | ------------------------------------------------------- | ---------------------- |
| Documentation  | `docs/*.md`                                             | on deploy              |
| Demo Report    | `assets/heatmap.html` via `scripts/inline_boot_json.py` | on deploy              |
| Pathogen Sheet | `assets/pathogen_sheet.csv`                             | **on every page load** |

The demo dist and the pinned-ref file are generated at build time and
gitignored (`docs/demo/`, `docs/javascripts/docs_ref.js`).

The pathogen sheet is not bundled at all — the browser fetches it from GitHub
when the page opens. `scripts/write_docs_ref.py` pins _which ref_ it reads, so
the `3.1` docs show the sheet that shipped with 3.1 while `latest` shows `main`.
That file is loaded before the table widget and sets `window.TAXTRIAGE_DOCS`.

## Local preview

```bash
pip install -r requirements-docs.txt

python scripts/write_docs_ref.py                       # defaults to current branch
python scripts/inline_boot_json.py                     # needs Node on PATH
mkdir -p docs/demo && cp _site/index.html docs/demo/index.html

mkdocs serve                                           # http://127.0.0.1:8000
```

`mkdocs serve` renders a single unversioned copy — the version selector only
appears on the deployed site. To preview the versioned layout:

```bash
mike deploy 0.0-test && mike serve
```

## Editing

1. Edit or add `docs/*.md`. Pages are flat; grouping comes from `nav:` in
   [`mkdocs.yml`](../mkdocs.yml).
2. Add new pages to `nav:`.
3. Link between pages with relative Markdown paths: `[Output](output.md)`.

## Checks

```bash
mkdocs build --strict                 # Markdown links, nav, anchors
python scripts/test_built_links.py    # raw-HTML src paths in the built site
node scripts/test_pathogen_table.js   # 60 checks, real widget in jsdom
```

`test_built_links.py` exists because `--strict` only validates Markdown links.
Raw HTML (`<img src>`, `<iframe src>`) and JS-built URLs are invisible to it,
and those resolve against the page's directory URL — a common source of 404s
that only appear once deployed.

`test_pathogen_table.js` builds its ground truth with a second, independent CSV
reader so a bug in the parser under test cannot make the assertions agree with
themselves. Override the input with `PATHOGEN_CSV=/path/to/sheet.csv`.

## Gotchas

- **Raw HTML asset paths need `../`.** Pages are served at directory URLs
  (`/output/`), so raw `<img src="images/x.png">` resolves to
  `/output/images/x.png`. Markdown image syntax is rewritten automatically.
- **Heading anchors match GitHub's**, via a `pymdownx.slugs` slugifier, so
  anchors bookmarked from the old wiki still resolve.
- **The table tolerates schema drift.** Any column the CSV lacks reads as blank
  and its facet is hidden. Two known upstream typos in `status`
  (`estbalished`, `etsablished`) are folded into `established` in the browser;
  the CSV is never modified.
- **`exclude_docs`** keeps `docs/novelty_sketch/` and reference PDFs in the repo
  but out of the published site.

## First-time setup

Settings → Pages → Source → **Deploy from a branch** → `gh-pages` / `/ (root)`.
The `gh-pages` branch appears after the first successful deploy, so run the
workflow once (or push to `main`) before setting this.
