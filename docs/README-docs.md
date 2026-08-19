# Documentation site

The contents of this directory are published to
<https://jhuapl-bio.github.io/taxtriage/> by
[`.github/workflows/deploy-pages.yml`](../.github/workflows/deploy-pages.yml),
using [MkDocs](https://www.mkdocs.org/) with the
[Material](https://squidfunk.github.io/mkdocs-material/) theme.

This replaces the old GitHub wiki. Documentation now lives in the repo and
changes go through normal pull requests, so docs and code review together.

## Tabs

| Tab | Source |
|---|---|
| Documentation | `docs/*.md` |
| Demo Report | generated from `assets/heatmap.html` by `scripts/inline_boot_json.py` |
| Pathogen Sheet | generated from `assets/pathogen_sheet.csv` by `scripts/build_pathogen_data.py` |

The last two are built during deploy and are gitignored (`docs/demo/`,
`docs/data/pathogen_sheet.json`), so they cannot drift from the repo.

## Local preview

```bash
pip install -r requirements-docs.txt

python scripts/inline_boot_json.py            # needs Node on PATH
mkdir -p docs/demo && cp _site/index.html docs/demo/index.html
python scripts/build_pathogen_data.py

mkdocs serve                                   # http://127.0.0.1:8000
```

## Editing

1. Edit or add a `docs/*.md` file. Pages are flat; grouping comes from `nav:`
   in [`mkdocs.yml`](../mkdocs.yml).
2. Add new pages to `nav:`.
3. Link between pages with relative Markdown paths: `[Output](output.md)`.

`mkdocs build --strict` fails on links to missing pages or anchors.

## Checks that run on deploy

```bash
mkdocs build --strict                 # Markdown links, nav, anchors
python scripts/test_built_links.py    # raw-HTML src and runtime fetch paths
node scripts/test_pathogen_table.js   # the table widget, in jsdom (needs: npm i jsdom)
```

`test_built_links.py` exists because `--strict` only validates Markdown links.
Raw HTML (`<img src=...>`, `<iframe src=...>`) and URLs built in JavaScript are
invisible to it, and those resolve against the page's directory URL — a common
source of 404s that only appear once deployed.

## Gotchas

- **Raw HTML asset paths need `../`.** Pages are served at directory URLs
  (`/output/`), so a raw `<img src="images/x.png">` resolves to
  `/output/images/x.png`. Markdown image syntax is rewritten automatically and
  does not need this. `scripts/migrate_wiki.py` applies the fix on import.
- **Heading anchors match GitHub's.** `mkdocs.yml` uses a `pymdownx.slugs`
  slugifier so anchors that people bookmarked from the wiki still resolve.
- **`exclude_docs`** keeps `docs/novelty_sketch/` and reference PDFs in the repo
  but out of the published site.
- **Two typos in `assets/pathogen_sheet.csv`** (`estbalished`, `etsablished` in
  the `status` column, one row each) are normalised by
  `scripts/build_pathogen_data.py`. Fixing them at the source makes that a no-op.
