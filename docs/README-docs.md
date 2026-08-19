# Documentation site

Source for <https://jhuapl-bio.github.io/taxtriage/>, built with
[MkDocs](https://www.mkdocs.org/) + [Material](https://squidfunk.github.io/mkdocs-material/)
and published with [mike](https://github.com/jimporter/mike), which keeps one
rendered copy per version on the `gh-pages` branch.

This replaces the old GitHub wiki. Docs live in the repo and change through
pull requests, so docs and code review together.

## Versions

| Version             | Built from      | Trigger                             |
| ------------------- | --------------- | ----------------------------------- |
| `latest`            | `main`          | every push to `main`                |
| `3.3.9`, `3.3.8`, … | the release tag | publishing **or editing** a Release |

`latest` tracks `main` and is the default, so the site root redirects there.
Releasing `v3.3.9` publishes version **`3.3.9`** — the exact patch, so the
dropdown shows `3.3.9` rather than `3.3`. There are no minor-series aliases.
Releases never move the default.

`latest` is a real version, not an alias. mike keeps versions and aliases in a
single namespace, so a release cannot also claim `latest` as an alias — it
fails with `alias 'latest' already specified as a version`. The workflow
guards against a tag that would resolve to `latest` and fails loudly.

Every deploy prunes entries left over from earlier versioning schemes:
`bleeding-edge`, and any bare `X.Y` minor entry. Matching on shape rather than
a hard-coded list means a stale entry cannot linger; `X.Y.Z` is never matched,
so real release docs are safe.

### Rebuilding a version

Retargeting a GitHub Release at a newer commit fires `release: edited`, which
the workflow listens for — it rebuilds that version from the tag's new commit
and overwrites the existing directory in place. No duplicate entry is created.

Force-moving a tag with `git push -f` alone does **not** fire a release event.
To rebuild in that case (or any other), run the Docs workflow manually:

| Input     | Value                                                  |
| --------- | ------------------------------------------------------ |
| `version` | `3.3.9` (blank rebuilds `latest` from `main`)          |
| `ref`     | blank uses `v<version>`; set it to pin a commit or ref |

Because `edited` covers every release edit — including title and body — a
cosmetic tweak also triggers a rebuild. That is harmless: the deploy is
idempotent.

`gh-pages` contains **only rendered HTML/CSS/JS** — no source, no `mkdocs.yml`,
no scripts. It is entirely machine-generated; never commit to it by hand.

## Branch model

Work happens on `dev` via pull requests, then `dev` merges into `main`.
`.github/workflows/docs-check.yml` builds and tests the site on PRs into `dev`
and `main` without deploying, so breakage is caught in review.
`.github/workflows/docs.yml` does the actual deploying.

## The three tabs

| Tab            | Source                                                  | Refreshes              | Published in  |
| -------------- | ------------------------------------------------------- | ---------------------- | ------------- |
| Documentation  | `docs/*.md`                                             | on deploy              | every version |
| Pathogen Sheet | `assets/pathogen_sheet.csv`                             | **on every page load** | every version |
| Demo Report    | `assets/heatmap.html` via `scripts/inline_boot_json.py` | on deploy              | `latest` only |

The **Pathogen Sheet** ships in every version. The CSV is fetched from GitHub
when the page opens rather than bundled, and `scripts/write_docs_ref.py` pins
which ref it reads — so the `3.3.9` docs show the sheet that shipped in v3.3.9
while `latest` shows `main`. Costs nothing to publish per version.

The **Demo Report** is `latest`-only: its ~4.5 MB dist would add megabytes to
`gh-pages` per release, and it illustrates the report UI rather than documenting
release-specific behaviour. In a versioned build
`scripts/docs_release_stub.py` replaces it with a redirect, so the tab stays in
the nav and old URLs keep working:

```
/taxtriage/3.3.9/demo-report/  ->  /taxtriage/latest/demo-report/
```

The redirect is relative (`../../latest/<slug>/`), so it survives a move to a
custom domain. A release directory is ~8 MB rather than ~13 MB.

The demo dist and the pinned-ref file are generated at build time and
gitignored (`docs/demo/`, `docs/javascripts/docs_ref.js`).

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
mkdocs build --strict                     # Markdown links, nav, anchors
python scripts/test_built_links.py        # raw-HTML src paths in the built site
node scripts/test_pathogen_table.js       # 60 checks, real widget in jsdom
python scripts/docs_release_stub.py --check   # release build can still redirect
```

CI builds both shapes (`latest` and `release`) on every docs PR. The release
shape is otherwise only exercised when a release is published, which is the
worst time to find out it is broken.

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
- **GitHub ````math` fences do not work here.** They are a GitHub-only
  extension; Python-Markdown prints them verbatim. Use `$$…$$`, which
  `pymdownx.arithmatex` renders via MathJax. `scripts/migrate_wiki.py`
  converts them on import.
- **`exclude_docs`** keeps `docs/novelty_sketch/` and reference PDFs in the repo
  but out of the published site.

## First-time setup

Settings → Pages → Source → **Deploy from a branch** → `gh-pages` / `/ (root)`.
The `gh-pages` branch appears after the first successful deploy, so run the
workflow once (or push to `main`) before setting this.

Only tags cut _after_ the docs landed in the repo can be published: older tags
predate `docs/` and `mkdocs.yml`, so there is nothing to build. The version
selector therefore fills in going forward rather than being backfilled.

If a version ever needs removing by hand:

```bash
pip install -r requirements-docs.txt
mike list   --branch gh-pages
mike delete --branch gh-pages --push <version>
```
