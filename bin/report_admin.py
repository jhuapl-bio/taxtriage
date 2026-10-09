#!/usr/bin/env python3
# ##############################################################################################
# # Copyright 2022 The Johns Hopkins University Applied Physics Laboratory LLC
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
"""Interactive-report configuration + Admin dialog payload.

Used by make_report.py; also runnable on its own to check a config file:

    python bin/report_admin.py --validate my_report_config.yml

Report config (--report_config, JSON or YAML)
---------------------------------------------
    title: "Outbreak run 42"          # optional banner title
    tabs:
      hidden:   [novelty, insilico]   # removed from the tab bar
      disabled: [explore]             # shown greyed out, not clickable
      show:     [summary, heatmap]    # whitelist (everything else hidden); optional
      default:  heatmap               # tab the report opens on
      order:    [summary, table]      # listed tabs first, the rest keep their order
      labels:   {table: Detections}   # rename tabs
      # A mapping form is accepted too:  tabs: {novelty: hidden, explore: disabled}
    admin:
      enabled: true                   # show the Admin button + dialog
      allow_tab_toggle: true          # viewers may change tab visibility in the dialog
      sections: [tabs, params, run, versions, config]
      redact: [paths, userName]       # blank these run-info fields / param values

Tab ids: see TABS below (aliases such as "vfamr", "metadata", "mapping" work).
"""

import argparse
import json
import os
import re
import sys

# (id, label) in tab-bar order. Must match the data-tab values in assets/heatmap.html.
TABS = [
    ("summary", "Summary"),
    ("heatmap", "Heatmap"),
    ("tass", "TASS Comparison"),
    ("sunburst", "Sunburst"),
    ("coverage", "Coverage"),
    ("proteins", "VF / AMR"),
    ("histogram", "Histograms"),
    ("novelty", "Novelty"),
    ("explore", "Explore"),
    ("table", "Table"),
    ("runmeta", "Metadata"),
    ("map", "Mapping"),
    ("trends", "Trends"),
    ("insilico", "In-Silico"),
]
TAB_IDS = [t for t, _ in TABS]

TAB_ALIASES = {
    "vfamr": "proteins", "vf_amr": "proteins", "vf/amr": "proteins", "amr": "proteins",
    "vf": "proteins", "protein": "proteins",
    "histograms": "histogram", "hist": "histogram",
    # Alignment Trends is a sub-tab of Trends, so its names address that tab.
    "alignment_trends": "trends", "align_trends": "trends", "aligntrends": "trends",
    "depth_trends": "trends",
    "tass_comparison": "tass",
    "metadata": "runmeta", "run_metadata": "runmeta", "meta": "runmeta",
    "mapping": "map", "geo": "map",
    "trend": "trends", "longitudinal": "trends",
    "in_silico": "insilico", "in-silico": "insilico", "simulation": "insilico",
    "detections": "table",
}

ADMIN_SECTIONS = ["tabs", "params", "run", "versions", "config"]
ADMIN_SECTION_ALIASES = {
    "parameters": "params", "param": "params",
    "run_info": "run", "runinfo": "run", "workflow": "run",
    "software": "versions", "software_versions": "versions", "version": "versions",
    "report_config": "config",
    "tab": "tabs",
}

# Run-info fields blanked by the "paths" redact token.
PATH_FIELDS = ["launchDir", "workDir", "projectDir", "outdir", "configFiles",
               "commandLine", "userName"]
REDACTED = "(redacted)"

_HIDE_WORDS = {"hide", "hidden", "off", "false", "no", "remove", "removed"}
_DISABLE_WORDS = {"disable", "disabled", "grey", "gray", "greyed", "inactive"}
_SHOW_WORDS = {"show", "shown", "visible", "on", "true", "yes", "enable", "enabled"}


def default_config():
    return {
        "title": None,
        "tabs": {"hidden": [], "disabled": [], "default": None, "order": [], "labels": {}},
        "admin": {
            "enabled": True,
            "allow_tab_toggle": True,
            "sections": {s: True for s in ADMIN_SECTIONS},
            "redact": [],
        },
        "source": None,
    }


def load_config_file(path):
    """Read a JSON or YAML report config. YAML needs PyYAML (the pipeline
    converts YAML to JSON in Nextflow first, so the container never does)."""
    with open(path, encoding="utf-8") as fh:
        text = fh.read()
    if not text.strip():
        return {}
    low = path.lower()
    if low.endswith(".json"):
        return json.loads(text)
    try:
        return json.loads(text)          # .yml that is really JSON
    except ValueError:
        pass
    try:
        import yaml  # noqa: WPS433
    except ImportError:
        raise SystemExit(
            f"[report_admin] ERROR: {path} looks like YAML but PyYAML is not installed. "
            "Install pyyaml or supply the config as JSON."
        )
    return yaml.safe_load(text)


def _tab_id(name, warnings, where):
    key = str(name).strip().lower().replace(" ", "_")
    key = TAB_ALIASES.get(key, key)
    if key in TAB_IDS:
        return key
    warnings.append(f"{where}: unknown tab '{name}' (known: {', '.join(TAB_IDS)})")
    return None


def _as_list(v):
    if v is None:
        return []
    if isinstance(v, (list, tuple, set)):
        return list(v)
    if isinstance(v, str):
        return [x for x in re.split(r"[,\s]+", v) if x]
    return [v]


def _uniq(seq):
    out = []
    for x in seq:
        if x and x not in out:
            out.append(x)
    return out


def _as_bool(v, default):
    if v is None:
        return default
    if isinstance(v, bool):
        return v
    s = str(v).strip().lower()
    if s in ("1", "true", "yes", "on", "y"):
        return True
    if s in ("0", "false", "no", "off", "n"):
        return False
    return default


def normalize_config(raw, source=None):
    """Validate a raw config mapping. Returns (config, warnings). Unknown keys
    and tab names are reported, never fatal, so a typo cannot sink a report."""
    cfg = default_config()
    cfg["source"] = source
    warnings = []
    if raw is None:
        return cfg, warnings
    if not isinstance(raw, dict):
        raise SystemExit("[report_admin] ERROR: report config must be a mapping at the top level")

    known_top = {"title", "tabs", "admin", "hide_tabs", "disable_tabs", "default_tab"}
    for k in raw:
        if k not in known_top:
            warnings.append(f"unknown top-level key '{k}' ignored")

    if raw.get("title"):
        cfg["title"] = str(raw["title"])

    tabs = raw.get("tabs") or {}
    hidden, disabled, show = [], [], None
    t = cfg["tabs"]

    if isinstance(tabs, list):           # bare list == whitelist
        show = tabs
        tabs = {}
    elif not isinstance(tabs, dict):
        warnings.append("'tabs' must be a mapping or a list; ignored")
        tabs = {}

    structural = {"hidden", "hide", "disabled", "disable", "show", "visible", "only",
                  "default", "order", "labels", "rename"}
    for k, v in tabs.items():
        if k in structural:
            continue
        tid = _tab_id(k, warnings, "tabs")
        if not tid:
            continue
        word = str(v).strip().lower() if not isinstance(v, bool) else ("show" if v else "hide")
        if isinstance(v, dict):          # {explore: {state: disabled, label: X}}
            word = str(v.get("state", "show")).strip().lower()
            if v.get("label"):
                t["labels"][tid] = str(v["label"])
        if word in _HIDE_WORDS:
            hidden.append(tid)
        elif word in _DISABLE_WORDS:
            disabled.append(tid)
        elif word not in _SHOW_WORDS:
            warnings.append(f"tabs.{k}: unrecognised state '{v}' (use show / hidden / disabled)")

    hidden += _as_list(tabs.get("hidden", tabs.get("hide"))) + _as_list(raw.get("hide_tabs"))
    disabled += _as_list(tabs.get("disabled", tabs.get("disable"))) + _as_list(raw.get("disable_tabs"))
    if show is None:
        show = tabs.get("show", tabs.get("visible", tabs.get("only")))

    t["hidden"] = _uniq(_tab_id(x, warnings, "tabs.hidden") for x in hidden)
    t["disabled"] = _uniq(_tab_id(x, warnings, "tabs.disabled") for x in disabled)
    if show is not None:
        keep = set(_uniq(_tab_id(x, warnings, "tabs.show") for x in _as_list(show)))
        t["hidden"] = _uniq(t["hidden"] + [x for x in TAB_IDS if x not in keep])
    t["disabled"] = [x for x in t["disabled"] if x not in t["hidden"]]

    dflt = tabs.get("default", raw.get("default_tab"))
    if dflt:
        tid = _tab_id(dflt, warnings, "tabs.default")
        if tid and (tid in t["hidden"] or tid in t["disabled"]):
            warnings.append(f"tabs.default '{dflt}' is hidden/disabled; the first available tab is used")
            tid = None
        t["default"] = tid
    t["order"] = _uniq(_tab_id(x, warnings, "tabs.order") for x in _as_list(tabs.get("order")))
    labels = tabs.get("labels", tabs.get("rename")) or {}
    if isinstance(labels, dict):
        for k, v in labels.items():
            tid = _tab_id(k, warnings, "tabs.labels")
            if tid and v not in (None, ""):
                t["labels"][tid] = str(v)
    else:
        warnings.append("tabs.labels must be a mapping of tab -> label; ignored")
    if len(t["hidden"]) >= len(TAB_IDS):
        warnings.append("every tab is hidden; un-hiding 'summary' so the report is usable")
        t["hidden"].remove("summary")

    admin = raw.get("admin")
    a = cfg["admin"]
    if isinstance(admin, bool):
        a["enabled"] = admin
    elif isinstance(admin, dict):
        known_admin = {"enabled", "allow_tab_toggle", "sections", "redact"}
        for k in admin:
            if k not in known_admin:
                warnings.append(f"unknown key 'admin.{k}' ignored")
        a["enabled"] = _as_bool(admin.get("enabled"), True)
        a["allow_tab_toggle"] = _as_bool(admin.get("allow_tab_toggle"), True)
        secs = admin.get("sections")
        if secs is not None:
            if isinstance(secs, dict):
                on = {k for k, v in secs.items() if _as_bool(v, True)}
                named = set(secs)
            else:
                on = set(_as_list(secs))
                named = None
            resolved_on = set()
            for s in (named if named is not None else on):
                sid = ADMIN_SECTION_ALIASES.get(str(s).lower(), str(s).lower())
                if sid not in ADMIN_SECTIONS:
                    warnings.append(f"admin.sections: unknown section '{s}' (known: {', '.join(ADMIN_SECTIONS)})")
                elif s in on:
                    resolved_on.add(sid)
            if named is None:
                a["sections"] = {s: (s in resolved_on) for s in ADMIN_SECTIONS}
            else:   # mapping: unnamed sections stay on
                named_ids = {ADMIN_SECTION_ALIASES.get(str(s).lower(), str(s).lower()) for s in named}
                a["sections"] = {s: (s in resolved_on) if s in named_ids else True for s in ADMIN_SECTIONS}
        a["redact"] = [str(x) for x in _as_list(admin.get("redact"))]
    elif admin is not None:
        warnings.append("'admin' must be a mapping or true/false; ignored")
    return cfg, warnings


# ──────────────────────────────────────────────────────────────────────────────
# Software versions (concatenated nf-core versions.yml files)
# ──────────────────────────────────────────────────────────────────────────────
_JUNK_VERSION = re.compile(r"unrecognized option|invalid option|illegal option|command not found|usage:", re.I)


def parse_versions(path):
    """Line parser for concatenated versions.yml files:

        "TAXTRIAGE:KRAKEN2":
            kraken2: 2.1.3
            pigz: 2.6

    No PyYAML needed, and tolerant of the duplicate keys concatenation creates.
    Returns a sorted, de-duplicated list of {process, tool, version}."""
    rows, proc = {}, None
    # key may be quoted and contain ':' ("NFCORE_TAXTRIAGE:TAXTRIAGE:KRAKEN2").
    # An unquoted key ends at the first ':' followed by whitespace / end of line.
    line_re = re.compile(r"""^\s*(?:"([^"]+)"|'([^']+)'|([^\s"'].*?))\s*:(?:\s+(.*?))?\s*$""")
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            m = line_re.match(line.rstrip("\n"))
            if not m:
                continue
            key = (m.group(1) or m.group(2) or m.group(3) or "").strip()
            val = (m.group(4) or "").strip().strip("'\"")
            # A key with no value is a process header -- whatever its indent,
            # since a heredoc that kept its tabs indents the whole file.
            if not val:
                proc = key
                continue
            # Tool names are simple identifiers; anything else is stray
            # `--version` output that leaked into the file.
            if proc is None or not re.match(r"^[\w.+-]+$", key) or key.lower() == "usage":
                continue
            if _JUNK_VERSION.search(val):
                continue
            short = proc.split(":")[-1]
            rows[(short, key, val)] = {"process": short, "process_full": proc, "tool": key, "version": val}
    return sorted(rows.values(), key=lambda r: (r["process"].lower(), r["tool"].lower()))


def load_run_info(path):
    with open(path, encoding="utf-8") as fh:
        data = json.load(fh)
    return data if isinstance(data, dict) else None


def _redact_names(cfg):
    names = set()
    for r in cfg["admin"]["redact"]:
        if r.lower() == "paths":
            names.update(PATH_FIELDS)
        else:
            names.add(r)
    return names


def build_admin_payload(cfg, run_info=None, versions=None, argv=None, disabled=False):
    """The `run_info` block embedded in the report. None when the Admin
    dialog is switched off (--no_admin or admin.enabled: false), so nothing
    about the run environment ends up in the HTML."""
    if disabled or not cfg["admin"]["enabled"]:
        return None
    secs = cfg["admin"]["sections"]
    redact = _redact_names(cfg)

    wf = dict((run_info or {}).get("workflow") or {})
    for k in list(wf):
        if k in redact:
            wf[k] = REDACTED

    groups = []
    if secs.get("params"):
        for g in (run_info or {}).get("param_groups") or []:
            g = dict(g)
            g["params"] = [dict(p, value=REDACTED) if p.get("name") in redact else p
                           for p in (g.get("params") or [])]
            groups.append(g)

    build = {
        "python": sys.version.split()[0],
        "argv": None if ("commandLine" in redact) else " ".join(argv or []),
    }
    return {
        "captured": bool(run_info),
        "captured_at": (run_info or {}).get("captured_at"),
        "workflow": wf if secs.get("run") else {},
        "param_groups": groups,
        "software_versions": (versions or []) if secs.get("versions") else [],
        "report_build": build if secs.get("run") else {},
        "redacted": sorted(redact),
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--validate", metavar="FILE", required=True, help="report config (JSON/YAML) to check")
    a = ap.parse_args(argv)
    cfg, warns = normalize_config(load_config_file(a.validate), os.path.basename(a.validate))
    for w in warns:
        print(f"WARNING: {w}", file=sys.stderr)
    print(json.dumps(cfg, indent=2))
    return 1 if warns else 0


if __name__ == "__main__":
    sys.exit(main())
