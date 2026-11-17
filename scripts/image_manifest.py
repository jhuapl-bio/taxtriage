#!/usr/bin/env python3
"""
image_manifest.py -- list (and optionally pre-pull) the container images TaxTriage
needs for a given run mode, resolved by Nextflow itself.

How it works
------------
For each mode it runs two Nextflow commands against the pipeline (nothing executes):

  1. `nextflow run <pipeline> -preview -with-dag dag.mmd <mode args>`
       -> the processes the workflow actually wires in for those params
          (every `if (params.x)` branch is honoured).
  2. `nextflow inspect <pipeline> -profile <engine> -format json`
       -> the container each process resolves to for that engine
          (docker / singularity / apptainer / podman).

Process aliases (`include { X as Y }`) are mapped back to their module, and the
per-alias `ext.singularity_pull_docker_container` overrides in conf/*.config are
applied, so a Singularity manifest matches what a real run pulls.

Processes are tagged:
  core            -> runs for every sample in this mode
  data-dependent  -> wired in, but only gets tasks if the samplesheet has matching
                     rows (trim=TRUE, SRA accessions, BAM inputs, controls, ...)

Usage
-----
  # from the repo root; default = the four preset modes, docker
  python3 scripts/image_manifest.py

  # Singularity, light modes only, core images only, and pre-pull them
  python3 scripts/image_manifest.py --engine singularity \
      --mode light-k2 --mode light-fasta --core-only \
      --pull --cache-dir /shared/nxf_singularity_cache

  # your own combination of switches as a named mode
  python3 scripts/image_manifest.py --custom "denovo=--db {DB} --use_denovo --annotate"

  # every image any process could ever use (no preview, just inspect)
  python3 scripts/image_manifest.py --mode all

Preset modes (placeholders {DB} {FASTA} {TAXDUMP} {ASSEMBLY} are dummy inputs
created in a temp dir, so you do not need real databases to build a manifest):

  full          default params, Kraken2                    --db {DB}
  light-k2      -profile light, Kraken2 picks references   --db {DB}
  light-fasta   -profile light, local FASTA only           --skip_kraken2 --reference_fasta {FASTA}
  light-both    -profile light, Kraken2 + local FASTA      --db {DB} --reference_fasta {FASTA}
  lite          -profile lite: MCE settings, Kraken2 standard (DOWNLOAD_DB), no FASTA
  all           every process defined in the pipeline (inspect only)

--offline adds `--taxdump {TAXDUMP} --assembly {ASSEMBLY}` to every mode, which
removes the DOWNLOAD_TAXDUMP / GET_ASSEMBLIES images (supply those files yourself).

Outputs (in --outdir, default ./image_manifest):
  manifest.tsv               one row per (mode, process): engine, tier, condition, container
  images_<mode>.txt          unique images for that mode (data-dependent ones marked)
  images_union.txt           unique images across all requested modes
  pull_<engine>.sh           the pull commands (also run directly with --pull)

Requires: python3 (stdlib), nextflow on PATH (or --nextflow), plus docker/podman/
singularity/apptainer only if you use --pull. Needs network for Nextflow's plugin
(nf-schema) on first use, same as a normal run.
"""
import argparse
import gzip
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from collections import OrderedDict, defaultdict

# ── Mode presets ──────────────────────────────────────────────────────────────
PRESETS = OrderedDict([
    ("full",        {"profiles": [],        "args": "--db {DB}"}),
    ("light-k2",    {"profiles": ["light"], "args": "--db {DB}"}),
    ("light-fasta", {"profiles": ["light"], "args": "--skip_kraken2 --reference_fasta {FASTA}"}),
    ("light-both",  {"profiles": ["light"], "args": "--db {DB} --reference_fasta {FASTA}"}),
    ("lite",        {"profiles": ["lite"],  "args": ""}),   # Kraken2 standard (downloaded), no ref FASTA
    ("all",         {"profiles": [],        "args": "--db {DB}", "inspect_only": True}),
])
DEFAULT_MODES = ["full", "light-k2", "light-fasta", "light-both"]

# ── Processes that are wired in but only get tasks for certain samplesheet rows ─
# (Nextflow's preview DAG is static, so this part is annotation. Anything not
# listed here is treated as `core` for the mode it appears in.)
DATA_DEPENDENT = {
    "GENERATE_SAMPLESHEET":        "--fastq_1 / --bam given instead of --input",
    "SRA_RESOLVE":                 "samplesheet rows with SRA/ENA accessions",
    "SRA_FETCH_ENA":               "samplesheet rows with SRA/ENA accessions",
    "SRA_FETCH_SRATOOLS":          "samplesheet rows with SRA/ENA accessions",
    "ARTIC_GUPPYPLEX":             "rows whose fastq_1 is a directory of ONT FASTQs",
    "PIGZ_COMPRESS":               "uncompressed FASTQ inputs",
    "TRIMGALORE":                  "ILLUMINA rows with trim=TRUE (off under -profile light)",
    "PORECHOP":                    "OXFORD/PACBIO rows with trim=TRUE (off under -profile light)",
    "PYCOQC":                      "OXFORD rows with a sequencing_summary",
    "PREPARE_BAM":                 "BAM/CRAM input rows",
    "BAM_CONSENSUS":               "BAM/CRAM input rows without a reference FASTA",
    "SAMTOOLS_COVERAGE_BAM":       "BAM/CRAM input rows",
    "SAMTOOLS_HIST_COVERAGE_BAM":  "BAM/CRAM input rows",
    "BEDTOOLS_GENOMECOVERAGE_BAM": "BAM/CRAM input rows with --fast",
    "DECOMPRESS_CRAM_REFERENCE":   "gzipped --cram_reference",
    "MAP_SAMPLE_FASTA_TO_ASSEMBLY": "rows with a per-sample `fasta` column",
    "BOWTIE2_BUILD_SAMPLE":        "rows with a per-sample `fasta` column + --use_bt2",
    "BUILD_SAMPLE_MASK":           "rows with mask_fasta/mask_bed columns, or --mask_products",
    "SAMTOOLS_MERGE":              "a sample aligned to >1 reference file (e.g. Kraken2 refs + --reference_fasta, or several FASTAs)",
    "ALIGNMENT_PER_SAMPLE_CONTROLS": "negative/positive control rows",
    "ALIGNMENT_PER_SAMPLE_INSILICO": "simulated datasets (--generate_iss/--generate_nanosim/--background_reads/--spikein_sheet)",
    "SINGLE_REPORT_INSILICO":      "simulated datasets",
    "INDEX_HOPPING":               "--index_hopping or an index_hopping samplesheet column",
    "FILTER_DB_DOWNLOAD":          "--filter_kraken2 with a named db",
}

DUMMY_K2 = ["hash.k2d", "opts.k2d", "taxo.k2d"]


def log(msg):
    print(msg, file=sys.stderr, flush=True)


def run(cmd, cwd, env=None):
    log("  $ " + " ".join(shlex.quote(c) for c in cmd))
    p = subprocess.run(cmd, cwd=cwd, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        log(p.stdout[-3000:])
        log(p.stderr[-3000:])
        tail = os.path.join(cwd, ".nextflow.log")
        raise SystemExit(f"command failed (rc={p.returncode}); see {tail}")
    return p.stdout


# ── Dummy inputs so the workflow can be parsed without real data/databases ─────
def make_dummy_inputs(tmp):
    d = os.path.join(tmp, "inputs")
    os.makedirs(d, exist_ok=True)
    fq = "@r1\nACGTACGTACGT\n+\nIIIIIIIIIIII\n"
    for n in ("r1.fastq.gz", "r2.fastq.gz", "ont.fastq.gz"):
        with gzip.open(os.path.join(d, n), "wt") as fh:
            fh.write(fq)
    ss = os.path.join(d, "samplesheet.csv")
    with open(ss, "w") as fh:
        fh.write("sample,platform,fastq_1,fastq_2,trim\n")
        fh.write(f"S1,ILLUMINA,{d}/r1.fastq.gz,{d}/r2.fastq.gz,FALSE\n")
        fh.write(f"S2,OXFORD,{d}/ont.fastq.gz,,FALSE\n")
    k2 = os.path.join(d, "k2db")
    os.makedirs(k2, exist_ok=True)
    for n in DUMMY_K2:
        open(os.path.join(k2, n), "a").close()
    fa = os.path.join(d, "ref.fasta")
    with open(fa, "w") as fh:
        fh.write(">NC_000000.1 dummy\nACGTACGTACGTACGT\n")
    tx = os.path.join(d, "taxdump")
    os.makedirs(tx, exist_ok=True)
    for n in ("nodes.dmp", "names.dmp"):
        open(os.path.join(tx, n), "a").close()
    asm = os.path.join(d, "assembly_summary_refseq.txt")
    open(asm, "a").close()
    return {"INPUT": ss, "DB": k2, "FASTA": fa, "TAXDUMP": tx, "ASSEMBLY": asm}


# ── Static parsing of the pipeline source ─────────────────────────────────────
def nf_files(pipeline):
    for sub in ("workflows", "subworkflows", "modules"):
        for root, _, files in os.walk(os.path.join(pipeline, sub)):
            for f in files:
                if f.endswith(".nf"):
                    yield os.path.join(root, f)
    yield os.path.join(pipeline, "main.nf")


def alias_map(pipeline):
    """alias -> base process/workflow name, from every `include { A as B; C }`."""
    amap = {}
    inc = re.compile(r"include\s*\{([^}]*)\}\s*from", re.S)
    for f in nf_files(pipeline):
        if not os.path.exists(f):
            continue
        for m in inc.finditer(open(f).read()):
            for part in m.group(1).split(";"):
                part = part.strip()
                mm = re.match(r"^([A-Za-z0-9_]+)(?:\s+as\s+([A-Za-z0-9_]+))?$", part)
                if mm and mm.group(2):
                    amap[mm.group(2)] = mm.group(1)
    return amap


def process_files(pipeline):
    """process name -> module file defining it."""
    out = {}
    for f in nf_files(pipeline):
        if not os.path.exists(f):
            continue
        for m in re.finditer(r"^\s*process\s+([A-Za-z0-9_]+)\s*\{", open(f).read(), re.M):
            out[m.group(1)] = f
    return out


def container_branches(module_file):
    """('singularity-branch', 'docker-branch') from the nf-core style ternary, or None."""
    src = open(module_file).read()
    # drop block comments so commented-out container lines are ignored
    src = re.sub(r"/\*.*?\*/", "", src, flags=re.S)
    m = re.search(r"container\s+\"\$\{.*?\?\s*'([^']+)'\s*:\s*'([^']+)'\s*\}\"", src, re.S)
    return (m.group(1), m.group(2)) if m else None


def pull_docker_selectors(pipeline):
    """withName selectors that set ext.singularity_pull_docker_container = true."""
    sels = []
    for cfg in ("conf/modules.config", "conf/base.config"):
        path = os.path.join(pipeline, cfg)
        if not os.path.exists(path):
            continue
        s = open(path).read()
        for m in re.finditer(r"withName\s*:\s*(['\"]?)([^'\"{\n]+?)\1\s*\{", s):
            i, depth = m.end(), 1
            while i < len(s) and depth:
                depth += {"{": 1, "}": -1}.get(s[i], 0)
                i += 1
            if re.search(r"singularity_pull_docker_container\s*=\s*true", s[m.end():i]):
                sels.append(m.group(2).strip())
    return sels


def selector_hits(sels, name):
    for sel in sels:
        try:
            if re.fullmatch(sel, name):
                return True
        except re.error:
            if sel == name:
                return True
    return False


# ── Nextflow calls ────────────────────────────────────────────────────────────
def mode_cmdline(mode_args, engine, profiles, ph, offline, extra):
    args = mode_args
    if offline:
        args += " --taxdump {TAXDUMP} --assembly {ASSEMBLY}"
    args = args.format(**ph)
    prof = ",".join(profiles + [engine])
    return ["-profile", prof, "--input", ph["INPUT"], "--outdir", "out"] + shlex.split(args) + shlex.split(extra)


def dag_processes(nxf, pipeline, cmdline, wdir, env):
    dag = os.path.join(wdir, "dag.mmd")
    if os.path.exists(dag):
        os.remove(dag)
    run([nxf, "-q", "run", pipeline, "-preview", "-with-dag", dag] + cmdline, wdir, env)
    names = re.findall(r'\(\["([A-Za-z0-9_]+)"\]\)', open(dag).read())
    return sorted(set(names))


def inspect_containers(nxf, pipeline, cmdline, wdir, env):
    out = run([nxf, "-q", "inspect", pipeline, "-format", "json"] + cmdline, wdir, env)
    out = out[out.index("{"):]
    return {p["name"]: p["container"] for p in json.loads(out)["processes"]}


def engine_registry(nxf, pipeline, engine, wdir, env):
    out = run([nxf, "config", "-flat", "-profile", engine, pipeline], wdir, env)
    m = re.search(rf"^{engine}\.registry\s*=\s*'([^']*)'", out, re.M)
    return m.group(1) if m else ""


# ── Engine-specific normalisation (mirrors Nextflow's ContainerHandler) ───────
def is_absolute_docker_name(img):
    p = img.find("/")
    if p == -1:
        return False
    head = img[:p]
    return "." in head or ":" in head or head == "localhost"


def normalise(img, engine, registry):
    if engine in ("docker", "podman"):
        if registry and not is_absolute_docker_name(img):
            return registry.rstrip("/") + "/" + img
        return img
    # singularity / apptainer
    if img.startswith("/") or re.match(r"^[^/:. ]+://", img):
        return img
    name = img
    if registry and not is_absolute_docker_name(img):
        name = registry.rstrip("/") + "/" + img
    return "docker://" + name


def sif_name(url):
    """SingularityCache.simpleName"""
    p = url.find("://")
    name = url[p + 3:] if p != -1 else url
    ext = ".img"
    if ".sif:" in name:
        ext, name = ".sif", name.replace(".sif:", "-")
    elif name.endswith(".sif"):
        ext, name = ".sif", name[:-4]
    return name.replace(":", "-").replace("/", "-") + ext


def pull_lines(engine, images, cache_dir):
    lines = ["#!/usr/bin/env bash", "set -uo pipefail", f"# {engine} pulls generated by scripts/image_manifest.py", ""]
    if engine in ("docker", "podman"):
        lines += [f"{engine} pull {shlex.quote(i)} || echo 'FAILED: {i}' >&2" for i in images]
    else:
        cache = cache_dir or f"${{NXF_{engine.upper()}_CACHEDIR:-./{engine}_cache}}"
        lines += [f'CACHE="{cache}"', 'mkdir -p "$CACHE"']
        for i in images:
            f = sif_name(i)
            lines.append(f'[ -s "$CACHE/{f}" ] || {engine} pull --name "$CACHE/{f}" {shlex.quote(i)} '
                         f"|| echo 'FAILED: {i}' >&2")
        lines += ["", f'echo "Point Nextflow at it with: export NXF_{engine.upper()}_CACHEDIR=$CACHE"']
    return "\n".join(lines) + "\n"


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pipeline", default=os.path.dirname(here), help="TaxTriage repo (default: this script's repo)")
    ap.add_argument("--engine", default="docker", choices=["docker", "singularity", "apptainer", "podman"])
    ap.add_argument("--mode", action="append", choices=list(PRESETS), help="preset mode(s); repeatable")
    ap.add_argument("--custom", action="append", default=[], metavar="NAME=ARGS",
                    help='extra mode, e.g. "denovo=--db {DB} --use_denovo". Add -profile light by starting ARGS with "light:"')
    ap.add_argument("--nf-args", default="", help="extra Nextflow args appended to every mode")
    ap.add_argument("--input", help="use a real samplesheet instead of the dummy one")
    ap.add_argument("--offline", action="store_true", help="add --taxdump/--assembly (drops the download images)")
    ap.add_argument("--core-only", action="store_true", help="leave data-dependent processes out of the image lists")
    ap.add_argument("--outdir", default="image_manifest")
    ap.add_argument("--pull", action="store_true", help="run the generated pull script")
    ap.add_argument("--cache-dir", help="singularity/apptainer image dir for --pull (default $NXF_<ENGINE>_CACHEDIR)")
    ap.add_argument("--nextflow", default=os.environ.get("NXF_BIN", "nextflow"))
    ap.add_argument("--keep-tmp", action="store_true")
    a = ap.parse_args()

    pipeline = os.path.abspath(a.pipeline)
    if not os.path.exists(os.path.join(pipeline, "main.nf")):
        sys.exit(f"no main.nf in {pipeline}")
    nxf = shutil.which(a.nextflow) or a.nextflow

    modes = OrderedDict()
    for m in (a.mode or DEFAULT_MODES):
        modes[m] = dict(PRESETS[m])
    for c in a.custom:
        name, _, args = c.partition("=")
        profiles = []
        if args.startswith("light:"):
            profiles, args = ["light"], args[len("light:"):]
        modes[name.strip()] = {"profiles": profiles, "args": args.strip()}

    tmp = tempfile.mkdtemp(prefix="tt_imgman_")
    env = dict(os.environ, NXF_ANSI_LOG="false")
    try:
        ph = make_dummy_inputs(tmp)
        if a.input:
            ph["INPUT"] = os.path.abspath(a.input)

        amap = alias_map(pipeline)
        pfiles = process_files(pipeline)
        pull_sels = pull_docker_selectors(pipeline)
        registry = engine_registry(nxf, pipeline, a.engine, tmp, env)
        log(f"[{a.engine}] registry={registry or '(none)'}  pipeline={pipeline}")

        rows = []
        for mode, spec in modes.items():
            log(f"== mode {mode}")
            cmd = mode_cmdline(spec["args"], a.engine, spec["profiles"], ph, a.offline, a.nf_args)
            base_containers = inspect_containers(nxf, pipeline, cmd, tmp, env)
            procs = sorted(base_containers) if spec.get("inspect_only") else dag_processes(nxf, pipeline, cmd, tmp, env)
            for proc in procs:
                base = amap.get(proc, proc)
                raw = base_containers.get(proc) or base_containers.get(base)
                if raw is None:
                    continue  # a sub-workflow or operator node, not a process
                # Per-alias ext.singularity_pull_docker_container (inspect only sees base names)
                if a.engine == "singularity" and base in pfiles:
                    want, have = selector_hits(pull_sels, proc), selector_hits(pull_sels, base)
                    br = container_branches(pfiles[base]) if want != have else None
                    if br:
                        raw = br[1] if want else br[0]
                img = normalise(raw, a.engine, registry)
                cond = DATA_DEPENDENT.get(proc, "")
                tier = "data-dependent" if cond else "core"
                if spec.get("inspect_only"):
                    tier, cond = "defined", "any process in the pipeline"
                rows.append(OrderedDict(engine=a.engine, mode=mode, process=proc, module=base,
                                        tier=tier, condition=cond, container=img))

        os.makedirs(a.outdir, exist_ok=True)
        with open(os.path.join(a.outdir, "manifest.tsv"), "w") as fh:
            fh.write("\t".join(rows[0].keys()) + "\n")
            for r in rows:
                fh.write("\t".join(r.values()) + "\n")

        union = OrderedDict()
        print(f"\n{'mode':<14}{'processes':>10}{'core':>6}{'images':>8}   ({a.engine})")
        for mode, spec in modes.items():
            mrows = [r for r in rows if r["mode"] == mode]
            imgs = defaultdict(lambda: {"core": False, "procs": []})
            for r in mrows:
                if a.core_only and r["tier"] == "data-dependent":
                    continue
                imgs[r["container"]]["procs"].append(r["process"] + ("" if r["tier"] != "data-dependent" else "*"))
                imgs[r["container"]]["core"] |= r["tier"] != "data-dependent"
            cmd = mode_cmdline(spec["args"], a.engine, spec["profiles"],
                               {k: f"<{k.lower()}>" for k in ph}, a.offline, a.nf_args)
            with open(os.path.join(a.outdir, f"images_{mode}.txt"), "w") as fh:
                fh.write(f"# mode: {mode}  engine: {a.engine}\n# nextflow run main.nf {' '.join(cmd)}\n")
                fh.write("# '*' = data-dependent process (only runs for matching samplesheet rows)\n")
                for img, info in sorted(imgs.items(), key=lambda kv: (not kv[1]["core"], kv[0])):
                    fh.write(f"{img}\t# {'core' if info['core'] else 'data-dependent'}: {', '.join(sorted(info['procs']))}\n")
                    union.setdefault(img, set()).add(mode)
            ncore = sum(1 for r in mrows if r["tier"] in ("core", "defined"))
            print(f"{mode:<14}{len(mrows):>10}{ncore:>6}{len(imgs):>8}")

        with open(os.path.join(a.outdir, "images_union.txt"), "w") as fh:
            for img, ms in union.items():
                fh.write(f"{img}\t# {','.join(sorted(ms))}\n")
        pull_sh = os.path.join(a.outdir, f"pull_{a.engine}.sh")
        with open(pull_sh, "w") as fh:
            fh.write(pull_lines(a.engine, list(union), a.cache_dir))
        os.chmod(pull_sh, 0o755)
        print(f"\nwrote {a.outdir}/manifest.tsv, images_<mode>.txt, images_union.txt ({len(union)} images), {os.path.basename(pull_sh)}")

        if a.pull:
            subprocess.run(["bash", pull_sh], check=False)
    finally:
        if a.keep_tmp:
            log(f"temp dir kept: {tmp}")
        else:
            shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    main()
