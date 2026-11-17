#!/usr/bin/env python3
"""
image_manifest.py -- which modules run, and which container images they need, for each
TaxTriage run mode, resolved by Nextflow itself.

How it works
------------
For every mode it runs (nothing executes, no data or databases needed):

  1. `nextflow run <pipeline> -preview -with-dag dag.mmd <mode params>`
       -> the processes the workflow wires in for those params (every `if (params.x)`
          branch is honoured).
  2. `nextflow inspect <pipeline> -profile <engine>[,<profile>] -format json`
       -> the container each process resolves to for that engine.

Process aliases (`include { X as Y }`) are mapped back to their module, and per-alias
`ext.singularity_pull_docker_container` overrides in conf/*.config are applied, so a
Singularity manifest matches what a real run pulls.

Modes (default: all of them except `defined`)
---------------------------------------------
  full         pipeline defaults, Kraken2 (--db)
  lite-k2      -profile lite: MCE settings, Kraken2 standard picks the references
  lite-fasta   -profile lite, local FASTA only (--skip_kraken2 --reference_fasta)
  lite-both    -profile lite, Kraken2 standard + local FASTA
  mce          -profile mce
  all          every optional module switched on. Some modules are alternatives to each
               other (Kraken2 vs Centrifuge, minimap2 vs bowtie2 vs HISAT2, the novelty
               backends, ...), so `all` is the union of several preview runs (ALL_VARIANTS).
  defined      every process defined in the code, including ones nothing calls
               (inspect only; not a run mode)

Processes are marked per mode:
  on    runs for every sample in that mode
  on*   wired in, but only gets tasks for matching samplesheet rows
        (trim=TRUE, SRA accessions, BAM inputs, controls, simulated datasets, ...)
  off   not part of that mode

Usage
-----
  python3 scripts/image_manifest.py                                   # all modes, docker
  python3 scripts/image_manifest.py --engine singularity --mode lite-k2 --core-only \
      --pull --cache-dir /shared/nxf_singularity_cache
  python3 scripts/image_manifest.py --custom "denovo=lite:--use_denovo --annotate"
  python3 scripts/image_manifest.py --offline       # add --taxdump/--assembly: no download images

Outputs (--outdir, default ./image_manifest)
--------------------------------------------
  manifest.xlsx       tab 1 Images  : one row per image, one column per mode
                      tab 2 Modules : process x mode grid, green = on, red = off
                      tab 3 Modes   : params and counts per mode, legend
  manifest.tsv        the Images tab as TSV
  modules.tsv         the Modules tab as TSV
  images_<mode>.txt   unique images for one mode (data-dependent ones marked)
  pull_<engine>.sh    pull commands for the selected modes (run with --pull)

Requires python3 (standard library only) and nextflow on PATH (or --nextflow); docker /
podman / singularity / apptainer only for --pull. Nextflow needs network for its nf-schema
plugin on first use, as for a normal run.
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
import zipfile
from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor
from xml.sax.saxutils import escape

# ── Modes ─────────────────────────────────────────────────────────────────────
# Each mode = profiles + one or more param sets ("variants"); a mode with several
# variants is the union of their previews. "{X}" placeholders are dummy inputs.
ALL_VARIANTS = [
    # minimap2 path with nearly every optional module on
    {"db": "{DB}", "reference_fasta": "{FASTA}", "organisms": "562", "remove_taxids": "9606",
     "genome": "aedes", "enable_fastp": True, "enable_plots": True, "downsample": True,
     "subsample": 1000, "use_denovo": True, "use_diamond": True, "annotate": True,
     "pathogenicity": True, "novelty": "mmseqs2", "microbert": "{MODEL}", "get_variants": True,
     "reference_assembly": True, "metaphlan": "{MPA}", "fcs": True, "mask_bed": "{BED}",
     "mask_products": True, "generate_iss": True, "generate_nanosim": True,
     "nanosim_training": "{NANOSIM}", "sim_subsample": True, "sim_series_counts": "100,500",
     "fast": True, "enable_genbank": True, "index_hopping": "{HOP}", "cram_reference": "{CRAMREF}"},
    # alternatives: Centrifuge, bowtie2, Deacon (prebuilt human index), kaiju novelty,
    # ONT spike-in, the curated pathogen set, MEGAHIT for long reads
    {"db": "{DB}", "centrifuge": True, "use_bt2": True, "use_deacon": True, "genome": "human",
     "novelty": "kaiju", "spikein_sheet": "{SPIKE}", "nanosim_training": "{NANOSIM}",
     "get_pathogens": True, "use_denovo": True, "use_megahit_longreads": True},
    # alternatives: HISAT2, bracken novelty (own db download), background dilution series,
    # Kraken2 host filter, Kraken2 db download
    {"db": "test", "download_db": True, "use_hisat2": True, "reference_fasta": "{FASTA}",
     "novelty": "bracken", "novelty_db": "viral", "background_reads": "{READS}",
     "sim_series_counts": "100,500", "filter_kraken2": "viral"},
    # alternatives: Deacon index built from host accessions, Illumina spike-in
    {"db": "{DB}", "use_deacon": True, "genome": "aedes", "spikein_sheet": "{SPIKE}",
     "background_platform": "ILLUMINA"},
    # pre-aligned input via --bam (samplesheet generated) with consensus reference
    {"bam": "{BAM}", "skip_kraken2": True, "organisms": "562"},
]

PRESETS = OrderedDict([
    ("full",       {"profiles": [],       "variants": [{"db": "{DB}"}],
                    "about": "Pipeline defaults, Kraken2"}),
    ("lite-k2",    {"profiles": ["lite"], "variants": [{}],
                    "about": "-profile lite: MCE settings, Kraken2 standard picks references"}),
    ("lite-fasta", {"profiles": ["lite"], "variants": [{"skip_kraken2": True, "reference_fasta": "{FASTA}"}],
                    "about": "-profile lite, local FASTA only"}),
    ("lite-both",  {"profiles": ["lite"], "variants": [{"reference_fasta": "{FASTA}"}],
                    "about": "-profile lite, Kraken2 standard + local FASTA"}),
    ("mce",        {"profiles": ["mce"],  "variants": [{}],
                    "about": "-profile mce"}),
    ("all",        {"profiles": [],       "variants": ALL_VARIANTS,
                    "about": "Every optional module on (union of %d preview runs)" % len(ALL_VARIANTS)}),
    ("defined",    {"profiles": [],       "variants": [{"db": "{DB}"}], "inspect_only": True,
                    "about": "Every process defined in the code (inspect only)"}),
])
DEFAULT_MODES = ["full", "lite-k2", "lite-fasta", "lite-both", "mce", "all"]

# ── Processes that are wired in but only get tasks for certain samplesheet rows ─
# (the preview DAG is static, so this part is annotation)
DATA_DEPENDENT = {
    "GENERATE_SAMPLESHEET":        "--fastq_1 / --bam given instead of --input",
    "SRA_RESOLVE":                 "rows with SRA/ENA accessions",
    "SRA_FETCH_ENA":               "rows with SRA/ENA accessions",
    "SRA_FETCH_SRATOOLS":          "rows with SRA/ENA accessions",
    "ARTIC_GUPPYPLEX":             "rows whose fastq_1 is a directory of ONT FASTQs",
    "PIGZ_COMPRESS":               "uncompressed FASTQ inputs",
    "TRIMGALORE":                  "ILLUMINA rows with trim=TRUE (forced off by lite/light)",
    "PORECHOP":                    "OXFORD/PACBIO rows with trim=TRUE (forced off by lite/light)",
    "PYCOQC":                      "OXFORD rows with a sequencing_summary",
    "PREPARE_BAM":                 "BAM/CRAM input rows",
    "BAM_CONSENSUS":               "BAM/CRAM rows without a reference FASTA",
    "SAMTOOLS_COVERAGE_BAM":       "BAM/CRAM input rows",
    "SAMTOOLS_HIST_COVERAGE_BAM":  "BAM/CRAM input rows",
    "BEDTOOLS_GENOMECOVERAGE_BAM": "BAM/CRAM input rows with --fast",
    "DECOMPRESS_CRAM_REFERENCE":   "gzipped --cram_reference",
    "MAP_SAMPLE_FASTA_TO_ASSEMBLY": "rows with a per-sample `fasta` column",
    "BOWTIE2_BUILD_SAMPLE":        "rows with a per-sample `fasta` column (+ --use_bt2)",
    "BUILD_SAMPLE_MASK":           "rows with mask_fasta/mask_bed columns, or --mask_products",
    "SAMTOOLS_MERGE":              "a sample aligned to >1 reference file (Kraken2 refs + --reference_fasta, several FASTAs)",
    "ALIGNMENT_PER_SAMPLE_CONTROLS": "negative/positive control rows",
    "ALIGNMENT_PER_SAMPLE_INSILICO": "simulated datasets",
    "SINGLE_REPORT_INSILICO":      "simulated datasets",
    "INDEX_HOPPING":               "--index_hopping or an index_hopping samplesheet column",
}

# ── Stage labels for the Modules grid (subworkflow, then per-process overrides) ─
STAGES = ["Input", "Databases & downloads", "Host removal", "Read QC & reduction",
          "Classification", "Pre-aligned (BAM) input", "References", "Simulation",
          "Alignment", "Variants & consensus", "De novo & annotation", "MicrobeRT",
          "Novelty", "Masking", "TASS", "Report", "Aggregate outputs", "Other"]
SUBWF_STAGE = {"INPUT_CHECK": "Input", "HOST_REMOVAL": "Host removal", "CLASSIFIER": "Classification",
               "BAM_PREP": "Pre-aligned (BAM) input", "BAM_INPUT": "Pre-aligned (BAM) input",
               "REFERENCE_PREP": "References", "INSILICO": "Simulation", "BACKGROUND": "Simulation",
               "SPIKEIN": "Simulation", "ALIGNMENT": "Alignment", "ASSEMBLY": "De novo & annotation",
               "PATHOGENICITY": "De novo & annotation", "PROTEINS": "De novo & annotation",
               "NOVELTY": "Novelty", "REPORT": "Report"}
PROC_STAGE = {}
for _stage, _procs in {
    "Input": "ARTIC_GUPPYPLEX PIGZ_COMPRESS COUNT_READS GENERATE_SAMPLESHEET",
    "Databases & downloads": "DOWNLOAD_TAXDUMP GET_ASSEMBLIES DOWNLOAD_DB MMSEQS_DOWNLOADDB KAIJU_DOWNLOADDB "
                             "KRAKEN2_DOWNLOADDB DOWNLOAD_PATHOGENS DECOMPRESS_CRAM_REFERENCE FCS_DOWNLOAD",
    "Read QC & reduction": "BBMAP_BBNORM SEQTK_SAMPLE TRIMGALORE PORECHOP PYCOQC FASTP FASTPLONG FASTQC NANOPLOT",
    "MicrobeRT": "CLUSTER_ALIGNMENT MMSEQS_EASYCLUSTER MICROBERT_PREDICT MICROBERT_PARSE",
    "Variants & consensus": "BCFTOOLS_MPILEUP BCFTOOLS_CONSENSUS BEDTOOLS_BAMTOBED BEDTOOLS_MASKFASTA",
    "Masking": "BUILD_MASK BUILD_SAMPLE_MASK MASK_REPORT",
    "TASS": "ALIGNMENT_PER_SAMPLE ALIGNMENT_PER_SAMPLE_CONTROLS ALIGNMENT_PER_SAMPLE_INSILICO INDEX_HOPPING",
    "Aggregate outputs": "MULTIQC CUSTOM_DUMPSOFTWAREVERSIONS",
}.items():
    for _p in _procs.split():
        PROC_STAGE[_p] = _stage


def log(msg):
    print(msg, file=sys.stderr, flush=True)


def run(cmd, cwd, env):
    log("  $ " + " ".join(shlex.quote(c) for c in cmd))
    p = subprocess.run(cmd, cwd=cwd, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        msg = "\n".join(l for l in (p.stdout + p.stderr).splitlines()[-25:] if "proxy" not in l)
        raise RuntimeError(f"command failed (rc={p.returncode}) in {cwd}:\n{msg}\nsee {cwd}/.nextflow.log")
    return p.stdout


# ── Dummy inputs so the workflow can be parsed without real data/databases ─────
def make_dummy_inputs(tmp):
    d = os.path.join(tmp, "inputs")
    os.makedirs(d, exist_ok=True)

    def touch(*parts):
        path = os.path.join(d, *parts)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        open(path, "a").close()
        return path

    fq = "@r1\nACGTACGTACGT\n+\nIIIIIIIIIIII\n"
    for n in ("r1.fastq.gz", "r2.fastq.gz", "ont.fastq.gz"):
        with gzip.open(os.path.join(d, n), "wt") as fh:
            fh.write(fq)
    ss = os.path.join(d, "samplesheet.csv")
    with open(ss, "w") as fh:
        fh.write("sample,platform,fastq_1,fastq_2,trim\n")
        fh.write(f"S1,ILLUMINA,{d}/r1.fastq.gz,{d}/r2.fastq.gz,FALSE\n")
        fh.write(f"S2,OXFORD,{d}/ont.fastq.gz,,FALSE\n")
    fa = os.path.join(d, "ref.fasta")
    with open(fa, "w") as fh:
        fh.write(">NC_000000.1 dummy\nACGTACGTACGTACGT\n")
    with gzip.open(os.path.join(d, "cram_ref.fa.gz"), "wt") as fh:
        fh.write(">x\nACGT\n")
    with open(os.path.join(d, "mask.bed"), "w") as fh:
        fh.write("NC_000000.1\t0\t4\n")
    with open(os.path.join(d, "spike.csv"), "w") as fh:
        fh.write("accession,count,replicates\nGCF_000005845.2,100,1\n")
    for n in ("hash.k2d", "opts.k2d", "taxo.k2d"):
        touch("k2db", n)
    for n in ("nodes.dmp", "names.dmp"):
        touch("taxdump", n)
    os.makedirs(os.path.join(d, "nanosim_model"), exist_ok=True)
    return {"INPUT": ss, "DB": os.path.join(d, "k2db"), "FASTA": fa,
            "TAXDUMP": os.path.join(d, "taxdump"), "ASSEMBLY": touch("assembly_summary_refseq.txt"),
            "MODEL": touch("microbert_model.pt"), "MPA": os.path.dirname(touch("mpa", "mpa_vJan21.pkl")),
            "BED": os.path.join(d, "mask.bed"), "HOP": touch("index_hopping.tsv"),
            "CRAMREF": os.path.join(d, "cram_ref.fa.gz"), "SPIKE": os.path.join(d, "spike.csv"),
            "NANOSIM": os.path.join(d, "nanosim_model"), "READS": os.path.join(d, "r1.fastq.gz"),
            "BAM": touch("sample.bam")}


def fill(v, ph):
    return v.format(**ph) if isinstance(v, str) else v


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
    for f in nf_files(pipeline):
        if os.path.exists(f):
            for m in re.finditer(r"include\s*\{([^}]*)\}\s*from", open(f).read(), re.S):
                for part in m.group(1).split(";"):
                    mm = re.match(r"^\s*([A-Za-z0-9_]+)(?:\s+as\s+([A-Za-z0-9_]+))?\s*$", part)
                    if mm and mm.group(2):
                        amap[mm.group(2)] = mm.group(1)
    return amap


def process_files(pipeline):
    out = {}
    for f in nf_files(pipeline):
        if os.path.exists(f):
            for m in re.finditer(r"^\s*process\s+([A-Za-z0-9_]+)\s*\{", open(f).read(), re.M):
                out[m.group(1)] = f
    return out


def container_branches(module_file):
    """('singularity-branch', 'docker-branch') of the nf-core style ternary, or None."""
    src = re.sub(r"/\*.*?\*/", "", open(module_file).read(), flags=re.S)
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
def parse_dag(path):
    """process -> (innermost subworkflow, vertex index) from a mermaid DAG."""
    out, stack = {}, []
    for line in open(path):
        line = line.strip()
        m = re.match(r'subgraph "([^"]*)"', line)
        if m:
            inner = re.search(r"\[([A-Za-z0-9_]+)\]$", m.group(1))
            stack.append(inner.group(1) if inner else "")
            continue
        if line == "end":
            if stack:
                stack.pop()
            continue
        m = re.match(r'v(\d+)\(\["([A-Za-z0-9_]+)"\]\)', line)
        if m:
            sub = next((s for s in reversed(stack) if s and s not in ("NFCORE_TAXTRIAGE", "TAXTRIAGE")), "")
            out[m.group(2)] = (sub, int(m.group(1)))
    return out


def preview(nxf, pipeline, cmd, wdir, env):
    os.makedirs(wdir, exist_ok=True)
    dag = os.path.join(wdir, "dag.mmd")
    run([nxf, "-q", "run", pipeline, "-preview", "-with-dag", dag] + cmd, wdir, env)
    return parse_dag(dag)


def inspect(nxf, pipeline, cmd, wdir, env):
    os.makedirs(wdir, exist_ok=True)
    out = run([nxf, "-q", "inspect", pipeline, "-format", "json"] + cmd, wdir, env)
    return {p["name"]: p["container"] for p in json.loads(out[out.index("{"):])["processes"]}


def engine_registry(nxf, pipeline, engine, wdir, env):
    os.makedirs(wdir, exist_ok=True)
    out = run([nxf, "config", "-flat", "-profile", engine, pipeline], wdir, env)
    m = re.search(rf"^{engine}\.registry\s*=\s*'([^']*)'", out, re.M)
    return m.group(1) if m else ""


# ── Engine-specific naming (mirrors Nextflow's ContainerHandler / SingularityCache) ─
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
    if img.startswith("/") or re.match(r"^[^/:. ]+://", img):
        return img
    if registry and not is_absolute_docker_name(img):
        img = registry.rstrip("/") + "/" + img
    return "docker://" + img


def sif_name(url):
    p = url.find("://")
    name = url[p + 3:] if p != -1 else url
    ext = ".img"
    if ".sif:" in name:
        ext, name = ".sif", name.replace(".sif:", "-")
    elif name.endswith(".sif"):
        ext, name = ".sif", name[:-4]
    return name.replace(":", "-").replace("/", "-") + ext


def pull_script(engine, images, cache_dir):
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


# ── Minimal XLSX writer (standard library only) ───────────────────────────────
# styles: 0 plain, 1 header, 2 on, 3 off, 4 on*, 5 text, 6 bold text, 7 note, 8 wrap text
_STYLES = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
<fonts count="5">
<font><sz val="10"/><name val="Arial"/></font>
<font><b/><sz val="10"/><name val="Arial"/></font>
<font><sz val="10"/><color rgb="FF006100"/><name val="Arial"/></font>
<font><sz val="10"/><color rgb="FF9C0006"/><name val="Arial"/></font>
<font><i/><sz val="10"/><color rgb="FF595959"/><name val="Arial"/></font>
</fonts>
<fills count="6">
<fill><patternFill patternType="none"/></fill>
<fill><patternFill patternType="gray125"/></fill>
<fill><patternFill patternType="solid"><fgColor rgb="FFC6EFCE"/><bgColor indexed="64"/></patternFill></fill>
<fill><patternFill patternType="solid"><fgColor rgb="FFFFC7CE"/><bgColor indexed="64"/></patternFill></fill>
<fill><patternFill patternType="solid"><fgColor rgb="FFD9E1F2"/><bgColor indexed="64"/></patternFill></fill>
<fill><patternFill patternType="solid"><fgColor rgb="FFE8F5E3"/><bgColor indexed="64"/></patternFill></fill>
</fills>
<borders count="2">
<border><left/><right/><top/><bottom/><diagonal/></border>
<border><left style="thin"><color rgb="FFD9D9D9"/></left><right style="thin"><color rgb="FFD9D9D9"/></right><top style="thin"><color rgb="FFD9D9D9"/></top><bottom style="thin"><color rgb="FFD9D9D9"/></bottom><diagonal/></border>
</borders>
<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>
<cellXfs count="9">
<xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>
<xf numFmtId="0" fontId="1" fillId="4" borderId="1" xfId="0" applyFont="1" applyFill="1" applyBorder="1" applyAlignment="1"><alignment horizontal="center" vertical="center" wrapText="1"/></xf>
<xf numFmtId="0" fontId="2" fillId="2" borderId="1" xfId="0" applyFont="1" applyFill="1" applyBorder="1" applyAlignment="1"><alignment horizontal="center"/></xf>
<xf numFmtId="0" fontId="3" fillId="3" borderId="1" xfId="0" applyFont="1" applyFill="1" applyBorder="1" applyAlignment="1"><alignment horizontal="center"/></xf>
<xf numFmtId="0" fontId="2" fillId="5" borderId="1" xfId="0" applyFont="1" applyFill="1" applyBorder="1" applyAlignment="1"><alignment horizontal="center"/></xf>
<xf numFmtId="0" fontId="0" fillId="0" borderId="1" xfId="0" applyBorder="1"/>
<xf numFmtId="0" fontId="1" fillId="0" borderId="1" xfId="0" applyFont="1" applyBorder="1"/>
<xf numFmtId="0" fontId="4" fillId="0" borderId="0" xfId="0" applyFont="1"/>
<xf numFmtId="0" fontId="0" fillId="0" borderId="1" xfId="0" applyBorder="1" applyAlignment="1"><alignment wrapText="1" vertical="top"/></xf>
</cellXfs>
<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>
</styleSheet>"""
ST = {"plain": 0, "header": 1, "on": 2, "off": 3, "on*": 4, "text": 5, "bold": 6, "note": 7, "wrap": 8}


def _col(n):
    s, n = "", n + 1
    while n:
        n, r = divmod(n - 1, 26)
        s = chr(65 + r) + s
    return s


def _sheet_xml(sheet):
    rows, widths = sheet["rows"], sheet["widths"]
    fx, fy = sheet.get("freeze", (0, 0))
    out = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
           '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
           'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">',
           '<sheetViews><sheetView workbookViewId="0"' + (' tabSelected="1"' if sheet.get("selected") else "") + ">"]
    if fx or fy:
        attrs = (f' xSplit="{fx}"' if fx else "") + (f' ySplit="{fy}"' if fy else "")
        pane = "bottomRight" if fx and fy else ("topRight" if fx else "bottomLeft")
        out.append(f'<pane{attrs} topLeftCell="{_col(fx)}{fy + 1}" activePane="{pane}" state="frozen"/>')
    out.append('</sheetView></sheetViews><sheetFormatPr defaultRowHeight="15"/><cols>')
    for i, w in enumerate(widths):
        out.append(f'<col min="{i + 1}" max="{i + 1}" width="{w}" customWidth="1"/>')
    out.append("</cols><sheetData>")
    for r, row in enumerate(rows):
        out.append(f'<row r="{r + 1}"' + (' ht="30" customHeight="1"' if r == 0 else "") + ">")
        for c, cell in enumerate(row):
            if cell is None:
                continue
            val, style = cell if isinstance(cell, tuple) else (cell, "text")
            ref, s = f"{_col(c)}{r + 1}", ST[style]
            if isinstance(val, (int, float)) and not isinstance(val, bool):
                out.append(f'<c r="{ref}" s="{s}"><v>{val}</v></c>')
            else:
                out.append(f'<c r="{ref}" s="{s}" t="inlineStr"><is><t xml:space="preserve">{escape(str(val))}</t></is></c>')
        out.append("</row>")
    out.append("</sheetData>")
    if sheet.get("autofilter"):
        out.append(f'<autoFilter ref="A1:{_col(sheet["autofilter"] - 1)}{sheet.get("table_rows", len(rows))}"/>')
    out.append("</worksheet>")
    return "".join(out)


def write_xlsx(path, sheets):
    n = len(sheets)
    ns_pkg = "http://schemas.openxmlformats.org/package/2006"
    ns_rel = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
    ct_sheet = "application/vnd.openxmlformats-officedocument.spreadsheetml"
    ct = (f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Types xmlns="{ns_pkg}/content-types">'
          '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
          '<Default Extension="xml" ContentType="application/xml"/>'
          f'<Override PartName="/xl/workbook.xml" ContentType="{ct_sheet}.sheet.main+xml"/>'
          f'<Override PartName="/xl/styles.xml" ContentType="{ct_sheet}.styles+xml"/>'
          + "".join(f'<Override PartName="/xl/worksheets/sheet{i + 1}.xml" ContentType="{ct_sheet}.worksheet+xml"/>'
                    for i in range(n)) + "</Types>")
    rels = (f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Relationships xmlns="{ns_pkg}/relationships">'
            f'<Relationship Id="rId1" Type="{ns_rel}/officeDocument" Target="xl/workbook.xml"/></Relationships>')
    defined = "".join(
        f'<definedName name="_xlnm._FilterDatabase" localSheetId="{i}" hidden="1">'
        f"'{escape(s['name'])}'!$A$1:${_col(s['autofilter'] - 1)}${s.get('table_rows', len(s['rows']))}</definedName>"
        for i, s in enumerate(sheets) if s.get("autofilter"))
    wb = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
          f'<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="{ns_rel}">'
          '<bookViews><workbookView activeTab="0"/></bookViews><sheets>'
          + "".join(f'<sheet name="{escape(s["name"])}" sheetId="{i + 1}" r:id="rId{i + 1}"/>' for i, s in enumerate(sheets))
          + "</sheets>" + (f"<definedNames>{defined}</definedNames>" if defined else "") + "</workbook>")
    wbrels = (f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Relationships xmlns="{ns_pkg}/relationships">'
              + "".join(f'<Relationship Id="rId{i + 1}" Type="{ns_rel}/worksheet" Target="worksheets/sheet{i + 1}.xml"/>'
                        for i in range(n))
              + f'<Relationship Id="rId{n + 1}" Type="{ns_rel}/styles" Target="styles.xml"/></Relationships>')
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml", ct)
        z.writestr("_rels/.rels", rels)
        z.writestr("xl/workbook.xml", wb)
        z.writestr("xl/_rels/workbook.xml.rels", wbrels)
        z.writestr("xl/styles.xml", _STYLES)
        for i, s in enumerate(sheets):
            z.writestr(f"xl/worksheets/sheet{i + 1}.xml", _sheet_xml(s))


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pipeline", default=os.path.dirname(here), help="TaxTriage repo (default: this script's repo)")
    ap.add_argument("--engine", default="docker", choices=["docker", "singularity", "apptainer", "podman"])
    ap.add_argument("--mode", action="append", choices=list(PRESETS), help="mode(s) to include; repeatable")
    ap.add_argument("--custom", action="append", default=[], metavar="NAME=[PROFILES:]ARGS",
                    help='extra mode, e.g. "denovo=lite:--use_denovo --annotate" ({DB} {FASTA} ... placeholders allowed)')
    ap.add_argument("--nf-args", default="", help="extra Nextflow args appended to every run")
    ap.add_argument("--input", help="use a real samplesheet instead of the dummy one")
    ap.add_argument("--offline", action="store_true", help="add --taxdump/--assembly (drops the download images)")
    ap.add_argument("--core-only", action="store_true", help="leave data-dependent processes out of image lists / pulls")
    ap.add_argument("--outdir", default="image_manifest")
    ap.add_argument("--pull", action="store_true", help="run the generated pull script")
    ap.add_argument("--cache-dir", help="singularity/apptainer image dir for --pull (default $NXF_<ENGINE>_CACHEDIR)")
    ap.add_argument("--nextflow", default=os.environ.get("NXF_BIN", "nextflow"))
    ap.add_argument("--jobs", type=int, default=4, help="parallel Nextflow runs (default 4)")
    ap.add_argument("--keep-tmp", action="store_true")
    a = ap.parse_args()

    pipeline = os.path.abspath(a.pipeline)
    if not os.path.exists(os.path.join(pipeline, "main.nf")):
        sys.exit(f"no main.nf in {pipeline}")
    nxf = shutil.which(a.nextflow) or a.nextflow

    modes = OrderedDict((m, dict(PRESETS[m])) for m in (a.mode or DEFAULT_MODES))
    for c in a.custom:
        name, _, args = c.partition("=")
        profiles = []
        m = re.match(r"^([a-z0-9_,]+):(.*)$", args.strip(), re.S)
        if m:
            profiles, args = m.group(1).split(","), m.group(2)
        modes[name.strip()] = {"profiles": profiles, "variants": [args.strip()], "about": "custom: " + c}

    tmp = tempfile.mkdtemp(prefix="tt_imgman_")
    env = dict(os.environ, NXF_ANSI_LOG="false")
    try:
        ph = make_dummy_inputs(tmp)
        if a.input:
            ph["INPUT"] = os.path.abspath(a.input)

        def shown(v):
            return f"<{os.path.basename(v)}>" if isinstance(v, str) and v.startswith(tmp) else v

        # One Nextflow command line per (mode, variant)
        jobs = []      # (mode, idx, cmd, display)
        for mode, spec in modes.items():
            prof = ",".join(spec["profiles"] + [a.engine])
            disp_prof = ",".join(spec["profiles"] + ["<engine>"])
            for i, var in enumerate(spec["variants"]):
                base = ["-profile", prof, "--input", ph["INPUT"], "--outdir", "out"]
                if isinstance(var, dict):
                    params = {k: fill(v, ph) for k, v in var.items()}
                    if a.offline:
                        params.update(taxdump=ph["TAXDUMP"], assembly=ph["ASSEMBLY"])
                    pf = os.path.join(tmp, f"params_{mode}_{i}.json")
                    with open(pf, "w") as fh:
                        json.dump(params, fh, indent=1)
                    cmd = base + ["-params-file", pf]
                    disp = f"-profile {disp_prof} " + " ".join(
                        f"--{k}" if v is True else f"--{k} {shown(v)}" for k, v in params.items())
                else:
                    extra = fill(var, ph) + (f" --taxdump {ph['TAXDUMP']} --assembly {ph['ASSEMBLY']}" if a.offline else "")
                    cmd = base + shlex.split(extra)
                    disp = f"-profile {disp_prof} {var}"
                jobs.append((mode, i, cmd + shlex.split(a.nf_args), disp.strip()))

        amap = alias_map(pipeline)
        pfiles = process_files(pipeline)
        pull_sels = pull_docker_selectors(pipeline)
        registry = engine_registry(nxf, pipeline, a.engine, os.path.join(tmp, "cfg"), env)
        log(f"[{a.engine}] registry={registry or '(none)'}  pipeline={pipeline}  preview runs={len(jobs)}")

        # inspect once per profile set; preview every (mode, variant); in parallel
        insp_cmds = OrderedDict()
        for mode, i, cmd, _ in jobs:
            insp_cmds.setdefault(tuple(modes[mode]["profiles"]), cmd)

        def do_inspect(item):
            key, cmd = item
            return key, inspect(nxf, pipeline, cmd, os.path.join(tmp, "insp_" + ("_".join(key) or "none")), env)

        def do_preview(job):
            mode, i, cmd, _ = job
            if modes[mode].get("inspect_only"):
                return (mode, i), None
            return (mode, i), preview(nxf, pipeline, cmd, os.path.join(tmp, f"run_{mode}_{i}"), env)

        with ThreadPoolExecutor(max_workers=max(1, a.jobs)) as ex:
            containers = dict(ex.map(do_inspect, insp_cmds.items()))
            dags = dict(ex.map(do_preview, jobs))
        default_cmap = containers[next(iter(containers))]

        # status: mode -> process -> 'on' | 'on*' ; where: process -> (stage, DAG order)
        status = OrderedDict((m, {}) for m in modes)
        where = {}
        for (mode, i), dag in dags.items():
            cmap = containers[tuple(modes[mode]["profiles"])]
            procs = {p: ("", 10 ** 6) for p in cmap} if dag is None else dag
            for proc, (sub, vid) in procs.items():
                base = amap.get(proc, proc)
                if proc not in cmap and base not in cmap:
                    continue      # not a process (e.g. a sub-workflow node)
                if status[mode].get(proc) != "on":
                    status[mode][proc] = "on*" if proc in DATA_DEPENDENT else "on"
                stage = PROC_STAGE.get(proc) or PROC_STAGE.get(base) or SUBWF_STAGE.get(sub) or "Other"
                old = where.get(proc)
                where[proc] = (stage, vid) if old is None else \
                    (old[0] if old[0] != "Other" else stage, min(old[1], vid))

        def resolve(proc):
            base = amap.get(proc, proc)
            raw = default_cmap.get(proc) or default_cmap.get(base)
            if a.engine == "singularity" and base in pfiles:
                want, have = selector_hits(pull_sels, proc), selector_hits(pull_sels, base)
                br = container_branches(pfiles[base]) if want != have else None
                if br:
                    raw = br[1] if want else br[0]
            return normalise(raw, a.engine, registry)

        all_procs = sorted(where, key=lambda p: (STAGES.index(where[p][0]), where[p][1], p))
        img_of = {p: resolve(p) for p in all_procs}
        mode_names = list(modes)

        # images: image -> processes + {mode: 'core' | 'data'}
        images = OrderedDict()
        for p in all_procs:
            info = images.setdefault(img_of[p], {"procs": [], "modes": {}})
            info["procs"].append(p)
            for m in mode_names:
                t = status[m].get(p)
                if t and info["modes"].get(m) != "core":
                    info["modes"][m] = "core" if t == "on" else "data"
        img_order = sorted(images, key=lambda i: (-sum(v == "core" for v in images[i]["modes"].values()),
                                                  -len(images[i]["modes"]), i))

        def proc_list(img, mode=None):
            return ", ".join(p + ("*" if p in DATA_DEPENDENT else "")
                             for p in images[img]["procs"] if mode is None or p in status[mode])

        os.makedirs(a.outdir, exist_ok=True)

        # ── TSVs ──
        with open(os.path.join(a.outdir, "manifest.tsv"), "w") as fh:
            fh.write("\t".join(["image", "processes"] + mode_names) + "\n")
            for img in img_order:
                fh.write("\t".join([img, proc_list(img)] + [images[img]["modes"].get(m, "") for m in mode_names]) + "\n")
        with open(os.path.join(a.outdir, "modules.tsv"), "w") as fh:
            fh.write("\t".join(["stage", "process", "module"] + mode_names + ["runs_only_when", "container"]) + "\n")
            for p in all_procs:
                fh.write("\t".join([where[p][0], p, amap.get(p, p)] + [status[m].get(p, "off") for m in mode_names]
                                   + [DATA_DEPENDENT.get(p, ""), img_of[p]]) + "\n")

        # ── per-mode image lists + pull script ──
        union, summary = OrderedDict(), []
        for m in mode_names:
            keep = [i for i in img_order if images[i]["modes"].get(m) == "core"
                    or (images[i]["modes"].get(m) == "data" and not a.core_only)]
            with open(os.path.join(a.outdir, f"images_{m}.txt"), "w") as fh:
                fh.write(f"# mode: {m}  engine: {a.engine}  ({modes[m]['about']})\n")
                fh.write("# '*' = data-dependent process (only runs for matching samplesheet rows)\n")
                for img in keep:
                    fh.write(f"{img}\t# {images[img]['modes'][m]}: {proc_list(img, m)}\n")
                    union.setdefault(img, None)
            summary.append((m,
                            sum(1 for t in status[m].values() if t == "on"),
                            sum(1 for t in status[m].values() if t == "on*"),
                            sum(1 for i in img_order if images[i]["modes"].get(m) == "core"),
                            sum(1 for i in img_order if m in images[i]["modes"])))
        pull_sh = os.path.join(a.outdir, f"pull_{a.engine}.sh")
        with open(pull_sh, "w") as fh:
            fh.write(pull_script(a.engine, list(union), a.cache_dir))
        os.chmod(pull_sh, 0o755)

        # ── XLSX: Images | Modules | Modes ──
        nm = len(mode_names)
        img_cell = {"core": ("core", "on"), "data": ("data*", "on*")}
        img_rows = [[("Image", "header"), ("Processes (* = data-dependent)", "header")]
                    + [(m, "header") for m in mode_names]]
        for img in img_order:
            img_rows.append([(img, "text"), (proc_list(img), "wrap")]
                            + [img_cell.get(images[img]["modes"].get(m), ("–", "off")) for m in mode_names])
        n_img_table = len(img_rows)
        img_rows += [[], [("core = needed for every sample in that mode;  data* = only needed for matching "
                           "samplesheet rows (see the Modules tab);  – = not used.", "note")]]

        mod_rows = [[("Stage", "header"), ("Process", "header"), ("Module", "header")]
                    + [(m, "header") for m in mode_names] + [("Runs only when", "header"), ("Container", "header")]]
        for p in all_procs:
            base = amap.get(p, p)
            mod_rows.append([(where[p][0], "text"), (p, "bold"), (base if base != p else "", "text")]
                            + [(status[m].get(p, "off"), status[m].get(p, "off")) for m in mode_names]
                            + [(DATA_DEPENDENT.get(p, ""), "text"), (img_of[p], "text")])
        n_mod_table = len(mod_rows)
        mod_rows += [[], [("on = runs for every sample in that mode;  on* = wired in, but only gets tasks for "
                           "matching samplesheet rows (see 'Runs only when');  off = not part of that mode.", "note")]]

        mode_rows = [[("Mode", "header"), ("What it is", "header"), ("Modules on", "header"),
                      ("Modules on*", "header"), ("Core images", "header"), ("All images", "header"),
                      ("Nextflow options (one line per preview run)", "header")]]
        for (m, n_on, n_cond, n_core_img, n_img) in summary:
            cmds = [d for (mm, _, _, d) in jobs if mm == m]
            mode_rows.append([(m, "bold"), (modes[m]["about"], "wrap"), (n_on, "text"), (n_cond, "text"),
                              (n_core_img, "text"), (n_img, "text"), ("\n".join(cmds), "wrap")])
        mode_rows += [[], [(f"Engine: {a.engine}. Generated by scripts/image_manifest.py from `nextflow run -preview` "
                            "+ `nextflow inspect`. <name> values are placeholder inputs; every run also gets "
                            "--input <samplesheet.csv> --outdir out.", "note")]]
        if "all" in modes:
            reached = {amap.get(p, p) for p in status["all"]}
            dead = sorted(b for b in default_cmap if b not in reached)
            if dead:
                mode_rows.append([(f"Defined in the code but not called in any mode ({len(dead)}): "
                                   + ", ".join(dead), "note")])

        write_xlsx(os.path.join(a.outdir, "manifest.xlsx"), [
            {"name": "Images", "rows": img_rows, "widths": [72, 60] + [11] * nm, "freeze": (1, 1),
             "autofilter": 2 + nm, "table_rows": n_img_table, "selected": True},
            {"name": "Modules", "rows": mod_rows, "widths": [22, 38, 30] + [11] * nm + [52, 72],
             "freeze": (2, 1), "autofilter": 5 + nm, "table_rows": n_mod_table},
            {"name": "Modes", "rows": mode_rows, "widths": [12, 46, 11, 11, 11, 11, 130], "freeze": (1, 1)},
        ])

        print(f"\n{'mode':<12}{'on':>5}{'on*':>6}{'core imgs':>11}{'all imgs':>10}   ({a.engine})")
        for (m, n_on, n_cond, n_core_img, n_img) in summary:
            print(f"{m:<12}{n_on:>5}{n_cond:>6}{n_core_img:>11}{n_img:>10}")
        print(f"\nwrote {a.outdir}/manifest.xlsx (Images, Modules, Modes), manifest.tsv, modules.tsv, "
              f"images_<mode>.txt, {os.path.basename(pull_sh)} ({len(union)} images)")
        if a.pull:
            subprocess.run(["bash", pull_sh], check=False)
    except RuntimeError as e:
        sys.exit(str(e))
    finally:
        if a.keep_tmp:
            log(f"temp dir kept: {tmp}")
        else:
            shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    main()
