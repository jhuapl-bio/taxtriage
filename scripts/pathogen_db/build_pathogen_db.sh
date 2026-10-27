#!/usr/bin/env bash
##############################################################################################
# Copyright 2026 The Johns Hopkins University Applied Physics Laboratory LLC
# (MIT-style license, see LICENSE)
##############################################################################################
# Build the TaxTriage pathogen reference from the pathogen sheet:
#
#   <prefix>.fasta(.gz)              one assembly per sheet taxid   -> --reference_fasta
#   <prefix>.mapping.taxid.tsv       Acc, Assembly, Organism_Name, Description, Mapped_Value
#   <prefix>.gcfmapping.tsv          Acc, Assembly, Organism_Name, Description (download map)
#   <prefix>.assemblies.txt          unique assembly accessions
#   <prefix>.mapping.taxid.validated.tsv   sheet name vs. downloaded organism name
#   pathogen_sheet.with_assembly.csv sheet + assembly_accession column
#   <prefix>.download_input.tsv      taxid<TAB>assembly_accession, reused by the next run
#   <prefix>.missing.txt             sheet entries no assembly could be found for
#
# This is the standalone form of the steps used to produce the pathogens.fasta.gz
# that --get_pathogens downloads. It calls bin/download_fastas.py, bin/append_taxid.py
# and the helpers next to this script, so run it from a TaxTriage checkout.
#
# Requirements: python3 with biopython + pandas; optional samtools (index), bgzip or pigz.
#
# Example:
#   scripts/pathogen_db/build_pathogen_db.sh -o pathogen_build/ -e you@org.org -t 8
#   scripts/pathogen_db/build_pathogen_db.sh -o pathogen_build/ -a assembly_summary_refseq.txt --skip-download
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
BIN="$REPO/bin"

SHEET="$REPO/assets/pathogen_sheet.csv"
OUTDIR="."
ASSEMBLY=""
EMAIL=""
THREADS=4
PREFIX="pathogens"
SKIP_DOWNLOAD=0
NO_COMPRESS=0
FRESH_INPUT=0

usage() {
    awk 'NR>5 && /^#/ {sub(/^# ?/, ""); print; next} NR>5 {exit}' "$0"
    cat <<EOF

Options:
  -s, --sheet FILE        pathogen sheet CSV with name + taxid columns
                          (default: $REPO/assets/pathogen_sheet.csv)
  -o, --outdir DIR        output directory (default: .)
  -a, --assembly FILE     assembly_summary_refseq.txt (default: <outdir>/assembly_summary_refseq.txt,
                          downloaded from NCBI if missing)
  -e, --email EMAIL       email for NCBI E-utilities (NCBI fallback + taxid backup)
  -t, --threads N         compression threads (default: 4)
  -p, --prefix NAME       output prefix (default: pathogens)
      --fresh-input       rebuild <prefix>.download_input.tsv from the sheet instead of reusing
                          the one written by the previous run
      --skip-download     reuse an existing <prefix>.fasta; only rebuild the mappings/sheet
      --no-compress       do not write <prefix>.fasta.gz
  -h, --help
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -s|--sheet)      SHEET=$2; shift 2 ;;
        -o|--outdir)     OUTDIR=$2; shift 2 ;;
        -a|--assembly)   ASSEMBLY=$2; shift 2 ;;
        -e|--email)      EMAIL=$2; shift 2 ;;
        -t|--threads)    THREADS=$2; shift 2 ;;
        -p|--prefix)     PREFIX=$2; shift 2 ;;
        --fresh-input)   FRESH_INPUT=1; shift ;;
        --skip-download) SKIP_DOWNLOAD=1; shift ;;
        --no-compress)   NO_COMPRESS=1; shift ;;
        -h|--help)       usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
done

[[ -f "$SHEET" ]] || { echo "ERROR: sheet not found: $SHEET" >&2; exit 1; }
mkdir -p "$OUTDIR"
SHEET=$(cd "$(dirname "$SHEET")" && pwd)/$(basename "$SHEET")
[[ -n "$ASSEMBLY" ]] && ASSEMBLY=$(cd "$(dirname "$ASSEMBLY")" && pwd)/$(basename "$ASSEMBLY")
cd "$OUTDIR"
LOG="$PREFIX.build.log"
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') $*" | tee -a "$LOG" >&2; }
run() { log "+ $*"; "$@" 2>&1 | tee -a "$LOG"; }

log "=== build_pathogen_db.sh  sheet=$SHEET  outdir=$(pwd)  prefix=$PREFIX"
EMAIL_ARGS=(); [[ -n "$EMAIL" ]] && EMAIL_ARGS=(-e "$EMAIL")

# ── 1. assembly summary ─────────────────────────────────────────────────────
if [[ -z "$ASSEMBLY" ]]; then
    ASSEMBLY="$(pwd)/assembly_summary_refseq.txt"
    if [[ ! -s "$ASSEMBLY" ]]; then
        URL=https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt
        log "Downloading $URL"
        if command -v curl >/dev/null; then curl -fsSL -o "$ASSEMBLY.part" "$URL"
        else wget -q -O "$ASSEMBLY.part" "$URL"; fi
        mv "$ASSEMBLY.part" "$ASSEMBLY"
    fi
fi
[[ -s "$ASSEMBLY" ]] || { log "ERROR: assembly summary not found: $ASSEMBLY"; exit 1; }

# ── 2. taxids + download input (taxid<TAB>assembly_accession) ───────────────
# Parsed with the csv module: names and references contain quoted commas.
python3 - "$SHEET" "$PREFIX" "$FRESH_INPUT" <<'PY'
import csv, os, sys
sheet, prefix, fresh = sys.argv[1], sys.argv[2], sys.argv[3] == "1"
with open(sheet, newline="", encoding="utf-8-sig") as fh:
    rows = list(csv.DictReader(fh))
taxids = sorted({r["taxid"].strip() for r in rows if (r.get("taxid") or "").strip().isdigit()}, key=int)
with open(f"{prefix}.taxids.txt", "w") as out:
    out.write("\n".join(taxids) + "\n")
di = f"{prefix}.download_input.tsv"
if fresh or not os.path.exists(di):
    seen = {}
    for r in rows:
        t = (r.get("taxid") or "").strip()
        if t.isdigit() and t not in seen:
            seen[t] = (r.get("assembly_accession") or "").strip()
    with open(di, "w") as out:
        for t, a in sorted(seen.items(), key=lambda x: int(x[0])):
            out.write(f"{t}\t{a}\n")
    print(f"Wrote {di} from the sheet ({len(seen)} taxids, "
          f"{sum(1 for a in seen.values() if a)} with an assembly accession)")
else:
    print(f"Reusing {di} from the previous run (use --fresh-input to rebuild it)")
print(f"{len(rows)} sheet rows, {len(taxids)} unique taxids")
PY

# ── 3. download references: accession first, then taxid, then NCBI ────────
if [[ $SKIP_DOWNLOAD -eq 1 && -s "$PREFIX.fasta" ]]; then
    log "Skipping download, reusing $PREFIX.fasta"
else
    run python3 "$BIN/download_fastas.py" \
        -i "$PREFIX.download_input.tsv" \
        -o "$PREFIX.fasta" \
        -m "$PREFIX.missing.txt" -f file -g "$PREFIX.gcfmapping.tsv" \
        -t "$ASSEMBLY" -k -c 1 --accession_col 2 -a 5 -y 7 \
        --ncbi_fallback ${EMAIL_ARGS[@]+"${EMAIL_ARGS[@]}"}
fi
[[ -s "$PREFIX.fasta" ]] || { log "ERROR: $PREFIX.fasta was not produced"; exit 1; }

cut -f 2 "$PREFIX.gcfmapping.tsv" | sort -u > "$PREFIX.assemblies.txt"
grep '^>' "$PREFIX.fasta" > "headers.$PREFIX.fasta" || true

# ── 4. index + compress ─────────────────────────────────────────────────────
# Always re-index: a .fai left over from an older FASTA silently points at the
# wrong bytes (scripts/pull_contigs.py checks for this, other tools do not).
if command -v samtools >/dev/null; then
    rm -f "$PREFIX.fasta.fai"
    run samtools faidx "$PREFIX.fasta"
else
    rm -f "$PREFIX.fasta.fai"
    log "samtools not found: no .fai written (any old one was removed)"
fi
if [[ $NO_COMPRESS -eq 0 ]]; then
    if command -v bgzip >/dev/null; then
        run bgzip -f -k -@ "$THREADS" "$PREFIX.fasta"
        command -v samtools >/dev/null && run samtools faidx "$PREFIX.fasta.gz"
    elif command -v pigz >/dev/null; then
        run pigz -f -k -p "$THREADS" "$PREFIX.fasta"
    else
        run gzip -f -k "$PREFIX.fasta"
    fi
fi

# ── 5. accession -> taxid map, name check, sheet with assembly column ───────
BACKUP=(); [[ -n "$EMAIL" ]] && BACKUP=(--ncbi-backup -e "$EMAIL")
run python3 "$BIN/append_taxid.py" \
    -i "$PREFIX.gcfmapping.tsv" \
    -r "$ASSEMBLY" \
    -o "$PREFIX.mapping.taxid.tsv" ${BACKUP[@]+"${BACKUP[@]}"}

run python3 "$SCRIPT_DIR/validate_names.py" \
    --mapsheet "$PREFIX.mapping.taxid.tsv" \
    --pathogens_sheet "$SHEET" \
    -o "$PREFIX.mapping.taxid.validated.tsv"

# Also regenerates the download input so the next run goes accession-first.
run python3 "$SCRIPT_DIR/add_assembly_column.py" \
    --pathogens_sheet "$SHEET" \
    --mapsheet "$PREFIX.mapping.taxid.tsv" \
    -o pathogen_sheet.with_assembly.csv \
    --download_input "$PREFIX.download_input.tsv"

# ── 6. summary ──────────────────────────────────────────────────────────────
n_rec=$(grep -c '^>' "$PREFIX.fasta" || true)
n_asm=$(grep -c . "$PREFIX.assemblies.txt" || true)
n_tax=$(grep -c . "$PREFIX.taxids.txt" || true)
n_miss=0; [[ -f "$PREFIX.missing.txt" ]] && n_miss=$(grep -c . "$PREFIX.missing.txt" || true)
log "Done: $n_tax sheet taxids -> $n_asm assemblies, $n_rec sequences; $n_miss missing (see $PREFIX.missing.txt)"
log "Use with: --reference_fasta $(pwd)/$PREFIX.fasta$([[ $NO_COMPRESS -eq 0 ]] && echo .gz)"
