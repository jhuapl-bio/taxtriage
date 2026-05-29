#!/usr/bin/env bash
# Copyright 2022 The Johns Hopkins University Applied Physics Laboratory LLC
# All rights reserved.
# Permission is hereby granted, free of charge, to any person obtaining a copy of this
# software and associated documentation files (the "Software"), to deal in the Software
# without restriction, including without limitation the rights to use, copy, modify,
# merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
# PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
# OR OTHER DEALINGS IN THE SOFTWARE.

set -euo pipefail

# ---------------------------------------------------------------------------
# remove_hostreads.sh
#
# Extract non-host reads from a host-alignment BAM.
#
# Single-end mode: FLAG-based (samtools view -f 4).
#
# Paired-end mode: tries the fast FLAG-based path first (-f 12), then falls
# back to a QNAME-level FASTQ scan if the fast path returns 0 reads while the
# BAM still contains unmapped reads.  The fallback handles:
#   • Illumina/Casava 1.8+ headers (description field after whitespace)
#   • Desynchronised R1/R2 files (orphan reads, different QNAME order)
#     where minimap2 emits singletons (FLAG 0x1 never set) so -f 12 fails
#
# USAGE
#   remove_hostreads.sh -b BAM -p PREFIX [options]
#
# REQUIRED
#   -b FILE    input BAM (host alignment, coordinate-sorted)
#   -p STR     sample ID / output prefix
#
# OPTIONAL
#   -1 FILE    R1 FASTQ (required for paired-end mode)
#   -2 FILE    R2 FASTQ (required for paired-end mode)
#   -q INT     minimum MAPQ; reads below this threshold are kept as non-host
#              even if they have a mapping position (default: 0 = disabled)
#   -s FILE    output path for non-host singleton / orphan reads
#              (default: /dev/null — orphans are discarded)
#   -e         single-end mode (default: paired-end when -1/-2 are supplied)
#   -i         include singletons in host removal (use -f 4 instead of -f 12
#              for the fast path; incompatible with -d)
#   -d         denovo / diamond mode: force pair-both-unmapped (-f 12) even
#              when -i is set (Megahit requires equal R1/R2 read counts)
#   -x STR     extra arguments forwarded to samtools fastq
# ---------------------------------------------------------------------------

usage() {
cat <<EOF
Usage: remove_hostreads.sh -b BAM -p PREFIX [options]

  -b FILE    input BAM (host alignment)
  -p STR     sample prefix for output filenames
  -1 FILE    R1 FASTQ  (paired-end)
  -2 FILE    R2 FASTQ  (paired-end)
  -q INT     min MAPQ threshold (default: 0)
  -s FILE    singleton/orphan output path (default: /dev/null)
  -e         single-end mode
  -i         include singletons (use -f 4 instead of -f 12)
  -d         denovo/diamond mode (forces -f 12 regardless of -i)
  -x STR     extra args for samtools fastq
  -h         show this help
EOF
}

# ── Defaults ────────────────────────────────────────────────────────────────
BAM=""
R1=""
R2=""
PREFIX=""
MIN_MAPQ=0
SINGLETON_OUT="/dev/null"
SINGLE_END=false
INCLUDE_SINGLETONS=false
DENOVO_MODE=false
EXTRA_ARGS=""

while getopts "b:1:2:p:q:s:eidx:h" OPT; do
    case "$OPT" in
        b) BAM="$OPTARG" ;;
        1) R1="$OPTARG" ;;
        2) R2="$OPTARG" ;;
        p) PREFIX="$OPTARG" ;;
        q) MIN_MAPQ="$OPTARG" ;;
        s) SINGLETON_OUT="$OPTARG" ;;
        e) SINGLE_END=true ;;
        i) INCLUDE_SINGLETONS=true ;;
        d) DENOVO_MODE=true ;;
        x) EXTRA_ARGS="$OPTARG" ;;
        h) usage; exit 0 ;;
        ?) usage; exit 1 ;;
    esac
done

# ── Validate ─────────────────────────────────────────────────────────────────
[[ -z "$BAM"    ]] && { echo "ERROR: -b BAM is required";    usage; exit 1; }
[[ -z "$PREFIX" ]] && { echo "ERROR: -p PREFIX is required"; usage; exit 1; }
[[ ! -f "$BAM"  ]] && { echo "ERROR: BAM file not found: $BAM"; exit 1; }

if ! $SINGLE_END; then
    [[ -z "$R1" || -z "$R2" ]] && {
        echo "ERROR: paired-end mode requires -1 R1 and -2 R2 (or pass -e for single-end)"
        exit 1
    }
fi

# ── Helpers ───────────────────────────────────────────────────────────────────

count_reads_in_fastqs() {
    # Sum read counts across all output fastq.gz files for this prefix.
    local total=0
    for fq in "${PREFIX}"*.hostremoved.fastq.gz "${PREFIX}"*.hostremoved.fq.gz; do
        [[ -f "$fq" ]] || continue
        local lines
        lines=$(gzip -dc "$fq" | wc -l)
        total=$(( total + lines / 4 ))
    done
    echo "$total"
}

# Filter a gzipped FASTQ, keeping only reads whose base QNAME is in ids_file.
# Handles:
#   - Illumina/Casava description field (strips everything after first space/tab)
#   - Legacy /1 /2 read-number suffixes
filter_fastq_by_ids() {
    local ids_file="$1"
    local fq_in="$2"
    gzip -dc "$fq_in" | awk -v ids_file="$ids_file" '
        BEGIN { while ((getline line < ids_file) > 0) keep[line] = 1 }
        FNR % 4 == 1 {
            id = $0
            sub(/^@/, "", id)
            split(id, a, /[ \t]/)
            id = a[1]
            sub(/\/[12]$/, "", id)
            in_keep = (id in keep)
        }
        in_keep { print }
    '
}

write_stats() {
    local total="$1"
    local retained="$2"
    local removed=$(( total - retained ))
    local pct_removed="0.00"
    local all_removed="NO"

    if [[ "$total" -gt 0 ]]; then
        pct_removed=$(awk "BEGIN {printf \"%.2f\", ($removed / $total) * 100}")
    fi

    echo "==========================================="
    echo "  HOST REMOVAL STATS: ${PREFIX}"
    echo "==========================================="
    echo "  Total input reads:    $total"
    echo "  Retained (non-host):  $retained"
    echo "  Removed (host):       $removed (${pct_removed}%)"
    if [[ "$retained" -eq 0 && "$total" -gt 0 ]]; then
        echo ""
        echo "  WARNING: ALL reads were classified as host for sample '${PREFIX}'."
        echo "  The sample will fall back to its original (unfiltered) reads."
        echo "  Consider checking the host reference or adjusting --min_mapq_host."
        echo ""
        all_removed="YES"
    fi
    echo "==========================================="

    printf "Sample\tTotal Reads\tRetained Reads\tRemoved (Host) Reads\tPercent Host\tAll Removed\n" \
        > "${PREFIX}.host_removal_stats_mqc.tsv"
    printf "%s\t%s\t%s\t%s\t%s%%\t%s\n" \
        "$PREFIX" "$total" "$retained" "$removed" "$pct_removed" "$all_removed" \
        >> "${PREFIX}.host_removal_stats_mqc.tsv"
}

# ── Single-end path ──────────────────────────────────────────────────────────
if $SINGLE_END; then
    total_reads=$(samtools view -c -F 0x900 "$BAM")

    if [[ "$MIN_MAPQ" -gt 0 ]]; then
        view_filter="-F 0x900 -e flag.unmap || mapq < ${MIN_MAPQ}"
        samtools view -b -F 0x900 -e "flag.unmap || mapq < ${MIN_MAPQ}" "$BAM" \
            | samtools fastq -n $EXTRA_ARGS \
              -0 "${PREFIX}.hostremoved.fastq.gz" -s /dev/null -
    else
        samtools view -b -f 4 -F 0x900 "$BAM" \
            | samtools fastq -n $EXTRA_ARGS \
              -0 "${PREFIX}.hostremoved.fastq.gz" -s /dev/null -
    fi

    retained_reads=$(count_reads_in_fastqs)
    write_stats "$total_reads" "$retained_reads"
    exit 0
fi

# ── Paired-end path ──────────────────────────────────────────────────────────
total_reads=$(samtools view -c -F 0x900 "$BAM")

# Determine the fast-path view filter (FLAG-based, operates on BAM directly).
# Denovo/diamond mode forces -f 12 to keep R1/R2 counts equal for Megahit.
if [[ "$MIN_MAPQ" -gt 0 ]]; then
    if $DENOVO_MODE || ! $INCLUDE_SINGLETONS; then
        FAST_FILTER="-F 0x900 -e (flag & 12) == 12 || mapq < ${MIN_MAPQ}"
    else
        FAST_FILTER="-F 0x900 -e flag.unmap || mapq < ${MIN_MAPQ}"
    fi
else
    if $DENOVO_MODE || ! $INCLUDE_SINGLETONS; then
        FAST_FILTER="-f 12 -F 0x900"
    else
        FAST_FILTER="-f 4 -F 0x900"
    fi
fi

# Fallback filter: per-read unmapped (no mate requirement).
if [[ "$MIN_MAPQ" -gt 0 ]]; then
    FALLBACK_FILTER="-F 0x900 -e flag.unmap || mapq < ${MIN_MAPQ}"
else
    FALLBACK_FILTER="-f 4 -F 0x900"
fi

# ── Step 1: fast FLAG-based extraction ──────────────────────────────────────
echo "INFO: Attempting fast FLAG-based host removal for '${PREFIX}'..."

# samtools view with expression filters needs the -e flag separately; for
# simple flag filters we use the combined string directly.
if [[ "$MIN_MAPQ" -gt 0 ]]; then
    if $DENOVO_MODE || ! $INCLUDE_SINGLETONS; then
        samtools view -b -F 0x900 -e "(flag & 12) == 12 || mapq < ${MIN_MAPQ}" "$BAM"
    else
        samtools view -b -F 0x900 -e "flag.unmap || mapq < ${MIN_MAPQ}" "$BAM"
    fi
else
    if $DENOVO_MODE || ! $INCLUDE_SINGLETONS; then
        samtools view -b -f 12 -F 0x900 "$BAM"
    else
        samtools view -b -f 4 -F 0x900 "$BAM"
    fi
fi | samtools fastq -n $EXTRA_ARGS \
      -1 "${PREFIX}_1.hostremoved.fastq.gz" \
      -2 "${PREFIX}_2.hostremoved.fastq.gz" \
      -0 "$SINGLETON_OUT" -

retained_reads=$(count_reads_in_fastqs)

# ── Step 2: QNAME fallback if fast path returned nothing ────────────────────
# Trigger condition: fast path got 0 reads AND the BAM has unmapped reads.
#   - 0 retained + unmapped reads present  → singletons/desynced BAM → fallback
#   - 0 retained + no unmapped reads       → all reads are host     → no fallback
if [[ "$retained_reads" -eq 0 && "$total_reads" -gt 0 ]]; then
    unmapped_in_bam=$(samtools view -c -f 4 -F 0x900 "$BAM")

    if [[ "$unmapped_in_bam" -gt 0 ]]; then
        echo "INFO: Fast path returned 0 reads; BAM has ${unmapped_in_bam} unmapped reads."
        echo "INFO: BAM likely lacks paired-end flags (desynchronised or Casava-format input)."
        echo "INFO: Retrying with QNAME-based FASTQ filtering..."

        # Clean up empty outputs from the fast attempt
        rm -f "${PREFIX}_1.hostremoved.fastq.gz" "${PREFIX}_2.hostremoved.fastq.gz"
        [[ "$SINGLETON_OUT" != "/dev/null" ]] && rm -f "$SINGLETON_OUT"

        # Extract per-read non-host QNAMEs from BAM
        if [[ "$MIN_MAPQ" -gt 0 ]]; then
            samtools view -F 0x900 -e "flag.unmap || mapq < ${MIN_MAPQ}" "$BAM" \
                | awk '{id=$1; sub(/\/[12]$/, "", id); print id}' \
                | sort -u > nonhost_bam_ids.txt
        else
            samtools view -f 4 -F 0x900 "$BAM" \
                | awk '{id=$1; sub(/\/[12]$/, "", id); print id}' \
                | sort -u > nonhost_bam_ids.txt
        fi

        # Filter R1 and R2 independently against the BAM ID list
        filter_fastq_by_ids nonhost_bam_ids.txt "$R1" | gzip -c > r1_nonhost.fastq.gz
        filter_fastq_by_ids nonhost_bam_ids.txt "$R2" | gzip -c > r2_nonhost.fastq.gz

        # Intersect R1 and R2 non-host IDs → equal-length paired output files
        gzip -dc r1_nonhost.fastq.gz \
            | awk 'NR%4==1{sub(/^@/,""); split($0,a," "); sub(/\/[12]$/,"",a[1]); print a[1]}' \
            | sort -u > r1_nonhost_ids.txt
        gzip -dc r2_nonhost.fastq.gz \
            | awk 'NR%4==1{sub(/^@/,""); split($0,a," "); sub(/\/[12]$/,"",a[1]); print a[1]}' \
            | sort -u > r2_nonhost_ids.txt

        comm -12 r1_nonhost_ids.txt r2_nonhost_ids.txt > paired_nonhost_ids.txt

        filter_fastq_by_ids paired_nonhost_ids.txt "$R1" \
            | gzip -c > "${PREFIX}_1.hostremoved.fastq.gz"
        filter_fastq_by_ids paired_nonhost_ids.txt "$R2" \
            | gzip -c > "${PREFIX}_2.hostremoved.fastq.gz"

        # Optional singletons: non-host reads whose mate mapped to host (or was absent)
        if [[ "$SINGLETON_OUT" != "/dev/null" ]]; then
            comm -23 r1_nonhost_ids.txt r2_nonhost_ids.txt > r1_orphan_ids.txt
            comm -13 r1_nonhost_ids.txt r2_nonhost_ids.txt > r2_orphan_ids.txt
            { filter_fastq_by_ids r1_orphan_ids.txt "$R1"
              filter_fastq_by_ids r2_orphan_ids.txt "$R2"; } \
                | gzip -c > "$SINGLETON_OUT"
        fi

        retained_reads=$(count_reads_in_fastqs)
    fi
fi

write_stats "$total_reads" "$retained_reads"
