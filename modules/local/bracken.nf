// ##############################################################################################
// # Copyright 2022 The Johns Hopkins University Applied Physics Laboratory LLC
// # (header trimmed for the sketch -- copy the full APL header from the other modules)
// ##############################################################################################
//
// BRACKEN novelty backend (closed-set, count-weighted).
//
// Bracken is NOT a translated/open-set search like mmseqs2 or kaiju -- it re-estimates
// abundance from a Kraken2 report. To fit the novelty model we run Kraken2 on the de novo
// contigs upstream (KRAKEN2_NOVELTY in the subworkflow) and hand its report to bracken here.
//
// We then flatten bracken's adjusted kraken-style report (`-w`) into the SAME lca.tsv schema
// NOVELTY_SCORE consumes, but COUNT-WEIGHTED: each taxon node emits ONE row carrying its
// taxon-level read count (NOVELTY_SCORE is told to read that count column instead of weighting
// each row as 1). Using taxon_reads (reads assigned AT that node, not the clade total) keeps
// every read counted once -- after bracken redistribution most land at species, and whatever
// remains at genus+ is exactly the "could only place above species" signal the score wants.
//
// Caveat carried into the score: the read-accounting denominators (total/ref/k2) are read-level
// from the BAM, while these counts are over CONTIGS, so dark_fraction here is approximate.
//
// Output contract (parity with MMSEQS_TAXONOMY / KAIJU):
//   *_lca.tsv      query, taxid, rank, name, count   (count consumed via --count-col)
//   *_report       bracken-adjusted kraken-style report
//   *_tophit_aln   EMPTY (no per-hit identity) -> low-identity tail component is 0
//
process BRACKEN {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bracken=2.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bracken:2.9--py39h1f90b4d_0' :
        'quay.io/biocontainers/bracken:2.9--py39h1f90b4d_0' }"

    input:
    tuple val(meta), path(kreport)   // Kraken2 report from KRAKEN2_NOVELTY (contigs)
    path  db                         // kraken2 db DIRECTORY (must hold databaseNmers.kmer_distrib)

    output:
    tuple val(meta), path("*_lca.tsv")    , emit: lca
    tuple val(meta), path("*_report")     , emit: report
    tuple val(meta), path("*_tophit_aln") , emit: tophit
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}.novelty"
    def readlen = params.novelty_bracken_readlen ?: 100
    def level   = params.novelty_bracken_level   ?: 'S'
    """
    set -eo pipefail

    # Bracken can exit non-zero when nothing was classified; don't let that kill the run --
    # we still emit empty outputs so NOVELTY_SCORE produces a (zeroed) row for the sample.
    bracken \\
        -d ${db} \\
        -i ${kreport} \\
        -o ${prefix}.bracken.tsv \\
        -w ${prefix}_report \\
        -r ${readlen} \\
        -l ${level} ${args} || echo "[novelty] bracken returned non-zero (likely nothing classified)" >&2

    [ -f ${prefix}_report ] || cp ${kreport} ${prefix}_report

    # Flatten the bracken-adjusted report into the count-weighted lca schema.
    python3 - ${prefix}_report ${prefix}_lca.tsv <<'PY'
import sys
rep_p, out_p = sys.argv[1:3]

# kraken/bracken rank code -> full rank word used by novelty_score's HIGH_RANKS set.
CODE = {'D': 'superkingdom', 'K': 'kingdom', 'P': 'phylum', 'C': 'class',
        'O': 'order', 'F': 'family', 'G': 'genus', 'S': 'species'}

rows = []
with open(rep_p) as fh:
    for line in fh:
        c = line.rstrip('\\n').split('\\t')
        if len(c) < 6:
            continue
        taxon_reads = c[2].strip()
        code = c[3].strip()
        taxid = c[4].strip()
        name = c[5].strip()
        if code in ('U', 'R') or taxid in ('', '0'):
            continue
        try:
            n = int(taxon_reads)
        except ValueError:
            continue
        if n <= 0:
            continue
        rank = CODE.get(code[0], 'no rank')   # 'S1','G1',... -> first char
        rows.append((taxid, rank, name, n))

with open(out_p, 'w') as out:
    for taxid, rank, name, n in rows:
        out.write('\\t'.join([taxid, taxid, rank, name, str(n)]) + '\\n')
PY

    : > ${prefix}_tophit_aln

    cat <<END_VERSIONS > versions.yml
"${task.process}":
    bracken: \$(bracken 2>&1 | grep -oP '(?i)(?<=version )[0-9.]+' | head -1 || echo "unknown")
END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.novelty"
    """
    touch ${prefix}_lca.tsv ${prefix}_report ${prefix}_tophit_aln
    echo '"${task.process}": {bracken: stub}' > versions.yml
    """
}
