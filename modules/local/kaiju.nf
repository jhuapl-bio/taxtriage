// ##############################################################################################
// # Copyright 2022 The Johns Hopkins University Applied Physics Laboratory LLC
// # (header trimmed for the sketch -- copy the full APL header from the other modules)
// ##############################################################################################
//
// KAIJU novelty backend (drop-in alternative to MMSEQS_TAXONOMY).
//
// Greedy translated search of the de novo contigs against a Kaiju protein DB (.fmi). Kaiju does
// its OWN 6-frame translation, so it consumes nucleotide contigs directly -- no Pyrodigal ORF
// prediction is needed (that is the whole point of routing this path around Pyrodigal).
//
// Output contract matches MMSEQS_TAXONOMY so it drops straight into NOVELTY_SCORE:
//   *_lca.tsv      per-contig: query, taxid, rank, name   (rank/name resolved from the db dmps)
//   *_report       kaiju2table summary (parity with the kraken-style report; not score-critical)
//   *_tophit_aln   EMPTY -- kaiju reports match length/score, not %identity, so the low-identity
//                  tail component of the novelty score is simply 0 on this backend.
//
// The db DIRECTORY is expected to contain: nodes.dmp, names.dmp, and exactly one *.fmi index.
//
process KAIJU {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::kaiju=1.10.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kaiju:1.10.1--h43eeafb_0' :
        'quay.io/biocontainers/kaiju:1.10.1--h43eeafb_0' }"

    input:
    tuple val(meta), path(contigs)   // de novo assembly (fasta / fasta.gz)
    path  db                         // kaiju db DIRECTORY (nodes.dmp, names.dmp, *.fmi)

    output:
    tuple val(meta), path("*_lca.tsv")    , emit: lca
    tuple val(meta), path("*_report")     , emit: report
    tuple val(meta), path("*_tophit_aln") , emit: tophit
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: '-a greedy -e 5'   // greedy mode, up to 5 mismatches
    def prefix = task.ext.prefix ?: "${meta.id}.novelty"
    def input_cmd = contigs.name.endsWith('.gz') ? "zcat ${contigs} > _in.fa" : "ln -sf ${contigs} _in.fa"
    """
    set -eo pipefail

    DBDIR='${db}'
    NODES="\${DBDIR}/nodes.dmp"
    NAMES="\${DBDIR}/names.dmp"
    FMI=\$(ls "\${DBDIR}"/*.fmi 2>/dev/null | head -n1)
    if [ -z "\$FMI" ] || [ ! -f "\$NODES" ] || [ ! -f "\$NAMES" ]; then
        echo "ERROR: kaiju db dir '\${DBDIR}' must contain nodes.dmp, names.dmp and a *.fmi index." >&2
        ls -la "\${DBDIR}" >&2
        exit 1
    fi
    echo "[novelty] kaiju db: nodes=\$NODES names=\$NAMES fmi=\$FMI" >&2

    ${input_cmd}

    # Per-contig classification. -v = verbose (taxid in col 3, score/length in col 4).
    kaiju \\
        -t "\$NODES" \\
        -f "\$FMI" \\
        -i _in.fa \\
        -o ${prefix}.kaiju.out \\
        -z ${task.cpus} \\
        -v ${args}

    # ------------------------------------------------------------------------------------
    # Build the NOVELTY_SCORE lca.tsv: query, taxid, rank, name.
    # Kaiju emits per-contig taxids only; resolve rank (nodes.dmp) and scientific name
    # (names.dmp) ourselves so the table matches the mmseqs LCA schema exactly.
    # ------------------------------------------------------------------------------------
    python3 - "\$NODES" "\$NAMES" ${prefix}.kaiju.out ${prefix}_lca.tsv <<'PY'
import sys
nodes_p, names_p, kaiju_p, out_p = sys.argv[1:5]

rank = {}
with open(nodes_p) as fh:
    for line in fh:
        f = [c.strip() for c in line.split('|')]
        if len(f) >= 3:
            rank[f[0]] = f[2]

name = {}
with open(names_p) as fh:
    for line in fh:
        f = [c.strip() for c in line.split('|')]
        if len(f) >= 4 and f[3] == 'scientific name':
            name[f[0]] = f[1]

with open(kaiju_p) as fh, open(out_p, 'w') as out:
    for line in fh:
        cols = line.rstrip('\\n').split('\\t')
        if len(cols) < 3:
            continue
        flag, query, taxid = cols[0], cols[1], cols[2]
        if flag != 'C' or taxid in ('', '0'):
            continue  # unclassified -> dark matter, counted via read accounting
        out.write('\\t'.join([query, taxid, rank.get(taxid, 'no rank'),
                              name.get(taxid, taxid)]) + '\\n')
PY

    # Kraken-style summary table for parity (not consumed by the score).
    kaiju2table -t "\$NODES" -n "\$NAMES" -r genus -o ${prefix}_report ${prefix}.kaiju.out || touch ${prefix}_report

    # Kaiju has no %identity -> empty tophit so the low-identity tail component is 0.
    : > ${prefix}_tophit_aln

    cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	    kaiju: \$(kaiju -h 2>&1 | head -n1 | sed 's/^Kaiju //')
	END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.novelty"
    """
    touch ${prefix}_lca.tsv ${prefix}_report ${prefix}_tophit_aln
    echo '"${task.process}": {kaiju: stub}' > versions.yml
    """
}
