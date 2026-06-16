// ##############################################################################################
// # Copyright 2022 The Johns Hopkins University Applied Physics Laboratory LLC
// # (header trimmed for the sketch -- copy the full APL header from the other modules)
// ##############################################################################################
//
// Translated-search taxonomic assignment (LCA) for reads/contigs that fell through
// Kraken2 + reference alignment. Protein space stays alignable across genus/family/order
// divergence, so this recovers a high-rank call for organisms with no nucleotide neighbor.
//
// Reuses the mmseqs2 container already vendored for MMSEQS_EASYCLUSTER.
// `easy-taxonomy` emits:
//   *_lca.tsv        per-query: query, taxid, rank, name, ...
//   *_report         Kraken-style report (so it drops into your existing kreport tooling)
//   *_tophit_aln     per-query best hit incl. % identity (pident) -> feeds the novelty score
//
process MMSEQS_TAXONOMY {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mmseqs2:18.8cc5c--hd6d6fdc_0' :
        'staphb/mmseqs2' }"

    input:
    tuple val(meta), path(query)          // unmapped/unclassified reads OR de novo contigs (fasta/fastq)
    path  db                              // prebuilt mmseqs seqTaxDB (e.g. UniRef50 + taxonomy), staged once

    output:
    tuple val(meta), path("*_lca.tsv")      , emit: lca
    tuple val(meta), path("*_report")       , emit: report
    tuple val(meta), path("*_tophit_aln")   , emit: tophit
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: '--lca-mode 3 --tax-lineage 1 -s 6.0'   // sensitivity knob
    def prefix = task.ext.prefix ?: "${meta.id}.novelty"
    """
    mmseqs createdb ${query} ${prefix}_qdb

    mmseqs taxonomy \\
        ${prefix}_qdb \\
        ${db}/targetDB \\
        ${prefix}_taxres \\
        tmp \\
        $args \\
        --threads ${task.cpus}

    # per-query LCA (taxid, rank, name, lineage)
    mmseqs createtsv ${prefix}_qdb ${prefix}_taxres ${prefix}_lca.tsv

    # Kraken-style report -> reuse your existing kreport parsers / krona steps
    mmseqs taxonomyreport ${db}/targetDB ${prefix}_taxres ${prefix}_report

    # best-hit alignment table incl. pident, for the identity-tail novelty signal
    mmseqs filtertaxdb ${db}/targetDB ${prefix}_taxres ${prefix}_tophit || true
    mmseqs convertalis ${prefix}_qdb ${db}/targetDB ${prefix}_taxres ${prefix}_tophit_aln \\
        --format-output "query,target,pident,evalue,bits,taxid,taxname" || touch ${prefix}_tophit_aln

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.novelty"
    """
    touch ${prefix}_lca.tsv ${prefix}_report ${prefix}_tophit_aln
    echo '"${task.process}": {mmseqs: stub}' > versions.yml
    """
}

// -----------------------------------------------------------------------------------------
// DROP-IN ALTERNATIVE: Kaiju instead of mmseqs (greedy translated search, very fast on reads).
// Same contract: emit a per-read taxon+rank table and a names-resolved table.
//
//   kaiju  -t nodes.dmp -f kaiju_db.fmi -i reads.fq -o ${prefix}.kaiju.out -z ${task.cpus} -a greedy -e 5
//   kaiju-addTaxonNames -t nodes.dmp -n names.dmp -i ${prefix}.kaiju.out -o ${prefix}_lca.tsv -r superkingdom,phylum,class,order,family,genus,species
//   kaiju2table -t nodes.dmp -n names.dmp -r genus -o ${prefix}_report ${prefix}.kaiju.out
//
// Kaiju gives match length/score rather than pident; in novelty_score.py use
// (match_aa_length / read_aa_length) as the divergence proxy instead of pident.
// -----------------------------------------------------------------------------------------
