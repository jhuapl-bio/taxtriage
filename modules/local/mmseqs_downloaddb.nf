// ##############################################################################################
// # Copyright 2022 The Johns Hopkins University Applied Physics Laboratory LLC
// # (header trimmed for the sketch -- copy the full APL header from modules/local/download_db.nf)
// ##############################################################################################
//
// Download a prebuilt mmseqs seqTaxDB via `mmseqs databases`, the same way DOWNLOAD_DB pulls
// the Kraken2 db. Mirrors the centrifuge/kraken2 behaviour: the WORKFLOW checks for a local
// path first and only calls this when nothing is found locally.
//
// Caching: `storeDir` writes the finished db to a persistent, per-db-name folder and SKIPS the
// process entirely on any later run where those files already exist. That covers both `-resume`
// (within a run) and fresh runs (across runs) -- the multi-GB download happens exactly once.
//
// Output contract: a directory `mmseqs_db/` containing the db under a fixed prefix
// (`params.novelty_db_prefix`, default `seqTaxDB`), so everything downstream references it as
// `<dir>/seqTaxDB` regardless of which source db was downloaded.
//
process MMSEQS_DOWNLOADDB {
    tag "$db_name"
    label 'process_medium'

    // Per-db-name cache dir -> switching --novelty_db doesn't clobber a previous download,
    // and each is reused independently. Closure form lets the directive see the input value.
    storeDir { "${params.novelty_db_cache}/${db_name.replaceAll('[^A-Za-z0-9._-]', '_')}" }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mmseqs2:18.8cc5c--hd6d6fdc_0' :
        'staphb/mmseqs2' }"

    input:
    val(db_name)        // e.g. UniProtKB, UniProtKB/Swiss-Prot, UniRef50, UniRef90, GTDB, NR ...

    output:
    path("mmseqs_db")        , emit: db
    path "versions.yml"      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = params.novelty_db_prefix ?: 'seqTaxDB'
    """
    mkdir -p mmseqs_db tmp_dl

    # `mmseqs databases <Name> <outDB> <tmp>` -- e.g. UniProtKB/Swiss-Prot -> mmseqs_db/seqTaxDB
    mmseqs databases \\
        '${db_name}' \\
        mmseqs_db/${prefix} \\
        tmp_dl \\
        --threads ${task.cpus} ${args}

    rm -rf tmp_dl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs version)
    END_VERSIONS
    """

    stub:
    def prefix = params.novelty_db_prefix ?: 'seqTaxDB'
    """
    mkdir -p mmseqs_db
    touch mmseqs_db/${prefix} mmseqs_db/${prefix}.dbtype mmseqs_db/${prefix}_taxonomy
    echo '"${task.process}": {mmseqs: stub}' > versions.yml
    """
}
