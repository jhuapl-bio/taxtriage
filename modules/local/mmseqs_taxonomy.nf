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
// Emits:
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
    path  db                              // staged seqTaxDB DIRECTORY (contains <prefix>* files)

    output:
    tuple val(meta), path("*_lca.tsv")      , emit: lca
    tuple val(meta), path("*_report")       , emit: report
    tuple val(meta), path("*_tophit_aln")   , emit: tophit
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: '--lca-mode 3 --tax-lineage 1 -s 6.0'   // sensitivity knob
    def prefix   = task.ext.prefix ?: "${meta.id}.novelty"
    def dbprefix = params.novelty_db_prefix ?: 'seqTaxDB'   // preferred base name inside the db dir

    // ----------------------------------------------------------------------------------------
    // Memory-aware target-DB splitting (mirrors the task.memory approach in MINIMAP2_ALIGN).
    // mmseqs holds the target k-mer index in RAM during the prefilter; when the DB doesn't fit
    // it must split it into chunks. A HARDCODED --split-memory-limit is the OOM trap: an 8G
    // limit on a 7G box lets the prefilter allocate past what's available and it gets killed.
    // Instead derive the limit from the Nextflow allocation -- 70% of task.memory, the same
    // headroom factor minimap2 uses for sort -- and let mmseqs auto-optimize the split count
    // (--split 0). Only inject these when the user hasn't already set them in ext.args, so an
    // explicit override still wins. \b is unsafe here (it matches inside --split-memory-limit),
    // so the split-count guard uses a negative lookahead for -memory.
    def ram_70pct_mb = (task.memory.toMega() * 0.70).longValue()
    def mem_flags = ''
    if (!(args =~ /--split-memory-limit/))  { mem_flags += " --split-memory-limit ${ram_70pct_mb}M" }
    if (!(args =~ /--split(?!-memory)/))     { mem_flags += " --split 0" }
    args = (args + mem_flags).trim()
    """
    set -eo pipefail

    # ------------------------------------------------------------------------------------
    # Locate the seqTaxDB base path inside the staged db directory.
    # Auto-downloaded DBs use the configured prefix (default 'seqTaxDB'); a user-supplied
    # local db dir may use any base name, so detect it: the sequence db that carries a
    # taxonomy companion (_mapping / _taxonomy) is the seqTaxDB we want.
    # ------------------------------------------------------------------------------------
    DBDIR='${db}'
    TARGET=""
    if [ -f "\${DBDIR}/${dbprefix}.dbtype" ]; then
        TARGET="\${DBDIR}/${dbprefix}"
    else
        for f in "\${DBDIR}"/*.dbtype; do
            [ -e "\$f" ] || continue
            b="\${f%.dbtype}"
            case "\$b" in *_h) continue ;; esac          # skip the header db
            if [ -e "\${b}_mapping" ] || [ -e "\${b}_taxonomy" ]; then
                TARGET="\$b"; break
            fi
        done
        # last resort: first non-header .dbtype in the dir
        if [ -z "\$TARGET" ]; then
            for f in "\${DBDIR}"/*.dbtype; do
                [ -e "\$f" ] || continue
                b="\${f%.dbtype}"
                case "\$b" in *_h) continue ;; esac
                TARGET="\$b"; break
            done
        fi
    fi
    if [ -z "\$TARGET" ] || [ ! -f "\${TARGET}.dbtype" ]; then
        echo "ERROR: no mmseqs seqTaxDB found in '\${DBDIR}' (looked for prefix '${dbprefix}')." >&2
        echo "Directory contents:" >&2; ls -la "\${DBDIR}" >&2
        exit 1
    fi
    echo "[novelty] using seqTaxDB target: \$TARGET" >&2

    mmseqs createdb ${query} ${prefix}_qdb

    mmseqs taxonomy \\
        ${prefix}_qdb \\
        "\$TARGET" \\
        ${prefix}_taxres \\
        tmp \\
        $args \\
        --threads ${task.cpus}

    # per-query LCA (taxid, rank, name, lineage)
    mmseqs createtsv ${prefix}_qdb ${prefix}_taxres ${prefix}_lca.tsv

    # Kraken-style report -> reuse your existing kreport parsers / krona steps
    mmseqs taxonomyreport "\$TARGET" ${prefix}_taxres ${prefix}_report

    # best-hit alignment table incl. pident, for the identity-tail novelty signal
    mmseqs convertalis ${prefix}_qdb "\$TARGET" ${prefix}_taxres ${prefix}_tophit_aln \\
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
