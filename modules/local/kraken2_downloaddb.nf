// ##############################################################################################
// # Copyright 2022 The Johns Hopkins University Applied Physics Laboratory LLC
// # (header trimmed for the sketch -- copy the full APL header from modules/local/download_db.nf)
// ##############################################################################################
//
// Download a prebuilt Kraken2 + Bracken database by a GENERIC NAME, the mmseqs-download way.
//
// Why this exists: mmseqs has `mmseqs databases <Name>` which resolves a name server-side and
// downloads. Kraken2/Bracken have NO such CLI -- the prebuilt indexes live as date-stamped
// tarballs on the genome-idx S3 bucket (Ben Langmead's aws-indexes). So we keep a small alias
// table here (name -> current tarball), wget + extract + validate, and present the result the
// same way MMSEQS_DOWNLOADDB does.
//
// These tarballs serve the BRACKEN backend completely: each one bundles BOTH the Kraken2 db
// (hash.k2d / opts.k2d / taxo.k2d ...) AND the Bracken k-mer distributions
// (database{100,150,200,...}mers.kmer_distrib). The NOVELTY subworkflow points KRAKEN2_NOVELTY
// and BRACKEN at the same directory, so one download covers both.
//
// Accepted db_name forms (case-insensitive for the aliases):
//   * an alias        -> 'standard', 'standard_8gb', 'standard_16gb', 'viral', 'minusb',
//                        'pluspf', 'pluspf_8gb', 'pluspf_16gb', 'pluspfp', 'pluspfp_8gb',
//                        'pluspfp_16gb', 'core_nt', 'gtdb', 'eupathdb', 'ncbi_reference'
//                        ('k2_' prefix optional, e.g. 'k2_standard' == 'standard')
//   * an exact file   -> 'k2_standard_20260226.tar.gz' (joined onto the genome-idx base URL)
//   * a full URL      -> 'https://.../something.tar.gz' (downloaded verbatim)
//
// Caching: `storeDir` writes the finished db to a persistent, per-name folder and SKIPS the
// process on any later run where it already exists -- the multi-GB download happens once.
//
// Output contract: a directory `kraken2_db/` holding the db files at the top level, so both
// KRAKEN2_NOVELTY and BRACKEN can reference `<dir>` directly.
//
process KRAKEN2_DOWNLOADDB {
    tag "$db_name"
    label 'process_medium'

    // Per-name cache dir -> switching --novelty_db doesn't clobber a previous download.
    storeDir { "${params.novelty_kraken2_db_cache}/${db_name.toString().replaceAll('[^A-Za-z0-9._-]', '_')}" }

    conda "conda-forge::gnu-wget=1.18"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h7132678_6' :
        'biocontainers/gnu-wget:1.18--h36e9172_9' }"

    input:
    val(db_name)

    output:
    path("kraken2_db")  , emit: db
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def base_url = params.novelty_kraken2_db_baseurl ?: 'https://genome-idx.s3.amazonaws.com/kraken'
    """
    set -eo pipefail

    name='${db_name}'
    base='${base_url}'

    # ---- resolve the generic name -> download URL --------------------------------------------
    # Aliases pin the latest tarball known at authoring time. To bump them, edit this case block
    # (or just pass an exact filename / full URL as --novelty_db). See genome-idx listing:
    #   https://benlangmead.github.io/aws-indexes/k2
    key=\$(echo "\$name" | tr '[:upper:]' '[:lower:]' | sed 's/^k2_//')
    case "\$name" in
        http://*|https://*)  url="\$name" ;;                 # full URL, verbatim
        *.tar.gz|*.tgz)      url="\$base/\$name" ;;          # exact filename on the base bucket
        *)
            case "\$key" in
                standard)                  file='k2_standard_20260226.tar.gz' ;;
                standard_8gb|standard_08gb)file='k2_standard_08_GB_20260226.tar.gz' ;;
                standard_16gb)             file='k2_standard_16_GB_20260226.tar.gz' ;;
                viral|viruses)             file='k2_viral_20260226.tar.gz' ;;
                minusb)                    file='k2_minusb_20260226.tar.gz' ;;
                pluspf)                    file='k2_pluspf_20260226.tar.gz' ;;
                pluspf_8gb|pluspf_08gb)    file='k2_pluspf_08_GB_20260226.tar.gz' ;;
                pluspf_16gb)               file='k2_pluspf_16_GB_20260226.tar.gz' ;;
                pluspfp)                   file='k2_pluspfp_20260226.tar.gz' ;;
                pluspfp_8gb|pluspfp_08gb)  file='k2_pluspfp_08_GB_20260226.tar.gz' ;;
                pluspfp_16gb)              file='k2_pluspfp_16_GB_20260226.tar.gz' ;;
                core_nt|corent|nt)         file='k2_core_nt_20251015.tar.gz' ;;
                gtdb|gtdb_genome_reps)     file='k2_gtdb_genome_reps_20250609.tar.gz' ;;
                eupathdb|eupathdb48)       file='k2_eupathdb48_20230407.tar.gz' ;;
                ncbi_reference|prackendb)  file='k2_NCBI_reference_20251007.tar.gz' ;;
                *)
                    echo "ERROR: unknown kraken2/bracken db name '\$name'." >&2
                    echo "Valid aliases: standard, standard_8gb, standard_16gb, viral, minusb," >&2
                    echo "  pluspf[_8gb|_16gb], pluspfp[_8gb|_16gb], core_nt, gtdb, eupathdb," >&2
                    echo "  ncbi_reference. Or pass an exact *.tar.gz filename or a full URL." >&2
                    exit 1 ;;
            esac
            url="\$base/\$file" ;;
    esac

    echo "[novelty] kraken2/bracken db '\$name' -> \$url" >&2

    mkdir -p kraken2_db
    # -nv quiet-ish; retries for flaky S3; stream straight into tar to avoid a 2nd disk copy.
    wget -nv --tries=3 --continue -O db.tar.gz "\$url" ${args}
    tar -xzf db.tar.gz -C kraken2_db
    rm -f db.tar.gz

    # Langmead tarballs extract files at the top level, but tolerate a single wrapper dir.
    if [ ! -f kraken2_db/hash.k2d ]; then
        inner=\$(find kraken2_db -maxdepth 2 -name hash.k2d -print -quit)
        if [ -n "\$inner" ]; then
            mv "\$(dirname "\$inner")"/* kraken2_db/ 2>/dev/null || true
        fi
    fi

    # ---- validate -----------------------------------------------------------------------------
    for f in hash.k2d opts.k2d taxo.k2d; do
        if [ ! -f "kraken2_db/\$f" ]; then
            echo "ERROR: extracted kraken2 db is missing \$f (got: \$(ls kraken2_db))" >&2
            exit 1
        fi
    done
    if ! ls kraken2_db/*.kmer_distrib >/dev/null 2>&1; then
        echo "[novelty] WARNING: no *.kmer_distrib in the db -- BRACKEN needs one matching" >&2
        echo "          --novelty_bracken_readlen (e.g. database100mers.kmer_distrib)." >&2
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version 2>/dev/null | head -n1 | grep -oE '[0-9.]+' | head -1 || echo "unknown")
    END_VERSIONS
    """

    stub:
    def readlen = params.novelty_bracken_readlen ?: 100
    """
    mkdir -p kraken2_db
    touch kraken2_db/hash.k2d kraken2_db/opts.k2d kraken2_db/taxo.k2d
    touch kraken2_db/seqid2taxid.map kraken2_db/database${readlen}mers.kmer_distrib
    echo '"${task.process}": {wget: stub}' > versions.yml
    """
}
