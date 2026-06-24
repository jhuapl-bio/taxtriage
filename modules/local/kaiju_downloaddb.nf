// ##############################################################################################
// # Copyright 2022 The Johns Hopkins University Applied Physics Laboratory LLC
// # (header trimmed for the sketch -- copy the full APL header from modules/local/download_db.nf)
// ##############################################################################################
//
// Download a prebuilt Kaiju protein index by a GENERIC NAME, the mmseqs-download way.
//
// Kaiju (like Kraken2/Bracken) has no `databases <Name>` CLI -- `kaiju-makedb` builds from
// source, which is huge and slow. The ready-made indexes live as date-stamped .tgz archives on
// the kaiju-idx S3 bucket, each containing the `.fmi` index plus `nodes.dmp` and `names.dmp`.
// We keep a small alias table here, wget + extract + validate, and present the result exactly
// like MMSEQS_DOWNLOADDB so the KAIJU module can consume it unchanged.
//
// Accepted db_name forms (case-insensitive for the aliases):
//   * 'test'     -> tiny git-lfs-hosted index for CI / -profile test (parity with kraken2 'test')
//   * an alias   -> 'viruses' (a.k.a. 'viral'), 'nr', 'nr_euk', 'refseq', 'refseq_nr',
//                   'refseq_ref', 'progenomes', 'fungi', 'plasmids', 'rvdb'
//   * an exact file -> 'kaiju_db_viruses_2024-08-15.tgz' (joined onto the base bucket URL)
//   * a full URL    -> 'https://.../kaiju_db_*.tgz' (downloaded verbatim)
//
// Output contract (parity with MMSEQS_DOWNLOADDB): a directory `kaiju_db/` containing exactly one
// `*.fmi`, plus `nodes.dmp` and `names.dmp` at the top level -- the layout the KAIJU module expects.
//
process KAIJU_DOWNLOADDB {
    tag "$db_name"
    label 'process_medium'

    // 'test' is a tiny CI index -> no persistent storeDir; it just lives in work/ like the
    // main kraken2 'test' download. Real indexes still persist to the per-name dbs/kaiju cache.
    storeDir { db_name.toString().equalsIgnoreCase('test') ? null :
               "${params.novelty_kaiju_db_cache}/${db_name.toString().replaceAll('[^A-Za-z0-9._-]', '_')}" }

    conda "conda-forge::gnu-wget=1.18"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h7132678_6' :
        'biocontainers/gnu-wget:1.18--h36e9172_9' }"

    input:
    val(db_name)

    output:
    path("kaiju_db")    , emit: db
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def base_url = params.novelty_kaiju_db_baseurl ?: 'https://kaiju-idx.s3.eu-central-1.amazonaws.com'
    """
    set -eo pipefail

    name='${db_name}'
    base='${base_url}'

    # ---- resolve the generic name -> download URL --------------------------------------------
    # Aliases pin the latest .tgz known at authoring time. Bump them here, or pass an exact
    # filename / full URL via --novelty_db. Listing:
    #   https://bioinformatics-centre.github.io/kaiju/downloads.html
    key=\$(echo "\$name" | tr '[:upper:]' '[:lower:]')
    case "\$name" in
        http://*|https://*)  url="\$name" ;;                 # full URL, verbatim
        *.tgz|*.tar.gz)      url="\$base/\$name" ;;          # exact filename on the base bucket
        [Tt][Ee][Ss][Tt])   # tiny git-lfs-hosted index for CI / -profile test (parity with kraken2 'test')
            url='https://github.com/Merritt-Brian/databases/raw/refs/heads/main/kaiju/test_kaiju.tar.gz' ;;
        *)
            case "\$key" in
                viruses|viral) rel='2024/kaiju_db_viruses_2024-08-15.tgz' ;;
                nr)            rel='2024/kaiju_db_nr_2024-08-25.tgz' ;;
                nr_euk)        rel='2023/kaiju_db_nr_euk_2023-05-10.tgz' ;;
                refseq)        rel='2024/kaiju_db_refseq_2024-08-14.tgz' ;;
                refseq_nr)     rel='2024/kaiju_db_refseq_nr_2024-08-13.tgz' ;;
                refseq_ref)    rel='2024/kaiju_db_refseq_ref_2024-08-14.tgz' ;;
                progenomes)    rel='2023/kaiju_db_progenomes_2023-05-25.tgz' ;;
                fungi)         rel='2024/kaiju_db_fungi_2024-08-16.tgz' ;;
                plasmids)      rel='2024/kaiju_db_plasmids_2024-08-15.tgz' ;;
                rvdb)          rel='2024/kaiju_db_rvdb_2024-12-20.tgz' ;;
                *)
                    echo "ERROR: unknown kaiju db name '\$name'." >&2
                    echo "Valid aliases: test, viruses (viral), nr, nr_euk, refseq, refseq_nr," >&2
                    echo "  refseq_ref, progenomes, fungi, plasmids, rvdb. Or pass an exact" >&2
                    echo "  *.tgz filename or a full URL." >&2
                    exit 1 ;;
            esac
            url="\$base/\$rel" ;;
    esac

    echo "[novelty] kaiju db '\$name' -> \$url" >&2

    mkdir -p kaiju_db
    wget -nv --tries=3 --continue -O db.tgz "\$url" ${args}
    tar -xzf db.tgz -C kaiju_db
    rm -f db.tgz

    # Tolerate a single wrapper dir inside the archive.
    if ! ls kaiju_db/*.fmi >/dev/null 2>&1; then
        inner=\$(find kaiju_db -maxdepth 2 -name '*.fmi' -print -quit)
        if [ -n "\$inner" ]; then
            mv "\$(dirname "\$inner")"/* kaiju_db/ 2>/dev/null || true
        fi
    fi

    # ---- validate -----------------------------------------------------------------------------
    if ! ls kaiju_db/*.fmi >/dev/null 2>&1; then
        echo "ERROR: no *.fmi index in extracted kaiju db (got: \$(ls kaiju_db))" >&2
        exit 1
    fi
    for f in nodes.dmp names.dmp; do
        if [ ! -f "kaiju_db/\$f" ]; then
            echo "ERROR: extracted kaiju db is missing \$f (got: \$(ls kaiju_db))" >&2
            exit 1
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version 2>/dev/null | head -n1 | grep -oE '[0-9.]+' | head -1 || echo "unknown")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p kaiju_db
    touch kaiju_db/kaiju_db.fmi kaiju_db/nodes.dmp kaiju_db/names.dmp
    echo '"${task.process}": {wget: stub}' > versions.yml
    """
}
