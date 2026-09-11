// Fetch (and concatenate) the reference FASTAs for a named host target.
//
// Takes the target name plus the space-separated RefSeq assembly accessions
// declared for it in conf/hosts.config, and resolves each accession in the
// order that is cheapest and most reproducible first:
//   1. the NCBI `datasets` CLI, when the image has it
//   2. the NCBI datasets REST API (no CLI needed)
//   3. the assembly_summary the pipeline already downloads -> direct FTP
// All resolved genomes are concatenated into one FASTA, which HOST_REMOVAL
// hands to minimap2 as the de-hosting reference.
//
// The output is written through `storeDir`, so the (large) host genomes are
// downloaded once and re-used by every later run pointing at the same
// --host_reference_dir.
process FETCH_HOST_REFS {
    tag "$target"
    label 'process_single'
    maxForks 1          // be polite to NCBI; these are whole eukaryotic genomes
    storeDir params.host_reference_dir ?: "${params.outdir}/host_references"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/jhuaplbio/taxtriage_confidence:2.1' :
        'jhuaplbio/taxtriage_confidence:2.1' }"

    input:
    tuple val(target), val(accessions), path(assembly_summary)

    output:
    tuple val(target), path("${target}.host.fasta"), emit: fasta
    path "${target}.host.accessions.tsv"           , emit: manifest

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -o pipefail
    mkdir -p parts
    : > ${target}.host.accessions.tsv

    for ACC in ${accessions}; do
        OUT=parts/\${ACC}.fasta

        # (0) A local FASTA path handed to us directly.
        if [ -f "\$ACC" ]; then
            case "\$ACC" in
                *.gz) zcat "\$ACC" > \$OUT ;;
                *)    cp "\$ACC" \$OUT ;;
            esac
        fi

        # (1) NCBI datasets CLI, when present.
        if [ ! -s \$OUT ] && command -v datasets >/dev/null 2>&1; then
            echo "[host-refs] \$ACC: trying datasets CLI" >&2
            rm -rf ds ds.zip
            datasets download genome accession "\$ACC" --include genome --filename ds.zip >/dev/null 2>&1 \\
                && unzip -o -q ds.zip -d ds >/dev/null 2>&1 \\
                && cat ds/ncbi_dataset/data/*/*.fna > \$OUT 2>/dev/null || true
        fi

        # (2) NCBI datasets REST API - same source, no CLI required.
        if [ ! -s \$OUT ]; then
            echo "[host-refs] \$ACC: trying datasets REST API" >&2
            rm -rf ds ds.zip
            curl -sSL --retry 3 --retry-delay 5 \\
                "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/\${ACC}/download?include_annotation_type=GENOME_FASTA" \\
                -o ds.zip \\
                && unzip -o -q ds.zip -d ds >/dev/null 2>&1 \\
                && cat ds/ncbi_dataset/data/*/*.fna > \$OUT 2>/dev/null || true
        fi

        # (3) assembly_summary -> FTP directory -> <basename>_genomic.fna.gz
        if [ ! -s \$OUT ] && [ -s "${assembly_summary}" ] && [ "${assembly_summary}" != "NO_FILE" ]; then
            FTP=\$(awk -F'\\t' -v acc="\$ACC" '\$1 == acc {print \$20; exit}' ${assembly_summary} || true)
            if [ -n "\$FTP" ] && [ "\$FTP" != "na" ]; then
                BASE=\$(basename "\$FTP")
                HTTP=\$(echo "\$FTP" | sed 's|^ftp://|https://|')
                echo "[host-refs] \$ACC: assembly_summary -> \$HTTP" >&2
                curl -sSL --retry 3 --retry-delay 5 "\$HTTP/\${BASE}_genomic.fna.gz" -o ref.fna.gz \\
                    && zcat ref.fna.gz > \$OUT || true
                rm -f ref.fna.gz
            fi
        fi

        if [ ! -s \$OUT ] || ! grep -q '^>' \$OUT; then
            echo "ERROR: could not fetch the host reference for accession '\$ACC' (target '${target}')." >&2
            echo "       Tried the NCBI datasets CLI, the datasets REST API and the pipeline's assembly_summary." >&2
            echo "       Check network access, or point --remove_reference_file at a local FASTA instead." >&2
            exit 1
        fi

        printf '%s\\t%s\\t%s\\n' "${target}" "\$ACC" "\$(grep -c '^>' \$OUT)" >> ${target}.host.accessions.tsv
    done

    cat parts/*.fasta > ${target}.host.fasta
    rm -rf parts ds ds.zip

    echo "[host-refs] ${target}: \$(grep -c '^>' ${target}.host.fasta) sequences from ${accessions}" >&2
    """

    stub:
    """
    echo ">stub_${target}" > ${target}.host.fasta
    echo "ACGT" >> ${target}.host.fasta
    printf '%s\\t%s\\t%s\\n' "${target}" "${accessions}" "1" > ${target}.host.accessions.tsv
    """
}
