// Fetch the reference FASTA for one spike-in accession.
//
// Accepts an assembly accession (GCF_/GCA_) or a nuccore accession, and resolves
// it in the order that is cheapest and most reproducible first:
//   1. a row in the assembly_summary the pipeline already downloads -> direct FTP
//   2. the NCBI `datasets` CLI, when the image has it
//   3. Entrez efetch (the only route for a bare nuccore accession)
// A local FASTA path in the accession column is passed straight through, so a
// user can spike in a genome that has no accession at all.
process FETCH_SPIKEIN_REFS {
    tag "$accession"
    label 'process_single'
    maxForks 3          // be polite to NCBI; matches the pipeline's other fetches

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/jhuaplbio/taxtriage_confidence:2.1' :
        'jhuaplbio/taxtriage_confidence:2.1' }"

    input:
    tuple val(accession), path(assembly_summary)

    output:
    tuple val(accession), path("refs/*.fasta"), emit: reference
    // accession -> taxid + organism name. The report needs the TAXID to tie a
    // spiked organism to a detection; matching on the sheet's optional free-text
    // name silently fails whenever that column is blank.
    path("taxids/*.taxid.tsv")                , emit: taxid
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def safe = accession.replaceAll(/[^A-Za-z0-9._-]/, '_')
    def api_key = params.ncbi_api_key ? "-api_key ${params.ncbi_api_key}" : ""
    """
    set -o pipefail
    mkdir -p refs taxids
    OUT=refs/${safe}.fasta
    TAXOUT=taxids/${safe}.taxid.tsv

    # (0) A path the user handed us directly.
    if [ -f "${accession}" ]; then
        case "${accession}" in
            *.gz) zcat "${accession}" > \$OUT ;;
            *)    cp "${accession}" \$OUT ;;
        esac
    fi

    # (1) assembly_summary -> FTP directory -> <basename>_genomic.fna.gz
    if [ ! -s \$OUT ] && [ -s "${assembly_summary}" ] && [ "${assembly_summary}" != "NO_FILE" ]; then
        FTP=\$(awk -F'\\t' -v acc="${accession}" '\$1 == acc {print \$20; exit}' ${assembly_summary} || true)
            # assembly_summary columns: 6 = taxid, 7 = species_taxid, 8 = organism_name
            awk -F'\\t' -v acc="${accession}" '\$1 == acc {printf "%s\\t%s\\t%s\\n", acc, \$6, \$8; exit}' \\
                ${assembly_summary} > \$TAXOUT || true
        if [ -n "\$FTP" ] && [ "\$FTP" != "na" ]; then
            BASE=\$(basename "\$FTP")
            HTTP=\$(echo "\$FTP" | sed 's|^ftp://|https://|')
            echo "[spikein-refs] ${accession}: assembly_summary -> \$HTTP" >&2
            curl -sSL --retry 3 --retry-delay 2 "\$HTTP/\${BASE}_genomic.fna.gz" -o ref.fna.gz \\
                && zcat ref.fna.gz > \$OUT || true
        fi
    fi

    # (2) NCBI datasets CLI, when present.
    if [ ! -s \$OUT ] && command -v datasets >/dev/null 2>&1; then
        case "${accession}" in
            GCF_*|GCA_*)
                echo "[spikein-refs] ${accession}: trying datasets CLI" >&2
                datasets download genome accession ${accession} --include genome --filename ds.zip >/dev/null 2>&1 \\
                    && unzip -o -q ds.zip -d ds >/dev/null 2>&1 \\
                    && cat ds/ncbi_dataset/data/*/*.fna > \$OUT 2>/dev/null || true
                ;;
        esac
    fi

    # (3) Entrez efetch — the route for a bare nuccore accession.
    if [ ! -s \$OUT ]; then
        echo "[spikein-refs] ${accession}: trying Entrez efetch" >&2
        if command -v efetch >/dev/null 2>&1; then
            efetch -db nuccore -id "${accession}" -format fasta ${api_key} > \$OUT 2>/dev/null || true
        else
            KEY=""
            if [ -n "${params.ncbi_api_key ?: ''}" ]; then KEY="&api_key=${params.ncbi_api_key ?: ''}"; fi
            curl -sSL --retry 3 --retry-delay 2 \\
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${accession}&rettype=fasta&retmode=text\${KEY}" \\
                -o \$OUT || true
        fi
    fi

    # Fall back to Entrez for the taxid when assembly_summary did not supply it.
    if [ ! -s \$TAXOUT ]; then
        TID=""
        ORG=""
        if command -v esearch >/dev/null 2>&1; then
            TID=\$(esearch -db nuccore -query "${accession}" 2>/dev/null \\
                   | esummary 2>/dev/null \\
                   | xtract -pattern DocumentSummary -element TaxId 2>/dev/null | head -1 || true)
            ORG=\$(esearch -db nuccore -query "${accession}" 2>/dev/null \\
                   | esummary 2>/dev/null \\
                   | xtract -pattern DocumentSummary -element Organism 2>/dev/null | head -1 || true)
        fi
        if [ -z "\$TID" ]; then
            TID=\$(curl -sSL --retry 2 \\
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&id=${accession}&retmode=json" 2>/dev/null \\
                | tr ',' '\\n' | grep -m1 '"taxid"' | tr -dc '0-9' || true)
        fi
        if [ -z "\$ORG" ]; then
            # Last resort: the description on the FASTA header we just fetched.
            ORG=\$(head -1 \$OUT 2>/dev/null | sed 's/^>[^ ]* //' | cut -d',' -f1 || true)
        fi
        printf '%s\\t%s\\t%s\\n' "${accession}" "\${TID:-}" "\${ORG:-}" > \$TAXOUT
    fi
    echo "[spikein-refs] ${accession}: taxid/org -> \$(cat \$TAXOUT)" >&2

    if [ ! -s \$OUT ] || ! grep -q '^>' \$OUT; then
        echo "ERROR: could not fetch a reference for spike-in accession '${accession}'." >&2
        echo "       Tried the pipeline's assembly_summary, the NCBI datasets CLI and Entrez." >&2
        echo "       Check the accession, or give a local FASTA path in that column instead." >&2
        exit 1
    fi

    echo "[spikein-refs] ${accession}: \$(grep -c '^>' \$OUT) sequence(s), \$(grep -v '^>' \$OUT | tr -d '\\n' | wc -c) bp" >&2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        curl: \$(curl --version | head -1 | sed 's/curl //' | cut -d' ' -f1)
    END_VERSIONS
    """
}
