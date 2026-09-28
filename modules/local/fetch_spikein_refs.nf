// Fetch the reference FASTA for one spike-in accession.
//
// Accepts an assembly accession (GCF_/GCA_) or a nuccore accession, and resolves
// it in the order that is cheapest and most reproducible first:
//   1. a row in any assembly_summary given (RefSeq, then GenBank) -> direct FTP
//   2. the NCBI `datasets` CLI, when the image has it
//   3. Entrez efetch (the only route for a bare nuccore accession)
// A local FASTA given in the sheet arrives staged as `local_ref` (NO_FILE when the
// row is an accession), so a user can spike in a genome with no accession at all.
// Its taxid comes from the sheet's `taxid` column, else from NCBI via the record
// accessions in its headers (resolve_fasta_taxids.py) — which also refuses a FASTA
// that holds several organisms, since each row is scored against ONE taxid.
process FETCH_SPIKEIN_REFS {
    tag "$accession"
    label 'process_single'
    maxForks 3          // be polite to NCBI; matches the pipeline's other fetches

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/jhuaplbio/taxtriage_confidence:2.1' :
        'jhuaplbio/taxtriage_confidence:2.1' }"

    input:
    // One or more assembly_summary tables (RefSeq, GenBank). Staged into numbered
    // dirs so two files with the same basename cannot collide.
    // accession = the spike-in id; fetch_acc = what to download when that differs
    // (a `record`-narrowed accession row); records = FASTA record ids to keep ('' = all)
    tuple val(accession), val(sheet_taxid), val(records), val(fetch_acc), path(local_ref, stageAs: 'local_ref/*'), path(assembly_summaries, stageAs: 'asm_summary_??/*')
    // --custom_accession_map (or NO_FILE): the same accession->taxid map the main
    // pipeline's MAP_TAXID_ASSEMBLY uses, so a spiked organism gets the SAME taxid
    // the report assigns to its detection.
    path(custom_map, stageAs: 'custom_map/*')

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
    def is_local = local_ref.name != 'NO_FILE'
    def query    = fetch_acc ?: accession
    def extract  = records ? "extract_fasta_records.py --fasta refs/full.fasta --records '${records}' --output \$OUT" : ''
    def cmap_arg = !custom_map.name.startsWith('NO_FILE') ? "--custom-map ${custom_map}" : ''
    """
    set -o pipefail
    mkdir -p refs taxids
    OUT=refs/${safe}.fasta
    TAXOUT=taxids/${safe}.taxid.tsv

    # (0) A local FASTA from the sheet: use it as-is, and resolve its taxid here
    #     (the Entrez fallbacks below key on the accession, which for a local file
    #     is just an id made from the file name).
    if [ "${is_local}" = "true" ]; then
        case "${local_ref}" in
            *.gz) zcat "${local_ref}" > \$OUT ;;
            *)    cp -L "${local_ref}" \$OUT ;;
        esac
        if [ -n "${records}" ]; then
            mv \$OUT refs/full.fasta
            ${extract}
            rm -f refs/full.fasta
        fi
        NCBI_API_KEY="${params.ncbi_api_key ?: ''}" resolve_fasta_taxids.py \\
            --fasta \$OUT \\
            --accession "${accession}" \\
            --taxid "${sheet_taxid ?: ''}" \\
            ${cmap_arg} \\
            --output \$TAXOUT
    fi

    # (1) assembly_summary (RefSeq, then GenBank) -> FTP directory -> <basename>_genomic.fna.gz
    #     GCA_ accessions only appear in the GenBank summary, so scan every table.
    for SUMMARY in ${assembly_summaries}; do
        [ -s \$OUT ] && break
        [ -s "\$SUMMARY" ] || continue
        [ "\$(basename \$SUMMARY)" = "NO_FILE" ] && continue
        FTP=\$(awk -F'\\t' -v acc="${query}" '\$1 == acc {print \$20; exit}' "\$SUMMARY" || true)
        [ -n "\$FTP" ] || continue
        echo "[spikein-refs] ${query}: found in \$(basename \$SUMMARY)" >&2
        # assembly_summary columns: 6 = taxid, 7 = species_taxid, 8 = organism_name
        awk -F'\\t' -v acc="${query}" -v id="${accession}" '\$1 == acc {printf "%s\\t%s\\t%s\\n", id, \$6, \$8; exit}' \\
            "\$SUMMARY" > \$TAXOUT || true
        if [ "\$FTP" != "na" ]; then
            BASE=\$(basename "\$FTP")
            HTTP=\$(echo "\$FTP" | sed 's|^ftp://|https://|')
            echo "[spikein-refs] ${query}: assembly_summary -> \$HTTP" >&2
            curl -sSL --retry 3 --retry-delay 2 "\$HTTP/\${BASE}_genomic.fna.gz" -o ref.fna.gz \\
                && zcat ref.fna.gz > \$OUT || true
        fi
    done

    # (2) NCBI datasets CLI, when present.
    if [ ! -s \$OUT ] && command -v datasets >/dev/null 2>&1; then
        case "${query}" in
            GCF_*|GCA_*)
                echo "[spikein-refs] ${query}: trying datasets CLI" >&2
                datasets download genome accession ${query} --include genome --filename ds.zip >/dev/null 2>&1 \\
                    && unzip -o -q ds.zip -d ds >/dev/null 2>&1 \\
                    && cat ds/ncbi_dataset/data/*/*.fna > \$OUT 2>/dev/null || true
                ;;
        esac
    fi

    # (3) Entrez efetch — the route for a bare nuccore accession.
    if [ ! -s \$OUT ]; then
        echo "[spikein-refs] ${query}: trying Entrez efetch" >&2
        if command -v efetch >/dev/null 2>&1; then
            efetch -db nuccore -id "${query}" -format fasta ${api_key} > \$OUT 2>/dev/null || true
        else
            KEY=""
            if [ -n "${params.ncbi_api_key ?: ''}" ]; then KEY="&api_key=${params.ncbi_api_key ?: ''}"; fi
            curl -sSL --retry 3 --retry-delay 2 \\
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${query}&rettype=fasta&retmode=text\${KEY}" \\
                -o \$OUT || true
        fi
    fi

    # A `record`-narrowed accession: keep just those records of what was fetched.
    if [ "${is_local}" != "true" ] && [ -n "${records}" ] && [ -s \$OUT ]; then
        mv \$OUT refs/full.fasta
        ${extract}
        rm -f refs/full.fasta
    fi

    # Fall back to Entrez for the taxid when assembly_summary did not supply it.
    if [ ! -s \$TAXOUT ]; then
        TID=""
        ORG=""
        if command -v esearch >/dev/null 2>&1; then
            TID=\$(esearch -db nuccore -query "${query}" 2>/dev/null \\
                   | esummary 2>/dev/null \\
                   | xtract -pattern DocumentSummary -element TaxId 2>/dev/null | head -1 || true)
            ORG=\$(esearch -db nuccore -query "${query}" 2>/dev/null \\
                   | esummary 2>/dev/null \\
                   | xtract -pattern DocumentSummary -element Organism 2>/dev/null | head -1 || true)
        fi
        if [ -z "\$TID" ]; then
            TID=\$(curl -sSL --retry 2 \\
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&id=${query}&retmode=json" 2>/dev/null \\
                | tr ',' '\\n' | grep -m1 '"taxid"' | tr -dc '0-9' || true)
        fi
        if [ -z "\$ORG" ]; then
            # Last resort: the description on the FASTA header we just fetched.
            ORG=\$(head -1 \$OUT 2>/dev/null | sed 's/^>[^ ]* //' | cut -d',' -f1 || true)
        fi
        printf '%s\\t%s\\t%s\\n' "${accession}" "\${TID:-}" "\${ORG:-}" > \$TAXOUT
    fi
    # --custom_accession_map overrides NCBI for accession rows too (matched on the
    # sheet accession or the fetched record ids). Local rows already used it above.
    if [ "${is_local}" != "true" ] && [ -n "${cmap_arg}" ] && [ -s \$OUT ]; then
        resolve_fasta_taxids.py --fasta \$OUT --accession "${accession}" \\
            ${cmap_arg} --no-ncbi --only-if-mapped --output \$TAXOUT
    fi

    # A taxid in the sheet always wins over what NCBI returned.
    if [ -n "${sheet_taxid ?: ''}" ]; then
        awk -F'\\t' -v OFS='\\t' -v t="${sheet_taxid}" '{ print \$1, t, \$3 }' \$TAXOUT > tx.tmp && mv tx.tmp \$TAXOUT
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
