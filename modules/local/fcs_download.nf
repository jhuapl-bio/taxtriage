// ##############################################################################################
// # Copyright 2022 The Johns Hopkins University Applied Physics Laboratory LLC
// # All rights reserved.
// # Permission is hereby granted, free of charge, to any person obtaining a copy of this
// # software and associated documentation files (the "Software"), to deal in the Software
// # without restriction, including without limitation the rights to use, copy, modify,
// # merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
// # permit persons to whom the Software is furnished to do so.
// #
// # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
// # INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
// # PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// # LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
// # TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
// # OR OTHER DEALINGS IN THE SOFTWARE.
// #

/**
 * FCS_DOWNLOAD
 *
 * Fetches NCBI's precomputed Foreign Contamination Screen (FCS-GX) reports from
 * https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/ :
 *   fcs_details_{genbank,refseq}.txt.gz  per-sequence contaminant spans (masked)
 *   fcs_summary_{genbank,refseq}.txt.gz  per-assembly contamination stats
 *
 * storeDir caches them (--fcs_db, default <outdir>/fcs_db) so later runs skip
 * the download. Only used under --fcs when no local --fcs_details is given.
 * Files stay gzipped: build_mask.py streams .gz directly.
 */
process FCS_DOWNLOAD {
    tag "ncbi_fcs"
    label 'process_low'
    storeDir "${ params.fcs_db ?: params.outdir + '/fcs_db' }"

    conda (params.enable_conda ? "conda-forge::curl" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/jhuaplbio/taxtriage_confidence:2.2' :
        'jhuaplbio/taxtriage_confidence:2.2' }"

    output:
    path "fcs_details_*.txt.gz" , emit: details
    path "fcs_summary_*.txt.gz" , emit: summary

    script:
    def base = params.fcs_url ?: 'https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS'
    def sets = (params.fcs_sets ?: 'genbank,refseq').tokenize(',')*.trim()
    def cmds = sets.collect { s ->
        "fetch ${base}/fcs_details_${s}.txt.gz fcs_details_${s}.txt.gz\n" +
        "fetch ${base}/fcs_summary_${s}.txt.gz fcs_summary_${s}.txt.gz"
    }.join('\n')
    """
    fetch() {
        if command -v curl >/dev/null 2>&1; then
            curl -fsSL --retry 3 -o "\$2" "\$1"
        else
            wget -q -O "\$2" "\$1"
        fi
        gzip -t "\$2"
    }
    ${cmds}
    """
}
