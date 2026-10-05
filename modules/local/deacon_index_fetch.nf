// Fetch a prebuilt Deacon minimizer index by name (e.g. panhuman-1, panmouse-1).
//
// This downloads the file directly with curl/wget rather than running
// `deacon index fetch`: the deacon biocontainer ships no CA certificates, so
// its built-in (rustls) downloader fails with "No CA certificates were loaded
// from the system".  The URL is the one `deacon index fetch` uses:
//   <deacon_index_url>/<name>.k<k>w<w>.idx
// where deacon_index_url ends in the index-format version ("/deacon/3" for
// deacon 0.13-0.15).  --deacon_insecure_download skips TLS certificate checks
// (curl -k / wget --no-check-certificate) for hosts with broken CA setups.
//
// The result is written through `storeDir`, so it is downloaded once and
// re-used by every later run pointing at the same --deacon_index_dir.
process DEACON_INDEX_FETCH {
    tag "$name"
    label 'process_single'
    maxForks 1
    storeDir params.deacon_index_dir ?: "${params.host_reference_dir ?: "${params.outdir}/host_references"}/deacon"

    // Same image FETCH_HOST_REFS uses for its NCBI downloads (curl + CA certs).
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/jhuaplbio/taxtriage_confidence:2.1' :
        'jhuaplbio/taxtriage_confidence:2.1' }"

    input:
    val(name)

    output:
    path("${name}.k31w15.idx"), emit: index

    when:
    task.ext.when == null || task.ext.when

    script:
    def base     = params.deacon_index_url.toString().replaceAll('/+$', '')
    def fname    = "${name}.k31w15.idx"
    def url      = "${base}/${fname}"
    def insecure = params.deacon_insecure_download
    def curl_k   = insecure ? '-k' : ''
    def wget_k   = insecure ? '--no-check-certificate' : ''
    """
    set -o pipefail
    echo "[deacon-fetch] ${url}${insecure ? ' (TLS certificate checks disabled)' : ''}" >&2

    if command -v curl >/dev/null 2>&1; then
        curl -fSL ${curl_k} --retry 5 --retry-delay 10 -C - -o ${fname}.tmp "${url}"
    elif command -v wget >/dev/null 2>&1; then
        wget ${wget_k} --tries=5 -c -O ${fname}.tmp "${url}"
    else
        echo "ERROR: neither curl nor wget is available to download ${url}" >&2
        exit 1
    fi

    # Guard against an HTML/XML error page or a truncated download (the
    # prebuilt indexes are several GB).
    SIZE=\$(wc -c < ${fname}.tmp)
    if [ "\$SIZE" -lt 1000000 ]; then
        echo "ERROR: ${url} returned only \$SIZE bytes - not a Deacon index:" >&2
        head -c 500 ${fname}.tmp >&2
        echo "" >&2
        echo "       Download it manually and pass --deacon_index /path/to/${fname}" >&2
        exit 1
    fi
    mv ${fname}.tmp ${fname}
    echo "[deacon-fetch] saved ${fname} (\$SIZE bytes)" >&2
    """

    stub:
    """
    touch ${name}.k31w15.idx
    """
}
