process COUNT_READS {
    tag { meta instanceof Map ? meta.id : meta }
    label 'process_low'

    conda ( params.enable_conda ? "conda-forge::gzip conda-forge::gawk" : null )
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:4.2.0' :
        'biocontainers/gawk:4.2.0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_count.txt"), emit: count

    script:
    def sample_id = meta instanceof Map ? meta.id : meta
    def read_list = reads instanceof List ? reads : [reads]
    def quoted_reads = read_list.collect { "\"${it}\"" }.join(' ')

    """
    set -euo pipefail

    total=0

    for f in ${quoted_reads}; do
        [[ -s "\$f" ]] || { echo "Missing or empty file: \$f" >&2; exit 1; }

        if [[ "\$f" == *.gz ]]; then
            cnt=\$(gzip -cd "\$f" | awk 'NR % 4 == 1 {count++} END {print count+0}')
        else
            cnt=\$(awk 'NR % 4 == 1 {count++} END {print count+0}' "\$f")
        fi

        total=\$(( total + cnt ))
    done

    echo "\$total" > "${sample_id}_count.txt"
    """
}
