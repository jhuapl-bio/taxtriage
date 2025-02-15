process COUNT_READS {
    // Use a closure for the tag so it works whether meta is a map or string
    tag { meta?.id ? meta.id : meta }

    label 'process_low'

    // Set up your container/conda as needed
    conda ( params.enable_conda ? "conda-forge::gawk" : null )
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'https://depot.galaxyproject.org/singularity/gawk:4.2.0' : 'biocontainers/gawk:4.2.0' }"

    input:
        // Expect a tuple: [ meta, reads ]
        tuple val(meta), path(reads)

    output:
        // Output the original metadata together with the count file
        tuple val(meta), path("count.txt"), path(reads), emit: count

    script:
        // Loop over each file in the 'reads' variable (works for one or multiple files)
        """
        total=0
        for f in $reads; do
            if [[ "\$f" =~ \\.gz\$ ]]; then
                cnt=\$(zcat "\$f" | awk 'NR % 4 == 1 {count++} END {print count}')
            else
                cnt=\$(awk 'NR % 4 == 1 {count++} END {print count}' "\$f")
            fi
            total=\$(( total + cnt ))
        done
        echo \$total > count.txt
        """
}
