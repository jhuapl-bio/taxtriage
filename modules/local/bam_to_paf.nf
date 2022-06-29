process BAM_TO_PAF {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconvert:0.5.2--pyhdfd78af_0' :
        'lbmc/bioconvert:0.4.0' }"

    input:
    tuple val(meta), path(bamfiles)

    output:
    path "versions.yml"           , emit: versions
    tuple val(meta), path("*.paf")     , optional:false, emit: paf



    when:
    task.ext.when == null || task.ext.when


    script: // This script is bundled with the pipeline, in nf-core/taxtriage/bin/
    """
    bioconvert bam2sam ${bamfiles} ${meta.id}.sam -a --force; 
    bioconvert sam2paf ${meta.id}.sam ${meta.id}.paf -a --force


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version)
    END_VERSIONS

    """
}


