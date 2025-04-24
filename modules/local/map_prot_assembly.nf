process MAP_PROT_ASSEMBLY {
    label 'process_medium'
    tag "$meta.id"

    conda (params.enable_conda ? "conda-forge::python=3.8.3 pysam" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.21.0--py39hcada746_1' :
        'biocontainers/pysam:0.21.0--py39hcada746_1' }"

    input:
    tuple val(meta), file(diamondoutput), file(features), file(mapping), file(map_names)
    file(assembly)

    output:
    tuple val(meta), path("*assembly_protein_map.txt"), optional:false, emit: promap

    script:

    def mapnames = map_names ? "--mapnames ${map_names}" : ""
    """

    map_prot_assembly.py \\
        -d $diamondoutput \\
        -f $features \\
        -a $assembly \\
        -n $mapping $mapnames \\
        -o ${meta.id}.assembly_protein_map.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python3: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
