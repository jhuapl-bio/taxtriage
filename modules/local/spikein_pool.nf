// Simulate ONE read pool per spike-in organism, from that organism's own
// reference. The pool is deliberately larger than the largest requested spike
// (params.spikein_pool_factor) so replicates can draw different reads from it
// without re-running the simulator per replicate — which is both slow and would
// make the draws harder to reproduce.
//
// Two processes rather than one with a branch, so each can use the SAME public
// biocontainer the pipeline's existing simulator modules use (insilicoseq.nf /
// nanosim.nf). Exact per-dataset counts are taken from these pools later by
// spike_into_background.py.
//
// n_reads travels inside the input tuple: each accession needs a pool sized to
// its own largest request, and a separate val channel would be zipped
// positionally (so every organism would silently get the first one's size).

process SPIKEIN_POOL_ISS {
    tag "$accession"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::insilicoseq=2.0.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/insilicoseq:2.0.1--pyh7cba7a3_0' :
        'biocontainers/insilicoseq:2.0.1--pyh7cba7a3_0' }"

    input:
    tuple val(accession), path(reference), val(n_reads)
    val(iss_model)

    output:
    tuple val(accession), path("pool/*.fastq.gz"), emit: pool
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def safe  = accession.replaceAll(/[^A-Za-z0-9._-]/, '_')
    def model = iss_model ?: 'miseq'
    def mode  = params.iss_mode ?: 'kde'
    """
    set -o pipefail
    mkdir -p pool

    # ISS keys abundance on the FASTA record id, so EVERY contig of a multi-contig
    # reference needs a row — listing only the first leaves the rest unsimulated.
    grep '^>' ${reference} | sed 's/^>//' | cut -d' ' -f1 > ids.txt
    N=\$(wc -l < ids.txt)
    if [ "\$N" -eq 0 ]; then
        echo "ERROR: ${reference} contains no FASTA records for ${accession}" >&2
        exit 1
    fi
    awk -v n="\$N" '{ printf "%s\\t%.10f\\n", \$1, 1.0/n }' ids.txt > abundance.tsv

    iss generate \\
        --genomes ${reference} \\
        --model ${model} \\
        --mode ${mode} \\
        --abundance_file abundance.tsv \\
        --output ${safe}.pool \\
        -n ${n_reads} \\
        --cpus ${task.cpus}

    for f in ${safe}.pool_R1.fastq ${safe}.pool_R2.fastq; do
        [ -f "\$f" ] || { echo "ERROR: InSilicoSeq did not produce \$f" >&2; exit 1; }
    done
    mv ${safe}.pool_R1.fastq pool/${safe}_R1.fastq
    mv ${safe}.pool_R2.fastq pool/${safe}_R2.fastq
    gzip pool/${safe}_R1.fastq pool/${safe}_R2.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        insilicoseq: \$(iss --version 2>&1 | sed 's/iss version //')
    END_VERSIONS
    """
}


process SPIKEIN_POOL_NANOSIM {
    tag "$accession"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::nanosim=3.2.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanosim:3.2.3--hdfd78af_0' :
        'biocontainers/nanosim:3.2.3--hdfd78af_2' }"

    input:
    tuple val(accession), path(reference), val(n_reads)
    path(nanosim_training)

    output:
    tuple val(accession), path("pool/*.fastq.gz"), emit: pool
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def safe   = accession.replaceAll(/[^A-Za-z0-9._-]/, '_')
    def nsbase = params.nanosim_base ?: 'training'
    """
    set -o pipefail
    mkdir -p pool

    # NanoSim's metagenome mode wants a genome list and an abundance file even for
    # a single genome.
    printf '${safe}\\t${reference}\\n' > genome_list.tsv
    printf '${safe}\\t100.0\\n' > abundance.tsv

    simulator.py metagenome \\
        -gl genome_list.tsv \\
        -a abundance.tsv \\
        -c ${nanosim_training}/${nsbase} \\
        -o ${safe}.pool \\
        --fastq \\
        --perfect \\
        -t ${task.cpus} \\
        -n ${n_reads}

    cat ${safe}.pool*.fastq > pool/${safe}.fastq
    gzip pool/${safe}.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanosim: \$(simulator.py --version 2>&1 | tail -1)
    END_VERSIONS
    """
}
