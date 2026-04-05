process CORRECTUMI {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*_namesorted.bam"), emit: bam
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    set -euo pipefail

    # Name-sort so mates are adjacent for deduplicate_bismark -p / --barcode
    samtools sort -n -@ ${task.cpus} -o ${prefix}_namesorted.bam ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | sed -n '1p' | sed 's/samtools //')
    END_VERSIONS
    """
}
