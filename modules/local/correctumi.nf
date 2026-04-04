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
    tuple val(meta), path("*umiCorrect.bam"), emit: bam
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    samtools view -h ${bam} | \\
    awk 'BEGIN{FS=OFS="\\t"}
         /^@/ { print; next }
         {
           # Extract UMI from a field like UMI_<seq>_1:N:0:<index>
           if (match(\$1, /:UMI_([^:_]+)/, m)) {
               # Keep everything before :UMI_ and append :<UMI>
               sub(/:UMI_.*/, ":" m[1], \$1)
           }
           print
         }' | \\
    samtools view -b -o ${prefix}_umiCorrect.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | sed -n '1p' | sed 's/samtools //')
    END_VERSIONS
    """
}
