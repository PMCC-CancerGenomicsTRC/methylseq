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
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
   
    """
    samtools \\
        view \\
        -@ $task.cpus \\
        $bam | awk -F"\t" \
        '{ split($1,name, ":"); split(name[8],umi,"_"); print name[1]":"name[2]":"name[3]":"name[4]":"name[5]":"name[6]":"name[7]":"umi[2]":"name[9]":"name[10]":"name[11]":"umi[1]"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15" }' | samtools view -b -p ${prefix}_umiCorrect.bam


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        correctumi: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        correctumi: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
