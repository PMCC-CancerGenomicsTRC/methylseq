process FASTQ_UMI_TO_BARCODE {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_umiheader_{1,2}.fastq.gz"), emit: reads
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    set -euo pipefail

    fix() {
      local in_fq="\$1"
      local out_fq="\$2"

      zcat "\$in_fq" | \\
      awk 'NR%4==1 {
             # Header line. Split into first token and the rest (keeps " 1:N:0:..." unchanged)
             split(\$0, a, " ")
             h=a[1]
             rest=""
             if (length(\$0) > length(h)) rest=substr(\$0, length(h)+1)

             # Convert trailing ":UMI_<SEQ>" into ":<SEQ>"
             # Example:
             #   @VH...:1019:UMI_CCCT... -> @VH...:1019:CCCT...
             sub(/:UMI_/, ":", h)

             print h rest
             next
           }
           { print }' \\
      | gzip -c > "\$out_fq"
    }

    # reads is a list: [R1, R2]
    fix "${reads[0]}" "${prefix}_umiheader_1.fastq.gz"
    fix "${reads[1]}" "${prefix}_umiheader_2.fastq.gz"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: "builtin"
    END_VERSIONS
    """
}
