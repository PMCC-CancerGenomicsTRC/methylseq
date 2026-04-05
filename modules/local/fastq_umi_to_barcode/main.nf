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
             # Header line. Keep only the first token (up to first space).
             # This prevents aligners from converting the trailing " 1:N:0:..." into "_1:N:0:..." in BAM QNAME,
             # which breaks deduplicate_bismark --barcode parsing.
             split(\$0, a, " ")
             h=a[1]

             # Convert ":UMI_<SEQ>" into ":<SEQ>"
             sub(/:UMI_/, ":", h)

             print h
             next
           }
           { print }' \\
      | gzip -c > "\$out_fq"
    }

    fix "${reads[0]}" "${prefix}_umiheader_1.fastq.gz"
    fix "${reads[1]}" "${prefix}_umiheader_2.fastq.gz"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: "builtin"
    END_VERSIONS
    """
}
