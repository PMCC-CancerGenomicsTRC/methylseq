process FASTQ_UMI_TO_BARCODE {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_umiheader_*.fastq.gz"), emit: reads
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
             # Split header into first token and the rest (keeps " 1:N:0:..." unchanged)
             split(\$0, a, " ")
             h=a[1]
             rest=""
             if (length(\$0) > length(h)) rest=substr(\$0, length(h)+1)

             # Replace ":UMI_<SEQ>" with ":<SEQ>"
             sub(/:UMI_([A-Za-z]+)/, ":\\\\1", h)

             # Uppercase the barcode and optionally restrict to ACGTN
             if (match(h, /:([A-Za-z]+)$/, m)) {
               b=toupper(m[1])
               gsub(/[^ACGTN]/, "N", b)
               sub(/:[A-Za-z]+$/, ":" b, h)
             }

             print h rest
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
