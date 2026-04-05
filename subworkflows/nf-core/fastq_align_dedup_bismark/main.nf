include { BISMARK_ALIGN                } from '../../../modules/nf-core/bismark/align/main'
include { CORRECTUMI                   } from '../../../modules/local/correctumi'
include { BISMARK_DEDUPLICATE          } from '../../../modules/nf-core/bismark/deduplicate/main'
include { SAMTOOLS_SORT                } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX               } from '../../../modules/nf-core/samtools/index/main'
include { BISMARK_METHYLATIONEXTRACTOR } from '../../../modules/nf-core/bismark/methylationextractor/main'
include { BISMARK_COVERAGE2CYTOSINE    } from '../../../modules/nf-core/bismark/coverage2cytosine/main'
include { BISMARK_REPORT               } from '../../../modules/nf-core/bismark/report/main'
include { BISMARK_SUMMARY              } from '../../../modules/nf-core/bismark/summary/main'

workflow FASTQ_ALIGN_DEDUP_BISMARK {

    take:
    ch_reads
    ch_fasta
    ch_bismark_index
    skip_deduplication
    cytosine_report

    main:
    ch_alignments                 = Channel.empty()
    ch_alignment_reports          = Channel.empty()
    ch_bam_final                  = Channel.empty()
    ch_bam_for_extractor          = Channel.empty()
    ch_methylation_bedgraph       = Channel.empty()
    ch_methylation_calls          = Channel.empty()
    ch_methylation_coverage       = Channel.empty()
    ch_methylation_report         = Channel.empty()
    ch_methylation_mbias          = Channel.empty()
    ch_coverage2cytosine_coverage = Channel.empty()
    ch_coverage2cytosine_report   = Channel.empty()
    ch_coverage2cytosine_summary  = Channel.empty()
    ch_bismark_report             = Channel.empty()
    ch_bismark_summary            = Channel.empty()
    ch_multiqc_files              = Channel.empty()
    ch_versions                   = Channel.empty()

    /*
     * Align with bismark
     */
    BISMARK_ALIGN (
        ch_reads,
        ch_fasta,
        ch_bismark_index
    )
    ch_alignments        = BISMARK_ALIGN.out.bam
    ch_alignment_reports = BISMARK_ALIGN.out.report.map{ meta, report -> [ meta, report, [] ] }
    ch_versions          = ch_versions.mix(BISMARK_ALIGN.out.versions)

    if (!skip_deduplication) {

        /*
         * Name-sort BAM so paired reads are adjacent for deduplicate_bismark -p/--barcode
         */
        CORRECTUMI (
            ch_alignments
        )
        ch_versions = ch_versions.mix(CORRECTUMI.out.versions)

        BISMARK_DEDUPLICATE (
            CORRECTUMI.out.bam
        )
        ch_bam_final = BISMARK_DEDUPLICATE.out.bam
        ch_versions  = ch_versions.mix(BISMARK_DEDUPLICATE.out.versions)

    } else {
        ch_bam_final = ch_alignments
    }

    /*
     * Name-sort (again) for methylation extractor.
     * Methylation extractor in paired-end mode expects read pairs adjacent (queryname sort).
     */
    CORRECTUMI (
        ch_bam_final
    )
    ch_bam_for_extractor = CORRECTUMI.out.bam
    ch_versions          = ch_versions.mix(CORRECTUMI.out.versions)

    /*
     * Coordinate-sort for BAM indexing + downstream tools
     */
    SAMTOOLS_SORT (
        ch_bam_final,
        [[:],[]]
    )
    ch_bam_final = SAMTOOLS_SORT.out.bam
    ch_versions  = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    /*
     * Index BAM (coordinate-sorted)
     */
    SAMTOOLS_INDEX ( ch_bam_final )
    ch_bai      = SAMTOOLS_INDEX.out.bai
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    /*
     * Methylation extractor (use NAME-SORTED BAM)
     */
    BISMARK_METHYLATIONEXTRACTOR ( ch_bam_for_extractor, ch_bismark_index )
    ch_methylation_bedgraph = BISMARK_METHYLATIONEXTRACTOR.out.bedgraph
    ch_methylation_calls    = BISMARK_METHYLATIONEXTRACTOR.out.methylation_calls
    ch_methylation_coverage = BISMARK_METHYLATIONEXTRACTOR.out.coverage
    ch_methylation_report   = BISMARK_METHYLATIONEXTRACTOR.out.report
    ch_methylation_mbias    = BISMARK_METHYLATIONEXTRACTOR.out.mbias
    ch_versions             = ch_versions.mix(BISMARK_METHYLATIONEXTRACTOR.out.versions)

    /*
     * coverage2cytosine
     */
    if (cytosine_report) {
        BISMARK_COVERAGE2CYTOSINE ( ch_methylation_coverage, ch_fasta, ch_bismark_index )
        ch_coverage2cytosine_coverage = BISMARK_COVERAGE2CYTOSINE.out.coverage
        ch_coverage2cytosine_report   = BISMARK_COVERAGE2CYTOSINE.out.report
        ch_coverage2cytosine_summary  = BISMARK_COVERAGE2CYTOSINE.out.summary
        ch_versions                   = ch_versions.mix(BISMARK_COVERAGE2CYTOSINE.out.versions)
    }

    /*
     * Reports + summary
     */
    BISMARK_REPORT (
        ch_alignment_reports
            .join(ch_methylation_report)
            .join(ch_methylation_mbias)
    )
    ch_bismark_report = BISMARK_REPORT.out.report
    ch_versions       = ch_versions.mix(BISMARK_REPORT.out.versions)

    BISMARK_SUMMARY (
        BISMARK_ALIGN.out.bam.collect{ meta, bam -> bam.name },
        ch_alignment_reports.collect{ meta, align_report, dedup_report -> align_report },
        ch_alignment_reports.collect{ meta, align_report, dedup_report -> dedup_report }.ifEmpty([]),
        ch_methylation_report.collect{ meta, report -> report },
        ch_methylation_mbias.collect{ meta, mbias -> mbias }
    )
    ch_bismark_summary = BISMARK_SUMMARY.out.summary
    ch_versions        = ch_versions.mix(BISMARK_SUMMARY.out.versions)

    /*
     * MultiQC inputs
     */
    ch_multiqc_files = ch_bismark_summary
        .mix(ch_alignment_reports.collect{ meta, align_report, dedup_report -> align_report })
        .mix(ch_alignment_reports.collect{ meta, align_report, dedup_report -> dedup_report })
        .mix(ch_methylation_report.collect{ meta, report -> report })
        .mix(ch_methylation_mbias.collect{ meta, mbias -> mbias })
        .mix(ch_bismark_report.collect{ meta, report -> report })

    emit:
    bam      = ch_bam_final
    bai      = ch_bai
    multiqc  = ch_multiqc_files
    versions = ch_versions
}
