#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/methylseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/methylseq
    Website: https://nf-co.re/methylseq
    Slack  : https://nfcore.slack.com/channels/methylseq
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTA_INDEX_BISMARK_BWAMETH } from './subworkflows/nf-core/fasta_index_bismark_bwameth/main'
include { PIPELINE_INITIALISATION     } from './subworkflows/local/utils_nfcore_methylseq_pipeline'
include { PIPELINE_COMPLETION         } from './subworkflows/local/utils_nfcore_methylseq_pipeline'
include { getGenomeAttribute          } from './subworkflows/local/utils_nfcore_methylseq_pipeline'
include { METHYLSEQ                   } from './workflows/methylseq/'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
params.fasta         = getGenomeAttribute('fasta')
params.fasta_index   = getGenomeAttribute('fasta_index')
params.bwameth_index = getGenomeAttribute('bwameth')
params.bismark_index = params.aligner == 'bismark_hisat' ? getGenomeAttribute('bismark_hisat2') : getGenomeAttribute('bismark')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INITIALISE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
ch_primer_fasta3 = params.primer_fasta3 ? Channel.fromPath(params.primer_fasta3).first(): Channel.value([])
ch_primer_fasta5 = params.primer_fasta5 ? Channel.fromPath(params.primer_fasta5).first(): Channel.value([])

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_METHYLSEQ {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Prepare any required reference genome indices
    //
    FASTA_INDEX_BISMARK_BWAMETH(
        params.fasta,
        params.fasta_index,
        params.bismark_index,
        params.bwameth_index,
    )
    ch_versions = ch_versions.mix(FASTA_INDEX_BISMARK_BWAMETH.out.versions)

    //
    // WORKFLOW: Run pipeline
    //

    METHYLSEQ (
        samplesheet,
        ch_versions,
        FASTA_INDEX_BISMARK_BWAMETH.out.fasta,
        FASTA_INDEX_BISMARK_BWAMETH.out.fasta_index,
        FASTA_INDEX_BISMARK_BWAMETH.out.bismark_index,
        FASTA_INDEX_BISMARK_BWAMETH.out.bwameth_index,
        ch_primer_fasta3,
        ch_primer_fasta5
    )
    ch_versions = ch_versions.mix(METHYLSEQ.out.versions)

    emit:
    multiqc_report = METHYLSEQ.out.multiqc_report // channel: [ path(multiqc_report.html )  ]
    versions       = ch_versions                  // channel: [ path(versions.yml) ]

}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_METHYLSEQ (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_METHYLSEQ.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
