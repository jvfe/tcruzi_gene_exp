/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowKrakenclassify.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.kraken2_db, params.fasta, params.gtf ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { UNTAR                       } from '../modules/nf-core/untar/main'
include { GUNZIP as GUNZIP_FASTA      } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GTF        } from '../modules/nf-core/gunzip/main'
include { FASTP                       } from '../modules/nf-core/fastp/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2  } from '../modules/nf-core/kraken2/kraken2/main'
include { HISAT2_ALIGN                } from '../modules/nf-core/hisat2/align/main'
include { HISAT2_BUILD                } from '../modules/nf-core/hisat2/build/main'
include { STAR_ALIGN                  } from '../modules/nf-core/star/align/main'
include { STAR_GENOMEGENERATE         } from '../modules/nf-core/star/genomegenerate/main'
include { SUBREAD_FEATURECOUNTS       } from '../modules/nf-core/subread/featurecounts/main'
include { GATHER_COUNTS_HISAT2        } from '../modules/local/gather_counts_hisat2'
include { GATHER_COUNTS_STAR          } from '../modules/local/gather_counts_star'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow KRAKENCLASSIFY {

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    kraken_targz = Channel.of([[id: 'krakendb'], file(params.kraken2_db)])

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    ch_reads = INPUT_CHECK.out.reads

    ch_reads.dump(tag: "raw_reads")

    FASTQC (
        ch_reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    FASTP(
        ch_reads
    )
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    UNTAR (
        kraken_targz
    )

    UNTAR.out.untar
        .map { meta, path -> path }
        .first()
        .set { krakendb }

    if (params.fasta.endsWith('.gz')) {
        ch_fasta    = GUNZIP_FASTA ( [ [:], params.fasta ] ).gunzip.first()
        ch_fasta_star = ch_fasta.map{it[1]}
    } else {
        ch_fasta = tuple([:], file(params.fasta))
        ch_fasta_star = ch_fasta[1]
    }

    if (params.gtf.endsWith('.gz')) {
            ch_gtf      = GUNZIP_GTF ( [ [:], params.gtf ] ).gunzip.first()
            ch_gtf_star = ch_gtf.map{it[1]}
            ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
    } else {
        ch_gtf = tuple([:], file(params.gtf))
        ch_gtf_star = ch_gtf[1]
    }

    KRAKEN2 (
        FASTP.out.reads,
        krakendb,
        true,
        true
    )
    ch_versions = ch_versions.mix(KRAKEN2.out.versions.first())

    classified_reads = KRAKEN2.out.classified_reads_fastq.map{ meta, classified_reads_fastq -> tuple( meta, classified_reads_fastq ) }

    HISAT2_BUILD( ch_fasta, ch_gtf )

    HISAT2_ALIGN(
        classified_reads,
        HISAT2_BUILD.out.index,
    )

    ch_versions = ch_versions.mix(HISAT2_ALIGN.out.versions.first())

    STAR_GENOMEGENERATE(
        ch_fasta_star, ch_gtf_star
    )

    STAR_ALIGN(
        classified_reads,
        STAR_GENOMEGENERATE.out.index, ch_gtf_star, false, '', ''
    )
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())

    SUBREAD_FEATURECOUNTS(
        HISAT2_ALIGN.out.bam.map{ meta, path -> tuple( meta, path,  file(params.gtf)  ) }, 'gene_id'
    )
    ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions)

    GATHER_COUNTS_HISAT2(
        SUBREAD_FEATURECOUNTS.out.counts.collect{it[1]}
    )

    GATHER_COUNTS_STAR(
        STAR_ALIGN.out.read_per_gene_tab.collect{it[1]}
    )

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowKrakenclassify.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowKrakenclassify.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.log.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(KRAKEN2.out.report.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(STAR_ALIGN.out.log_final.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(HISAT2_ALIGN.out.summary.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(SUBREAD_FEATURECOUNTS.out.summary.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
