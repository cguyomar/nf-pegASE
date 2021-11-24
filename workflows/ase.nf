/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowAse.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta, params.star_align ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )
include { REMOVE_MULTIMAP } from '../modules/local/remove_multimap' addParams( options: [:] )
include { FILTER_PROPERLY_PAIRED } from '../modules/local/filter_properly_paired'
include { FILTER_JUNCTIONS } from '../modules/local/filter_junctions'
include { SELECT_VARIANTS } from '../modules/local/select_variants'
include { PHASER_MAIN } from '../modules/local/phaser/phaser'
include { PHASER_GENE_AE } from '../modules/local/phaser/gene_ae'
include { VARIANT_FILTRATION } from '../modules/local/variant_filtration'
//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check' addParams( options: [:] )

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

def star_align_options            = modules['star_align']
star_align_options.args          += params.save_unaligned ? Utils.joinModuleArgs(['--outReadsUnmapped Fastx']) : ''
if (params.save_align_intermeds)  { star_align_options.publish_files.put('bam','') }
if (params.save_unaligned)        { star_align_options.publish_files.put('fastq.gz','unmapped') }

def picard_markduplicates_options            = modules['picard_markduplicates']
picard_markduplicates_options.args          += " --READ_NAME_REGEX " + params.read_name_regex

def variantfiltration_options     = modules['gatk4_variantfiltration']

def samtools_view_options     = modules['samtools_view']
//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC  } from '../modules/nf-core/modules/fastqc/main'  addParams( options: modules['fastqc'] )
include { TRIMGALORE } from '../modules/nf-core/modules/trimgalore/main'
include { MULTIQC } from '../modules/nf-core/modules/multiqc/main' addParams( options: multiqc_options   )
include { STAR_ALIGN } from '../modules/nf-core/modules/star/align/main' addParams( options: star_align_options )
include { STAR_ALIGN as STAR_ALIGN_2 } from '../modules/nf-core/modules/star/align/main' addParams( options: star_align_options )
include { SAMTOOLS_MERGE } from '../modules/nf-core/modules/samtools/merge/main'
include { SAMTOOLS_FAIDX } from '../modules/nf-core/modules/samtools/faidx/main'
include { SAMTOOLS_VIEW } from '../modules/nf-core/modules/samtools/view/main' addParams( options: samtools_view_options)
include { STAR_GENOMEGENERATE } from '../modules/nf-core/modules/star/genomegenerate/main'
include { PICARD_MARKDUPLICATES } from '../modules/nf-core/modules/picard/markduplicates/main' addParams( options: picard_markduplicates_options )
include { GATK4_CREATESEQUENCEDICTIONARY } from '../modules/nf-core/modules/gatk4/createsequencedictionary/main'
include { GATK4_SPLITNCIGARREADS } from '../modules/nf-core/modules/gatk4/splitncigarreads/main'
include { GATK4_VARIANTFILTRATION } from '../modules/nf-core/modules/gatk4/variantfiltration/main' addParams( options: variantfiltration_options )






/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow ASE {

    ch_software_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    ).set { ch_fastq }

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_fastq
    )
    ch_software_versions = ch_software_versions.mix(FASTQC.out.version.first().ifEmpty(null))

    //
    // MODULE: Run Trim Galore
    //
    TRIMGALORE (
        ch_fastq
    )
    ch_software_versions = ch_software_versions.mix(TRIMGALORE.out.version.first().ifEmpty(null))

    //
    // Prepare index
    //
    STAR_GENOMEGENERATE (
      params.fasta,
      params.gtf
      )
    fai = SAMTOOLS_FAIDX(params.fasta).fai
    dict = GATK4_CREATESEQUENCEDICTIONARY(params.fasta).dict

    //
    // MODULE: Run first STAR alignment
    //
    STAR_ALIGN (
        TRIMGALORE.out.reads,
        STAR_GENOMEGENERATE.out.index,
        params.gtf,
        []
      )

    //
    // MODULE: FILTER_JUNCTIONS
    //
    FILTER_JUNCTIONS (
        STAR_ALIGN.out.bam_sorted
      )

    //
    // MODULE: Run second STAR alignment
    //
    STAR_ALIGN_2 (
        TRIMGALORE.out.reads,
        STAR_GENOMEGENERATE.out.index,
        params.gtf,
        FILTER_JUNCTIONS.out.sjout
      )

    // Merge BAMs
    STAR_ALIGN_2.out.bam_sorted.set { bams }

    bams.map{
        meta, bam ->
            meta.id = meta.id.split('_')[0..-2].join('_')
            [ meta, bam ] }
        .groupTuple(by: [0])
        .branch {
            meta, bam ->
               single  : bam.size() == 1
                   return [ meta, bam.flatten() ]
               multiple: bam.size() > 1
                   return [ meta, bam.flatten() ]
        }
        .set { bam_to_merge }

    SAMTOOLS_MERGE (
        bam_to_merge.multiple
    ).bam
    .mix(bam_to_merge.single)
    .set { bam_merged }

    // Filter properly paired alignments
    bam_merged.branch {
        meta, bam ->
            unpaired : meta.single_end
                return [ meta, bam]
            paired : !meta.single_end
                return [ meta, bam]
    }
    .set { bam_merged }
    SAMTOOLS_VIEW (
        bam_merged.paired,
        []
    )

    // Merge paired and unpaired bams
    bam_properly_paired = bam_merged.unpaired.mix(SAMTOOLS_VIEW.out.bam)

    //
    // MODULE: Run second STAR alignment
    //
    PICARD_MARKDUPLICATES ( 
        bam_properly_paired
    )
    
    //
    // MODULE: Run second STAR alignment
    //
    REMOVE_MULTIMAP (
        PICARD_MARKDUPLICATES.out.bam
    )

    fasta = [ params.fasta,
              fai,
              dict.collect()
    ]

    //
    // MODULE: Run Split N Cigar reads
    //
    GATK4_SPLITNCIGARREADS (
        REMOVE_MULTIMAP.out.bam,
        params.fasta,
        fai,
        GATK4_CREATESEQUENCEDICTIONARY.out.dict
    )

    //
    // MODULE: Run SELECT_VARIANTS
    //
    SELECT_VARIANTS([ ["id":"vcf"],  params.vcf])

    //
    // MODULE: Run GATK4_VARIANTFILTRATION
    //
    GATK4_VARIANTFILTRATION(
        SELECT_VARIANTS.out.vcf,
        params.fasta,
        fai,
        GATK4_CREATESEQUENCEDICTIONARY.out.dict
    )

    // Filter by genotyping proportion
    VARIANT_FILTRATION(
        GATK4_VARIANTFILTRATION.out.vcf,
        params.gtf,
        params.fasta
    )

    //
    // MODULE: Run PHASER
    //
    PHASER_MAIN(
        REMOVE_MULTIMAP.out.bam,
        VARIANT_FILTRATION.out.gt_vcf.collect()
    )

    PHASER_GENE_AE(
        PHASER_MAIN.out.hap_count,
        params.features
    )




    //
    // MODULE: Pipeline reporting
    //
    ch_software_versions
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { ch_software_versions }

    GET_SOFTWARE_VERSIONS (
        ch_software_versions.map { it }.collect()
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowAse.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report       = MULTIQC.out.report.toList()
    ch_software_versions = ch_software_versions.mix(MULTIQC.out.version.ifEmpty(null))
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
