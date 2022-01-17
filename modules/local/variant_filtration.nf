// Import generic module functions
include { saveFiles } from './functions'

params.options = [:]

process VARIANT_FILTRATION {

    label 'process_high'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'pipeline_info', meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::pyfaidx" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "library://cguyomar/default/variantfiltration:sha256.a2a9c3ff8ecf924899e195ede4d9fddf93359c099af2b0a182aafe01cf9a7a34"
    } else {
        container "library://cguyomar/default/variantfiltration:sha256.a2a9c3ff8ecf924899e195ede4d9fddf93359c099af2b0a182aafe01cf9a7a34"
    }

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path(gtf)
    path(fasta)

    output:
    tuple val(meta), path('*.annotated.vcf') , emit:annotated_vcf
    path('*.gt.vcf') , emit:gt_vcf
    tuple val(meta), path('*.gt.tsv') , emit:gt_tsv
    tuple val(meta), path('*.dp.tsv') , emit:dp_tsv
    tuple val(meta), path('*.ad.tsv') , emit:ad_tsv
    tuple val(meta), path('*.filtration_summary.txt') , emit:summary

    script:
    """
    samtools faidx $fasta -o ${fasta}.fai

    variantFiltration.py \
    --nb-cpus $task.cpus \
    --in-vcf $vcf \
    --in-gtf $gtf \
    --in-fasta $fasta \
    --out-vcf ${meta.id}.annotated.vcf \
    --out-vcf-gt ${meta.id}.gt.vcf \
    --out-gt ${meta.id}.gt.tsv \
    --out-dp ${meta.id}.dp.tsv \
    --out-ad ${meta.id}.ad.tsv \
    --summary ${meta.id}.filtration_summary.txt
    
    """
}
