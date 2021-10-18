// Import generic module functions
include { saveFiles } from './functions'

params.options = [:]

process PHASER {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'phaser', meta:[:], publish_by_meta:[]) }

    if (workflow.containerEngine == 'singularity' ) {
        container "library://cguyomar/phaser/phaser:latest"
    }

    input:
    tuple val(meta_bam), path(bam)
    tuple val(meta_vcf), path(vcf), path(vcf_index)

    output:
    path("*allelic_counts.txt"), emit: csv
    path("*haplotypic_counts.txt"), emit: hap_count


    script:
    """
    bgzip $vcf
    tabix ${vcf}.gz
    samtools index $bam
    phaser.sh \
        --vcf ${vcf}.gz \
        --bam $bam \
        --sample $meta_bam.id \
        --threads 1 \
        --o $meta_bam.id \
        --paired_end ${meta_bam.single_end ? "0" : "1"} \
        --unphased_vars 1 \
        --gw_phase_vcf 1 \
        --baseq $params.baseQuality \
        --mapq $params.mappingQuality #\
        #--id_separator $params.idSeparator
    """
}



