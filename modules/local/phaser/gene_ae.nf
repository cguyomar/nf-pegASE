// Import generic module functions
include { saveFiles } from '../functions'

params.options = [:]

process PHASER_GENE_AE {
    tag "$meta_bam.id"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'phaser', meta:[:], publish_by_meta:[]) }

    if (workflow.containerEngine == 'singularity' ) {
        container "library://cguyomar/phaser/phaser:0.1"
    }

    input:
    tuple val(meta_bam), path(hap_count)
    path(features)

    output:
    tuple val(meta_bam), path("*_phASER.gene_ae.txt"), emit: gene_ae


    script:
    def idSeparator   = params.idSeparator ? params.idSeparator : "_"
    """
    phaser_gene_ae.sh \
        --haplotypic_counts $hap_count \
        --features $features \
        --o ${meta_bam.id}_phASER.gene_ae.txt \
        --id_separator $idSeparator
    """
}



