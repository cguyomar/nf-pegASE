params.options = [:]

process ASE_READ_COUNTER {
    //tag "$meta_bam.id"

    conda (params.enable_conda ? "bioconda::gatk4=4.2.0.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gatk4:4.2.0.0--0"
    } else {
        container "quay.io/biocontainers/gatk4:4.2.0.0--0"
    }

    input:
    tuple val(meta_bam), path(bam)
    tuple val(meta_vcf), path(vcf), path(vcf_index)
    tuple path(ref), path(ref_index), path(dict)
    
    output:
    path '*_Aligned.sortedByCoord_rg.bam' , emit:bam_rg

    script:
    """
    gatk ASEReadCounter \
        -I $bam \
        -V $vcf \
        -R $ref \
        -O $meta_bam.id \
        --min-mapping-quality  $params.mappingQuality
        --min-base-quality  $params.baseQuality
    """
}
