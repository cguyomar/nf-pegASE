params.options = [:]

process FILTER_PROPERLY_PAIRED {
    tag "$meta.id"
 
    conda (params.enable_conda ? 'bioconda::samtools=1.13' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/samtools:1.13--h8c37831_0"
    } else {
        container "quay.io/biocontainers/samtools:1.13--h8c37831_0"
    }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${bam.baseName}_pp.bam"), emit: bam_pp

    script:
    """
    samtools view -h -f 2 $bam -bo ${bam.baseName}_pp.bam
    """
}
