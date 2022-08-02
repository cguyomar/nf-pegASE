params.options = [:]

process REMOVE_MULTIMAP {
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
    tuple val(meta), path("${bam.baseName}_uniq.bam"), emit: bam

    script:
    """
    samtools view -h -q 255 $bam -bo ${bam.baseName}_uniq.bam
    """
}
