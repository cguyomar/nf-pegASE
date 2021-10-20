// Import generic module functions
include { saveFiles } from './functions'

params.options = [:]

process FILTER_JUNCTIONS {
    tag "$meta.id"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'pipeline_info', meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::gatk4=4.2.0.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gatk4:4.2.0.0--0"
    } else {
        container "quay.io/biocontainers/gatk4:4.2.0.0--0"
    }

    input:
    tuple val(meta), path(bam)

    output:
    path '*_Aligned.sortedByCoord_rg.bam' , emit:bam_rg

    script:
    """
    gatk --java-options "-Xmx{params.mem}" AddOrReplaceReadGroups -I {input.bam} \
        -O {output.bam} --RGID {params.idx} --RGLB {params.name}  \
        --RGPL {params.sequencer} --RGPU "-" --RGSM {params.name} 
    """
}
