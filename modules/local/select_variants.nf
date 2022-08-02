params.options = [:]

process SELECT_VARIANTS {
    // publishDir "${params.outdir}",
    //     mode: params.publish_dir_mode,
        //saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'pipeline_info', meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::gatk4=4.2.0.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gatk4:4.2.0.0--0"
    } else {
        container "quay.io/biocontainers/gatk4:4.2.0.0--0"
    }

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path('*_SNP.vcf') , emit:vcf

    script:
    """

    gatk SelectVariants  \\
        -V $vcf  \\
        -O ${vcf.baseName}_SNP.vcf \\
        --select-type-to-include SNP \\
        --restrict-alleles-to BIALLELIC

    """
}
