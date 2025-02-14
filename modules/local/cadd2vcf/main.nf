process CADD2VCF {
    tag "cadd2vcf"
    label 'process_low'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), val(download_path)

    output:
    path('*.vcf'), emit: vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cadd2vcf.py $download_path cadd.vcf
    """

    stub:
    """
    touch ${meta.id}.vcf
    """

}
