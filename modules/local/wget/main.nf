process WGET {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::gnu-wget=1.18"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h5bf99c6_5' :
        'quay.io/biocontainers/gnu-wget:1.18--h5bf99c6_5' }"

    input:
    tuple val(meta), val(download_path)

    output:
    path('*'), emit: download

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    wget $download_path
    """

    stub:
    """
    touch $path
    """

}
