process ECHTVAR_ENCODE {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/fellen31/echtvar:0.2.1"

    input:
    tuple val(meta), path(bcf)
    tuple val(meta2), path(json)

    output:
    tuple val(meta), path("${meta.id}.zip"), emit: db
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    prefix      = task.ext.prefix ?: "${meta.id}"

    """
    echtvar encode ${prefix}.zip ${json} ${bcf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        echtvar: \$(echo \$(echtvar -V) | sed 's/echtvar //' )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.zip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        echtvar: \$(echo \$(echtvar -V) | sed 's/echtvar //' )
    END_VERSIONS
    """

}

