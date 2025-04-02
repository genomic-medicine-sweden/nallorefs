process VCFEXPRESS_FILTER {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/fellen31/vcfexpress:0.3.4"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(pre_lua)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    prefix      = task.ext.prefix ?: "${meta.id}"
    expression = args.contains(' --expression ') || args.contains(' -e ') ? '' : '--expression "return true"'
    """
    vcfexpress filter \\
        $args \\
        $expression \\
        --lua-prelude $pre_lua \\
        --output ${prefix}.vcf.gz \\
        $vcf
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcfexpress: \$(echo \$(vcfexpress -V) | sed 's/vcfexpress //' )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcfexpress: \$(echo \$(vcfexpress -V) | sed 's/vcfexpress //' )
    END_VERSIONS
    """

}

