process DOWNLOADGNOMAD {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52bef205b624ee5656d0932e38bab39148ec36fbcc6fdff08f4ae99d49e91b34/data':
        'community.wave.seqera.io/library/bcftools_curl:423214203ff8c6b0' }"

    input:
    tuple val(meta), path(download_path)

    output:
    tuple val(meta), path("*.vcf.gz")    , emit: vcf, optional: true
    tuple val(meta), path("*.bcf.gz")    , emit: bcf, optional: true
    tuple val(meta), path("*.bcf.gz.csi"), emit: csi, optional: true
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi, optional: true
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("--output-type b") || args.contains("-Ob") ? "bcf.gz" :
                    args.contains("--output-type u") || args.contains("-Ou") ? "bcf" :
                    args.contains("--output-type z") || args.contains("-Oz") ? "vcf.gz" :
                    args.contains("--output-type v") || args.contains("-Ov") ? "vcf" :
                    "vcf"
    // TODO: Perhaps remove annotate header line here so that md5sum can be snapshot
    """
        bcftools \\
            annotate \\
            $args \\
            --threads ${task.cpus-1} \\
            --include "FILTER='PASS'" \\
            --remove ^INFO/AF,INFO/AF_grpmax \\
            --output ${prefix}.${extension} \\
            ${download_path}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version |& sed '1!d ; s/bcftools //')
        curl: \$(curl --version | grep "curl" | sed -E 's/^curl ([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version |& sed '1!d ; s/bcftools //')
        curl: \$(curl --version | grep "curl" | sed -E 's/^curl ([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/')
    END_VERSIONS
    """
}
