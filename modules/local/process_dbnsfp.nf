process PROCESS_DBNSFP {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/htslib:1.20--h5efdd21_2' :
        'biocontainers/htslib:1.20--h5efdd21_2' }"

    input:
    tuple val(meta), path(files)
    val version

    output:
    tuple val(meta), path("*.gz") , emit: gz
    tuple val(meta), path("*.tbi"), emit: tbi
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // TODO: Instead of sorting, gunzip separately, then cat them together in a predetermined order without sorting,
    // and add header
    // https://www.ensembl.org/info/docs/tools/vep/script/vep_example.html#dbNSFP
    """
    zcat dbNSFP${version}_variant.chr1.gz | head -n1 > h || true

    zcat dbNSFP${version}_variant.chr* | grep -v "^#chr" | sort -k1,1 -k2,2n - | cat h - | bgzip -c > grch38_dbNSFP${version}.gz
    tabix -s 1 -b 2 -e 2 grch38_dbNSFP${version}.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch grch38_dbNSFP${version}.gz
    touch grch38_dbNSFP${version}.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
