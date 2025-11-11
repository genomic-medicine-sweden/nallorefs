process CAT_BGZIP_TABIX_DBNSFP {
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
    // https://www.ensembl.org/info/docs/tools/vep/script/vep_example.html#dbNSFP
    """
    cat \\
        header \\
        dbNSFP${version}_variant.chr1 \\
        dbNSFP${version}_variant.chr10 \\
        dbNSFP${version}_variant.chr11 \\
        dbNSFP${version}_variant.chr12 \\
        dbNSFP${version}_variant.chr13 \\
        dbNSFP${version}_variant.chr14 \\
        dbNSFP${version}_variant.chr15 \\
        dbNSFP${version}_variant.chr16 \\
        dbNSFP${version}_variant.chr17 \\
        dbNSFP${version}_variant.chr18 \\
        dbNSFP${version}_variant.chr19 \\
        dbNSFP${version}_variant.chr2 \\
        dbNSFP${version}_variant.chr20 \\
        dbNSFP${version}_variant.chr21 \\
        dbNSFP${version}_variant.chr22 \\
        dbNSFP${version}_variant.chr3 \\
        dbNSFP${version}_variant.chr4 \\
        dbNSFP${version}_variant.chr5 \\
        dbNSFP${version}_variant.chr6 \\
        dbNSFP${version}_variant.chr7 \\
        dbNSFP${version}_variant.chr8 \\
        dbNSFP${version}_variant.chr9 \\
        dbNSFP${version}_variant.chrX \\
        dbNSFP${version}_variant.chrY |\\
        sort -k1,1 -k2,2n |\\
        bgzip -c > grch38_dbNSFP${version}.gz
    
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
