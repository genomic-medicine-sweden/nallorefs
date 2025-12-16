process GUNZIP_REMOVE_HEADER_SORT_DBNSFP {
    tag "${archive}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

    input:
    tuple val(meta), path(archive)

    output:
    tuple val(meta), path("${gunzip}"), emit: gunzip
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def extension = (archive.toString() - '.gz').tokenize('.')[-1]
    def name = archive.toString() - '.gz' - ".${extension}"
    def prefix = task.ext.prefix ?: name
    gunzip = prefix + ".${extension}"

    // Removes header (tail -n +2), selects specific columns (cut -f ...),
    // 1 #chr
    // 2 pos(1-based)
    // 3 ref
    // 4 alt
    // 5 aaref
    // 6 aaalt
    // 7 rs_dbSNP
    // 15 Ensembl_transcriptid
    // 83 REVEL_score
    // 84 REVEL_rankscore
    // 184 GERP++_NR
    // 185 GERP++_RS
    // 187 phyloP100way_vertebrate
    // 193 phastCons100way_vertebrate
    // sorts by chrom and pos (sort -k1,1 -k2,2n)
    // to match expected order for dbNSFP usage with VEP
    // Needs to match the fields set in EXTRACT_HEADER_DBNSFP
    """
    # Not calling gunzip itself because it creates files
    # with the original group ownership rather than the
    # default one for that user / the work directory

    gzip \\
        -cd \\
        ${args} \\
        ${archive} | \\
        tail -n +2 | \\
        cut -f 1,2,3,4,5,6,7,15,83,84,184,185,187,193 | \\
        sort -k1,1 -k2,2n \\
        > ${gunzip}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gunzip: \$(echo \$(gunzip --version 2>&1) | sed 's/^.*(gzip) //; s/ Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def extension = (archive.toString() - '.gz').tokenize('.')[-1]
    def name = archive.toString() - '.gz' - ".${extension}"
    def prefix = task.ext.prefix ?: name
    gunzip = prefix + ".${extension}"
    """
    touch ${gunzip}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gunzip: \$(echo \$(gunzip --version 2>&1) | sed 's/^.*(gzip) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
