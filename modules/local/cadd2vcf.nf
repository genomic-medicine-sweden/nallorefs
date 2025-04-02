process CADD2VCF {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52bef205b624ee5656d0932e38bab39148ec36fbcc6fdff08f4ae99d49e91b34/data':
        'community.wave.seqera.io/library/bcftools_curl:423214203ff8c6b0' }"

    input:
    tuple val(meta), path(download_path)

    output:
    tuple val(meta), path("*.bcf.gz"), emit: bcf
    tuple val(meta), path("*.csi")   , emit: csi
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO: Perhaps remove annotate header line here so that md5sum can be snapshot
    """
    gunzip -c ${download_path} | \\
        awk '
            BEGIN {
                # Print the VCF headere
                print "##fileformat=VCFv4.2";
                print "##INFO=<ID=raw,Number=1,Type=Float,Description=\"raw cadd score\">";
                print "##INFO=<ID=phred,Number=1,Type=Float,Description=\"phred-scaled cadd score\">";
                for (i = 1; i <= 22; i++) print "##contig=<ID=" i ">";
                print "##contig=<ID=X>";
                print "##contig=<ID=Y>";
            }
            NR==1 {
                sub(/^# /, "", \$0);
                print "##CADDCOMMENT=<ID=comment,comment=\\"" \$0 "\\">";# Skip the second line
                print "#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO";
                next;
            }
            NR > 2 {
                # Process data lines
                printf "%s\\t%s\\t.\\t%s\\t%s\\t1\\tPASS\\traw=%.2f;phred=%.2f\\n", \$1, \$2, \$3, \$4, \$5, \$6;
            }' | \\
        bcftools \\
            view \\
            --threads ${task.cpus-1} \\
            --output-type b \\
            --write-index \\
            --output ${prefix}.bcf.gz

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
