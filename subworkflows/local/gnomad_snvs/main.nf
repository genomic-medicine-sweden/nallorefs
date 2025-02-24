include { MD5SUM as MD5SUM_GNOMAD_SNVS                       } from '../../../modules/nf-core/md5sum/main'
include { BCFTOOLS_ANNOTATE as BCFTOOLS_ANNOTATE_GNOMAD_SNVS } from '../../../modules/nf-core/bcftools/annotate/main'
include { BCFTOOLS_CONCAT                                    } from '../../../modules/nf-core/bcftools/concat/main' 
include { WGET as WGET_GNOMAD_SNVS                           } from '../../../modules/local/wget/'

workflow GNOMAD_SNVS {

    take:
    ch_gnomad_snvs_local  // channel: [ val(meta), path(vcf) ]
    ch_gnomad_snvs_remote // channel: [ val(meta), val(vcf) ]

    main:
    ch_versions = Channel.empty()

    // Should be stored in a file ...
    gnomad_checksums = [
        'gnomad.genomes.v4.1.sites.chr1.vcf.bgz'  : '328b4578212afec2cde394a1b02d544f',
        'gnomad.genomes.v4.1.sites.chr2.vcf.bgz'  : '518ca01e6757a68bc0abe76f85af644d',
        'gnomad.genomes.v4.1.sites.chr3.vcf.bgz'  : '716c181431a3c11a2eb18c5f50a3542d',
        'gnomad.genomes.v4.1.sites.chr4.vcf.bgz'  : 'd9b913f3e30c8f410f9ce7dee5dee6d4',
        'gnomad.genomes.v4.1.sites.chr5.vcf.bgz'  : 'a2a3b9014af5c8f9bbaaf743968e48f8',
        'gnomad.genomes.v4.1.sites.chr6.vcf.bgz'  : 'e65c2aa321c5e272548fb3d901bae382',
        'gnomad.genomes.v4.1.sites.chr7.vcf.bgz'  : '58ee22cf3dcc8cb8b493d218e19432cc',
        'gnomad.genomes.v4.1.sites.chr8.vcf.bgz'  : '9854d9df22977cf5bac0f9fc05f4e8f5',
        'gnomad.genomes.v4.1.sites.chr9.vcf.bgz'  : '6adfc9c47000cf66d1305392051b391d',
        'gnomad.genomes.v4.1.sites.chr10.vcf.bgz' : 'c2cd760130d2339f7135fc70700db1e1',
        'gnomad.genomes.v4.1.sites.chr11.vcf.bgz' : 'd1e7a4dcf3ff62eeffca57afefc5a33a',
        'gnomad.genomes.v4.1.sites.chr12.vcf.bgz' : '644bdbc5c53d9112edbacad401a28d1a',
        'gnomad.genomes.v4.1.sites.chr13.vcf.bgz' : '84b12f299210d2a7e390c56a234b7b68',
        'gnomad.genomes.v4.1.sites.chr14.vcf.bgz' : '0e524b414faf5a51b74d153939b3ddcd',
        'gnomad.genomes.v4.1.sites.chr15.vcf.bgz' : '00d1386eadbfcb653eb810a6c08ed250',
        'gnomad.genomes.v4.1.sites.chr16.vcf.bgz' : '2b106c12ca8cca9bd58e98d2c248ef4d',
        'gnomad.genomes.v4.1.sites.chr17.vcf.bgz' : '9d58b459e75312b487c660666ea540c4',
        'gnomad.genomes.v4.1.sites.chr18.vcf.bgz' : 'b95d30a6f5a45242d834fc7351afd760',
        'gnomad.genomes.v4.1.sites.chr19.vcf.bgz' : '73fa1c7c09072e8ca5e26a1443f0af2a',
        'gnomad.genomes.v4.1.sites.chr20.vcf.bgz' : '6c6f326c67b288ec99932905da72c1e6',
        'gnomad.genomes.v4.1.sites.chr21.vcf.bgz' : '9b6584cfe62f6e8bd0c3cea90a4ce56a',
        'gnomad.genomes.v4.1.sites.chr22.vcf.bgz' : 'a7bcf712a6b8d29e690468bf1dd8913d',
        'gnomad.genomes.v4.1.sites.chrX.vcf.bgz'  : '8b91766906865b0795c653af51cb73b8',
        'gnomad.genomes.v4.1.sites.chrY.vcf.bgz'  : '1ffb9c683674f41ff7cf524e5bb56bb8'
    ]
    
    // Only executed if remote path 
    WGET_GNOMAD_SNVS (
        ch_gnomad_snvs_remote
    )
    
    ch_gnomad_snvs_local
        .mix(WGET_GNOMAD_SNVS.out.download)
        .set { ch_gnomad_snvs }
    
    MD5SUM_GNOMAD_SNVS (
        ch_gnomad_snvs, 
        false
    )
    ch_versions = ch_versions.mix(MD5SUM_GNOMAD_SNVS.out.versions)

    // Get expected - local checksum pairs for all gnomAD files
    Channel.from(gnomad_checksums.collect { key, value -> [ ['id': key ], value ] } )
        .join(MD5SUM_GNOMAD_SNVS.out.checksum) // Failonmismatch!!!
        .set { ch_md5sum_check_gnomad }

    assertChecksumChannel( ch_md5sum_check_gnomad )

    BCFTOOLS_ANNOTATE_GNOMAD_SNVS (
        ch_gnomad_snvs.map { meta, vcf -> [ meta, vcf, [], [] , [] ] },
        [],
        []
    )
    ch_versions = ch_versions.mix(BCFTOOLS_ANNOTATE_GNOMAD_SNVS.out.versions)
    
    BCFTOOLS_ANNOTATE_GNOMAD_SNVS.out.vcf
        .map { _meta, bcf -> [ [ 'id': 'gnomad_snvs' ], bcf ] }
        .groupTuple()
        .map { meta, bcfs -> [ meta, bcfs, [] ] }
        .set { ch_bcftools_concat_gnomad }

    BCFTOOLS_CONCAT (
        ch_bcftools_concat_gnomad
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)
 
    emit:
    vcf      = BCFTOOLS_CONCAT.out.vcf // channel: [ val(meta), path(vcf) ]
    tbi      = BCFTOOLS_CONCAT.out.tbi // channel: [ val(meta), path(vcf) ]
    versions = ch_versions             // channel: [ versions.yml ]
}

def assertChecksumChannel(channel) {
    channel
        .map { _meta, expected_checksum, file ->
            def checksum = file.readLines().collect { line -> line.split('  ')[0] }[0]
            def filename = file.readLines().collect { line -> line.split('  ')[1] }[0]
            assert checksum == expected_checksum : "Checksum for ${filename} was '${checksum}' but expected '${expected_checksum}'"
        }
}
