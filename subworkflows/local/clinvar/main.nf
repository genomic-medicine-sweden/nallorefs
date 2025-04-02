include { BCFTOOLS_ANNOTATE as BCFTOOLS_ANNOTATE_RENAME_CLINVAR_CHRS } from '../../../modules/nf-core/bcftools/annotate/main' 
include { MD5SUM as MD5SUM_CLINVAR                                   } from '../../../modules/nf-core/md5sum/main'
include { TABIX_TABIX                                                } from '../../../modules/nf-core/tabix/tabix/main' 
include { VCFEXPRESS_FILTER as VCFEXPRESS_ADD_CLNVID                 } from '../../../modules/local/vcfexpress/filter/'
include { assertChecksum                                             } from '../../../subworkflows/local/utils_nfcore_nallorefs_pipeline/main'
workflow CLINVAR {

    take:
    ch_clinvar_vcf        // channel: [ val(meta), path(vcf) ]
    md5sum                //  string: md5sum
    ch_vcfexpress_prelude // channel: [ val(meta), path(lua) ]
    ch_rename_chrs        // channel: [ path(rename_chrs) ]

    main:
    ch_versions = Channel.empty()
    // Could perhaps be in another subworkflow
    MD5SUM_CLINVAR (
        ch_clinvar_vcf,
        false
    )
    ch_versions = ch_versions.mix(MD5SUM_CLINVAR.out.versions)
    
    // Could perhaps be checked in another subworkflow
    assertChecksum (MD5SUM_CLINVAR.out.checksum, md5sum)
    
    BCFTOOLS_ANNOTATE_RENAME_CLINVAR_CHRS (
        ch_clinvar_vcf.map { meta, vcf -> [ meta, vcf, [], [], [] ] },
        [],
        ch_rename_chrs
    )
    ch_versions = ch_versions.mix(BCFTOOLS_ANNOTATE_RENAME_CLINVAR_CHRS.out.versions)

    VCFEXPRESS_ADD_CLNVID (
        BCFTOOLS_ANNOTATE_RENAME_CLINVAR_CHRS.out.vcf,
        ch_vcfexpress_prelude
    )
    ch_versions = ch_versions.mix(VCFEXPRESS_ADD_CLNVID.out.versions)

    TABIX_TABIX (
        VCFEXPRESS_ADD_CLNVID.out.vcf
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)

    emit:
    vcf      = VCFEXPRESS_ADD_CLNVID.out.vcf // channel: [ val(meta), path(vcf) ]
    tbi      = TABIX_TABIX.out.tbi           // channel: [ val(meta), path(vcf) ]
    versions = ch_versions                   // channel: [ versions.yml ]
}
