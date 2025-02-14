/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CADD2VCF                                 } from '../modules/local/cadd2vcf/'
include { GUNZIP                                   } from '../modules/nf-core/gunzip/'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_CADD_SNVS } from '../modules/nf-core/bcftools/view/'
include { WGET as WGET_CADD_INDELS                 } from '../modules/local/wget/'
include { WGET as WGET_CADD_SNVS                   } from '../modules/local/wget/'
include { WGET as WGET_GENERAL                     } from '../modules/local/wget/'
include { WGET as WGET_VEP_PLUGIN_FILES            } from '../modules/local/wget/'
include { paramsSummaryMap                         } from 'plugin/nf-schema'
include { softwareVersionsToYAML                   } from '../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// This workflow should:
// - Download and unzip all reference-files
// - Download and prepare CADD-indels
// - Download and prepare CADD-SNVs
// - Download and prepare SNV-databases

workflow NALLOREFS {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:
    // General files stored in reference-files
    base_reference_dir           = 'https://github.com/Clinical-Genomics/reference-files/raw/cc2cfdf643094b821c7cfc47f35aa295fe9c3592/'
    ch_genmod_reduced_penetrance = downloadChannelOf(base_reference_dir + 'nallo/annotation/grch38_reduced_penetrance_-v1.0-.tsv')
    ch_target_pathogenic_repeats = downloadChannelOf(base_reference_dir + 'nallo/annotation/grch38_trgt_pathogenic_repeats_-v1.0-.txt')
    ch_variant_consequences      = downloadChannelOf(base_reference_dir + 'nallo/annotation/grch38_variant_consequences_-v1.0-.txt')
    ch_somalier_sites            = downloadChannelOf(base_reference_dir + 'nallo/annotation/sites.hg38.vcf.gz')
    ch_stranger_repeat_catalog   = downloadChannelOf('https://github.com/Clinical-Genomics/stranger/raw/84b97cf23ce522fa94d360612e7ed930789bb277/stranger/resources/variant_catalog_grch38.json')
    ch_genmod_snvs_rank_model    = downloadChannelOf(base_reference_dir + 'nallo/rank_model/grch38_rank_model_snvs_-v1.0-.ini')
    ch_genmod_svs_rank_model     = downloadChannelOf(base_reference_dir + 'nallo/rank_model/grch38_rank_model_svs_-v1.0-.ini')
    ch_target_regions            = downloadChannelOf(base_reference_dir + 'nallo/region/grch38_chromosomes_split_at_centromeres_-v1.0-.bed') // Not really necessary
    ch_hificnv_excluded_regions  = downloadChannelOf(base_reference_dir + 'nallo/region/grch38_hificnv_excluded_regions_common_50_-v1.0-.bed.gz')
    ch_hificnv_expected_xx_cn    = downloadChannelOf(base_reference_dir + 'nallo/region/grch38_hificnv_expected_copynumer_xx_-v1.0-.bed')
    ch_hificnv_expected_xy_cn    = downloadChannelOf(base_reference_dir + 'nallo/region/grch38_hificnv_expected_copynumer_xy_-v1.0-.bed')
    ch_par_regions               = downloadChannelOf(base_reference_dir + 'nallo/region/grch38_par.bed')
    ch_vep_loftool_scores        = downloadChannelOf(base_reference_dir + 'nallo/annotation/grch38_vep_112_loftool_scores_-v1.0-.txt') // Should perhaps VEP download process
    ch_vep_pli_values            = downloadChannelOf(base_reference_dir + 'nallo/annotation/grch38_vep_112_pli_values_-v1.0-.txt') // Should perhaps VEP download process
    // CADD indels
    ch_cadd_indels               = downloadChannelOf('https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/gnomad.genomes.r3.0.indel.tsv.gz')
    ch_cadd_indels_tbi           = downloadChannelOf('https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/gnomad.genomes.r3.0.indel.tsv.gz.tbi')

    // Initialise channels
    ch_versions = Channel.empty()
    ch_general_files_to_download = Channel.empty()
    ch_files_to_unzip = Channel.empty()

    //
    // Files that needs unzipping
    //
    ch_reference_genome = Channel.fromPath(
        'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz'
    ).map { it -> [ [ it.simpleName ], it ] }

    //
    // Download general files
    //
    ch_general_files_to_download = ch_general_files_to_download
        .mix(ch_genmod_reduced_penetrance)
        .mix(ch_target_pathogenic_repeats)
        .mix(ch_variant_consequences)
        .mix(ch_somalier_sites)
        .mix(ch_stranger_repeat_catalog)
        .mix(ch_genmod_snvs_rank_model)
        .mix(ch_genmod_svs_rank_model)
        .mix(ch_target_regions)
        .mix(ch_hificnv_excluded_regions)
        .mix(ch_hificnv_expected_xx_cn)
        .mix(ch_hificnv_expected_xy_cn)
        .mix(ch_par_regions)
        .mix(ch_vep_loftool_scores)
        .mix(ch_vep_pli_values)

    WGET_GENERAL ( ch_general_files_to_download )

    //
    // Unzip downloaded files
    // 
    ch_files_to_unzip = ch_files_to_unzip
        .mix(ch_reference_genome)

    GUNZIP(ch_files_to_unzip)


    if(!params.skip_cadd) {
        //
        // Download CADD indels
        //
        WGET_CADD_INDELS (
            ch_cadd_indels.concat(ch_cadd_indels_tbi)
        )
        
        //
        // CADD SNVs to echtvar
        //
        ch_cadd_snvs = downloadChannelOf(
            'https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz'
        )

        WGET_CADD_SNVS (
            ch_cadd_snvs
        )

        CADD2VCF (
            WGET_CADD_SNVS.out.download
        )

        //BCFTOOLS_VIEW_CADD_SNVS (
        //    CADD2VCF.out.vcf.map{ it -> [ [ id:'cadd2vcf', ], it ] }
        //)
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'nallorefs_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def downloadChannelOf(path) {
    Channel.of(path)
        .map { it -> [ [ 'id': file(it).simpleName ], it ] }
}