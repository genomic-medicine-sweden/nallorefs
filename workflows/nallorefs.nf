include { samplesheetToList } from 'plugin/nf-schema'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CADD2VCF                                                   } from '../modules/local/cadd2vcf.nf'
include { CLINVAR                                                    } from '../subworkflows/local/clinvar/main'
include { GNOMAD_SNVS                                                } from '../subworkflows/local/gnomad_snvs/main'
include { GUNZIP                                                     } from '../modules/nf-core/gunzip/'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_CADD_SNVS                   } from '../modules/nf-core/bcftools/view/'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_GNOMAD_SVS                  } from '../modules/nf-core/bcftools/view/'
include { ECHTVAR_ENCODE                                             } from '../modules/local/echtvar/encode/'
include { MD5SUM as MD5SUM_CADD_ANNOTATIONS                          } from '../modules/nf-core/md5sum/main'
include { MD5SUM as MD5SUM_CADD_SNVS                                 } from '../modules/nf-core/md5sum/main'
include { MD5SUM as MD5SUM_LOCAL_ECTHVAR_DATABASES                   } from '../modules/nf-core/md5sum/main'
include { MD5SUM as MD5SUM_LOCAL_SVDB_DATABASES                      } from '../modules/nf-core/md5sum/main'
include { UNTAR as UNTAR_VEP_CACHE                                   } from '../modules/nf-core/untar/main'
include { UNTAR as UNTAR_CADD_ANNOTATIONS                            } from '../modules/nf-core/untar/main'
include { WGET as WGET_CADD_ANNOTATIONS                              } from '../modules/local/wget/'
include { WGET as WGET_CADD_INDELS                                   } from '../modules/local/wget/'
include { WGET as WGET_CADD_SNVS                                     } from '../modules/local/wget/'
include { WGET as WGET_GENERAL                                       } from '../modules/local/wget/'
include { WGET as WGET_VEP_PLUGIN_FILES                              } from '../modules/local/wget/'
include { WGET as WGET_LOCAL_SVDB_DATABASES                          } from '../modules/local/wget/'
include { paramsSummaryMap                                           } from 'plugin/nf-schema'
include { softwareVersionsToYAML                                     } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { assertChecksum; assertChecksumChannel                      } from '../subworkflows/local/utils_nfcore_nallorefs_pipeline/main'

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
    ch_samplesheet // channel: mock samplesheet read in from --input

    main:
    // All files that can be downloaded without post-processing are created from downloadChannelOf
    // and downloaded with a WGET process to control publishDir. Files that can be put in directly into the
    // pipeline outdir are downloaded though WGET_GENERAL. Files that needs to be put into specific directory
    // structures are downloaded with their respective WGET_* process.
    //
    // Most files that requries some form of post-processing are created with Channel.fromPath('...') and are
    // staged by the Nextflow headjob, input into a process and published by the final process, except for some of
    // the bigger ones e.g. the CADD resources, that may be better downloaded in a process.

    // General files stored in reference-files
    ch_genmod_reduced_penetrance = downloadChannelOf(params.base_reference_dir + 'nallo/annotation/grch38_reduced_penetrance_-v1.0-.tsv')
    ch_trgt_pathogenic_repeats   = downloadChannelOf(params.base_reference_dir + 'nallo/annotation/grch38_trgt_pathogenic_repeats.bed')
    ch_variant_consequences      = downloadChannelOf(params.base_reference_dir + 'nallo/annotation/grch38_variant_consequences_-v1.0-.txt')
    ch_somalier_sites            = downloadChannelOf(params.base_reference_dir + 'nallo/annotation/sites.hg38.vcf.gz')
    ch_stranger_repeat_catalog   = downloadChannelOf(params.base_reference_dir + 'nallo/annotation/grch38_variant_catalog_stranger.json')
    ch_genmod_snvs_rank_model    = downloadChannelOf(params.base_reference_dir + 'nallo/rank_model/grch38_rank_model_snvs_-v1.0-.ini')
    ch_genmod_svs_rank_model     = downloadChannelOf(params.base_reference_dir + 'nallo/rank_model/grch38_rank_model_svs_-v1.0-.ini')
    ch_target_regions            = downloadChannelOf(params.base_reference_dir + 'nallo/region/grch38_chromosomes_split_at_centromeres_-v1.0-.bed') // Not really necessary
    ch_hificnv_excluded_regions  = downloadChannelOf(params.base_reference_dir + 'nallo/region/grch38_hificnv_excluded_regions_common_50_-v1.0-.bed.gz')
    ch_hificnv_expected_xx_cn    = downloadChannelOf(params.base_reference_dir + 'nallo/region/grch38_hificnv_expected_copynumer_xx_-v1.0-.bed')
    ch_hificnv_expected_xy_cn    = downloadChannelOf(params.base_reference_dir + 'nallo/region/grch38_hificnv_expected_copynumer_xy_-v1.0-.bed')
    ch_par_regions               = downloadChannelOf(params.base_reference_dir + 'nallo/region/grch38_par_-v1.0-.bed')
    ch_vep_loftool_scores        = downloadChannelOf(params.base_reference_dir + 'nallo/annotation/grch38_vep_112_loftool_scores_-v1.0-.txt') // Should perhaps VEP download process
    ch_vep_pli_values            = downloadChannelOf(params.base_reference_dir + 'nallo/annotation/grch38_vep_112_pli_values_-v1.0-.txt') // Should perhaps VEP download process
    // Files that needs unzipping
    ch_reference_genome          = fileChannelOf('https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz')
    // CADD resources. Either download and checksum remote file, or just checksum local.
    ch_cadd_annotations          = remoteLocalBranchedChannelOf(params.cadd_annotations)
    // CADD prescored indels
    ch_cadd_indels               = downloadChannelOf('https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/gnomad.genomes.r3.0.indel.tsv.gz')
    ch_cadd_indels_tbi           = downloadChannelOf('https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/gnomad.genomes.r3.0.indel.tsv.gz.tbi')
    // CADD SNVs. Either download and checksum remote file, or just checksum local.
    ch_cadd_snvs                     = remoteLocalBranchedChannelOf(params.cadd_snvs)
    ch_echtvar_encode_cadd_snvs_json = fileChannelOf("${projectDir}/assets/ecthvar_encode_cadd_snvs.json") // TODO: move to reference-files?
    // ClinVar
    ch_clinvar_vcf                       = fileChannelOf(params.clinvar)
    ch_clinvar_rename_chrs               = fileChannelOf(params.base_reference_dir + 'nallo/annotation/grch38_clinvar_rename_chromosomes_-v1.0-.txt').map { _meta, file -> file }
    ch_vcfexpress_add_clnvid_prelude     = fileChannelOf("${projectDir}/assets/vcfexpress_add_clnvid_prelude.lua")
    // CoLoRSdb SNVs - path changes every new Zenodo release
    ch_colorsdb_snvs                     = fileChannelOf("https://zenodo.org/records/14814308/files/CoLoRSdb.GRCh38.v1.2.0.deepvariant.glnexus.vcf.gz")
    ch_echtvar_encode_colorsdb_snvs_json = fileChannelOf("${projectDir}/assets/echtvar_encode_colorsdb_snvs.json") // TODO: move to reference-files?
    // CoLoRSdb SVs - could think of renaming these, or storing in a separate location
    ch_colorsdb_svs_vcf                  = downloadChannelOf('https://zenodo.org/records/14814308/files/CoLoRSdb.GRCh38.v1.2.0.pbsv.jasmine.vcf.gz')
    // GnomAD SVs
    ch_gnomad_svs_vcf                    = fileChannelOf("https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz")
    // GnomAD SNVs
    ch_echtvar_encode_gnomad_snvs_json   = fileChannelOf("${projectDir}/assets/gnomad.json") // TODO: move to reference-files?
    chromosomes = (1..22) + ['X', 'Y']
    ch_gnomad_snvs_per_chr = Channel.fromList(chromosomes)
        .map { chr -> "${params.gnomad_base_path}gnomad.genomes.v${params.gnomad_version}.sites.chr${chr}.vcf.bgz" }
        .map { it -> [ [ 'id': file(it).name ], it ] }
        .branch { meta, file_path ->
            remote: file_path.startsWith('https://')
                return [ meta, file_path ] // return file_path as string
            local: !file_path.startsWith('https://')
                return [ meta, file(file_path) ]
        }

    // Local SNV databases
    ch_echtvar_local_databases = params.local_echtvar_databases ? Channel.fromList(samplesheetToList(params.local_echtvar_databases, 'assets/schema_local_echtvar_databases.json'))
        .map { meta, vcf, json ->
            [ meta + [ 'id': file(vcf).name ], vcf, json ] } : Channel.empty()

    // VEP cache - 26 Gb
    vep_cache_type_string = params.vep_cache_type == "merged" ? '_merged' : ''
    ch_vep_cache = fileChannelOf("https://ftp.ensembl.org/pub/release-${params.vep_cache_version}/variation/indexed_vep_cache/homo_sapiens${vep_cache_type_string}_vep_${params.vep_cache_version}_GRCh38.tar.gz")

    // Initialise channels
    ch_versions                  = Channel.empty()
    ch_general_files_to_download = Channel.empty()
    ch_echtvar_encode_files      = Channel.empty()

    //
    // Download general files
    //
    if (!params.skip_general_files) {
        ch_general_files_to_download = ch_general_files_to_download
            .mix(ch_genmod_reduced_penetrance)
            .mix(ch_trgt_pathogenic_repeats)
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
            .mix(ch_colorsdb_svs_vcf) // TODO: think about renaming the output files of these - and check md5
        WGET_GENERAL ( ch_general_files_to_download )

    }

    //
    // Unzip/untar downloaded files
    //

    GUNZIP (
        ch_reference_genome
    )

    if(!params.skip_vep_cache) {
        UNTAR_VEP_CACHE (
            ch_vep_cache
        )
    }

    //
    // Download CADD resources (200GB)
    //
    if(!params.skip_cadd_annotations) {

        // Only executed if input file is remote
        WGET_CADD_ANNOTATIONS (
            ch_cadd_annotations.remote
        )

        ch_cadd_annotations.local
            .mix(WGET_CADD_ANNOTATIONS.out.download)
            .set { ch_cadd_annotations_to_untar }

        MD5SUM_CADD_ANNOTATIONS (
            ch_cadd_annotations_to_untar,
            false
        )

        // Might be okay to just checksum 'final files' sometimes?
        // But here we would just untar the files, and that's a lot of files to checksum...without having a 'truth'
        assertChecksum (MD5SUM_CADD_ANNOTATIONS.out.checksum, params.cadd_annotations_md5sum)

        UNTAR_CADD_ANNOTATIONS (
            ch_cadd_annotations_to_untar
        )
    }

    //
    // Download CADD indels
    //
    if (!params.skip_cadd_indels) {
        WGET_CADD_INDELS (
            ch_cadd_indels.concat(ch_cadd_indels_tbi)
        )
    }

    //
    // Local echtvar SNV databases
    //
    if (!params.skip_local_echtvar_databases) {

        ch_echtvar_encode_files = ch_echtvar_encode_files
            .mix(
                ch_echtvar_local_databases
            )

        // To ensure local files are correct and not corrupted.
        MD5SUM_LOCAL_ECTHVAR_DATABASES (
            ch_echtvar_encode_files.map { meta, vcf, json -> [ meta, vcf ] },
            false
        )

        assertChecksumChannel (
            MD5SUM_LOCAL_ECTHVAR_DATABASES.out.checksum,
            MD5SUM_LOCAL_ECTHVAR_DATABASES.out.checksum.map { meta, _checksum -> meta.md5sum }
        )
    }


    //
    // Download CADD SNVs and convert to echtvar zip
    //
    if (!params.skip_cadd_snvs) {

        if(!params.cadd_snvs_echtvar_zip) {
            WGET_CADD_SNVS (
                ch_cadd_snvs.remote
            )

            ch_cadd_snvs.local
                .mix(WGET_CADD_SNVS.out.download)
                .set { ch_cadd_snvs_to_convert }

            MD5SUM_CADD_SNVS (
                ch_cadd_snvs_to_convert,
                false
            )

            assertChecksum (MD5SUM_CADD_SNVS.out.checksum, params.cadd_snvs_md5sum)

            // TODO: if echtvar-archive is provided, then just md5-sumcheck the echtvar file.
            // This is because the CADD-file contains all possible variants in the genome,
            // therefore takes a really long time to go through...

            CADD2VCF (
                ch_cadd_snvs_to_convert
            )

        }

        ch_echtvar_encode_files = ch_echtvar_encode_files
            .mix(
                CADD2VCF.out.bcf
                    .combine(ch_echtvar_encode_cadd_snvs_json)
                    .map { vcf_meta, vcf, _json_meta, json -> [ vcf_meta, vcf, json ] }
            )

    }

    //
    // Reformat ClinVar to reference genome and Scout compatibility
    // (1 -> chr1 & ID -> INFO/CLNVID, annotated via VEP)
    //
    if(!params.skip_clinvar) {
        CLINVAR (
            ch_clinvar_vcf,
            params.clinvar_md5sum,
            ch_vcfexpress_add_clnvid_prelude,
            ch_clinvar_rename_chrs
        )
    }

    //
    // Download/Get GnomAD local files and strip info/convert
    //
    if(!params.skip_gnomad_snvs) {

        GNOMAD_SNVS (
            ch_gnomad_snvs_per_chr.local,
            ch_gnomad_snvs_per_chr.remote
        )

        // Add to echtvar files
        ch_echtvar_encode_files = ch_echtvar_encode_files
        .mix(
            GNOMAD_SNVS.out.vcf
                .combine(ch_echtvar_encode_gnomad_snvs_json)
                .map { vcf_meta, vcf, _json_meta, json -> [ vcf_meta, vcf, json ] }
            )
    }

    //
    // Remove `CNV` SV-type, else SVDB will fail, e.g. https://github.com/nf-core/raredisease/issues/615
    //
    BCFTOOLS_VIEW_GNOMAD_SVS (
        ch_gnomad_svs_vcf.map { meta, vcf -> [ meta, vcf, [] ] },
        [],
        [],
        []
    )

    //
    // Echtvar SNV databases - GnomAD, CADD, CoLoRsDB
    //

    ch_echtvar_encode_files = ch_echtvar_encode_files
        .mix(
            ch_colorsdb_snvs
                .combine(ch_echtvar_encode_colorsdb_snvs_json)
                .map { vcf_meta, vcf, _json_meta, json -> [ vcf_meta, vcf, json ] }
            )

    ch_echtvar_encode_files
        .multiMap { vcf_meta, vcf, json ->
            vcf:  [ vcf_meta, vcf ]
            json: [ vcf_meta, json ]
        }
        .set { ch_echtvar_encode_in }

    ECHTVAR_ENCODE (
        ch_echtvar_encode_in.vcf,
        ch_echtvar_encode_in.json
    )

    // TODO: Could do a MD5-sum check of echtvar encode,
    // especially CADD SNVs if params.cadd_snvs_ecthvar_zip has been set

    //
    // TODO: Could/should do a md5sumcheck of all final files and write to output directories
    //

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'nallorefs_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Creates a value channel with meta, input to download processes
def downloadChannelOf(path) {
    // Channel.of(file(path, checkifexists: true)) ?
    Channel.of(path)
        .map { it -> [ [ 'id': file(it).name ], it ] }
}

def fileChannelOf(path) {
    Channel.fromPath(path)
        .map { it -> [ [ 'id': it.name ], it ] }
}


def remoteLocalBranchedChannelOf(path) {
    Channel.of(path)
        .map { it -> [ [ 'id': file(it).name ], it ] }
        .branch { meta, file_path ->
            remote: file_path.startsWith('https://')
                return [ meta, file_path ] // return file_path as string
            local: !file_path.startsWith('https://')
                return [ meta, file(file_path) ]
        }
}
