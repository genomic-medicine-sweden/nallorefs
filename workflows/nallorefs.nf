/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CADD2VCF                                                   } from '../modules/local/cadd2vcf.nf'
include { GUNZIP                                                     } from '../modules/nf-core/gunzip/'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_CADD_SNVS                   } from '../modules/nf-core/bcftools/view/'
include { BCFTOOLS_ANNOTATE as BCFTOOLS_ANNOTATE_RENAME_CLINVAR_CHRS } from '../modules/nf-core/bcftools/annotate/main' 
include { BCFTOOLS_CONCAT                                            } from '../modules/nf-core/bcftools/concat/main' 
include { DOWNLOADGNOMAD as DOWNLOAD_GNOMAD_SNVS                     } from '../modules/local/downloadgnomad.nf'
include { DOWNLOADGNOMAD as DOWNLOAD_GNOMAD_SVS                      } from '../modules/local/downloadgnomad.nf'
include { ECHTVAR_ENCODE                                             } from '../modules/local/echtvar/encode/'
include { MD5SUM as MD5SUM_CADD_ANNOTATIONS                          } from '../modules/nf-core/md5sum/main'
include { VCFEXPRESS_FILTER as VCFEXPRESS_ADD_CLNVID                 } from '../modules/local/vcfexpress/filter/'
include { TABIX_TABIX                                                } from '../modules/nf-core/tabix/tabix/main' 
include { UNTAR as UNTAR_VEP_CACHE                                   } from '../modules/nf-core/untar/main'
include { UNTAR as UNTAR_CADD_ANNOTATIONS                            } from '../modules/nf-core/untar/main'
include { WGET as WGET_CADD_ANNOTATIONS                              } from '../modules/local/wget/'
include { WGET as WGET_CADD_INDELS                                   } from '../modules/local/wget/'
include { WGET as WGET_CADD_SNVS                                     } from '../modules/local/wget/'
include { WGET as WGET_GENERAL                                       } from '../modules/local/wget/'
include { WGET as WGET_VEP_PLUGIN_FILES                              } from '../modules/local/wget/'
include { paramsSummaryMap                                           } from 'plugin/nf-schema'
include { softwareVersionsToYAML                                     } from '../subworkflows/nf-core/utils_nfcore_pipeline'

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
    //
    // All files that can be downloaded without post-processing are created from downloadChannelOf 
    // and downloaded with a WGET process to control publishDir. Files that can be put in directly into the
    // pipeline outdir are downloaded though WGET_GENERAL. Files that needs to be put into specific directory
    // structures are downloaded with their respective WGET_* process.
    //
    // All files that requries some form of post-processing are created with Channel.fromPath('...') and are
    // staged by the Nextflow headjob, input into a process and published by the final process. Some of these
    // are quite large, e.g. the CADD resources, and may be better downloaded in a process rather than by the head job.
    // Then they would each require a download process, or meta to branch later. Unclear what is best here.
    // If it's one download at a time maybe head job staging is not the way to go.
    //

    // General files stored in reference-files
    base_reference_dir           = 'https://github.com/Clinical-Genomics/reference-files/raw/cc2cfdf643094b821c7cfc47f35aa295fe9c3592/'
    ch_genmod_reduced_penetrance = downloadChannelOf(base_reference_dir + 'nallo/annotation/grch38_reduced_penetrance_-v1.0-.tsv')
    ch_trgt_pathogenic_repeats   = downloadChannelOf(base_reference_dir + 'nallo/annotation/grch38_trgt_pathogenic_repeats.bed')
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
    // Files that needs unzipping
    ch_reference_genome          = Channel.fromPath('https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz').map { it -> [ [ 'id': file(it).name ], it ] }
    // CADD resources - needs to be untarred. 
    // Since this file is so big, I think this pipeline _should_ perhaps be able to take a local file here,
    // that has been checked with md5. TODO: Perhaps the pipeline could validate then instead of download?
    ch_cadd_annotations          = remoteLocalBranchedChannelOf(params.cadd_annotations)
    // CADD prescored indels
    ch_cadd_indels               = downloadChannelOf('https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/gnomad.genomes.r3.0.indel.tsv.gz')
    ch_cadd_indels_tbi           = downloadChannelOf('https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/gnomad.genomes.r3.0.indel.tsv.gz.tbi')
    // CADD SNVs - 200 GB
    // Same for this file, it's so big so it might be okay to have a local file and then the pipeline validates against the md5
    // 
    ch_cadd_snvs = downloadChannelOf('https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz') 
    // echvtar jsons
    ch_echtvar_encode_cadd_snvs_json = Channel.fromPath("${projectDir}/assets/echtvar_encode_cadd_snvs.json").map { it -> [ [ 'id': it.name ], it ] } // TODO: move to reference-files?
    // ClinVar
    ch_clinvar_vcf                   = Channel.fromPath("https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar_${params.clinvar_version}.vcf.gz").map { it -> [ [ 'id': file(it).name ], it ] }
    ch_clinvar_rename_chrs           = Channel.fromPath(base_reference_dir + 'nallo/annotation/grch38_clinvar_rename_chromosomes_-v1.0-.txt')
    ch_vcfexpress_add_clnvid_prelude = Channel.fromPath("${projectDir}/assets/vcfexpress_add_clnvid_prelude.lua").map { it -> [ [ 'id': it.name ], it ] }
    // CoLoRSdb SNVs - path changes every new Zenodo release
    ch_colorsdb_snvs                     = Channel.fromPath("https://zenodo.org/records/14814308/files/CoLoRSdb.GRCh38.v1.2.0.deepvariant.glnexus.vcf.gz").map { it -> [ [ 'id': it.name ], it ] } 
    ch_echtvar_encode_colorsdb_snvs_json = Channel.fromPath("${projectDir}/assets/echtvar_encode_colorsdb_snvs.json").map { it -> [ [ 'id': it.name ], it ] } // TODO: move to reference-files?
    // CoLoRSdb SVs - TODO: make sure nothing needs to be done with these?
    ch_colorsdb_svs_vcf = downloadChannelOf('https://zenodo.org/records/14814308/files/CoLoRSdb.GRCh38.v1.2.0.pbsv.jasmine.vcf.gz')
    ch_colorsdb_svs_tbi = downloadChannelOf('https://zenodo.org/records/14814308/files/CoLoRSdb.GRCh38.v1.2.0.pbsv.jasmine.vcf.gz.tbi')
    // VEP cache
    vep_cache_type_string = params.vep_cache_type == "merged" ? '_merged' : ''
    ch_vep_cache = Channel.fromPath("https://ftp.ensembl.org/pub/release-${params.vep_cache_version}/variation/indexed_vep_cache/homo_sapiens${vep_cache_type_string}_vep_${params.vep_cache_version}_GRCh38.tar.gz").map { it -> [ [ 'id': it.name ], it ] } 
    // GnomAD SNVs - these should not be staged by head job - so need to be passed as val
    chromosomes = (1..22) + ['X', 'Y']
    gnomad_base_path = "https://storage.googleapis.com/gcp-public-data--gnomad/release/${params.gnomad_version}/vcf/genomes/gnomad.genomes.v${params.gnomad_version}.sites.chr"
    ch_gnomad_snvs_per_chr = Channel.fromList(chromosomes)
        .map { chr -> 
            def gnomad_full_path = "${gnomad_base_path}${chr}.vcf.bgz"
            return [ [ 'id': file(gnomad_full_path).name ], gnomad_full_path ]
        }
    // GnomAD SVs
    ch_gnomad_svs = downloadChannelOf("https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz")
    // Initialise channels
    ch_versions = Channel.empty()
    ch_general_files_to_download = Channel.empty()
    ch_files_to_unzip = Channel.empty()
    ch_echtvar_encode_files = Channel.empty() 
    
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
            .mix(ch_colorsdb_svs_vcf) // TODO: think about renaming the output files of these
            .mix(ch_colorsdb_svs_tbi) // TODO: think about renaming the output files of these 

        WGET_GENERAL ( ch_general_files_to_download )

    }
    
    //
    // Unzip downloaded files
    // 
    ch_files_to_unzip = ch_files_to_unzip
        .mix(ch_reference_genome)

    GUNZIP (
        ch_files_to_unzip
    )
    
    //
    // Untar VEP cache - 26GB download
    //
    /* TODO: Fix error with untar
    if(!params.skip_vep_cache) {
        UNTAR_VEP_CACHE (
            ch_vep_cache
        )
    }
    */
    //
    // Untar CADD resources - 200GB download
    //
    if(!params.skip_cadd_annotations) {
        
        // Only executed if input file is remote
        WGET_CADD_ANNOTATIONS (
            ch_cadd_annotations.remote
        )

        ch_cadd_annotations.local
            .mix(WGET_CADD_ANNOTATIONS.out.download)
            .set { ch_cadd_annotations_to_untar }
        
        // Should take ~20 min for 200 GB file
        MD5SUM_CADD_ANNOTATIONS (
            ch_cadd_annotations_to_untar,
            false
        )

        // Might be okay to just checksum 'final files' sometimes?
        // But here we would just untar the files, and that's a lot of files to checksum...without having a 'truth'
        assertChecksum (MD5SUM_CADD_ANNOTATIONS.out.checksum, 'd4eafcf7beeab093bea441e6f7cab73a')

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
    // CADD SNVs to echtvar
    //
    if (!params.skip_cadd_snvs) {
        
        if(!params.cadd_snvs_echtvar_zip) {
            // Would convert on the fly...so we can could possibly  check the checksum after?
            // ooor...we check them before.
            WGET_CADD_SNVS (
                ch_cadd_snvs.remote
            )

            ch_cadd_snvs.local
                .mix(WGET_CADD_SNVS.out.remote)
                .set { ch_cadd_snvs_to_convert }
            
            // TODO: Check checksum here...
            // Best solution I think is to download whole file (pretty quick for ~80 GB)
            // Do a MD5SUM-check
            // Convert to VCF and ecthvar archive
            // md5sum-check the ecthvar archive
            // if echtvar-archive is provided, then just md5-sumcheck the echtvar file
            
            // TODO: Rename this
            DOWNLOADCADD (
                ch_cadd_snvs_to_convert
            )

        }

        ch_echtvar_encode_files = ch_echtvar_encode_files
            .mix(
                DOWNLOADCADD.out.bcf
                    .combine(ch_echtvar_encode_cadd_snvs_json)
                    .map { vcf_meta, vcf, _json_meta, json -> [ vcf_meta, vcf, json ] }
            )

    }
    
    //
    // Reformat ClinVar to reference genome and Scout compatibility (1 -> chr1 & ID -> INFO/CLNVID, annotated via VEP)
    //
    if(!params.skip_clinvar) {
        BCFTOOLS_ANNOTATE_RENAME_CLINVAR_CHRS (
            ch_clinvar_vcf.map { meta, vcf -> [ meta, vcf, [], [], [] ] },
            [],
            ch_clinvar_rename_chrs
        )

        VCFEXPRESS_ADD_CLNVID (
            BCFTOOLS_ANNOTATE_RENAME_CLINVAR_CHRS.out.vcf,
            ch_vcfexpress_add_clnvid_prelude
        )
        
        TABIX_TABIX (
            VCFEXPRESS_ADD_CLNVID.out.vcf
        )
    }

    // CoLoRSdb SNVs
    ch_echtvar_encode_files = ch_echtvar_encode_files
        .mix(
            ch_colorsdb_snvs
                .combine(ch_echtvar_encode_colorsdb_snvs_json)
                .map { vcf_meta, vcf, _json_meta, json -> [ vcf_meta, vcf, json ] }
            )
    
    //
    // Echtvar SNV databases
    //  
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
    // TODO: Do a MD5-sum check of echtvar encode, especially CADD SNVs if params.cadd_snvs_ecthvar_zip has been set

        
    
    //
    // Download GnomAD and strip info/convert on the fly
    //
    if(!params.skip_gnomad_snvs) {
        DOWNLOAD_GNOMAD_SNVS (
            ch_gnomad_snvs_per_chr
        )

        DOWNLOAD_GNOMAD_SNVS.out.bcf
            .join(DOWNLOAD_GNOMAD_SNVS.out.csi)
            .map { _meta, bcf, csi -> [ [ 'id': 'gnomad_snvs' ], bcf, csi ] }
            .groupTuple()
            .set { ch_bcftools_concat_gnomad }

        BCFTOOLS_CONCAT (
            ch_bcftools_concat_gnomad
        )
    }

    // Download GnomAD SVS and strip info/convert on the fly
    if(!params.skip_gnomad_svs) {
        DOWNLOAD_GNOMAD_SVS (
            ch_gnomad_svs
        )
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

// Creates a value channel with meta, input to download processes
def downloadChannelOf(path) {
    // Channel.of(file(path, checkifexists: true)) ?
    Channel.of(path)
        .map { it -> [ [ 'id': file(it).name ], it ] }
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

def assertChecksum(channel, expected_checksum) {
    channel
        .map { _meta, file -> 
            def checksum = file.readLines().collect { line -> line.split('  ')[0] }[0]
            def filename = file.readLines().collect { line -> line.split('  ')[1] }[0]
            assert checksum == expected_checksum : "Checksum for ${filename} was '${checksum}' but expected '${expected_checksum}'"
        }
}