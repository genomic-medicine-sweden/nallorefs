# genomic-medicine-sweden/nallorefs pipeline parameters

A Nexflow pipeline to download and create genomic-medicine-sweden/nallo references

## Input/output options

Define where the pipeline should find input data and save output data.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `input` | Path to comma-separated file containing information about the samples in the experiment. <details><summary>Help</summary><small>You will need to create a design file with information about the samples in your experiment before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row.</small></details>| `string` | /home/felix/projects/genomic-medicine-sweden/nallorefs/assets/samplesheet.csv | True |  |
| `outdir` | The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure. | `string` |  | True |  |

## Institutional config options

Parameters used to describe centralised config profiles. These should not be edited.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `custom_config_version` | Git commit id for Institutional configs. | `string` | master |  | True |
| `custom_config_base` | Base directory for Institutional configs. <details><summary>Help</summary><small>If you're running offline, Nextflow will not be able to fetch the institutional config files from the internet. If you don't need them, then this is not a problem. If you do need them, you should download the files from the repo and tell Nextflow where to find them with this parameter.</small></details>| `string` | https://raw.githubusercontent.com/nf-core/configs/master |  | True |
| `config_profile_name` | Institutional config name. | `string` |  |  | True |
| `config_profile_description` | Institutional config description. | `string` |  |  | True |
| `config_profile_contact` | Institutional config contact information. | `string` |  |  | True |
| `config_profile_url` | Institutional config URL link. | `string` |  |  | True |

## Generic options

Less common options for the pipeline, typically set in a config file.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `version` | Display version and exit. | `boolean` |  |  | True |
| `publish_dir_mode` | Method used to save pipeline results to output directory. (accepted: `symlink`\|`rellink`\|`link`\|`copy`\|`copyNoFollow`\|`move`) <details><summary>Help</summary><small>The Nextflow `publishDir` option specifies which intermediate files should be saved to the output directory. This option tells the pipeline what method should be used to move these files. See [Nextflow docs](https://www.nextflow.io/docs/latest/process.html#publishdir) for details.</small></details>| `string` | copy |  | True |
| `monochrome_logs` | Do not use coloured log outputs. | `boolean` |  |  | True |
| `validate_params` | Boolean whether to validate parameters against the schema at runtime | `boolean` | True |  | True |
| `pipelines_testdata_base_path` | Base URL or local path to location of pipeline test dataset files | `string` | https://raw.githubusercontent.com/nf-core/test-datasets/ |  | True |
| `trace_report_suffix` | Suffix to add to the trace report filename. Default is the date and time in the format yyyy-MM-dd_HH-mm-ss. | `string` |  |  | True |

## Other parameters

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `base_reference_version` |  | `string` | 2f03fb600de11fc769fff3dce68167556786ce27 |  |  |
| `base_reference_dir` |  | `string` | https://github.com/Clinical-Genomics/reference-files/raw/2f03fb600de11fc769fff3dce68167556786ce27/ |  |  |
| `clinvar_version` |  | `integer` | 20251215 |  |  |
| `clinvar` |  | `string` | https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar_20251215.vcf.gz |  |  |
| `clinvar_md5sum` |  | `string` | 8569c7ded09f17d4704caee5ebb263cd |  |  |
| `vep_cache_version` |  | `integer` | 110 |  |  |
| `vep_cache_type` |  | `string` | merged |  |  |
| `cadd_annotations` |  | `string` | https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh38/GRCh38_v1.7.tar.gz |  |  |
| `cadd_annotations_md5sum` |  | `string` | 205d3e702df3565efb424e2ca80c9d25 |  |  |
| `cadd_snvs` |  | `string` | https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz |  |  |
| `cadd_snvs_md5sum` |  | `string` | 88577a55f1cd519d44e0f415ba248eb9 |  |  |
| `gnomad_version` |  | `number` | 4.1 |  |  |
| `gnomad_base_path` |  | `string` | https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/ |  |  |
| `local_echtvar_databases` |  | `string` |  |  |  |
| `cadd_snvs_echtvar_zip` |  | `string` |  |  |  |
| `skip_general_files` |  | `boolean` |  |  |  |
| `skip_gens_baf_positions` |  | `boolean` |  |  |  |
| `skip_vep_cache` |  | `boolean` |  |  |  |
| `skip_cadd_annotations` |  | `boolean` |  |  |  |
| `skip_cadd_indels` |  | `boolean` |  |  |  |
| `skip_local_echtvar_databases` |  | `boolean` |  |  |  |
| `skip_cadd_snvs` |  | `boolean` |  |  |  |
| `skip_clinvar` |  | `boolean` |  |  |  |
| `skip_gnomad_snvs` |  | `boolean` |  |  |  |
| `dbnsfp` |  | `string` | https://usf.box.com/shared/static/2hzcx5s6p1xui7oen16xqzndfrkt8l9l |  |  |
| `dbnsfp_md5sum` |  | `string` | 7e01ebb7a4dcb679b9b93b30c7ba131e |  |  |
| `dbnsfp_version` |  | `string` | 4.5a |  |  |
| `strdrop_training_set_version` |  | `integer` | 260120 |  |  |
| `skip_dbnsfp` |  | `boolean` |  |  |  |
| `cadd_prescored_indels_tsv` |  | `string` | https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz |  |  |
| `cadd_prescored_indels_tbi` |  | `string` | https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz.tbi |  |  |
