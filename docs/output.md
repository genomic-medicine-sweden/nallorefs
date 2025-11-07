# genomic-medicine-sweden/nallorefs: Output

## Introduction

This document describes the output produced by the pipeline.

The directories listed and files below will be created in the results directory after the pipeline has finished.

```
├── CADD-scripts
│   └── data
│       └── annotations
├── CoLoRSdb.GRCh38.v1.2.0.deepvariant.glnexus.vcf.gz.zip
├── CoLoRSdb.GRCh38.v1.2.0.pbsv.jasmine.vcf.gz
├── CoLoRSdb.GRCh38.v1.2.0.pbsv.jasmine.vcf.gz.tbi

├── gnomad_snvs.zip
├── gnomad.v4.1.sv.sites.no_cnv.vcf.gz
├── grch38_chromosomes_split_at_centromeres_-v1.0-.bed
├── grch38_clinvar_20250217_renamed_reformatted.vcf.gz
├── grch38_clinvar_20250217_renamed_reformatted.vcf.gz.tbi
├── GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta
├── grch38_hificnv_excluded_regions_common_50_-v1.0-.bed.gz
├── grch38_hificnv_expected_copynumer_xx_-v1.0-.bed
├── grch38_hificnv_expected_copynumer_xy_-v1.0-.bed
├── grch38_par_-v1.0-.bed
├── grch38_rank_model_snvs_-v1.0-.ini
├── grch38_rank_model_svs_-v1.0-.ini
├── grch38_reduced_penetrance_-v1.0-.tsv
├── grch38_trgt_pathogenic_repeats.bed
├── grch38_variant_consequences_-v1.0-.txt
├── grch38_vep_112_loftool_scores_-v1.0-.txt
├── grch38_vep_112_pli_values_-v1.0-.txt
├── md5sum
│   ├── annotationsGRCh38_v1.6.tar.gz.md5
│   ├── clinvar_20250217.vcf.gz.md5
│   └── whole_genome_SNVs.tsv.gz.md5
├── pipeline_info
├── prescored
│   └── GRCh38_v1.6
│       └── no_anno
├── sites.hg38.vcf.gz
├── variant_catalog_grch38.json
├── vep_cache
│   └── homo_sapiens_merged
│       └── 110_GRCh38
└── whole_genome_SNVs.tsv.gz.zip
```

### CADD-scripts

To calcuate CADD scores for small indels, input with `--cadd_resources` to nallo. Equivalent of `data/annotations/` folder described [here](https://github.com/kircherlab/CADD-scripts/#manual-installation).

### Echtvar databases

```
CoLoRSdb.GRCh38.v1.2.0.deepvariant.glnexus.vcf.gz.zip
gnomad_snvs.zip
whole_genome_SNVs.tsv.gz.zip
```

May be input with `--echtvar_snv_databases` to nallo, as described [here](https://genomic-medicine-sweden.github.io/nallo/latest/usage/#snv-annotation).

### SVDB databases

```
CoLoRSdb.GRCh38.v1.2.0.pbsv.jasmine.vcf.gz
gnomad.v4.1.sv.sites.no_cnv.vcf.gz
```

May be input with `--svdb_sv_databases` to nallo, as described [here](https://genomic-medicine-sweden.github.io/nallo/latest/usage/#sv-annotation).

### Target regions

```
grch38_chromosomes_split_at_centromeres_-v1.0-.bed
```

May be input with `--target_regions`.

### VEP files

```
grch38_clinvar_20250217_renamed_reformatted.vcf.gz
grch38_clinvar_20250217_renamed_reformatted.vcf.gz.tbi
grch38_vep_112_loftool_scores_-v1.0-.txt
grch38_vep_112_pli_values_-v1.0-.txt
```

May be input to VEP with `--vep_plugin_files`, as described [here](https://genomic-medicine-sweden.github.io/nallo/latest/usage/#snv-annotation).
If using `grch38_vep_112_loftool_scores_-v1.0-.txt` and `grch38_vep_112_pli_values_-v1.0-.txt` you probably need to update your config file for `ENSEMBLVEP_SNV` and `ÈNSEMBLVEP_SV`, becuase Nallo by default expects them to be named `pLI_values.txt` and `LoFtool_scores.txt`. You can find an example of this under the `vep_cache` section below.

### Reference fasta

```
GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta
```

May be used as reference genome with `--fasta`.

### CNV-calling files

```
grch38_hificnv_excluded_regions_common_50_-v1.0-.bed.gz
grch38_hificnv_expected_copynumer_xx_-v1.0-.bed
grch38_hificnv_expected_copynumer_xy_-v1.0-.bed
```

May be input to HifiCNV as described [here](https://genomic-medicine-sweden.github.io/nallo/latest/usage/#cnv-calling).

### PARs

```
grch38_par_-v1.0-.bed
```

May be input with `--par_regions` as decribed [here](https://genomic-medicine-sweden.github.io/nallo/latest/usage/#snv-calling)

### Ranking files

```
grch38_rank_model_snvs_-v1.0-.ini
grch38_rank_model_svs_-v1.0-.ini
grch38_reduced_penetrance_-v1.0-.tsv
```

May be input to genmod with `--genmod_score_config_snvs`, `--genmod_score_config_svs` and `--genmod_reduced_penetrance`, as described [here](https://genomic-medicine-sweden.github.io/nallo/latest/usage/#rank-snvs-and-indels) and [here](https://genomic-medicine-sweden.github.io/nallo/latest/usage/#sv-annotation)

> [!NOTE]
The rank models downloaded with this pipeline currently assumes you are using a local [loqusdb](https://github.com/Clinical-Genomics/loqusdb) database. To run without this, remove the `[loqusdb]` entries in the rank models.

### Repeat calling files

```
grch38_trgt_pathogenic_repeats.bed
```
May be input to TRGT with `--trgt_repeats` as described [here](https://genomic-medicine-sweden.github.io/nallo/latest/usage/#repeat-calling).

### Annotation files

```
grch38_variant_consequences_-v1.0-.txt
```

May be input to both `--variant_consequences_snvs` and `--variant_consequences_svs` as described [here](https://genomic-medicine-sweden.github.io/nallo/latest/usage/#snv-annotation) and [here](https://genomic-medicine-sweden.github.io/nallo/latest/usage/#sv-annotation).

### md5sum

```
md5sum
```

Contains md5sums for the input files of CADD-annotations, CADD SNVs and Clinvar.

### CADD prescored indels

```
prescored
```

To give prescored CADD scores for small indels, input with `--cadd_prescored_indels` to nallo. Equivalent of `data/prescored/` folder described [here](https://github.com/kircherlab/CADD-scripts/#manual-installation).

### Somalier sites

```
sites.hg38.vcf.gz
```

May be given to somalier with `--somalier_sites` as described [here](https://genomic-medicine-sweden.github.io/nallo/latest/usage/#alignment)

### Repeat annotation

```
variant_catalog_grch38.json
```

May be given to stranger with `--stranger_repeat_catalog`, as described [here](https://genomic-medicine-sweden.github.io/nallo/latest/usage/#repeat-annotation).

### vep_cache

```
vep_cache
```

Contains merged vep cache input with `--vep_cache`, as described [here](https://genomic-medicine-sweden.github.io/nallo/latest/usage/#snv-annotation). By default this pipeline downloads the merged cache, so you may need to update your config file for `ENSEMBLVEP_SNV` and `ENSEMBLVEP_SV`, becuase Nallo by default does not expect the merged cache.

The config could for example look like this:

```
process {
    withName: '.*:ANNOTATE_SVS:ENSEMBLVEP_SV' {
        ext.args = [
          '--dir_plugins .',
          '--plugin pLI,grch38_vep_112_pli_values_-v1.0-.txt',
          '--distance 5000',
          '--buffer_size 5000',
          '--format vcf --max_sv_size 999999999',
          '--appris --biotype --cache --canonical --ccds --compress_output bgzip',
          '--domains --exclude_predicted --force_overwrite',
          '--hgvs --humdiv --no_progress --numbers',
          '--per_gene',
          '--polyphen p --protein --offline --sift p --regulatory --symbol --tsl',
          '--uniprot --vcf',
          '--no_stats',
          '--merged'
        ].join(' ')
    }
    withName: '.*:ANNOTATE_SNVS:ENSEMBLVEP_SNV' {
        ext.args = [
            '--dir_plugins .',
            '--plugin LoFtool,grch38_vep_112_loftool_scores_-v1.0-.txt',
            '--plugin pLI,grch38_vep_112_pli_values_-v1.0-.txt',
            '--plugin SpliceAI,snv=spliceai_scores.raw.snv.hg38.vcf.gz,indel=spliceai_scores.raw.indel.hg38.vcf.gz',
            '--plugin dbNSFP,grch38_dbNSFP4.5a.gz,GERP++_RS,GERP++_NR,phyloP100way_vertebrate,phastCons100way_vertebrate,REVEL_rankscore,REVEL_score,rs_dbSNP150',
            '--distance 5000',
            '--buffer_size 20000',
            '--format vcf --max_sv_size 999999999',
            '--appris --biotype --cache --canonical --ccds --compress_output bgzip',
            '--domains --exclude_predicted --force_overwrite',
            '--hgvs --humdiv --no_progress --numbers',
            '--polyphen p --protein --offline --regulatory --sift p --symbol --tsl',
            '--uniprot --vcf',
            '--no_stats',
            '--merged',
            '--custom grch38_clinvar_20250217_renamed_reformatted.vcf.gz,CLINVAR,vcf,exact,0,CLNSIG,CLNVID,CLNREVSTAT'
        ].join(' ')
    }
}
```

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
