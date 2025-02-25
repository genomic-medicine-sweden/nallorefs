# genomic-medicine-sweden/nallorefs

[![GitHub Actions CI Status](https://github.com/genomic-medicine-sweden/nallorefs/actions/workflows/ci.yml/badge.svg)](https://github.com/genomic-medicine-sweden/nallorefs/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/genomic-medicine-sweden/nallorefs/actions/workflows/linting.yml/badge.svg)](https://github.com/genomic-medicine-sweden/nallorefs/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/genomic-medicine-sweden/nallorefs)

## Introduction

**genomic-medicine-sweden/nallorefs** is a bioinformatics pipeline that downloads _most_ of a set of references that [genomic-medicine-sweden/nallo](https://github.com/genomic-medicine-sweden/nallo) can use (but is not limited to).

For CADD-resources, CADD SNVs and GnomAD SNVs, these files are so large that while the pipeline may be able to download them, it might be better to download them manually. They can then be input into the pipeline with `--cadd_annotations`, `--cadd_snvs` and `--gnomad_base_path` (which should point to the directory where the gnomad chromosome VCFs are stored). The md5sum will be cheched, and they will be processed by the pipeline.

These files can be downloaded with:
```
# Download CADD annotations 
wget -c https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/annotationsGRCh38_v1.6.tar.gz
# Download CADD SNVs
wget -c https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz
# Download gnomad
mkdir -p gnomad
cd gnomad
echo {1..22} X Y | xargs -d ' ' -n 1 -P 10 -I {} wget -c https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr{}.vcf.bgz
```

##### Obtaining non-public reference files

The rank model downloaded with this pipeline also assumes you will input SpliceAI scores to VEP via `--vep_plugin_files` in Nallo. The precomputed spliceAI scores are available on [basespace](https://basespace.illumina.com/s/otSPW8hnhaZR). Download the raw SNV and INDEL score files, rename them and place them in the reference directory.

```
SpliceAI/spliceai_scores.raw.indel.hg38.vcf.gz
SpliceAI/spliceai_scores.raw.indel.hg38.vcf.gz.tbi
SpliceAI/spliceai_scores.raw.snv.hg38.vcf.gz
SpliceAI/spliceai_scores.raw.snv.hg38.vcf.gz.tbi
```

<!-- TODO nf-core:
   Complete this sentence with a 2-3 sentence summary of what types of data the pipeline ingests, a brief overview of the
   major pipeline sections and the types of output it produces. You're giving an overview to someone new
   to nf-core here, in 15-20 seconds. For an example, see https://github.com/nf-core/rnaseq/blob/master/README.md#introduction
-->

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/contributing/design_guidelines#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

<!-- TODO nf-core: Describe the minimum required steps to execute the pipeline, e.g. how to prepare samplesheets.
     Explain what rows and columns represent. For instance (please edit as appropriate):

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

-->

Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run genomic-medicine-sweden/nallorefs \
   -profile <docker/singularity/.../institute> \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Credits

genomic-medicine-sweden/nallorefs was originally written by Felix Lenner.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use genomic-medicine-sweden/nallorefs for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->



This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
