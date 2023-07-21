# SNP genotyping
This turorial contians my pipeline and will help you perform read alignment, genotyping, and varient filtering for WGS data. Going from high throughput Illimuna reads (fq.gz) to a ready to use SNP dataset (VCF) for population genomic analyses. This tutorial is for people who are just beginning with WGS analyses, but also for people who want to learn how to improve their analyses.
This tutorial assumes you have/know te following:
* Reference genome for your species
* Paired-End (PE) Illumina short-read data for multiple individuals
* Location information for each of your individuals
* Access to a high performance cluster
* basic UNIX/bash scripting to submit jobs to SLURM (for VUW users on raapoi please see their [documentation](https://vuw-research-computing.github.io/raapoi-docs/))

This tutorial includes the following steps (program)
## 1) Read alignment
* Quality control of you sequencing data (fastqc)
* Adapter trimming (adapter removal)
* Alignment of reads to a reference genome (bwa)
* Prosessing of alignment files (picard)
## 2) genotyping
* Genotyping (bcftools)
## 3) varient filtering
* Varient filtering (bcftools, vcftools, R)

It's imporant to understand each of the diffent files produced by these analyses
* [fastq.gz](https://en.wikipedia.org/wiki/FASTQ_format)
* [SAM/BAM](https://en.wikipedia.org/wiki/SAM_(file_format))
* [VCF/BCF](https://en.wikipedia.org/wiki/Variant_Call_Format)

# STEP 1: Read alignment
Here we take the the sequencing output you've received from your sequencing provider and map all your data to the reference for each individuals seperately. Sequencing output for each individual consists of two files, read1 and read2. These are likely 150 bp sequences, each sequenced from one end of a DNA fragment. NOTE, your DNA fragment is normally larger than the total 300 bp of the two fragment, and so do not overlap. the read files wil look something like this:
```
read1: CA118001_R1.gq.gz
read2: CA118001_R2.gq.gz
```
## Quality checking (fastqc)

