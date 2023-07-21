# SNP genotyping
This turorial contians my pipeline and will help you perform read alignment, genotyping, and varient filtering for WGS data. Going from high throughput Illimuna reads (fq.gz) to a ready to use SNP dataset (VCF) for population genomic analyses. This tutorial is for people who are just beginning with WGS analyses, but also for people who want to learn how to improve their analyses. 
This tutorial assumes you have/know te following:
* Reference genome for your species
* Paired-End (PE) Illumina short-read data for multiple individuals
* Location information for each of your individuals
* Access to a high performance cluster
* basic UNIX/bash scripting to submit jobs to SLURM (for VUW users on raapoi please see their [documentation](https://vuw-research-computing.github.io/raapoi-docs/))

Shell scripts are included in the tutorial but will need to be modified to suite your own projects. Knowing how array jobs work with SBATCH will make your life a lot easier! 

You will need the following programs installed on your cluster (and load via module load), locally installed, or as an executable jar file. The version mentioned are the ones I used, ands would recommend using the same or newer. This tutorial uses custom Rscript which you can find in scipts directory.
1. fastqc
2. multiqc
2. AdapterRemoval
3. BWA-kit
4. SAMtools
5. BCFtools
6. VCFtools
7. R
8. htslib

This tutorial includes the following steps (program)
### 1) Read alignment
* Quality control of you sequencing data (fastqc)
* Adapter trimming (adapter removal)
* Alignment of reads to a reference genome (bwa)
* Prosessing of alignment files (picard)
### 2) genotyping
* Genotyping (bcftools)
### 3) varient filtering
* Varient filtering (bcftools, vcftools, R)

It's imporant to understand each of the diffent files produced by these analyses
1. [fastq.gz](https://en.wikipedia.org/wiki/FASTQ_format)
2. [SAM/BAM](https://en.wikipedia.org/wiki/SAM_(file_format))
3. [VCF/BCF](https://en.wikipedia.org/wiki/Variant_Call_Format)

# STEP 1: Read alignment
Here we take the the sequencing output you've received from your sequencing provider and map all your data to the reference for each individuals seperately. Sequencing output for each individual consists of two files, read1 and read2. These are likely 150 bp sequences, each sequenced from one end of a DNA fragment. NOTE, your DNA fragment is normally larger than the total 300 bp of the two fragment, and so do not overlap. the read files wil look something like this:
```
read1: CA118001_R1.gq.gz
read2: CA118001_R2.gq.gz
```
## Quality checking ([fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)/[multiqc](https://multiqc.info/))
First lets check whether the data we have is good to use. fastqc analyses the fastq file and generates a HTML summary file. Here is a summary file you could get.
```
fastqc READFILES -t 12 --noextract -o OUTPUTPATH
```
In the example, some quality metrics were flagged for potential issues.
1. Per base sequence content. I've found this to be common in all my reads and it's caused by non-random priming of hexamers. for more information see [Hansen et al., 2010](https://academic.oup.com/nar/article/38/12/e131/2409775)
2. per tile sequence quality. local reduction in sequence quality on the flow cell. Some local reductions is not bad, as long as your "Per base sequence quality" is still good.
3. Overrepresented sequences. Most likely low quality reads that will be filtered out in coming stages

if you have a lot of samples and wnat to summerise the output from your fastqc reports we can run multiqc which creates a single summary file
```
cd /folder/containing/fastqc/reports
multiqc .
```
## Adapter removal ([AdapterRemoval](https://adapterremoval.readthedocs.io/en/stable/))