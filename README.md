# SNP genotyping - **REPO STILL UNDER CONSTRUCTION**
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
3. paleomix
4. SAMtools
5. BCFtools
6. VCFtools
7. R
8. htslib
9. Plink

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

Also note that you can check that the Adapter Content. These we will remove in the following section

if you have a lot of samples and wnat to summerise the output from your fastqc reports we can run multiqc which creates a single summary file
```
cd /folder/containing/fastqc/reports
multiqc .
```
## Read alignment with [paleomix](https://paleomix.readthedocs.io/en/stable/bam_pipeline/index.html) - [BAM pipeline](https://paleomix.readthedocs.io/en/stable/bam_pipeline/index.html)
I highly recommand PALEOMIX to make your life easier and streamline your analyses when you have many individuals. It's been designed to work with ancietn DNA but also works great for modern samples. It caries out mutiple analytical steps and generates useful summary statistics. Another useful feature is that it allows you align reads to multiple genomes, and create seperate bam files. This allows you to obtain seperate bam files file the nulclear and mitochodnrial genomes and prevent unwanted misaligned reads to your nuclear genome. Paleomix uses a range of well known bioinformatic tools (i.e. AdapterRemoval, BWA, picard, samtools). For more information please have a look at their page which containes all you need to know.
Paleomix inclues the following steps:
* create prefixes for your reference
* quality trimming and removing adapters
* Read alignment 
* remove duplicates
* local realignment
To run paleomix you first need to create a makefile (YAML) which contains all the parameters for the different programs and provide the paths to your reference genome(s) and read files. I've provided an example makefile (example.yaml) and a detailed description of component can be foudn [here](https://paleomix.readthedocs.io/en/stable/bam_pipeline/makefile.html). To trim off any adapter sequences you need to which protocol was used to create the library. Usesually this will from Illumina and you can find the correct adapter sequences for each library preperation protol [here](https://knowledge.illumina.com/library-preparation/general/library-preparation-general-reference_material-list/000001314).
Ones you have your makefile you run it like:
```
paleomix bam_pipeline run PATH/TO/MAKEFILE.yaml
```
I also prefer to remove all the "clipped" from my alignments. These are reads where only part of the sequence mapped to the genome and the rest is masked (soft-clipped), or removed (hard-clipped) from the BAM file. I remove these reads using samtools
```
# create list of all clipped reads
samtools view paleomix_output.bam |awk '$6 ~ /H|S/ {print $1}' |sort -u > clipped_reads_list.txt
# filter clipped reads
samtools view -h paleomix_output.bam | fgrep -wvf clipped_reads_list.txt | samtools view -b - -o clipped_alignment.bam
# create index for new bam file
samtools index clipped_alignment.bam
```
Paleomix creates A LOT of temporary files you no longer need ones you're happy the analyses was performed correctly. I would recommend keeping the final alignment (BAM) files for the nuclear and if provided the mitochondrial genome, the makefile (YAML) that contains all the parameters used to generated the alignmetns, and the summary.txt file.
Ones you have generated alignmeng (BAM) files for all you individuals you can move forward to genotyping.

# STEP 2: Genotyping
With genotyping you take all the alignments and determine where in the genome in mutation have arisen, and save these in a "variant call format" or VCF file. I can't stress this enough, take the time to [understand](https://en.wikipedia.org/wiki/Variant_Call_Format) how that data is structured in a VCF file. To save (a lot!) of time you can run an array script to perform genotyping for each scafoold/chromosome in parallel. If you havce a good assembly I would only genoytpe for the scaffolds that represent the largest fraction of your genome. For example, if you have 6000 scaffolds in your reference genome but 95% of it is contained in the 25 largest scaffolds, I would only perform genotyping on those first 25 scaffolds. After genotyping with bcftools we will use the same program to generate tags for multiple fields.
**!!!NOTE!!!**
You will see me refer to PROJECT and SET in the following scripts. This is so that you can easily perform the same analyses for multiple species (PROJECT), and groups of indiviauls (SET).
**!!!NOTE!!!**
```
#!/bin/bash
#SBATCH --array 1-25
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH --partition=parallel
#SBATCH --time=2-0:00
#SBATCH --job-name=bcftools_mpileup

###!!!###                                       ###!!!####
# determine the number of scaffolds you want to genotype #
#               change --array accordingly               #
###!!!###                                       ###!!!####

#load modules
module load htslib/1.9
module load bcftools/1.10.1

#variables
SCAFFOLD=${SLURM_ARRAY_TASK_ID}
PROJECT=$1
SET_NEW=$PROJECT'_'$2
TMP_DIR=$SCRATCH/projects/$PROJECT/data/snp/$SET_NEW/tmp
mkdir -p $TMP_DIR

# Path to your reference genome
REF=$SCRATCH/projects/$PROJECT/resources/reference_genomes/nuclear/Chrysophrys_auratus.v.1.0.all.assembly.units.fasta

# List with paths to all the bam (individual) you want to analyse.
BAMLIST=$SCRATCH/projects/$PROJECT/resources/bam_lists/$SET_NEW'_bam.list'

# Obtain the scaffold name for the ith scaffold from your reference. A fai file should have been created when paoleomix indexed your genome.
REGION=$( head -n $SCAFFOLD $REF.fai | tail -n 1 | cut -f 1 )

#genotype
bcftools mpileup -Ov                                                              \
                 -a 'FORMAT/AD,FORMAT/DP,FORMAT/SP,FORMAT/ADF,FORMAT/ADR,INFO/AD' \
                 -f $REF                                                          \
                 -r $REGION                                                       \
                 -b $BAMLIST                                                      |
bcftools call -Ov -mv > $TMP_DIR/$REGION'_'$SET_NEW'_raw_tmp1.vcf'

#update INFO fields
bcftools +fill-tags $TMP_DIR/$REGION'_'$SET_NEW'_raw_tmp1.vcf'    \
        -Oz -o      $TMP_DIR/$REGION'_'$SET_NEW'_raw_tmp2.vcf.gz' \
         -- -t AC,AF,AN,MAF,NS,AC_Hom,AC_Het			
#compress output to reduce file size
bgzip --reindex $TMP_DIR/$REGION'_'$SET_NEW'_raw_tmp2.vcf.gz'

#clean up
rm $TMP_DIR/$REGION'_'$SET_NEW'_raw_tmp1.vcf'
```
Now you generated compressed vcf (vcf.gz) files for each scaffold which we want to merge into a single vcf.gz file containing all "raw" SNPS. With raw I mean that we have not perfermed any filtering of low quality genotypes.
```
#!/bin/bash
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --partition=bigmem
#SBATCH --time=1-0:00
#SBATCH --job-name=bcf_concat

#load modules
module load htslib/1.9
module load vcftools/0.1.16
module load bcftools/1.10.1

#variables
PROJECT=$1
SET_NEW=$PROJECT'_'$2

# Path where all vcf files are located
TMP_DIR=$SCRATCH/projects/$PROJECT/data/snp/$SET_NEW/tmp

# Extension to new VCF file
VCF_NEW=$SCRATCH/projects/$PROJECT/data/snp/$SET_NEW/$SET_NEW'_raw'

#merge vcf files in tmp dir
bcftools concat -Oz -o $VCF_NEW'.vcf.gz' $( ls -v $TMP_DIR/*'_raw_tmp2.vcf.gz' ) --threads 10

#create indexes for raw vcf file
bgzip --reindex $VCF_NEW'.vcf.gz'
tabix -p vcf $VCF_NEW'.vcf.gz'
```
Now that we have a file containing all genotypes we can start filtering and work toward vcf files that we can use to perform population genomic analyses.

# STEP 3: SNP filtering 
We will perform multiple filter steps and generate multiple VCF files. Different SNPs data sets can be used for different types of analyses.
After all the filter steps we will have the following data sets:
1. raw - vcf file we already have containing all genotypes
2. qc - vcf file containing all high-quality SNPs that have passed innitial quality check (qc)
3. outlier - vcf file containing outlier SNPs that show significnat signs of selection
4. neutral -  vcf file containing independantly segregating SNPs that show no sign of selection

For most population genomic analyses you will be using the neutral SNP dataset but the qc and outlier data sets can also provide unique insights. We also first need to identify high-quality SNPs and outlier SNPs before we can get our nuetral dataset, so there's really no reason not to generate them.
Depending on how many SNPs you have in your raw.vcf, the initial filter steps can take a lot of time if you filter the chromosomes in series. Just like we did for genotyping we will run an array script that performs these innitial filter steps per scaffold/linkage group/chromosome of your data. It may seem a little redundant to first merge all the seperate raw.vcf files into a sigle vcf and then split them up again for SNPs filtereing. But I like that after genotyping we and up with a single file and we delate all other intermediate files, making your data management much easier!

#### QC filter
Here we remove low quality genotypes and and select only biallilic SNPs, and involves the following steps:
1. subsample vcf
2. set individual genotypes with low coverage to missing
3. remove low quality SNPs 
4. test for allelic imbalance (custom Rscript)
5. remove SNPs with allelic imbalance
```
#!/bin/bash
#SBATCH --array 1-25
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=parallel
#SBATCH --time=2-0:00
#SBATCH --job-name=QC_filtering

###!!!###                                       ###!!!####
# determine the number of scaffolds you want to genotype #
#               change --array accordingly               #
###!!!###                                       ###!!!####

###run input
PROJECT=$1
SET=$PROJECT'_'$2
N=${SLURM_ARRAY_TASK_ID}

###load packages
module load htslib/1.9
module load vcftools/0.1.16
module load bcftools/1.10.1
module load R/4.0.2

#Rscript paths
AB_script=PATH/TO/allelelic_imbalance_4.0.R

###resrouces
REF=$SCRATCH/projects/$PROJECT/resources/reference_genomes/nuclear/Chrysophrys_auratus.v.1.0.all.assembly.units.fasta
AB_exclude=$SCRATCH/projects/$PROJECT/resources/sample_info/high_coverage_samples.list #only if innital tests for allelic imbalance suggest removal of certain individuals.

###Obtain the scaffold name for the ith scaffold from your reference. A fai file should have been created when paoleomix indexed your genome.
LG=$( head -n $N $REF.fai | tail -n 1 | cut -f 1 )

###set paths
VCF=$SCRATCH/projects/$PROJECT/data/snp/$SET/$SET
DIR=$SCRATCH/projects/$PROJECT/data/snp/$SET/tmp
mkdir -p $DIR
TMP=$DIR/$LG'_'$SET

##1## subsample vcf
bcftools view -Oz -o $TMP'_tmp1.vcf.gz' $VCF'_raw.vcf.gz' -r $LG

##2## filter genotyeps with DP < 3
vcftools --gzvcf $TMP'_tmp1.vcf.gz' \
         --out 	$TMP'_tmp2'         \
         --minDP 3                  \
         --remove-indels            \
         --recode-INFO-all --recode
mv $TMP'_tmp2.recode.vcf' $TMP'_tmp2.vcf'
bgzip $TMP'_tmp2.vcf'

##3## basic filter parameters
vcftools --gzvcf $TMP'_tmp2.vcf.gz'	          \
         --out $TMP'_tmp3' --max-missing 0.95 \
         --min-alleles 2   --max-alleles 2    \
         --min-meanDP  8   --max-meanDP 25    \
         --minQ 600        --maf 0.01         \
         --recode-INFO-all 	  --recode
mv $TMP'_tmp3.recode.vcf' $TMP'_tmp3.vcf'
bgzip $TMP'_tmp3.vcf' 

##4## test for allelic imbalance
#Select output from VCF (genotypes)
vcftools --gzvcf $TMP'_tmp3.vcf.gz' \
         --out 	$TMP'_tmp4'         \
         --extract-FORMAT-info GT 	
#Select output from VCF (allelic depth)
vcftools --gzvcf $TMP'_tmp3.vcf.gz' \
         --out 	$TMP'_tmp4'         \
         --extract-FORMAT-info AD
#run binomial test to filter sites with allelic imbalance - could require high mem when many SNPs are to be analysed
Rscript $AB_script --GT_file  $TMP'_tmp4.GT.FORMAT' \
                   --AD_file  $TMP'_tmp4.AD.FORMAT' \
                   --out_file $TMP'_tmp4_qc'        \
                   --conf.level 0.99                \
                   --plots TRUE                     \
                   --remove $AB_exclude	### make a list of sample names you want to exclude
##5## remove sites suffering of allelic imbalance
vcftools --gzvcf $TMP'_tmp3.vcf.gz' \
         --out 	$TMP'_tmp5_qc'      \
         --recode-INFO-all --recode \
         --exclude-positions $TMP'_tmp4_qc.exclude_pval0.01.list'
mv $TMP'_tmp5_qc.recode.vcf' $TMP'_tmp5_qc.vcf'
bgzip -fi $TMP'_tmp5_qc.vcf'
tabix -fp vcf $TMP'_tmp5_qc.vcf.gz'
```
Now we just have to merge all the qc vcf files in to a single vcf file. For a lot of downstream analyses you will want to load your SNPs data into R. I've become a big fan of the package [SNPRelate](https://www.bioconductor.org/packages/release/bioc/html/SNPRelate.html) from bioconductor for analysing SNP data in R. I wrote a small Rscript that converts a vcf.gz file to a gds file required to use SNPRelate. Large data sets can be loaded into R really quick with SNPRelate and the gds format. My population genomics repository also uses this package extensively. Similarly, I wrote a small Rscript to convert a gds file to plink (bed) format. To do this you can also supply a tab-delimited population file with two columns 1) Sample name 2) population name.     

```
#!/bin/bash
#SBATCH --cpus-per-task=10
#SBATCH --mem=40G
#SBATCH --partition=bigmem
#SBATCH --time=0-5:00
#SBATCH --job-name=merge_vcf_bcftools

###run input
PROJECT=$1
SET=$PROJECT'_'$2

###load packages
module load htslib/1.9
module load bcftools/1.10.1
module load R/4.0.2

#custom scripts
vcf2R=PATH/TO/vcf2Rinput.R
gds2plink=PATH/TO/gds2plink.R

#tab delimited popualtion file
POP=PATH/TO/$PROJECT'_pop_info.tsv'

###set paths
VCF=$SCRATCH/projects/$PROJECT/data/snp/$SET/$SET'_qc'
TMP=$SCRATCH/projects/$PROJECT/data/snp/$SET/tmp

#merge vcf files in tmp dir
bcftools concat -Oz -o $VCF.vcf.gz $( ls -v $TMP/*'_tmp5_qc.vcf.gz' ) --threads 10
bgzip --reindex $VCF.vcf.gz
tabix -p vcf $VCF.vcf.gz

#create gds
Rscript $vcf2R --gzvcf $VCF.vcf.gz  \
               --snprelate_out $VCF	

#create plink bed file
Rscript $gds2plink --gds_file $VCF.gds \
                   --out      $VCF     \
                   --pop_file $POP

#clean up - remove the # if you're sure you dont need the temporary output
#rm -r $TMP
```
#### Identify outlier SNPs
Now that we have a vcf file containing high-quality SNPs we can identify outlier loci. There are many programs that that can be used to identify outliers. Here we'll use [pcadapt](https://bcm-uga.github.io/pcadapt/articles/pcadapt.html). I like using pcadapt because it doesn't make any assumptions regarding population structure, and tends to be the least conservative. At this stage I like removing anything that could potentially be an outlier to obtain a neutral dataset. 
If you want to use outlier identification to make inferences regarding selection, I recommand combining resutls from multiple outlier analyses. For example the program [OutFlank](http://rstudio-pubs-static.s3.amazonaws.com/305384_9aee1c1046394fb9bd8e449453d72847.html) which is an FST based approach and tends to be much more conservative. Pcadapt (and OutFlank) also runs well with large WGS datasets. More "oldschool" applications inteded for smaller datasets simply can't handle large datasets or take ages to run. The next piece of code shows how to run the shell script which contains a custom R script to run pcadapt. 
See the pcadapt directory for explenation what the script does and which R packages need to be isntalled before you can run this script. This input for pcadapt in plink (bed) format which we have generated in the previous script, and if you provided a --pop_file you can already get some insights into possible population structure. If you prefer to perform outlier identification on your desktop I have also provided a markdown (.Rmd) version of this analyses.
```
#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH --partition=bigmem
#SBATCH --time=2-0:00
#SBATCH --job-name=outlier_selection

###run input
PROJECT=$1
SET=$PROJECT'_'$2

###load packages
module load htslib/1.9
module load vcftools/0.1.16
module load bcftools/1.10.1
module load R/4.0.2
module load plink/1.90

#Rscript paths
R=$SCRATCH/scripts/population_genomics/4_variant_filtering/R
pcadapt=$R/pcadapt_shell.R
vcf2R=$R/vcf2Rinput.R

###resrouces
pop_file=$SCRATCH/projects/$PROJECT/resources/sample_info/$PROJECT'_pop_info.tsv'

###set paths
VCF=$SCRATCH/projects/$PROJECT/data/snp/$SET/$SET
PCADAPT=$SCRATCH/projects/$PROJECT/output/$SET/outlier_analyses/pcadapt

### create directories
#mkdir $TMP
mkdir -p $PCADAPT

##################################### filter outlier loci  #############################
QVAL=0.05
K=1
Rscript $pcadapt --plink   $VCF'_qc'      \
                 --out     $PCADAPT/$SET  \
                 --K0      5              \
                 --K       $K             \
                 --maf     0.05           \
                 --q       $QVAL          \
                 --slw     50000          \
                 --minNsnp 5              \
                 --mode    full

#filter vcf file for outliers
if [ -e $PCADAPT/$SET*'_K'$K'_q'$QVAL'_all_outliers.tsv' ]
then
 #extract SNP info from outlier list
 tail -n +2 $PCADAPT/$SET*'_K'$K'_q'$QVAL'_all_outliers.tsv' |  awk '{print $2}' | sed -e 's/:/\t/p' > $PCADAPT/$SET'_all_outlier_LOC.tsv'
 #outputting outlier vcf file containing all snps
 vcftools --gzvcf $VCF'_qc.vcf.gz'                        \
          --positions $PCADAPT/$SET'_all_outlier_LOC.tsv' \
          --out $VCF'_outliers'                           \
          --recode-INFO-all                               \
          --recode
 mv $VCF'_outliers.recode.vcf' $VCF'_outliers.vcf'
 bgzip -i $VCF'_outliers.vcf' 

 Rscript $vcf2R --gzvcf $VCF'_outliers.vcf.gz' \
		--snprelate_out $VCF'_outliers'	

fi			
```
Almost there! Now that we have identified potential outlier loci we can remove them from our high-quality SNP data set. Finally, we will remove SNPs not in Hardy-Weinberg equilibrium (HWA), low frequency alleles (maf), and perform thinning and linkage disequilibrium (LD) pruning to obtain independently segregating loci. If you don't have any information on how quick linkage disequilibrium decays you have a look at the LD directory and run the LDdecay script. Thinning drastically reduces the dataset but could still contain a number of linked loci. To weed out the last SNPs in linkage we run plink to identify any remove these. The final step is to convert the vcf.gz file to .gds .genlight and .bed (plink) formats which will allow you to perform a range of population genomics analsyes.
```
#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --partition=quicktest
#SBATCH --time=0-5:00
#SBATCH --job-name=filter_neutral

###run input
PROJECT=$1
SET=$PROJECT'_'$2

###load packages
module load htslib/1.9
module load vcftools/0.1.16
module load bcftools/1.10.1
module load R/4.0.2
module load plink/1.90

#Rscript paths
R=$SCRATCH/scripts/population_genomics/4_variant_filtering/R
vcf2R=$R/vcf2Rinput.R
gds2plink=$R/gds2plink.R

###resrouces
pop_file=$SCRATCH/projects/$PROJECT/resources/sample_info/$PROJECT'_pop_info.tsv'

###set paths
VCF=$SCRATCH/projects/$PROJECT/data/snp/$SET/$SET
TMP=$SCRATCH/projects/$PROJECT/data/snp/$SET/tmp
PCADAPT=$SCRATCH/projects/$PROJECT/output/$SET/outlier_analyses/pcadapt

### create directories
mkdir $TMP
################################# filter neutral loci  ###################################
#parameters you used to for outlier detection using pcadapt to remove the correct outliers
QVAL=0.05
K=1
#remove outlier loci
if [ -e $PCADAPT/$SET*'_K'$K'_q'$QVAL'_all_outliers.tsv' ]
then
 tail -n +2 $PCADAPT/$SET*'_K'$K'_q'$QVAL'_all_outliers.tsv' |  awk '{print $2}' | sed -e 's/:/\t/p' > $PCADAPT/$SET'_all_outlier_LOC.tsv'
 vcftools --gzvcf $VCF'_qc.vcf.gz'                                \
          --out $TMP/$SET'_tmp6'                                  \
          --exclude-positions $PCADAPT/$SET'_all_outlier_LOC.tsv' \
          --hwe 0.05 --maf 0.05 --thin 10000                      \
          --recode-INFO-all --recode
else
 vcftools --gzvcf $VCF'_qc.vcf.gz'          \
          --out $TMP/$SET'_tmp6'            \
          --hwe 0.05 --maf 0.05 --thin 5000 \
          --recode-INFO-all --recode		
fi
bgzip -ci $TMP/$SET'_tmp6.recode.vcf' > $TMP/$SET'_tmp6.vcf.gz'
get vcf.gz
Rscript $vcf2R --gzvcf $TMP/$SET'_tmp6.vcf.gz' 	\
               --snprelate_out $TMP/$SET'_tmp6'
#get gds
Rscript $gds2plink --gds_file $TMP/$SET'_tmp6.gds' \
                   --out $TMP/$SET'_tmp6'          \
                   --pop_file $pop_file
#filter linked sites
plink --bfile $TMP/$SET'_tmp6'   \
      --indep-pairwise 50 5 0.2  \
      --out $TMP/$SET'_tmp6'
sed -i 's/:/\t/g' $TMP/$SET'_tmp6.prune.in'
# extract independant SNPs
vcftools --gzvcf $TMP/$SET'_tmp6.vcf.gz'       \
         --out $VCF'_neutral'                  \
         --positions $TMP/$SET'_tmp6.prune.in' \
         --recode-INFO-all --recode
#convert to vcf.gz
mv $VCF'_neutral.recode.vcf' $VCF'_neutral.vcf'
bgzip -i $VCF'_neutral.vcf'
#convert to gds
Rscript $vcf2R --gzvcf $VCF'_neutral.vcf.gz' 	\
               --snprelate_out $VCF'_neutral'	\
               --genlight_out  $VCF'_neutral'
#convert to plink					
Rscript $gds2plink --gds_file $VCF'_neutral.gds' \
                   --out $VCF'_neutral'          \
                   --pop_file $pop_file
#clean up
rm -r $TMP/
```