#######################################################################
# Project: Quercus macrocarpa HybSeq
# Topic: Chloroplast Variant Calling
# Title: GATK pipeline notes
# Author: Kasey Pham
# Notes: GATK documentation found here:
# https://software.broadinstitute.org/gatk/documentation/tooldocs/current/
# https://www.broadinstitute.org/partnerships/education/broade/best-practices-variant-calling-gatk-1
# If you notice a documentation discrepancy, it's probably because
# I used GATK 3.8 and the documentation is for 4.0.
#######################################################################

############
# RAW DATA #
############
# fastq R1 (forwards) & R2 (backwards) reads from HybSeq for each sample
# Reference chloroplast genome without inverted repeats: Q_macrocarpa_MOR-672.fas

########################
# PROGRAM DEPENDENCIES #
########################
# Python 3
# Trimmomatic
# FastQC
# MultiQC
# GATK + Picard (usually bundled together)
# BWA-MEM
# Docker
# htslib
# bcftools

#################
# DATA CLEANING #
#################
# Rationale: Raw reads from sequencing must have some form of quality control imposed
# so that sequencing errors and biases don't skew results.

# Cleaned reads using Trimmomatic (?). Mira did this step, so I'm not 100% sure.
# The cleaning program tosses reads beneath a certain quality threshhold. If only one read
# in a paired set is tossed, the leftover read is an unpaired single-ended read, which I
# call "orphaned" in most of my naming schemes below.

# Assessed quality of trimmed reads in FastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
# Used MultiQC (https://multiqc.info/) to summarize FastQC results for all samples and did a visual
# spotcheck on quality.

# Compiled all samples which passed quality check into the file gatk_samplelist.txt. This is just
# a simple txt file listing the ID of all samples I wanted to include in analysis separated by a carriage return.
# This file is how I automate analyses for all samples and easily control which are advanced from one step to
# the next.

# (If you need more information about this step, I'm happy to send over code. I just have it
# saved in a different file and figured it was possible that you did all of this already.)

#################
# MAPPING READS #
#################
# Rationale: Reads should be mapped to genome for two reasons:
# 1. To weed out off-target reads from bacterial contamination, plastid, or nuclear genes, depending
# on your targeted sequence.
# 2. To align on-target reads to your reference genome for future steps in GATK.

# dir: /mnt/USERS/kpham/oaks-hybseq/gatk/
# Index reference genome for mapping
bwa index Q_macrocarpa_MOR-672.fas
# Use BWA-MEM to map each sample's reads to reference genome
# You can also use bowtie/SNAP/etc. for this step, but BWA-MEM is recommended by GATK Best Practices.
# BWA documentation: http://bio-bwa.sourceforge.net/bwa.shtml
# general command form: bwa mem [reference genome] [forward paired reads] [reverse paired reads] > [mapping alignment file]
### could use any aligner --- you just need a SAM file out

# I use a bash loop to automate the command over multiple samples, which requires that sample list file from above.
# Documentation on bash loops here: https://www.cyberciti.biz/faq/bash-loop-over-file/
# I map paired and unpaired reads (the orphaned ones) separately, since BWA can't map both at once.
### was done for a run of ca 250-300 individuals
while read NAME
    # paired reads
    do /home/kpham/bin/bwa/bwa mem /mnt/USERS/kpham/oaks-hybseq/gatk/Q_macrocarpa_MOR-672.fas /mnt/USERS/kpham/oaks-hybseq/gatk/qced_reads/paired/"$NAME"_R1_paired.fastq /mnt/USERS/kpham/oaks-hybseq/gatk/qced_reads/paired/"$NAME"_R2_paired.fastq > /mnt/USERS/kpham/oaks-hybseq/gatk/bwa/paired/"$NAME"_paired.sam
    # unpaired reads
    /home/kpham/bin/bwa/bwa mem /mnt/USERS/kpham/oaks-hybseq/gatk/Q_macrocarpa_MOR-672.fas /mnt/USERS/kpham/oaks-hybseq/gatk/qced_reads/unpaired/"$NAME"_unpaired.fastq > /mnt/USERS/kpham/oaks-hybseq/gatk/bwa/"$NAME"_unpaired.sam
done < /mnt/USERS/kpham/oaks-hybseq/gatk/qced_reads/gatk_samplelist.txt

# merge paired and unpaired SAM files together
mkdir /mnt/USERS/kpham/oaks-hybseq/gatk/bwa/merged
while read NAME
    do /home/kpham/bin/gatk-4.1.2.0/gatk MergeSamFiles -I /mnt/USERS/kpham/oaks-hybseq/gatk/bwa/paired/"$NAME"_paired.sam -I /mnt/USERS/kpham/oaks-hybseq/gatk/bwa/unpaired/"$NAME"_unpaired.sam -O /mnt/USERS/kpham/oaks-hybseq/gatk/bwa/merged/"$NAME".sam
done < /mnt/USERS/kpham/oaks-hybseq/gatk/qced_reads/gatk_samplelist.txt

#######################
# GATK PRE-PROCESSING #
#######################
# Rationale: Mapped reads are almost ready for calling variants, but we need to take care of some housekeeping
# and reorganizing of the reads before GATK can be run on them efficiently.
# Steps taken from ShinHan Shiu's Bioinformatics course at MSU. See 1108_variant_calling_and_app_exercise_v18.docx

# clean sam files
while read NAME; do /home/kpham/bin/gatk-4.1.2.0/gatk CleanSam -I /mnt/USERS/kpham/oaks-hybseq/gatk/bwa/merged/"$NAME".sam -O /mnt/USERS/kpham/oaks-hybseq/gatk/preproc/1_clean/"$NAME"_clean.sam; done < /mnt/USERS/kpham/oaks-hybseq/gatk/qced_reads/gatk_samplelist.txt

# convert to bam
while read NAME; do /home/kpham/bin/gatk-4.1.2.0/gatk SamFormatConverter -I /mnt/USERS/kpham/oaks-hybseq/gatk/preproc/1_clean/"$NAME"_clean.sam -O /mnt/USERS/kpham/oaks-hybseq/gatk/preproc/2_bam/"$NAME"_clean.bam; done < /mnt/USERS/kpham/oaks-hybseq/gatk/qced_reads/gatk_samplelist.txt

# add read groups
### both Kasey and Nisa have had difficulties getting the headers fixed up for readgroup... because we are doing
# These are identifiers stored for each read in the SAM file. GATK may throw an error if you don't fill these in.
# I wrote custom python scripts to extract read groups from my reads. Your data may be in a different format and
# may not work with my scripts, but feel free to copy and tweak my them if they work for you.
# Read group documentation: https://software.broadinstitute.org/gatk/documentation/article.php?id=6472
python /mnt/USERS/kpham/oaks-hybseq/scripts/get_readgroup.py /mnt/USERS/kpham/oaks-hybseq/gatk/qced_reads/gatk_samplelist.txt /mnt/USERS/kpham/oaks-hybseq/gatk/qced_reads/unpaired/ _unpaired.fastq | python /mnt/USERS/kpham/oaks-hybseq/scripts/make_rg_job.py /mnt/USERS/kpham/oaks-hybseq/gatk/preproc/add_readgroup.job
chmod +755 /mnt/USERS/kpham/oaks-hybseq/gatk/preproc/add_readgroup.job
/mnt/USERS/kpham/oaks-hybseq/gatk/preproc/add_readgroup.job

# sort bam files by genomic position coordinate
while read NAME; do /home/kpham/bin/gatk-4.1.2.0/gatk SortSam -I /mnt/USERS/kpham/oaks-hybseq/gatk/preproc/3_readgroup/"$NAME"_wRG.bam -O /mnt/USERS/kpham/oaks-hybseq/gatk/preproc/4_sort/"$NAME"_sort.bam -SO coordinate; done < /mnt/USERS/kpham/oaks-hybseq/gatk/qced_reads/gatk_samplelist.txt

# mark duplicate reads
while read NAME; do /home/kpham/bin/gatk-4.1.2.0/gatk MarkDuplicates -I /mnt/USERS/kpham/oaks-hybseq/gatk/preproc/4_sort/"$NAME"_sort.bam -O /mnt/USERS/kpham/oaks-hybseq/gatk/preproc/5_markdupes/"$NAME"_dupe.bam -M /mnt/USERS/kpham/oaks-hybseq/gatk/preproc/5_markdupes/metrics/"$NAME"_dupe_metrics.txt; done < /mnt/USERS/kpham/oaks-hybseq/gatk/qced_reads/gatk_samplelist.txt

# index bam files
while read NAME; do /home/kpham/bin/gatk-4.1.2.0/gatk BuildBamIndex -I /mnt/USERS/kpham/oaks-hybseq/gatk/preproc/5_markdupes/"$NAME"_dupe.bam -O /mnt/USERS/kpham/oaks-hybseq/gatk/preproc/6_index/"$NAME".bai; ln -s /mnt/USERS/kpham/oaks-hybseq/gatk/preproc/5_markdupes/"$NAME"_dupe.bam /mnt/USERS/kpham/oaks-hybseq/gatk/preproc/6_index/"$NAME".bam; done < /mnt/USERS/kpham/oaks-hybseq/gatk/qced_reads/gatk_samplelist.txt

# validate bam files (this step is optional but I usually find it wise.)
while read NAME; do /home/kpham/bin/gatk-4.1.2.0/gatk ValidateSamFile -I /mnt/USERS/kpham/oaks-hybseq/gatk/preproc/6_index/"$NAME".bam -M SUMMARY -O /mnt/USERS/kpham/oaks-hybseq/gatk/preproc/7_validation/"$NAME"_summary.txt; done < /mnt/USERS/kpham/oaks-hybseq/gatk/qced_reads/gatk_samplelist.txt
python /mnt/USERS/kpham/oaks-hybseq/scripts/validation_summary.py /mnt/USERS/kpham/oaks-hybseq/gatk/preproc/7_validation/file_list.txt

# GATK Best Practices suggests doing indel realignment, but I don't do it
# because I'm using HaplotypeCaller instead of UnifiedGenotyper
# source: https://gatkforums.broadinstitute.org/gatk/discussion/3151/

####################
# CALLING VARIANTS #
####################
# Rationale: This is the step where GATK actually identifies sequence variants among mapped reads.
# We're actually going to do this in two parts -- first I use the HaplotypeCaller tool to call
# bases and variants within each sample, then I use GenotypeGVCFs to estimate likely SNPs among
# *all* samples.
# This is better than consolidating reads from all samples into one SAM file and trying to make
# the HaplotypeCaller guess which variants are heterozygosity/alignment error and which are
# inter-sample variation.
# It's important to note that my samples were all haploid, so the program tossed any intrasample
# variants. You'll need to specify ploidy level, and I think HaplotypeCaller requires they all be
# the same, though I might be wrong on that.


# index reference genome for GATK
/home/kpham/bin/gatk-4.1.2.0/gatk CreateSequenceDictionary -R /mnt/USERS/kpham/oaks-hybseq/gatk/Q_macrocarpa_MOR-672.fasta -O /mnt/USERS/kpham/oaks-hybseq/gatk/Q_macrocarpa_MOR-672.dict
/home/kpham/bin/samtools-1.9/samtools faidx /mnt/USERS/kpham/oaks-hybseq/gatk/Q_macrocarpa_MOR-672.fasta

# At this point, I realized that GATK doesn't work on our lab computer because we run a newer
# version of Java than the one it was built in. (Yeah, it doesn't have forward compatibility.)
# Andrew installed Docker onto the computer, which allows you to load a virtual image with the right
# version of Java and GATK in which to run your samples.
# See https://software.broadinstitute.org/gatk/documentation/article?id=11090
# Docker shortcuts:
#     - attach: docker attach <container id>
#     - detach: cntrl+p+q
# I imported my working files into the Docker image using the -v flag on the Docker run command.
# If you don't do that, you won't be able to read your own files because that virtual image is
# basically treated like its own isolated Linux system with its own separate directory structure.
# I did not account for all my symbolic links breaking because of the change in directory structure,
# so I had to remake all them them down below. (oops.) You might be able to avoid this by just copying
# files instead of making symlinks, but it depends on your volume of data.

# Make your personal Docker image
docker pull broadinstitute/gatk:4.1.2.0
docker run --name gatk -v /mnt/USERS/kpham/oaks-hybseq/gatk:/gatk/my_data -it broadinstitute/gatk:4.1.2.0
# --name option allows for easy specification of container id

# Make symbolic links to preprocessed files which can be used in Docker.
# The reason I'm making symlinks in the first place is because the SAM index files need to be in the same directory as the
# SAM files and need to have the same naming scheme. This is my fault; a smarter preprocessing naming and storing scheme
# could have prevented this.
while read NAME; do ln -s /gatk/my_data/preproc/5_markdupes/"$NAME"_dupe.bam ./"$NAME".bam; ln -s /gatk/my_data/preproc/6_index/"$NAME".bai .; done < /mnt/USERS/kpham/oaks-hybseq/gatk/qced_reads/gatk_samplelist.txt

# Run HaplotypeCaller for all samples using different sets of parameters
# I wasn't sure which parameters would be the best for my dataset, so I did two different sets:
# 1. Stringent calling using parameters from Robin Buell's Genomics course at MSU
# 2. Non-stringent default calling in GATK
#
# Stringent parameters
while read NAME
    do /gatk/gatk HaplotypeCaller -I /gatk/my_data/preproc/8_docker/"$NAME".bam -O /gatk/my_data/hapcall/stringent/gvcfs/"$NAME".g.vcf -R /gatk/my_data/Q_macrocarpa_MOR-672.fasta -ERC GVCF --standard-min-confidence-threshold-for-calling 30.0 --sample-ploidy 1 --heterozygosity .05 --indel-heterozygosity .001 --min-base-quality-score 15 --pcr-indel-model AGGRESSIVE
done < /gatk/my_data/qced_reads/gatk_samplelist.txt

# Default parameters
while read NAME
    do /gatk/gatk HaplotypeCaller -I /gatk/my_data/preproc/8_docker/"$NAME".bam -O /gatk/my_data/hapcall/default/gvcfs/"$NAME".g.vcf -R /gatk/my_data/Q_macrocarpa_MOR-672.fasta -ERC GVCF --sample-ploidy 1
done < /gatk/my_data/qced_reads/gatk_samplelist.txt

# Spoiler: I ended up using the stringent parameters to be more conservative, but I did not notice any huge differences between
# the two parameter sets in terms of SNPs called. The stringent set called slightly less.

# Make sample list for GenotypeGVCFs tool
# Needs to be in Docker environment so the addresses are right in future analyses!
ls -d /gatk/my_data/hapcall/stringent/gvcfs/*.vcf > /gatk/my_data/hapcall/stringent/stringent_samples.list
ls -d /gatk/my_data/hapcall/default/gvcfs/*.vcf > /gatk/my_data/hapcall/default/default_samples.list

# DETACH DOCKER: cntrl+p+q

# Create job file to run CombineGVCFs on all samples based on sample list
# I made a custom python python script to construct the CombineGVCFs command because
# otherwise I would have had to type every sample name into the command manually. You
# should be able to use this one as is. I document usage of all my scripts in their
# headers, so check that if you use them.
# CombineGVCFs documentation: https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_CombineGVCFs.php
python /mnt/USERS/kpham/oaks-hybseq/scripts/generate_combine_command.py /gatk/gatk /gatk/my_data/Q_macrocarpa_MOR-672.fasta /mnt/USERS/kpham/oaks-hybseq/gatk/hapcall/stringent/stringent_samples.list /gatk/my_data/hapcall/stringent/stringent_cohort.g.vcf > /mnt/USERS/kpham/oaks-hybseq/gatk/hapcall/stringent/run_combine_stringent.job
python /mnt/USERS/kpham/oaks-hybseq/scripts/generate_combine_command.py /gatk/gatk /gatk/my_data/Q_macrocarpa_MOR-672.fasta /mnt/USERS/kpham/oaks-hybseq/gatk/hapcall/default/default_samples.list /gatk/my_data/hapcall/default/default_cohort.g.vcf > /mnt/USERS/kpham/oaks-hybseq/gatk/hapcall/default/run_combine_default.job

# REATTACH DOCKER
docker attach gatk

# Combine sample GVCFs
/gatk/my_data/hapcall/stringent/run_combine_stringent.job
/gatk/my_data/hapcall/default/run_combine_default.job

# Joint genotyping
# Note that you need to specify ploidy here.
/gatk/gatk GenotypeGVCFs -R /gatk/my_data/Q_macrocarpa_MOR-672.fasta -V /gatk/my_data/hapcall/stringent/stringent_cohort.g.vcf --sample-ploidy 1 -O /gatk/my_data/hapcall/stringent/stringent_joint.vcf
/gatk/gatk GenotypeGVCFs -R /gatk/my_data/Q_macrocarpa_MOR-672.fasta -V /gatk/my_data/hapcall/default/default_cohort.g.vcf --sample-ploidy 1 -O /gatk/my_data/hapcall/default/default_joint.vcf

# validate formatting
/gatk/gatk ValidateVariants -R /gatk/my_data/Q_macrocarpa_MOR-672.fasta -V /gatk/my_data/hapcall/stringent/stringent_joint.vcf --validation-type-to-exclude ALL
/gatk/gatk ValidateVariants -R /gatk/my_data/Q_macrocarpa_MOR-672.fasta -V /gatk/my_data/hapcall/default/default_joint.vcf --validation-type-to-exclude ALL

# filter variants based on GenotypeGVCFs-assigned quality scores
# I have notes on what each of these filters mean, but they're in my office.
# Please remind me to send them if I haven't by afternoon 9/23.
#
# For SNPs called with stringent parameters
/gatk/gatk VariantFiltration -R /gatk/my_data/Q_macrocarpa_MOR-672.fasta -V /gatk/my_data/hapcall/stringent/stringent_joint.vcf -O /gatk/my_data/hapcall/stringent/stringent_joint_filtered.vcf -filter "MQ < 40.0" --filter-name "MQ40" -filter "FS > 60.0" --filter-name "FS60" -filter "MQRankSum < -12.5" --filter-name "MQRankSum" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum"
# For SNPs called with default parameters
/gatk/gatk VariantFiltration -R /gatk/my_data/Q_macrocarpa_MOR-672.fasta -V /gatk/my_data/hapcall/default/default_joint.vcf -O /gatk/my_data/hapcall/default/default_joint_filtered.vcf -filter "MQ < 40.0" --filter-name "MQ40" -filter "FS > 60.0" --filter-name "FS60" -filter "MQRankSum < -12.5" --filter-name "MQRankSum" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum"

# extract SNPs (exclude indels) from filtered variants
/gatk/gatk SelectVariants -R /gatk/my_data/Q_macrocarpa_MOR-672.fasta -V /gatk/my_data/hapcall/stringent/stringent_joint_filtered.vcf -O /gatk/my_data/hapcall/stringent/stringent_joint_snps.vcf --select-type-to-include SNP --exclude-filtered true
/gatk/gatk SelectVariants -R /gatk/my_data/Q_macrocarpa_MOR-672.fasta -V /gatk/my_data/hapcall/default/default_joint_filtered.vcf -O /gatk/my_data/hapcall/default/default_joint_snps.vcf --select-type-to-include SNP --exclude-filtered true

# validate variants (Optional. Checks that variant formatting is correct.)
/gatk/gatk ValidateVariants -R /gatk/my_data/Q_macrocarpa_MOR-672.fasta -V /gatk/my_data/hapcall/stringent/stringent_joint_snps.vcf --validation-type-to-exclude ALL
/gatk/gatk ValidateVariants -R /gatk/my_data/Q_macrocarpa_MOR-672.fasta -V /gatk/my_data/hapcall/default/default_joint_snps.vcf --validation-type-to-exclude ALL

# I decided here to just go with the stringent parameter SNP set.

# filter SNPs by % missing data (meaning the sample did not have sufficient reads to make a call at the position)
# <50% missing data
/gatk/gatk SelectVariants -R /gatk/my_data/Q_macrocarpa_MOR-672.fasta -V /gatk/my_data/hapcall/stringent/stringent_joint_filtered.vcf -O /gatk/my_data/hapcall/stringent/snpsets/stringent_joint_snps_50.vcf --select-type-to-include SNP --max-nocall-fraction 0.5 --exclude-filtered true
# <75% missing data
/gatk/gatk SelectVariants -R /gatk/my_data/Q_macrocarpa_MOR-672.fasta -V /gatk/my_data/hapcall/stringent/stringent_joint_filtered.vcf -O /gatk/my_data/hapcall/stringent/snpsets/stringent_joint_snps_75.vcf --select-type-to-include SNP --max-nocall-fraction 0.75 --exclude-filtered true
# <20% missing data
/gatk/gatk SelectVariants -R /gatk/my_data/Q_macrocarpa_MOR-672.fasta -V /gatk/my_data/hapcall/stringent/stringent_joint_filtered.vcf -O /gatk/my_data/hapcall/stringent/snpsets/stringent_joint_snps_20.vcf --select-type-to-include SNP --max-nocall-fraction 0.20 --exclude-filtered true

# EXIT DOCKER: cntrl+p+q
# is there any difference in the number of SNPs retained?
# This works because each line starting with a hashtag in a VCF file delineates a new SNP
grep -v "^#" /mnt/USERS/kpham/oaks-hybseq/gatk/hapcall/stringent/stringent_joint_snps.vcf | wc -l # 881
grep -v "^#" /mnt/USERS/kpham/oaks-hybseq/gatk/hapcall/stringent/snpsets/stringent_joint_snps_50.vcf | wc -l # 867
grep -v "^#" /mnt/USERS/kpham/oaks-hybseq/gatk/hapcall/stringent/snpsets/stringent_joint_snps_75.vcf | wc -l # 879
grep -v "^#" /mnt/USERS/kpham/oaks-hybseq/gatk/hapcall/stringent/snpsets/stringent_joint_snps_20.vcf | wc -l # 445

# Decided to go for the most conservative filtering criterion (<20% missing data) because 445 SNPs seemed
# good enough to work with.

# compress final VCF file because BCFtools is finicky and can't work with the uncompressed versions
/home/kpham/bin/htslib/bgzip -c /mnt/USERS/kpham/oaks-hybseq/gatk/hapcall/stringent/snpsets/stringent_joint_snps_20.vcf > /mnt/USERS/kpham/oaks-hybseq/gatk/hapcall/stringent/snpsets/stringent_joint_snps_20.vcf.gz
# Index final VCF file for BCFtools
/home/kpham/bin/bcftools/bcftools index /mnt/USERS/kpham/oaks-hybseq/gatk/hapcall/stringent/snpsets/stringent_joint_snps_20.vcf.gz

# Convert VCF file to FASTA file using BCFtools consensus.
# This applies the SNPs for each sample to the "background" of the reference genome.
# Later I remove the "background" sequence, leaving only the SNPs. This seems roundabout,
# but the alternative was custom-coding a python script to do this, and I was too lazy.
while read NAME; do /home/kpham/bin/bcftools/bcftools consensus -f /mnt/USERS/kpham/oaks-hybseq/gatk/Q_macrocarpa_MOR-672.fasta -I -M "N" -o /mnt/USERS/kpham/oaks-hybseq/gatk/hapcall/stringent/snpsets/ind_fastas/"$NAME".fas  -s "$NAME" /mnt/USERS/kpham/oaks-hybseq/gatk/hapcall/stringent/snpsets/stringent_joint_snps_20.vcf.gz; done < /mnt/USERS/kpham/oaks-hybseq/gatk/qced_reads/gatk_samplelist.txt

# Rename headers of individual FASTAs
# I needed to do this because BCFtools consensus just applies the reference genome's name
# to the FASTA file, and I want each individual sample to have its own name.
while read NAME; do sed "s/>macrocarpa_|_MOR672_|_USANM/>${NAME}/" /mnt/USERS/kpham/oaks-hybseq/gatk/hapcall/stringent/snpsets/ind_fastas/"$NAME".fas > /mnt/USERS/kpham/oaks-hybseq/gatk/hapcall/stringent/snpsets/ind_fastas_renamed/"$NAME".fas; done < /mnt/USERS/kpham/oaks-hybseq/gatk/qced_reads/gatk_samplelist.txt

# Extract position of SNPs in genomic reference from column 2 of VCF
cut -f 2 /mnt/USERS/kpham/oaks-hybseq/gatk/hapcall/stringent/snpsets/stringent_joint_snps_20.vcf > /mnt/USERS/kpham/oaks-hybseq/gatk/hapcall/stringent/snpsets/stringent_joint_snps_20_pos.txt
# manually removed VCF header lines from stringent_joint_snps_20_pos.txt

# create list of FASTAs to extract from
echo /mnt/USERS/kpham/oaks-hybseq/refs/plastomes/rich_plastomes/individuals/Q_chrysolepis_MOR-415.fas > /mnt/USERS/kpham/oaks-hybseq/gatk/hapcall/stringent/snpsets/ind_fastas_renamed/ind_fastas_list.txt
find /mnt/USERS/kpham/oaks-hybseq/gatk/hapcall/stringent/snpsets/ind_fastas_renamed/*.fas >> /mnt/USERS/kpham/oaks-hybseq/gatk/hapcall/stringent/snpsets/ind_fastas_renamed/ind_fastas_list.txt

# extract SNPs to multi-sample FASTA file
# I used a custom python script for this. It should be applicable as-is.
# Documentation on how to use it is in the body of the code.
python /mnt/USERS/kpham/oaks-hybseq/scripts/make_mult_align.py /mnt/USERS/kpham/oaks-hybseq/gatk/hapcall/stringent/snpsets/stringent_joint_snps_20_pos.txt /mnt/USERS/kpham/oaks-hybseq/gatk/hapcall/stringent/snpsets/ind_fastas_renamed/ind_fastas_list.txt /mnt/USERS/kpham/oaks-hybseq/gatk/hapcall/stringent/stringent_snps.fas
# manually substituted * for Ns using cntrl+f in Notepad++

# Maximum likelihood phylogenetic analysis
# RAXML manual recommends using the ASC versions of all models for SNP-only data, as it will try
# to correct for ascertainment bias in the model from not having any invariable sites.
#
# homogeneous rate of evolution
raxmlHPC-PTHREADS-AVX -T 16 -f a -p 33395 -s /mnt/USERS/kpham/oaks-hybseq/gatk/hapcall/stringent/stringent_snps.fas -x 33395 -# 100 -m ASC_GTRCAT -V -w /mnt/USERS/kpham/oaks-hybseq/gatk/phylo/raxml_20190819_1 -n stringent_snps.tre --asc-corr=lewis
# heterogeneous GAMMA model of evolution
raxmlHPC-PTHREADS-AVX -T 16 -f a -p 93742 -s /mnt/USERS/kpham/oaks-hybseq/gatk/hapcall/stringent/stringent_snps.fas -x 93742 -# 100 -m ASC_GTRGAMMA -w /mnt/USERS/kpham/oaks-hybseq/gatk/phylo/raxml_20190819_2 -n stringent_snps.tre --asc-corr=lewis
