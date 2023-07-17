#!/bin/bash

mkdir /Project/Planaria/Mutations/ploidy
cd /Project/Planaria/Mutations
BCFTOOLS="/bin/bcftools/bcftools"

# filter for high-confidence SNPs present in all samples

${BCFTOOLS} filter -i 'INFO/AN=24 && TYPE="snp" && QUAL > 1000 && COUNT(ALT)=1 &&  ABS(INFO/MQRankSum) < 3 && INFO/FS < 10' control_samples.raw.vcf.gz |\
cut -f1,2,4,5,10-17 | tr ":" "\t" | cut -f1,2,3,4,6,11,16,21,26,31,36,41 > ploidy/control_samples_for_ploidy.txt

${BCFTOOLS} filter -i 'INFO/AN=24 && TYPE="snp" && QUAL > 1000 && COUNT(ALT)=1 &&  ABS(INFO/MQRankSum) < 3 && INFO/FS < 10' regenerated_samples.raw.vcf.gz |\
cut -f1,2,4,5,10-21 |  tr ":" "\t" | cut -f1,2,3,4,6,11,16,21,26,31,36,41,46,51,56,61 > ploidy/regenerated_samples_for_ploidy.txt

# draw histograms of allele frequency distributions

Rscript ploidy_histogram.R ploidy/control_samples_for_ploidy.txt
Rscript ploidy_histogram.R ploidy/regenerated_samples_for_ploidy.txt
