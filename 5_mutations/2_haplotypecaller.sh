#!/bin/bash

REF="/Projects/Planaria/scaffolding/lrna/bm_Spol_g1.fasta"
GATK="/bin/GenomeAnalysisTK.jar"
PICARD="/bin/picard.jar"
OUTDIR="/Projects/Planaria/Mutations/hc/"

# running HaplotypeCaller for each sample in GVCF mode
for SAMPLE in $(/Projects/Planaria/Mutations/aligned/*bam); do
	java -jar ${GATK} -T HaplotypeCaller -R ${REF} -I ${TARGET} \
	-o ${OUTDIR}/$(basename $SAMPLE .bam).vcf.gz -ERC GVCF -ploidy 3
done

# joint genotyping for regenerated, control and parent samples samples

ls -1 ${OUTDIR}/*g.vcf.gz | grep -e "Pc" -e "Fc" | awk '{print "-V "$1}' | tr "\n" " " > control.list
java -jar ${GATK} -T CombineGVCFs -R ${REF} --arg_file control.list -o control_samples.g.vcf.gz
java -jar ${GATK} -T GenotypeGVCFs -R ${REF} -V control_samples.g.vcf.gz -o control_samples.raw.vcf.gz -ploidy 3

ls -1 ${OUTDIR}/*g.vcf.gz | grep -e "Pr" -e "Fr" -r "R" | awk '{print "-V "$1}' | tr "\n" " " > regenerated.list
java -jar ${GATK} -T CombineGVCFs -R ${REF} --arg_file regenerated.list -o regenerated_samples.g.vcf.gz
java -jar ${GATK} -T GenotypeGVCFs -R ${REF} -V regenerated_samples.g.vcf.gz -o regenerated_samples.raw.vcf.gz -ploidy 3

ls -1 ${OUTDIR}/*g.vcf.gz | grep -e "Pr" -e "Pc" | awk '{print "-V "$1}' | tr "\n" " " > parent.list
java -jar ${GATK} -T CombineGVCFs -R ${REF} --arg_file parent.list -o parent_samples.g.vcf.gz
java -jar ${GATK} -T GenotypeGVCFs -R ${REF} -V parent_samples.g.vcf.gz -o parent_samples.raw.vcf.gz -ploidy 3
