#!/bin/bash

# alignment of each Illumina fastq dataset for the final bm_Spol_g1 assembly

BWA="/bin/bwa/bwa"
SAMTOOLS="/bin/samtools/samtools"
SAMBLASTER="/bin/samblaster/samblaster"

REF="/Projects/Planaria/scaffolding/lrna/bm_Spol_g1.fasta"
OUTDIR="/Projects/Planaria/Mutations/aligned/"

${BWA} index ${REF}

for SAMPLE in $(ls -d /Projects/Planaria/Illumina/fastq/P*fastq); do
	ls ${SAMPLE}/*_1.fq.gz | sort | xargs cat > ${WHICH}_forw.fq.gz
	ls ${SAMPLE}/*_2.fq.gz | sort | xargs cat > ${WHICH}_rev.fq.gz

	${BWA} mem -t 18 -R "@RG\tID:group1\tSM:"$(basename $SAMPLE _fastq)"\tPL:illumina\tLB:lib1\tPU:unit1" \
	${REF} ${WHICH}_forw.fq.gz ${WHICH}_rev.fq.gz |\
	${SAMBLASTER} -r | ${SAMTOOLS} view -@18 -u - |\
	${SAMTOOLS} sort -@ 18 -m 500M - -o ${OUTDIR}/$(basename $SAMPLE _fastq)_RMdup.bam
	${SAMTOOLS} index ${OUTDIR}/$(basename $SAMPLE _fastq)_RMdup.bam

	rm {WHICH}_forw.fq.gz ${WHICH}_rev.fq.gz
done


