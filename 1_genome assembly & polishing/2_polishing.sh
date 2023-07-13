#!/bin/bash

# --------------------------------------------------
# arrow was tun with three iterations
#
# in each case, the PacBio reads (same as those used 
# for the assembly) had to be aligned against the 
# most recent assembly version
# --------------------------------------------------

BAMDIR="/Projects/Planaria/PacBio/aligned/"
POLDIR="/Projects/Planaria/PacBio/polishing/"

PB_BAM="/Projects/Planaria/PacBio/H201SC19061971/raw_data/m54220_190904_082545/m54220_190904_082545.subreads.bam"
FASTQ="/Projects/Planaria/PacBio/fastq/Spol_pacbio.fastq.gz"

MINIMAP="/bin/minimap2/minimap2"
SAMTOOLS="/bin/samtools/samtools"
PILON="/bin/pilon-1.23.jar"
BWA="/bin/bwa/bwa"
SAMBLASTER="/bin/samblaster/samblaster"
TRIMMOMATIC="/bin/Trimmomatic-0.39/trimmomatic-0.39.jar"

mkdir -p ${BAMDIR} ${POLDIR}

# ARROW ITERATION 1
# -----------------

REF="/Projects/Planaria/PacBio/canu_assembly/Spol.contigs.fasta"

${MINIMAP} -t 16 -ax map-pb ${REF} ${FASTQ} | ${SAMTOOLS} view -buS - -o ${OUTDIR}/Spol_long.bam

/bin/pacbio/bin/pbbamify--input ${BAMDIR}/Spol_long3.bam \
--output ${BAMDIR}/Spol_long3_pbbamified.bam ${REF} ${PB_BAM}

${SAMTOOLS} sort -m 2G -@ 16 ${BAMDIR}/Spol_long3_pbbamified.bam -o ${BAMDIR}/Spol_long3_pbbamified.sorted.bam

/bin/pacbio/bin/pbindex ${BAMDIR}/Spol_long3_pbbamified.sorted.bam

/bin/pacbio/bin/arrow ${BAMDIR}/Spol_long3_pbbamified.sorted.bam \
-o ${POLDIR}/Spol_ar1.fa \
-o ${POLDIR}/Spol_ar1.fq \
-o ${POLDIR}/Spol_ar1.gff

rm ${BAMDIR}/*

# ARROW ITERATION 2
# -----------------

REF="${POLDIR}/Spol_ar1.fa"

${MINIMAP} -t 16 -ax map-pb ${REF} ${FASTQ} | ${SAMTOOLS} view -buS - -o ${OUTDIR}/Spol_long.bam

/bin/pacbio/bin/pbbamify--input ${BAMDIR}/Spol_long.bam \
--output ${BAMDIR}/Spol_long_pbbamified.bam ${REF} ${PB_BAM}

${SAMTOOLS} sort -m 2G -@ 16 ${BAMDIR}/Spol_long_pbbamified.bam -o ${BAMDIR}/Spol_long_pbbamified.sorted.bam

/bin/pacbio/bin/pbindex ${BAMDIR}/Spol_long_pbbamified.sorted.bam

/bin/pacbio/bin/arrow ${BAMDIR}/Spol_long_pbbamified.sorted.bam \
-o ${POLDIR}/Spol_ar2.fa \
-o ${POLDIR}/Spol_ar2.fq \
-o ${POLDIR}/Spol_ar2.gff

rm ${BAMDIR}/*

# ARROW ITERATION 3
# -----------------

REF="${POLDIR}/Spol_ar2.fa"

${MINIMAP} -t 16 -ax map-pb ${REF} ${FASTQ} | ${SAMTOOLS} view -buS - -o ${OUTDIR}/Spol_long.bam

/bin/pacbio/bin/pbbamify--input ${BAMDIR}/Spol_long3.bam \
--output ${BAMDIR}/Spol_long3_pbbamified.bam ${REF} ${PB_BAM}

${SAMTOOLS} sort -m 2G -@ 16 ${BAMDIR}/Spol_long3_pbbamified.bam -o ${BAMDIR}/Spol_long3_pbbamified.sorted.bam

/bin/pacbio/bin/pbindex ${BAMDIR}/Spol_long3_pbbamified.sorted.bam

/bin/pacbio/bin/arrow ${BAMDIR}/Spol_long3_pbbamified.sorted.bam \
-o ${POLDIR}/Spol_ar3.fa \
-o ${POLDIR}/Spol_ar3.fq \
-o ${POLDIR}/Spol_ar3.gff

rm ${BAMDIR}/*

# PILON
# -----

REF="${POLDIR}/Spol_ar3.fa"

# trim adapters from Illumina reads, and do k-mer based error correction
F=/Projects/Planaria/Illumina/fastq/SM001_1.fq.gz
R=/Projects/Planaria/Illumina/fastq/SM001_2.fq.gz

java jar ${TRIMMOMATIC} PE -threads 16 -phred33 $F $R \
-baseout /Projects/Planaria/Illumina/fastq/SM001_trimmed.fq.gz \
ILLUMINACLIP:/bin/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10:2:TRUE \
MINLEN:100 &> /Projects/Planaria/Illumina/fastq/SM001.log

ls /Projects/Planaria/Illumina/fastq/SM001_trimmed_*fq.gz | tr "\n" " " > quake_filelist
/szutsnas/bin/Quake/bin/quake.py --no_jelly -f ${FILELIST} -k 19 -p 16

FREADS="/Projects/Planaria/Illumina/fastq/SM001_trimmed_1P.cor.fq.gz"
RREADS="/Projects/Planaria/Illumina/fastq/SM001_trimmed_2P.cor.fq.gz"

# align Illumina reads with bwa
${BWA} index ${REF}

${BWA} mem -t 16 ${REF} ${FREADS} ${RREADS} | ${SAMBLASTER} -r | ${SAMTOOLS} view -buS - \
| ${SAMTOOLS} sort -m 5G -@ 16 - -o ${BAMDIR}/Spol_illumina.bam

${SAMTOOLS} index ${BAMDIR}/Spol_illumina.bam

# run Pilon
java -Xmx240G -jar ${PILON} --genome ${REF} --bam ${BAMDIR}/Spol_illumina.bam --changes \
--fix all,breaks,circles --threads 16 --output ${POLDIR}/Spol_ar3_pi.fa

# -------------------------------------------------------------------
# the final output is the polished fasta file polished/Spol_ar3_pi.fa
# -------------------------------------------------------------------

