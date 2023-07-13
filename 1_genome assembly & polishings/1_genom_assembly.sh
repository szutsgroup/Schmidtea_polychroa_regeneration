#!/bin/bash


SAMTOOLS="/szutsnas/bin/samtools/samtools"
CANU="/bin/canu/Linux-amd64/bin/canu"

PB_BAM="/Projects/Planaria/PacBio/H201SC19061971/raw_data/m54220_190904_082545/m54220_190904_082545.subreads.bam"
PB_FASTQ="/Projects/Planaria/PacBio/fastq/Spol_pacbio.fastq.gz"

# generate fastq from PacBio bam
${SAMTOOLS} fastq -0 ${PB_FASTQ} -@ 16 ${PB_BAM}

# Canu is run with parameters biased for low GC genomes and somewhat low coverage, according to the Canu documentation

# Canu. correction
${CANU} -correct -p Spol -d canu_assembly genomeSize=700000000 \
useGrid=false maxThreads=32 correctedErrorRate=0.105 corMinCoverage=0 \
corMhapSensitivity="high" corOutCoverage=10000 corMaxEvidenceErate=0.15 \
 -pacbio-raw ${PB_FASTQ}

COR_FASTQ="/Projects/Planaria/PacBio/canu_assembly/Spol.correctedReads.fasta.gz""

# Canu trimming
${CANU} -trim -p Spol -d canu_assembly genomeSize=700000000 \
stopOnLowCoverage=0.1 useGrid=false maxThreads=48 \
correctedErrorRate=0.105 -pacbio ${COR_FASTQ}

TRIM_FASTQ="/Projects/Planaria/PacBio/canu_assembly/Spol.trimmedReads.fasta.gz"

# Canu assemby
${CANU} -assemble -p Spol -d canu_assembly \
genomeSize=700000000 useGrid=false correctedErrorRate=0.105 \
corOutCoverage=10000 corMhapSensitivity=high \
maxThreads=48 -pacbio ${TRIM_FASTQ}

# these commands generate the raw Canu output file canu_assembly/Spol.contigs.fasta


