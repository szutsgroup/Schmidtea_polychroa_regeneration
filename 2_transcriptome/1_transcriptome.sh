#!/bin/bash

# transcriptome assembly with Trinity
# -----------------------------------

TRINITY="/bin/trinityrnaseq-v2.13.2/Trinity"
LEFT="/Projects/Planaria/Illumina/fastq/Spol_RNA_1.fq.gz"
RIGHT="/Projects/Planaria/Illumina/fastq/Spol_RNA_2.fq.gz"
OUTDIR="LEFT="/Projects/Planaria/Transcriptome/trinity_output"

mkdir -p $OUTDIR

${TRINITY} --seqType fq --max_memory 100G --left ${LEFT} --right ${RIGHT} \
--CPU 16 --SS_lib_type RF --output ${OUTDIR}

# running all tools required by Trinotate
# ---------------------------------------

QUERY="/Projects/Planaria/Transcriptome/trinity_output/trinity_output.Trinity.fasta"
TRINOTATE="/bin/Trinotate-Trinotate-v4.0.0/Trinotate"

mkdir -p /Projects/Planaria/Transcriptome/Trinotate/data

${TRINOTATE} --create -db /Projects/Planaria/Transcriptome/Trinotate/Trinotate.sqlite \
--trinotate_data_dir /Projects/Planaria/Transcriptome/Trinotate/data/

# transdecoder

mkdir -p  /Projects/Planaria/Transcriptome/transdecoder_output

/bin/TransDecoder-TransDecoder-v5.5.0/TransDecoder.LongOrfs \
-t /Projects/Planaria/Transcriptome/trinity_output/trinity_output.Trinity.fasta \
--gene_trans_map /Projects/Planaria/Transcriptome/trinity_output/trinity_output.Trinity.fasta.gene_trans_map \
-O /Projects/Planaria/Transcriptome/transdecoder_output

PEP_QUERY="/Projects/Planaria/Transcriptome/transdecoder_output/longest_orfs.pep"

# BlastX and BlastP

mkdir -p /Projects/Planaria/Transcriptome/BLAST_output

/bin/ncbi-blast-2.10.0+/bin/blastx -query $QUERY -db /Projects/Planaria/Transcriptome/Trinotate/data/uniprot_sprot.pep \
-num_threads 12 -max_target_seqs 1 -max_hsps 1 -outfmt 6 \
-evalue 1e-3 > /Projects/Planaria/Transcriptome/BLAST_output/blastx.outfmt6

/bin/ncbi-blast-2.10.0+/bin/blastp -query $PEP_QUERY -db /Projects/Planaria/Transcriptome/Trinotate/data/uniprot_sprot.pep \
-num_threads 12 -max_target_seqs 1 -max_hsps 1 -outfmt 6 \
-evalue 1e-3 > /Projects/Planaria/Transcriptome/BLAST_output/blastp.outfmt6

# HMMScan

mkdir -p /Projects/Planaria/Transcriptome/HMMScan_output

/bin/hmmscan --cpu 12 --domtblout /Projects/Planaria/Transcriptome/HMMScan_output/pfam.out \
/Projects/Planaria/Transcriptome/Trinotate/data/Pfam-A.hmm ${PEP_QUERY}

# INFERNAL

mkdir -p /Projects/Planaria/Transcriptome/INFERNAL_output

/bin/infernal-1.1.4/src/cmsearch /Projects/Planaria/Transcriptome/Trinotate/data/Rfam.cm \
${QUERY} > /Projects/Planaria/Transcriptome/INFERNAL_output/infernal.out

# tmhmm

mkdir -p /Projects/Planaria/Transcriptome/tmhmm_output

/bin/tmhmm-2.0c/bin/tmhmm --short < ${PEP_QUERY} > /Projects/Planaria/Transcriptome/tmhmm_output/tmhmm.out

# SignalP

mkdir -p /Projects/Planaria/Transcriptome/signalP_output

/bin/signalp6 -fasta ${PEP_QUERY} -od /Projects/Planaria/Transcriptome/signalP_output

# running Trinotate
# -----------------

${TRINOTATE} --db /Projects/Planaria/Transcriptome/Trinotate/Trinotate.sqlite --init \
--gene_trans_map /Projects/Planaria/Transcriptome/trinity_output/trinity_output.Trinity.fasta.gene_trans_map \
--transcript_fasta ${QUERY}
--transdecoder_pep ${PEP_QUERY}

${TRINOTATE} --db /Projects/Planaria/Transcriptome/Trinotate/Trinotate.sqlite \
--LOAD_signalp /Projects/Planaria/Transcriptome/signalP_output/output.gff3

${TRINOTATE} --db /Projects/Planaria/Transcriptome/Trinotate/Trinotate.sqlite \
--LOAD_tmhmmv2 /Projects/Planaria/Transcriptome/tmhmm_output/tmhmm.out

${TRINOTATE} --db /Projects/Planaria/Transcriptome/Trinotate/Trinotate.sqlite \
--LOAD_swissprot_blastp /Projects/Planaria/Transcriptome/BLAST_output/blastp.outfmt6

${TRINOTATE} --db /Projects/Planaria/Transcriptome/Trinotate/Trinotate.sqlite \
--LOAD_pfam /Projects/Planaria/Transcriptome/HMMScan_output/pfam.out

${TRINOTATE} --db /Projects/Planaria/Transcriptome/Trinotate/Trinotate.sqlite \
--LOAD_infernal /Projects/Planaria/Transcriptome/INFERNAL_output/infernal.out

${TRINOTATE} --db /Projects/Planaria/Transcriptome/Trinotate/Trinotate.sqlite \
--LOAD_swissprot_blastx /Projects/Planaria/Transcriptome/BLAST_output/blastx.outfmt6


{TRINOTATE} --db /Projects/Planaria/Transcriptome/Trinotate/Trinotate.sqlite --report

# generate PCF and PCFL sets using the attached R script
# ------------------------------------------------------

Rscript pcfl_filtering.R --report /Projects/Planaria/Transcriptome/Trinotate/Trinotate.report.tsv --fasta ${QUERY}
