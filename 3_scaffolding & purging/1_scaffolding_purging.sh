#!/bin/bash

MINIMAP="/bin/minimap2/minimap2"
SAMTOOLS="/bin/samtools/samtools"

ASSEMBLY="/Projects/Planaria/PacBio/polishing/Spol_ar3_pi.fa"
OUTDIR="/Projects/Planaria/PacBio/scaffolding/"

mkdir -p $OUTDIR

# run purge_dups for the first round of purging
# ---------------------------------------------

mkdir /Projects/Planaria/PacBio/scaffolding/purge_dups
cd /Projects/Planaria/PacBio/scaffolding/purge_dups
export PATH="/bin/purge_dups/bin/":$PATH 

ngscstat ${BAMDIR}/Spol_illumina.bam

calcuts -l 7 -m 45 -u 70 TX.stat > cutoffs

split_fa ${ASSEMBLY} > $(basename ${ASSEMBLY} .fa).split.fasta

SPLIT=$(basename ${ASSEMBLY} .fa).split.fasta

${MINIMAP} -xasm5 -t 24 -DP ${SPLIT} ${SPLIT} \
| gzip -c - > Spol_ar3_pi.split.self.paf.gz

purge_dups -2 -T cutoffs -c TX.base.cov Spol_ar3_pi.split.self.paf.gz \
 > dups.bed

get_seqs dups.bed ${ASSEMBLY}

# this results in  /Projects/Planaria/PacBio/scaffolding/purge_dups/purged.fa

# run redundans.py for the second round of purging and scaffolding
# ----------------------------------------------------------------

# for the scaffolding, we use the same PacBio long read set as 
# for the assembly, and also 2 Illumina datasets, which are 
# different than the one used in Pilon polishing

mkdir -p /Projects/Planaria/PacBio/scaffolding/redundans/
cd /Projects/Planaria/PacBio/scaffolding/redundans/

/bin/redundans/redundans.py \
-i /Projects/Planaria/Illumina/fastq/SM00[23]*P.fq.gz \
-f ../purge_dups/purged.fa \
-l /Projects/Planaria/PacBio/fastq/Spol_pacbio.fastq.gz \
--overlap 0.50 -v -t 24 -o Spol_ar3_pi_pdill045_redundans_OL50

# this results in 
# /Projects/Planaria/PacBio/scaffolding/redundans/Spol_ar3_pi_pdill045_redundans_OL50/scaffolds.fa


# run the attached custom_purging.R script
# ----------------------------------------

mkdir -p /Projects/Planaria/PacBio/scaffolding/purge3/
cd /Projects/Planaria/PacBio/scaffolding/purge3/

# create a self-alignment with minimap2 for the assembly
$ASSEMBLY="/Projects/Planaria/PacBio/scaffolding/redundans/Spol_ar3_pi_pdill045_redundans_OL50/scaffolds.fa"

/bin/minimap2/minimap2 -x asm5 -t 12 -D -n 15 --secondary=no \ 
$ASSEMBLY $ASSEMBLY > assembly_self.paf

# split the assembly by sequences
mkdir split split_mod

/szutsnas/bin/seqkit split -i -o split/ ${ASSEMBLY}

# run the R script
Rscript custom_purging.R --dir split --outdir split_mod --paf assembly_self.paf \
--re ...

# replace the modified scaffolds with their purged variants
ls split_mod/*fa | sed -e s,split_mod/,, -e s/_deleted[0-9]*.fa// > changed_contigs
for i in $(ls split/*fasta | grep -v -f changed_contigs); do cp -v $i split_mod/; done
cat split_mod/* > Spol_ar3_pi_pdill045_redundans_OL50_custom.fa

rm -r split split_mod





