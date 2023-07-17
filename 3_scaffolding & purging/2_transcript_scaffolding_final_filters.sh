#!/bin/bash

# L_RNA_Scaffolder: scaffolding with transcriptome
# ------------------------------------------------

ASSEMBLY="/Projects/Planaria/PacBio/scaffolding/purge3/Spol_ar3_pi_pdill045_redundans_OL50_custom.fa"

PCFL="/Projects/Planaria/Transcriptome/trinity_output/trinity_output.Trinity.pcfl.fasta"

OUTDIR="/Projects/Planaria/PacBio/scaffolding/lrna/"

/bin/blat/blat -noHead $ASSEMBLY $PCFL ${OUTDIR}/$(basename $ASSEMBLY .fasta).blatout

sh /szutsnas/bin/L_RNA_scaffolder/L_RNA_scaffolder.sh \
-d /szutsnas/bin/L_RNA_scaffolder/ \
-i $(basename $ASSEMBLY .fasta).blatout \
-j $ASSEMBLY -o $OUTDIR

# this will generate /Projects/Planaria/PacBio/scaffolding/lrna/L_RNA_Scaffolder.fasta


# final filtering by removing scaffolds shorter than 1kb
# and those that have a GC ratio over 40% (see Fig S1A)
# ------------------------------------------------------

cd /Projects/Planaria/PacBio/scaffolding/

ASSEMBLY="/Projects/Planaria/PacBio/scaffolding/lrna/L_RNA_Scaffolder.fasta"
mkdir split

/szutsnas/bin/seqkit split -i -O split/ $ASSEMBLY

# shorten the names of the split scaffolds


# get short scaffolds
/bin/samtools/samtools faidx $ASSEMBLY

awk '$2 <= 1000 {print $1}' > short_scaffolds
for i in split/*fasta
do
    n=$(echo $i | sed 's/Spol_long3ar3pi1_pdill045_redundans_OL50_custom_lrna.id_//')
    mv $i $n
    echo "renaming " $i "-->" $n
done

# get GC ratios
/bin/EMBOSS-6.6.0/emboss/geecee -sequence $ASSEMBLY \
-outfile $(basename $ASSEMBLY .fasta).geecee

awk '$2 >= 0.4 {print $1}' > high_gc

while read x
do
    n=$(echo $x | sed 's/:/__/g')
    mv split/${n}.fasta split/${n}.not
done < <(cat short_scaffolds high gc)

cat split/${n}.fasta > bm_Spol_v1.fasta
