#!/bin/bash

GMAP_DIR="/bin/gmap-2019-06-10/bin/"
PCFL="/Projects/Planaria/Transcriptome/trinity_output/trinity_output.Trinity.pcfl.fasta"
ASSEMBLY="/Projects/Planaria/PacBio/scaffolding/bm_Spol_g1.fasta"

${GMAP_DIR}/gmap_build -d bm_Spol_g1 ${ASSEMBLY}
${GMAP_DIR}/gmap -d bm_Spol_g1 -S -t 8 ${PCFL} > /Projects/Planaria/Synteny/gmap_output/ddSpolv4_on_bmSpolg1.gmap

PCFL="/Projects/Planaria/Other "
ASSEMBLY="/Projects/Planaria/PacBio/scaffolding/bm_Spol_g1.fasta"

${GMAP_DIR}/gmap_build -d bm_Spol_g1 ${ASSEMBLY}
${GMAP_DIR}/gmap -d bm_Spol_g1 -S -t 8 ${PCFL} > /Projects/Planaria/Synteny/gmap_output/ddSpolv4_on_bmSpolg1.gmap
