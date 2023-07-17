#!/bin/bash

# re-check allele frequencies for all three raw vcf files (regenerated, control and parent) with the getaf.py script

BAMDIR="/Projects/Planaria/Mutations/aligned/"

cd /Projects/Planaria/Mutations/hc/

for i in regenerated control parent; do
    gzip -cd ${i}_samples.raw.vcf.gz | grep -v "#" | cut -f 1,2,4,5 > ${i}_samples_position.txt

    mkdir ${i}_symlinks
    case $i in
    
        regenerated)
            for j in $(${BAMDIR}/* | grep -e "Pr" -e "R" -e "Fr"); do
                ln -s $j ${i}_symlinks/
            done
            ;;
            
        control)
            for j in $(${BAMDIR}/* | grep -e "Pc" -e "Fc"); do
                ln -s $j ${i}_symlinks/
            done
            ;;
            
        parent)
            for j in $(${BAMDIR}/* | grep -e "Pr" -e "Pc"); do
                ln -s $j ${i}_symlinks/
            done
            ;;
    esac

    python3 getaf.py ${i}_samples_position.txt ${i}_symlinks/ ${i}_samples_getaf.out
done

# filter the getaf output for rows that:
# a) mean coverage across the considered samples is between 25 and 130
# b) the minimal mapping quality in the sample where the mutation is detected is over 40
# c) the mutation is "present" if it has at least 7 supporting reads
# d) the mutation is only present in 1-2 samples (to include shared mutations in siblings),
# also one or two extra sample with some support is allowed (because of possible subclonal presence in regenerants/parents)

Rscript getaf_filter.R --input control_samples_getaf.out \
--zeromin 5 --zeromax 7 --covmin 25 --covmax 130 --maxsuppmin 7 --minnonzerosupp 1

Rscript getaf_filter.R --input regenerated_samples_getaf.out \
--zeromin 9 --zeromax 11 --covmin 25 --covmax 130 --maxsuppmin 7 --minnonzerosupp 1

Rscript getaf_filter.R --input parent_samples_getaf.out \
--zeromin 5 --zeromax 7 --covmin 25 --covmax 130 --maxsuppmin 7 --minnonzerosupp 1

# these commands result in three sets of filtered positions, which were manually verified in IGV,
# and were used to find lineage-specific and de novo mutation sets


    
  
