#!/usr/bin/env python2.7

import sys,os,glob
os.environ["PATH"] += os.pathsep + /bin/isomut/ +'/src'
from isomut_wrappers_mod import run_isomut

# Isomut finds only unique mutations in a sample set
# It means that as we know that the Fr4A-Fr4B and Fc4A-Fc4B pairs are siblings, we ran two instances of the program

# defining running parameters

params=dict()
params['n_min_block']=200
params['n_conc_blocks']=12
params['ref_fasta']="/Projects/Planaria/PacBio/scaffolding/bm_Spol_g1.fasta"
params['input_dir']='/Projects/Planaria/Mutations/aligned/'

# defining mutation calling parameters

params['min_sample_freq']=0.15
params['min_other_ref_freq']=0.9
params['cov_limit']=20
params['base_quality_limit']=30
params['min_gap_dist_snv']=0
params['min_gap_dist_indel']=20

# Isomut run 1: without Fr4A and Fc4A

params['output_dir']='/Projects/Planaria/Mutations/IsoMut/1/'
params['bam_filenames']=list(map(os.path.basename, glob.glob(params['input_dir'] + "/*bam")))
params['bam_filenames'].remove("Fr4A_RMdup.bam")
params['bam_filenames'].remove("Fc4A_RMdup.bam")
run_isomut(params)

# Isomut run 1: without Fr4B and Fc4B

params['output_dir']='/Projects/Planaria/Mutations/IsoMut/2/'
params['bam_filenames']=list(map(os.path.basename, glob.glob(params['input_dir'] + "/*bam")))
params['bam_filenames'].remove("Fr4B_RMdup.bam")
params['bam_filenames'].remove("Fc4B_RMdup.bam")
run_isomut(params)

# this results in two-two SNV and indel out tables in /Projects/Planaria/Mutations/IsoMut/1/ and /Projects/Planaria/Mutations/IsoMut/2/
