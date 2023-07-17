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
params['ref_fasta']="/zfs2/szuts/poti/202303_planaria_review/assembly_versions/bm_Spol_v1.fasta"
#input dir output dir
params['input_dir']='/zfs2/szuts/poti/202303_planaria_review/aligned/'
