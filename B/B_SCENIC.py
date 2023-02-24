# Ref
### https://www.nature.com/articles/s41596-020-0336-2
### https://scenic.aertslab.org/tutorials/
### https://github.com/aertslab/SCENICprotocol/blob/master/notebooks/PBMC10k_SCENIC-protocol-CLI.ipynb

# import dependencies
#import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
import sklearn
import os
#from MulticoreTSNE import MulticoreTSNE as TSNE

import session_info
session_info.show()

#!~/.conda/envs/yuan2/bin/pyscenic



# set a working directory
wdir = "/home/zjr863/Proj2_VLMC/FGF_2022/2022-02/Double/network_FGF8/"
os.chdir( wdir )


# # path to loom file with basic filtering applied (this will be created in the "initial filtering" step below). Optional.
f_loom_path_scenic = wdir+"B_SCENIC/output/filtered_scenic.loom"
output_adjaceny = wdir+"B_SCENIC/output/adjacencies.csv"
output_regulon = wdir+"B_SCENIC/output/regulons.csv"
output_auc_mtx = wdir+"B_SCENIC/output/auc_mtx.csv"
THREADS=20



# 1. SCENIC analysis: predict activated regulons

### 1. Gene regulatory network inference and generation of co-expressiom modules from expression matrix

# Phase Ia GRN inference using the GRNBoost2 algorithm

#For this step the CLI version of SCENIC is used. This step can be deployed on an High performance computing system. We use the counts matrix (without log transformation or further processing) from the loom file we wrote earier. OUTPUT: list of adjancencies between a TF and its targets stored in ADACENCIES_FNAME

# transcription factors list
f_tfs = "/scratch/yuan/pyscenicdata/hs_hgnc_curated_tfs.txt" # human
# f_tfs = "/ddn1/vol1/staging/leuven/stg_00002/lcb/cflerin/resources/allTFs_dmel.txt" # drosophila
# f_tfs = "/ddn1/vol1/staging/leuven/stg_00002/lcb/cflerin/resources/allTFs_mm.txt"   # mouse
# tf_names = load_tf_names( f_tfs )

# f_loom_path_scenic: This dataset has been filtered, but not processed
!~/.conda/envs/yuan2/bin/pyscenic grn {f_loom_path_scenic} {f_tfs} -o {output_adjaceny} \
    --num_workers { THREADS }





### 2. Regulon prediction aka cisTarget from CLI (Find enriched motifs for a gene signature and opitionally prune targets from this signature based on cis-regulatory cues)
#For this step the CLI version of SCENIC is used. This step can be deployed on an High Performance Computing system.
#Output: List of adjacencies between a TF and its targets stored in MOTIFS_FNAME
import glob
# ranking databases
f_db_glob = "/scratch/yuan/pyscenicdata/*feather"
f_db_names = ' '.join( glob.glob(f_db_glob) )
print(f_db_names)

# motif databases (motifs of TFs)
f_motif_path = "/scratch/yuan/pyscenicdata/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"

#f_loom_path_scenic: this loom file is like counts.tsv
!~/.conda/envs/yuan2/bin/pyscenic ctx {output_adjaceny} \
  {f_db_names} \
  --annotations_fname {f_motif_path} \
  --expression_mtx_fname {f_loom_path_scenic} \
  --mode "dask_multiprocessing" \
  --output {output_regulon} \
  --mask_dropouts \
  --num_workers
#reg.csv: is the direct target output (regulon)
#f_loom_path_scenic: This dataset has been filtered, but not processed





### 3. aucell: quantify activity of gene signatures/regulons across single cells
#f_loom_path_scenic: this loom file is like counts.tsv
!~/.conda/envs/yuan2/bin/pyscenic aucell \
  {f_loom_path_scenic} \
  {output_regulon} \
  --output {output_auc_mtx} \
  --num_workers { THREADS }




