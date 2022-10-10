#!/home/zjr863/.conda/envs/yuan2 python
# From count matrix to loom file

import sys
#import os
import numpy as np
#import pandas as pd
import scanpy as sc
import loompy as lp
import session_info
session_info.show()



#The count data here has been filtered in Seurat
wdir = sys.argv[1]
f_exprMat = wdir+'/A/output/counts.tsv'
adata = sc.read_text( f_exprMat, delimiter='\t', first_column_names=True )

#Write to an loom file
row_attrs = { 
    "Gene": np.array(adata.var.index) ,
}

col_attrs = { 
    "CellID":  np.array(adata.obs.index) ,
    "nGene": np.array( np.sum(adata.X.transpose()>0 , axis=0)).flatten() ,
    "nUMI": np.array( np.sum(adata.X.transpose() , axis=0)).flatten() ,
}

f_loom_path_scenic = wdir+"/B/output/filtered_scenic.loom"
lp.create( f_loom_path_scenic, adata.X.transpose(), row_attrs, col_attrs )

