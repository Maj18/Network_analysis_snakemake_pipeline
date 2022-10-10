#!/home/zjr863/.conda/envs/yuan2 python
import sys
import pandas as pd
from sklearn.preprocessing import binarize

DIR=sys.argv[1]
output_auc_mtx = DIR+"/B/output/auc_mtx.csv"
output_auc_mtx_binary = DIR+"/B/output/auc_mtx_binary.csv"

auc_mtx = pd.read_csv(output_auc_mtx) # relative position
auc_mtx = pd.DataFrame(auc_mtx)
auc_mtx.loc[0:3, :]

indices = auc_mtx.loc[:, "Cell"]
indices

regs = list(auc_mtx.columns.values)[1:auc_mtx.shape[1]]
auc_mtx = auc_mtx.loc[:, regs]
#auc_mtx = pd.DataFrame(auc_mtx, index=indices)
auc_mtx=auc_mtx.rename(index=indices)

auc_mtx_binary = binarize(auc_mtx)
auc_mtx_binary=pd.DataFrame(auc_mtx_binary, columns=list(auc_mtx.columns), index=list(auc_mtx.index))

auc_mtx_binary.to_csv(output_auc_mtx_binary)

