#!/home/zjr863/.conda/envs/yuan2 R
args = commandArgs(trailingOnly=TRUE)
DIR= args[1]
print(DIR)
ind <- args[2] #ARC_neuron

library(Seurat)
name=readRDS(paste0(DIR, "/A/data/", ind, "_processed.RDS"))
write.table(t(as.matrix(name@assays$SCT@counts)), 
		paste0(paste0(DIR, '/A/output/counts.tsv')), 
		sep = '\t', col.names = NA)

sink(paste0(DIR, "/A/notebooks/sessionInfo.txt"))
sessionInfo()
sink()


