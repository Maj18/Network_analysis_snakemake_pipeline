# Based on target genes to look for TF

library(Seurat)
library(dplyr)

DIR <- "/home/zjr863/Proj3_hypothalamus_ref/gw10_new/network/"
OUT_DIR <- paste0(DIR, "C/output/target2TF/")
NAME <- "SCENIC_target2TF"
CLUSTERING <- "SCT_snn_res.1.3_modi"
targets <- "LHX9"
clr <- "15"




# 1. Define functions
target2TF <- function(regulons, target, DERs, OUT_DIR, clr) {
	# Regulons that has "target" as target genes
	index <- grepl("\\'", target, "\\'", regulons$Enrichment.6)
	regulon_all <- regulons$X
	regs <- regulon_all[index] %>% unique()

	sink(paste0(OUT_DIR, "From_", target, "2TF_in_cluster_", clr, ".txt"))
	#Diff enriched regulons in cluster clr
	print(paste0("Differentially enriched regulons in cluster ", clr, ":"))
	DER <- DERs$gene[DERs$cluster %in% clr]
	if (length(DER)>0) {
		print(DER)
	} else {
		print("None")
	}

	# See whether Regulons that has "target" as target genes are among the Diff enriched regulons in cluster clr
	print(paste0("Which regulons that have ", target, " as target genes are among the differentially enriched regulons in cluster ", clr, ":"))
	shared <- intersect(regs, DER)
	print(shared)

	for (r in regs) {
  		sub <- subset(regulons, X%in%r)
  		motifs = stringr::str_split(sub$Enrichment.6, "[']")
  		targets <- lapply(motifs, function(motif) motif[(1:abs(length(motif)/2)) *2]) %>% unlist() %>% unique()
  		print(r)
  		print(paste0(target, " is among the targets of ", r, ": "))
  		print(target %in% targets)
  		print(paste0("The targets of ", r, " are :"))
  		print(targets)
	}
	sink()
}





# 2. Load data
name <- readRDS(paste0(DIR, "gw10_processed.RDS"))
DERs <- read.csv(paste0(DIR, "C/output/DER/SCENIC_DER_marker_table_", CLUSTERING,".csv"),
                 header=TRUE, sep=";")
aucell <- read.csv(paste0(DIR, "B/output/auc_mtx.csv"), sep=",")
head(aucell)
rownames(aucell) <- aucell$Cell
aucell=aucell[, 2:ncol(aucell)]
colnames(aucell) <- gsub("[...]", "", colnames(aucell) )
head(aucell)

aucell_bin <- read.csv(paste0(DIR, "B/output/auc_mtx_binary.csv"), sep=",")
head(aucell_bin)
rownames(aucell_bin) <- aucell_bin$X
aucell_bin=aucell_bin[, 2:ncol(aucell_bin)]
colnames(aucell_bin) <- gsub("[...]", "", colnames(aucell_bin) )

regulons <- read.csv(paste0(DIR, "B/output/regulons.csv"), sep=",")


AUC <- CreateSeuratObject(t(aucell))
AUC[['AUCBinary']] <- CreateAssayObject(data = as.matrix(t(aucell_bin)))
AUC = AddMetaData(AUC, metadata=name@meta.data[colnames(AUC), ])
AUC@reductions=name@reductions
AUC@graphs=name@graphs
AUC@neighbors=name@neighbors

adj <- read.csv(paste0(DIR, "B/output/adjacencies.csv"), sep=",")


# 3. Sear regulons for target
for (target in targets) {
	target2TF(regulons, target, DERs, OUT_DIR, clr)
}





sink(paste0(OUT_DIR, "sessionInfo.txt"))
sessionInfo()
sink()


