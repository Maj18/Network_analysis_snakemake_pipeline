#!/home/zjr863/.conda/envs/yuan2 R
# Differentially enriched regulons
args=commandArgs(trailingOnly=TRUE)
DIR=args[1]
ind=args[2] #gw15_neuron

library(Seurat)
library(tidyverse)
library(magrittr)
library(viridis)
#library(SCENIC)
library(network)
library(igraph)


IN_DIR <- paste0(DIR, "/A/data/")
OUT_DIR <- paste0(DIR, "/C/output/DER/")
NAME <- "SCENIC_DER"
CLUSTERING <- "SCT_snn_res.0.5"
res <- "0.5"




# 2. Load data
name <- readRDS(paste0(IN_DIR, ind, "_processed.RDS"))
degs <- read.csv(paste0(DIR, "/", ind, "_marker_table_res.", res, ".csv"),
                 header=TRUE, sep=";")
aucell <- read.csv(paste0(DIR, "/B/output/auc_mtx.csv"), sep=",")
#head(aucell)
rownames(aucell) <- aucell[,"Cell"]
aucell=aucell[, 2:ncol(aucell)]
colnames(aucell) <- gsub("[...]", "", colnames(aucell) )
#head(aucell)

aucell_bin <- read.csv(paste0(DIR, "/B/output/auc_mtx_binary.csv"), sep=",")
#head(aucell_bin)
rownames(aucell_bin) <- aucell_bin[, "X"]
aucell_bin=aucell_bin[, 2:ncol(aucell_bin)]
colnames(aucell_bin) <- gsub("[...]", "", colnames(aucell_bin) )

regulons <- read.csv(paste0(DIR, "/B/output/regulons.csv"), sep=",")

AUC <- CreateSeuratObject(t(aucell))
AUC[['AUCBinary']] <- CreateAssayObject(data = as.matrix(t(aucell_bin)))
AUC = AddMetaData(AUC, metadata=name@meta.data[colnames(AUC), ])
AUC@reductions=name@reductions
AUC@graphs=name@graphs
AUC@neighbors=name@neighbors

adj <- read.csv(paste0(DIR, "/B/output/adjacencies.csv"), sep=",")


# name.markers.ARC_neuron <- fread("/home/zjr863/Proj3_hypothalamus_ref/gw15_neuron/network_analysis/C/output/DER/SCENIC_DER_marker_table_SCT_snn_res.0.5.csv")
# rownames(name.markers.ARC_neuron) <- name.markers.gw15_neuron$V1
# name.markers.ARC_neuron <- name.markers.ARC_neuron[, 2:ncol(name.markers.gw15_neuron)]
# colnames(name.markers.ARC_neuron)
# head(name.markers.ARC_neuron)

pdf(paste0(OUT_DIR, NAME, "_ARC_clr7_markers.pdf"), w=16, h=9)
FeaturePlot(AUC, 
            slot="counts",
            reduction = "umap",
            features = c("CREB5", "ISL1", "ONECUT1", "ONECUT2", "DLX5", "LHX8", "ESR2", "MAFG", "RARB"),
            pt.size = 0.6, keep.scale = "feature",
            label=F, ncol=4) #& NoLegend()
dev.off()


pdf(paste0(OUT_DIR, NAME, "_gw15_neuron_clr10_12_markers.pdf"), w=16, h=6)
FeaturePlot(AUC, 
            slot="counts",
            reduction = "umap",
            features = c("ZNF668", "EGR4", "CREBZF", "STAT5B", "TBX3", "ILS2", "MEF2C"),
            pt.size = 0.6, keep.scale = "feature",
            label=F, ncol=4) #& NoLegend()
dev.off()

# 3. Differential regulon identification based on regulon activity score
Idents(AUC) <- AUC@meta.data[, CLUSTERING]
name.markers <- FindAllMarkers(AUC, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.005)
name.markers <- subset(name.markers, subset=p_val_adj<0.05)
write.csv2(name.markers, paste0(OUT_DIR, NAME, "_marker_table_", CLUSTERING, ".csv"))

AUC <- ScaleData(AUC, features=rownames(AUC))

top_N = 10
tops <- name.markers %>% group_by(cluster) %>% top_n(n=top_N, wt=avg_log2FC)
pdf(paste0(OUT_DIR, NAME, "_Regulon_AUC_heatmap_", CLUSTERING, ".pdf"), h=13, w=19)
DoHeatmap(AUC, features = tops$gene, raster=F) +scale_fill_gradient2(
  low=rev(c("#d1e5f0", "#67a9cf", "#2166ac")),
  mid="white",
  high=rev(c("#b2182b", "#ef8a62", "#fddbc7")),
  midpoint=0,
  guide="colourbar",
  aesthetics="fill"
  )
dev.off()


top_N = 5
tops <- name.markers %>% group_by(cluster) %>% top_n(n=top_N, wt=avg_log2FC)
pdf(paste0(OUT_DIR, NAME, "_Regulon_AUC_heatmap_", CLUSTERING, "_small.pdf"), h=8, w=19)
DoHeatmap(AUC, features = tops$gene, raster=F) +scale_fill_gradient2(#+ scale_fill_viridis()
  low=rev(c("#d1e5f0", "#67a9cf", "#2166ac")),
  mid="white",
  high=rev(c("#b2182b", "#ef8a62", "#fddbc7")),
  midpoint=0,
  guide="colourbar",
  aesthetics="fill"
  )
dev.off()


top_N = 10
pt.size=0.6
tops <- name.markers %>% group_by(cluster) %>% top_n(n=top_N, wt=avg_log2FC)
pdf(paste0(OUT_DIR, NAME, "_cluster_top_regulonAUC_", CLUSTERING, ".pdf"), w=24, h=36)
FeaturePlot(AUC, 
      slot="counts",
      reduction = "umap",
      features = unique(tops$gene),
      pt.size = pt.size, keep.scale = "feature",
      label=F, ncol=7) #& NoLegend()
dev.off()

top_N = 30
pt.size=0.6
tops <- name.markers %>% group_by(cluster) %>% top_n(n=top_N, wt=avg_log2FC)
pdf(paste0(OUT_DIR, NAME, "_cluster_top_regulonAUC_", CLUSTERING, "_2.pdf"), w=24, h=36)
FeaturePlot(AUC, 
            slot="counts",
            reduction = "umap",
            features = unique(tops$gene)[113:length(unique(tops$gene))],
            pt.size = pt.size, keep.scale = "feature",
            label=F, ncol=8) #& NoLegend()
dev.off()



top_N = 2
pt.size=0.6
tops <- name.markers %>% group_by(cluster) %>% top_n(n=top_N, wt=avg_log2FC)
pdf(paste0(OUT_DIR, NAME, "_cluster_top_geneExpr_", CLUSTERING, ".pdf"), w=24, h=9)
FeaturePlot(name, 
      reduction = "umap",
      features = unique(tops$gene),
      pt.size = pt.size, keep.scale = "feature",
      label=F, ncol=6) #& NoLegend()
dev.off()




# 4. Binary regulon activity score
DefaultAssay(AUC) <- 'AUCBinary'
top_N = 5
tops <- name.markers %>% group_by(cluster) %>% top_n(n=top_N, wt=avg_log2FC)
pdf(paste0(OUT_DIR, NAME, "_Regulon_AUCbinary_heatmap_", CLUSTERING, ".pdf"), h=7, w=15)
DoHeatmap(AUC, features = tops$gene, slot = 'data') + scale_fill_gradientn(colors = c( "white", "pink"))
dev.off()





sink(paste0(OUT_DIR, "sessionInfo.txt"))
sessionInfo()
sink()


