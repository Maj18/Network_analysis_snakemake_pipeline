# Based on regulons to look for targets


DIR <- "/home/zjr863/Proj3_hypothalamus_ref/gw10_new/network/"
OUT_DIR <- paste0(DIR, "C/output/TF2target/")
NAME <- "TF2targets"
CLUSTERING <- "SCT_snn_res.1.3_modi"
regulon.names <- c("LHX2")
clr_top_DERs <- c("15")
clr_DEGs_for_regulon_targets <- c("15")


# 1. Define functions
visuNetwork <- function(adj, degs, regulons, regulon.name, clusters, top=50) {
  adj.ls <- regulon.name %>% map(~{
    tmp <- adj[which(adj$TF == .x), ]
    tmp <- tmp[order(tmp$importance, decreasing=T), ]
    loci <- unique(c(1:1000, grep(.x, tmp$target)))
    tmp <- tmp[loci, ]
    return(tmp)
  })
  adj.sub <- Reduce('rbind', adj.ls)
  ### generate network
  edge.df <- adj.sub
  colnames(edge.df) <- c("from", "to", "weights")
  edge.df$from <- as.character(edge.df$from)
  edge.df$to <- as.character(edge.df$to)
  
  reg <- subset(regulons, X%in%regulon.name)
  motifs = stringr::str_split(reg$Enrichment.6, "[']")
  targets <- lapply(motifs, function(motif) motif[(1:abs(length(motif)/2)) *2]) %>% unlist() %>% unique()
  edge.df <- subset(edge.df, to%in%targets)
  
  df <- edge.df
  
  
  # all top targes
  if (nrow(edge.df)>0) {
    top_n <- min(top, nrow(edge.df))
    edge.df <- edge.df[1:top_n, ]
    #net1 <- network(edge.df, directed=T, loops=T)
    #plot(net1)
    
    vertex <- unique(c(edge.df$from, edge.df$to))
    net <- graph_from_data_frame(d=edge.df, vertices=vertex, directed=T)
    
    ## vertex size scaled to weights
    vsize <- scale(edge.df$weights, center=min(quantile(edge.df$weights)), scale=max(quantile(edge.df$weights))-min(quantile(edge.df$weights)))
    
    plot(
      net,
      #vertex.label=vlabel,
      vertex.label.cex=1,
      vertex.label.dist=1,
      vertex.label.color="black",
      vertex.size=c(1, vsize+0.1)*10,
      vertex.frame.color=NA,
      edge.arrow.size=0,
      main=regulon.name
    )    
  }
  
  
  # Only for targets that are also DEGs of clusters
  degs_sub <- subset(degs, cluster%in%clusters)
  degs_targets <- intersect(df$to, degs_sub$gene)
  if (length(degs_targets) >1) {
    edge.degs <- subset(df, to%in%degs_targets[1:top])
    vertex.degs <- unique(c(edge.degs$from, edge.degs$to))
    net.degs <- graph_from_data_frame(d=edge.degs, vertices=vertex.degs, directed=T)
    
    ## vertex size scaled to weights
    vsize.degs <- scale(edge.degs$weights, center=min(quantile(edge.degs$weights)), 
                        scale=max(quantile(edge.degs$weights))-min(quantile(edge.degs$weights)))
    
    plot(
      net.degs,
      #vertex.label=vlabel,
      vertex.label.cex=1,
      vertex.label.dist=1,
      vertex.label.color="black",
      vertex.size=c(1, vsize.degs+0.1)*10,
      vertex.frame.color=NA,
      edge.arrow.size=0,
      main=regulon.name
    )    
  }
  return(df)
}



library(Seurat)
library(tidyverse)
library(magrittr)
library(viridis)
#library(SCENIC)
library(network)
library(igraph)


# 2. Load data
name <- readRDS(paste0(DIR, "gw10_processed.RDS"))
degs <- read.csv(paste0(DIR, "gw10_marker_table_cellType.res.1.3_modi.csv"),
                 header=TRUE, sep=";")
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


# 3. Look for targets of "reugon.name" in cluster clr, only top 50 targets shown, for each regulon, two plots will be generated, one is top 50 among all targets, and the other is for top50 among targets that are also differentially expressed genes in clr
if (length(clr_top_DERs)>0) {
	top_regulons <- DERs[DERs$cluster %in% clr_top_DERs, ]$gene
	pdf(paste0(OUT_DIR, NAME, "_Regulon_target_network_for_top_DERs_in", clr_top_DERs, "_clr_DEGs_for_regulon_targets", clr_DEGs_for_regulon_targets, ".pdf"), h=8, w=8)
	B=lapply(top_regulons, function(t) {
  		df <- visuNetwork(adj, degs, regulons, regulon.name=t, clusters=clr_DEGs_for_regulon_targets, top=50)
  	})
	dev.off()
}


if (length(regulon.names) > 0) {
	pdf(paste0(OUT_DIR, NAME, "_Regulon_target_network_for_", regulon.names, "_clr_DEGs_for_regulon_targets", clr_DEGs_for_regulon_targets, ".pdf"), h=8, w=8)
	B=lapply(regulon.names, function(t) {
  		df <- visuNetwork(adj, degs, regulons, regulon.name=t, clusters=clr_DEGs_for_regulon_targets, top=50)
  	})
	dev.off()
}




sink(paste0(OUT_DIR, "sessionInfo.txt"))
sessionInfo()
sink()




