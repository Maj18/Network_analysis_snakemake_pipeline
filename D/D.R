#' Create a heatmap of features organized by cluster
#' 
#' Creates a heatmap of feature expression (typically transcription factor
#' activation scores) by cells organized by cluster.
#' 
#' @param dom A domino object with network built (build_domino)
#' @param bool A boolean indicating whether the heatmap should be continuous or boolean. If boolean then bool_thresh will be used to determine how to define activity as positive or negative.
#' @param bool_thresh A numeric indicating the threshold separating 'on' or 'off' for feature activity if making a boolean heatmap.
#' @param title Either a string to use as the title or a boolean describing whether to include a title. In order to pass the 'main' parameter to NMF::aheatmap you must set title to FALSE.
#' @param norm Boolean indicating whether or not to normalize the transcrption factors to their max value.
#' @param feats Either a vector of features to include in the heatmap or 'all' for all features. If left NULL then the features selected for the signaling network will be shown.
#' @param ann_cols Boolean indicating whether to include cell cluster as a column annotation. Colors can be defined with cols. If FALSE then custom annotations can be passed to NMF.
#' @param cols A named vector of colors to annotate cells by cluster color. Values are taken as colors and names as cluster. If left as NULL then default ggplot colors will be generated.
#' @param min_thresh Minimum threshold for color scaling if not a boolean heatmap
#' @param max_thresh Maximum threshold for color scaling if not a boolean heatmap
#' @param ... Other parameters to pass to NMF::aheatmap. Note that to use the 'main' parameter of NMF::aheatmap you must set title = FALSE and to use 'annCol' or 'annColors' ann_cols must be FALSE.
#' @export
#' 
feat_heatmap = function(dom, feats = NULL, bool = TRUE, bool_thresh = .2, 
                        title = TRUE, norm = FALSE, cols = NULL, ann_cols = TRUE, min_thresh = NULL, 
                        max_thresh = NULL, CLUSTER_ORDER=NULL, top_n=5, ...){
  if(!length(dom@clusters)){
    warning("This domino object wasn't build with clusters. Cells will not be ordered.")
    ann_cols = FALSE
  }
  mat = dom@features
  cl = dom@clusters
  if (!is.null(CLUSTER_ORDER)) {
    cl = cl[order(match(cl, CLUSTER_ORDER))]
  } else {cl=sort(cl)}
  
  
  if(norm & (!is.null(min_thresh) | !is.null(max_thresh))){
    warning('You are using norm with min_thresh and max_thresh. Note that values will be thresholded AFTER normalization.')
  }
  
  if(norm){
    mat = domino:::do_norm(mat, 'row')
  }
  
  if(!is.null(min_thresh)){
    mat[which(mat < min_thresh)] = min_thresh
  }
  if(!is.null(max_thresh)){
    mat[which(mat > max_thresh)] = max_thresh
  }
  
  if(bool){
    cp = mat
    cp[which(mat >= bool_thresh)] = 1
    cp[which(mat < bool_thresh)] = 0
    mat = cp
  }
  
  if(title == TRUE){
    title = 'Feature expression by cluster'
  }
  
  if(is.null(feats)){
    feats = c()
    links = dom@linkages$clust_tf
    for(i in CLUSTER_ORDER){
      feats = c(feats, links[[i]][1:top_n])
    }
    feats = unique(feats)
  } else if(feats[1] != 'all'){
    mid = match(feats, rownames(dom@features))
    na = which(is.na(mid))
    na_feats = paste(feats[na], collapse = ' ')
    if(length(na) != 0){
      print(paste('Unable to find', na_feats))
      feats = feats[-na]
    } 
  } else if(feats == 'all'){
    feats = rownames(mat)
  }
  
  if(length(cl)){
    mat = mat[feats, names(cl)]
  }
  
  if(ann_cols){
    ac = list('Cluster' = cl)
    names(ac[[1]]) = c()
    if(is.null(cols)){
      cols = domino:::ggplot_col_gen(length(levels(cl)))
      names(cols) = levels(cl)
    }
    cols = list('Cluster' = cols)
  }
  
  if(title != FALSE & ann_cols != FALSE){
    NMF::aheatmap(mat, Colv = NA, annCol = ac, annColors = cols, main = title, Rowv = NA, ...)
  } else if(title == FALSE & ann_cols != FALSE){
    NMF::aheatmap(mat, Colv = NA, annCol = ac, annColors = cols, ...)
  } else if(title != FALSE & ann_cols == FALSE){
    NMF::aheatmap(mat, Colv = NA, main = title, ...)
  } else if(title == FALSE & ann_cols == FALSE){
    NMF::aheatmap(mat, Colv = NA, ...)
  }
  return(feats)
}



# https://github.com/chris-cherry/domino
DIR <- "Proj2_VLMC/FGF_2022/2022-02/Double/network_FGF17/"
OUT_DIR <- paste0(DIR, "D_domino/output/")
clr_col <- "SCT_snn_res.0.1"

#CLUSTER_ORDER=c("27", "30", "14", "17", "25", "10", "3", "22", "18", "20", "9", "4", "26", "13", "16", "2", "28", "5", "29", "11", "24", "8", "1", "12", "19", "7", "6", "0", "15", "23", "21")
top_N=10




library(Seurat)
library(domino)
library(dplyr)


# 1. Import data
name=readRDS("/home/zjr863/Proj2_VLMC/FGF_2022/2022-02/Double/seurat_basics/output/FGFprocessed.RDS")
ser <- subset(name, sample%in%"FGF17")
DefaultAssay(ser) <- "RNA"
ser <- NormalizeData(ser)
ser <- ScaleData(ser)


# 2. Create the Domino project
### build and visualize the global signaling network
ser@meta.data[, clr_col] <- droplevels(ser@meta.data[, clr_col])
Idents(ser) <- ser@meta.data[, clr_col]
z_scores = ser@assays$RNA@scale.data
counts = ser@assays$RNA@counts[rownames(z_scores),]
ser@meta.data[ , clr_col] <- ser@active.ident
clusters = ser@meta.data[ , clr_col] #active.ident

auc = t(read.table(paste0(DIR, 'B_SCENIC/output/auc_mtx.csv'), header = TRUE, row.names = 1, 
                   stringsAsFactors = FALSE, sep = ','))
rownames(auc) <- gsub("[...]", "", rownames(auc))

dom_pre = create_domino(signaling_db = '/scratch/yuan/accessory/cellphonedb/', 
                    features = auc, counts = counts, z_scores = z_scores, clusters = clusters, 
                    df = paste0(DIR, 'B_SCENIC/output/regulons.csv'))

saveRDS(dom_pre, paste0(OUT_DIR, "Global/dom_pre.RDS"))

#There are four major parameters you can play with when buidling the signaling network, two for selecting TFs for each cluster and two for connecting TFs to recs. 
#min_tf_pval is the minimum p val for a TF to be assigned to a cluster and max_tf_per_clust is the maximum number of transcription factors allowed in each cluster. 
#The same patter is true for rec_tf_cor_threshold and max_rec_per_tf except the thresholding is on the Pearon correlation between receptor and transcription factor. 
#Building the signaling network takes very little time so feel free to play around with these values. 
#In general we prefer to select a pval and cor threshold such that a good portion of the TFs and recs are not being trimmed by the maximum number thresholds.
#Here, the TFs were chosen based on DEG test, only DEGs with p<min_tf_pval are kept as cluster-specific TF when building subnetwork for each cluster separately
#rec_tf_cor_threshold: the receptors were selected based on their correlation to regulon activity, only receptors with rec_tf_cor>0.15 were kept here.
#
dom = build_domino(dom_pre, max_tf_per_clust = 10, 
                   min_tf_pval = .001, max_rec_per_tf = 10, rec_tf_cor_threshold = .23)

#

### Create a global tf-r-l network
pdf(paste0(OUT_DIR, "Global/TF_receptor_ligand_network.pdf"), h=10, w=15)
info = gene_network(dom, clust = levels(dom@clusters), 
                    lig_scale = FALSE, layout = 'fr')
plot(info$graph, layout = info$layout, vertex.size = 3, edge.color = 'grey', 
     vertex.frame.color = 'black') #, vertex.label = NA)
dev.off()

print(paste0(OUT_DIR, "Global/Intercellular_network_matrix.csv"))
write.csv2(dom@signaling, paste0(OUT_DIR, "Global/Intercellular_network_matrix.csv")) 



### Visualize the global intercellular signaling network
#We are using a max threshold here (max_thresh) because the Mk cluster's expression of ITGA2B is so much higher than the other clusters/ligands 
#it drowns other relevant signaling out. Give it a shot without the max_thresh and you'll see what I mean.
pdf(paste0(OUT_DIR, "Global/signaling_network.pdf"), h=4, w=4)
signaling_network(dom, edge_weight = .5, max_thresh = 2.5)
dev.off()
#The edges are weighted based on the strength of signaling between two specific clusters. 
#The color of the edge will match the color of the ligand cluster. 

### Generate a heatmap of the TF activation scores
#This heatmap is based on regulon heatmap scores
names(dom@clusters) <- colnames(dom@features)
#tops <- DERs %>% group_by(cluster) %>% top_n(n=top_N, wt=avg_log2FC) %>% unique()
pdf(paste0(OUT_DIR, "Global/heatmap_TFscores.pdf"), h=10, w=12)
features=feat_heatmap(dom, norm = FALSE, bool = FALSE, top_n=10, CLUSTER_ORDER = levels(dom@clusters))#, bool_thresh=0.2) #median(dom@features))
dev.off()

### Create a correlation heatmap between TFs and receptors
dom@features=dom@features[features,]
pdf(paste0(OUT_DIR, "Global/Correlation_TF_receptor.pdf"), h=4.5, w=10)
features2=cor_heatmap(dom, bool = FALSE, mark_connections = TRUE, feats=features, Rowv = NA, Colv = NA)
dev.off()



# 3. Check cluster wise signaling network for each cluster, as well as the incoming ligand expr heatmap
#degs <- read.csv2(paste0(DIR, "gw10_marker_table_cellType.res.1.3_modi.csv"))

for (clust in levels(dom@clusters)) {
  print("Now, we are looking at cluster: ")
  print(clust)
  print(paste0(OUT_DIR, "Network_single_cluster/Network_", clust, ".pdf"))
  pdf(paste0(OUT_DIR, "Network_single_cluster/Network_", clust, ".pdf"), h=20, w=30)
  gene_network(dom, clust = clust, layout = 'fr')
  dev.off()
  
  # Diff expr genes
  #degs_clust <- degs$gene[degs$cluster%in%clust]
  #degs_clust <- intersect(degs_clust, rownames(ser))
  #tf_clust are selected based on regulons
  tf_clust <- dom@linkages$clust_tf[[clust]]
  tf_clust <- intersect(tf_clust, rownames(ser))
  rec_clust <- dom@linkages$tf_rec[tf_clust] %>% unlist()
  rec_clust <- intersect(rec_clust, rownames(ser))
  lig_clust <- dom@linkages$rec_lig[rec_clust] %>% unlist()
  lig_clust <- intersect(lig_clust, rownames(ser))
  print(paste0(OUT_DIR, "FeaturePlots_single_cluster/TFs_", clust, ".pdf"))
  pdf(paste0(OUT_DIR, "FeaturePlots_single_cluster//TFs_", clust, ".pdf"), h=4, w=6)
  for (tf in tf_clust) {
    print(tf)
    p=FeaturePlot(ser, features=tf)
    print(p)
  }  
  dev.off()
  
  print(paste0(OUT_DIR, "FeaturePlots_single_cluster/Receptors_", clust, ".pdf"))
  pdf(paste0(OUT_DIR, "FeaturePlots_single_cluster/Receptors_", clust, ".pdf"), h=4, w=6)
  for (rec in rec_clust) {
    print(rec)
    print(FeaturePlot(ser, features=rec))
  }
  dev.off()
  
  print(paste0(OUT_DIR, "FeaturePlots_single_cluster/Ligands_", clust, ".pdf"))
  pdf(paste0(OUT_DIR, "FeaturePlots_single_cluster/Ligands_", clust, ".pdf"), h=4, w=6)
  for (lig in lig_clust) {
    print(lig)
    print(FeaturePlot(ser, features=lig))
  }
  dev.off()
  
  
  print(paste0(OUT_DIR, "Incoming_ligands_single_cluster/Incoming_ligand_heatmap_", clust, ".pdf"))
  pdf(paste0(OUT_DIR, "Incoming_ligands_single_cluster/Incoming_ligand_heatmap_", clust, ".pdf"), h=6, w=18)
  incoming_signaling_heatmap(dom, rec_clust = clust, max_thresh = 2.5)
  dev.off() 
  
  print(paste0(OUT_DIR, "Incoming_ligands_single_cluster/Incoming_ligands_matrix_", clust, ".csv"))
  write.csv2(dom@cl_signaling_matrices[[clust]], paste0(OUT_DIR, "Incoming_ligands_single_cluster/Incoming_ligands_matrix_", clust, ".csv")) 
}




sink(paste0(OUT_DIR, "C_sessionInfo.txt"))
sessionInfo()
sink()
