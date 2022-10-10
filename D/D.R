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

#' Create a heatmap of correlation between receptors and transcription factors
#' 
#' Creates a heatmap of correlation values between receptors and transcription 
#' factors.
#' 
#' @param dom A domino object with network built (build_domino)
#' @param bool A boolean indicating whether the heatmap should be continuous or boolean. If boolean then bool_thresh will be used to determine how to define activity as positive or negative.
#' @param bool_thresh A numeric indicating the threshold separating 'on' or 'off' for feature activity if making a boolean heatmap.
#' @param title Either a string to use as the title or a boolean describing whether to include a title. In order to pass the 'main' parameter to NMF::aheatmap you must set title to FALSE.
#' @param feats Either a vector of features to include in the heatmap or 'all' for all features. If left NULL then the features selected for the signaling network will be shown.
#' @param recs Either a vector of receptors to include in the heatmap or 'all' for all receptors. If left NULL then the receptors selected in the signaling network connected to the features plotted will be shown.
#' @param mark_connections A boolean indicating whether to add an 'x' in cells where there is a connected receptor or TF. Default FALSE.
#' @param ... Other parameters to pass to NMF::aheatmap. Note that to use the 'main' parameter of NMF::aheatmap you must set title = FALSE and to use 'annCol' or 'annColors' ann_cols must be FALSE.
#' @export
#' 
cor_heatmap = function(dom, bool = TRUE, bool_thresh = .15, title = TRUE, 
    feats = NULL, recs = NULL, mark_connections = FALSE, ...){
    mat = dom@cor

    if(bool){
        cp = mat
        cp[which(mat >= bool_thresh)] = 1
        cp[which(mat < bool_thresh)] = 0
        mat = cp
    }

    if(title == TRUE){
        title = 'Correlation of features and receptors'
    }

    if(is.null(feats)){
        feats = c()
        links = dom@linkages$clust_tf
        for(i in links){
            feats = c(feats, i)
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

    if(is.null(recs)){
        recs = c()
        links = dom@linkages$tf_rec
        for(feat in feats){
            feat_recs = links[[feat]]
            if(length(feat_recs) > 0){
                recs = c(recs, feat_recs)
            }
        }
        recs = unique(recs)
    } else if(recs == 'all'){
        recs = rownames(mat)
    }

    recs=c("B2M", "NTRK2", "NOTCH1", "ERBB4", "FGFR3", "CXCR4", "NRP2", "NRP1", "EPHA4", "NCAM1", "CD47", "EPHA7", "GFRA1", "INSR", "IFNGR1", "EFNB2", "GFRA2", "NTRK3", "EPHA5",  recs) %>% unique()
    mat = mat[recs, feats]

    if(mark_connections){
        cons = mat
        cons[] = ''
        for(feat in feats){
            feat_recs = dom@linkages$tf_rec[[feat]]
            if(length(feat_recs)){
                cons[feat_recs, feat] = 'X'
            }
        }
    }

    if(title != FALSE & mark_connections){
        NMF::aheatmap(mat, main = title, txt = cons, Rowv = NA,
  Colv = NA, ...)
    } else {
        NMF::aheatmap(mat, ...)
    }
    return(list(feats=feats, mat=mat))
}

#NMF::aheatmap(mat, Colv = NA, annCol = ac, annColors = cols, main = title, Rowv = NA, ...)


# BASE="/home/zjr863/Proj3_hypothalamus_ref/gw10_new/network/"
# WORK_DIR=${BASE}"D/output/"
# mkdir -p ${BASE}D/data ${BASE}D/notebooks ${BASE}D/output
# mkdir -p ${WORK_DIR}Global
# mkdir -p ${WORK_DIR}Network_single_cluster
# CLUSTERS=("0"  "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29" "30")
# for clust in ${CLUSTERS[@]}; do
# 	echo ${clust}
#  	mkdir -p ${WORK_DIR}FeaturePlots_single_cluster/${clust}
# 	mkdir -p ${WORK_DIR}Incoming_ligands_single_cluster/${clust};
#  done

# 

# https://github.com/chris-cherry/domino
DIR <- "/home/zjr863/Proj3_hypothalamus_ref/gw10_new/network/"
OUT_DIR <- paste0(DIR, "D/output/")
clr_col <- "SCT_snn_res.1.3_modi"
#CLUSTER_ORDER=c("27", "30", "14", "17", "25", "10", "3", "29", "1", "8", "19", "6", "7", "2", "4", "5", "9", "11", "12", "13", "15", "16", "18", "20", "21", "22", "23", "24", "26", "28", "0")
#top_N <- 5
CLUSTER_ORDER=c("27", "30", "14", "17", "25", "10", "3", "22", "18", "20", "9", "4", "26", "13", "16", "2", "28", "5", "29", "11", "24", "8", "1", "12", "19", "7", "6", "0", "15", "23", "21")
top_N=1000

#INFILES
print(paste0(DIR, "gw10_processed.RDS"))
print(paste0(DIR, 'B/output/auc_mtx.csv'))

#OUTFILES
print(paste0(OUT_DIR, "Global/Domino.obj.RDS"))
print(paste0(OUT_DIR, "Global/signaling_network.pdf"))
print(paste0(OUT_DIR, "Global/heatmap_TFscores.pdf"))
print(paste0(OUT_DIR, "Global/Correlation_TF_receptor.pdf"))
print(paste0(OUT_DIR, "Global/TF_receptor_ligand_network.pdf"))
#CLUSTERS=("0"  "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29" "30")
#for clust in ${CLUSTERS[@]}; do echo $clust; done
#print(paste0(OUT_DIR, "Network_single_cluster/Network_", clust, ".pdf"))
#print(paste0(OUT_DIR, "FeaturePlots_single_cluster/", clust, "/TFs_", clust, ".pdf"))
#print(paste0(OUT_DIR, "FeaturePlots_single_cluster/", clust, "/Receptors_", clust, ".pdf"))
#print(paste0(OUT_DIR, "FeaturePlots_single_cluster/", clust, "/Ligands_", clust, ".pdf"))
#print(paste0(OUT_DIR, "Incoming_ligands_single_cluster/", clust, "/Incoming_ligand_heatmap_", clust, ".pdf"))
#print(paste0(OUT_DIR, "Incoming_ligands_single_cluster/", clust, "/Incoming_ligands_matrix_", clust, ".csv"))
paste0(OUT_DIR, "C_sessionInfo.txt")



library(Seurat)
library(domino)
library(dplyr)


# 1. Import data
ser = readRDS(paste0(DIR, "gw10_processed.RDS"))


# 2. Create the Domino project

### build and visualize the global signaling network
ser@meta.data[, clr_col] <- droplevels(ser@meta.data[, clr_col])
Idents(ser) <- ser@meta.data[, clr_col]
z_scores = ser@assays$SCT@scale.data
counts = ser@assays$SCT@counts[rownames(z_scores),]
ser@meta.data[ , clr_col] <- ser@active.ident
clusters = ser@meta.data[ , clr_col] #active.ident

auc = t(read.table(paste0(DIR, 'B/output/auc_mtx.csv'), header = TRUE, row.names = 1, 
    stringsAsFactors = FALSE, sep = ','))
rownames(auc) <- gsub("[...]", "", rownames(auc))

DERs <- read.csv(paste0(DIR, "C/output/DER/SCENIC_DER_marker_table_", clr_col,".csv"),
                 header=TRUE, sep=";")

dom = create_domino(signaling_db = '/scratch/yuan/accessory/cellphonedb/', 
    features = auc, counts = counts, z_scores = z_scores, clusters = clusters, 
    df = paste0(DIR, 'B/output/regulons.csv'))

#There are four major parameters you can play with when buidling the signaling network, two for selecting TFs for each cluster and two for connecting TFs to recs. 
#min_tf_pval is the minimum p val for a TF to be assigned to a cluster and max_tf_per_clust is the maximum number of transcription factors allowed in each cluster. 
#The same patter is true for rec_tf_cor_threshold and max_rec_per_tf except the thresholding is on the Pearon correlation between receptor and transcription factor. 
#Building the signaling network takes very little time so feel free to play around with these values. 
#In general we prefer to select a pval and cor threshold such that a good portion of the TFs and recs are not being trimmed by the maximum number thresholds.
#Here, the TFs were chosen based on DEG test, only DEGs with p<min_tf_pval are kept as cluster-specific TF when building subnetwork for each cluster separately
#rec_tf_cor_threshold: the receptors were selected based on their correlation to regulon activity, only receptors with rec_tf_cor>0.15 were kept here.
#
dom = build_domino(dom, max_tf_per_clust = 1000, 
    min_tf_pval = .00001, max_rec_per_tf = 1000, rec_tf_cor_threshold = .15)

saveRDS(dom, paste0(OUT_DIR, "Global/Domino.obj.RDS"))


print("Check the dom object:")
print("dom@db_info:")
dom@db_info

print("dom@z_scores[1:3, 1:5]:")
dom@z_scores[1:3, 1:5]
dim(dom@z_scores)

print("dom@counts[1:3, 1:5]:")
dom@counts[1:3, 1:5]
dim(dom@counts)

print("dom@clusters: ")
dom@clusters

print("dom@features[1:3, 1:5]: ")
dom@features[1:3, 1:5]
dim(dom@features)

print("dom@cor[1:3, 1:5] (col: TF, row: receptors)")
dom@cor[1:3, 1:5]
dim(dom@cor)

print("dom@linkages$rec_lig[1:3]***")
dom@linkages$rec_lig[1:3]
length(dom@linkages$rec_lig)

print("dom@linkages$tf_targets[1:3]***")
dom@linkages$tf_targets[1:3]
length(dom@linkages$tf_targets)

print("dom@linkages$clust_tf[1:2]***, Here are TFs specific to each cluster, they are selected based on dom@clust_de with p<0.001 here")
dom@linkages$clust_tf[1:2]
length(dom@linkages$clust_tf)

print("dom@linkages$tf_rec[1:3]**")
dom@linkages$tf_rec[1:3]
length(dom@linkages$tf_rec)

print("dom@clust_de[1:3, 1:5]")
dom@clust_de[1:3, 1:5]
dim(dom@clust_de)

print("dom@misc$tar_lr_cols")
dom@misc$tar_lr_cols

print("dom@misc$create")
dom@misc$create

print("dom@misc$build")
dom@misc$build

print("head(dom@misc$rl_map)")
head(dom@misc$rl_map)

print("dom@misc$build_vars")
dom@misc$build_vars

print("dom@cl_signaling_matrices$`0`[1:3, 1:5] (This is incoming ligand signaling heatmap for each cluster based on)")
dom@cl_signaling_matrices$`0`[1:3, 1:5]

print("dom@signaling[1:3, 1:5] (this is globle intercellular network based on)")
dom@signaling[1:3, 1:5]
dim(dom@signaling)




### Visualize the global intercellular signaling network
#We are using a max threshold here (max_thresh) because the Mk cluster's expression of ITGA2B is so much higher than the other clusters/ligands 
#it drowns other relevant signaling out. Give it a shot without the max_thresh and you'll see what I mean.
pdf(paste0(OUT_DIR, "Global/signaling_network.pdf"), h=12, w=15)
	signaling_network(dom, edge_weight = .5, max_thresh = 3.5)
dev.off()
#The edges are weighted based on the strength of signaling between two specific clusters. 
#The color of the edge will match the color of the ligand cluster. 

### Generate a heatmap of the TF activation scores
#This heatmap is based on regulon heatmap scores
names(dom@clusters) <- colnames(dom@features)
#tops <- DERs %>% group_by(cluster) %>% top_n(n=top_N, wt=avg_log2FC) %>% unique()
pdf(paste0(OUT_DIR, "Global/heatmap_TFscores.pdf"), h=16, w=16)
	features=feat_heatmap(dom, norm = TRUE, bool = FALSE, CLUSTER_ORDER=CLUSTER_ORDER, top_n=5)#, bool_thresh=0.2) #median(dom@features))
dev.off()

### Create a correlation heatmap between TFs and receptors
dom@features=dom@features[features,]
pdf(paste0(OUT_DIR, "Global/Correlation_TF_receptor.pdf"), h=10, w=29)
	features2=cor_heatmap(dom, bool = FALSE, mark_connections = TRUE, feats=features)
dev.off()

### Create a global tf-r-l network
pdf(paste0(OUT_DIR, "Global/TF_receptor_ligand_network.pdf"), h=40, w=50)
	info = gene_network(dom, clust = levels(dom@clusters), 
    lig_scale = FALSE, layout = 'fr')
	plot(info$graph, layout = info$layout, vertex.size = 3, edge.color = 'grey', 
    	vertex.frame.color = 'black') #, vertex.label = NA)
dev.off()

print(paste0(OUT_DIR, "Global/Intercellular_network_matrix.csv"))
write.csv2(dom@signaling, paste0(OUT_DIR, "Global/Intercellular_network_matrix.csv")) 




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
  print(paste0(OUT_DIR, "FeaturePlots_single_cluster/", clust, "/TFs_", clust, ".pdf"))
  pdf(paste0(OUT_DIR, "FeaturePlots_single_cluster/", clust, "/TFs_", clust, ".pdf"), h=4, w=6)
  	for (tf in tf_clust) {
  		print(tf)
  		p=FeaturePlot(ser, features=tf)
  		print(p)
  	}  
  dev.off()

  print(paste0(OUT_DIR, "FeaturePlots_single_cluster/", clust, "/Receptors_", clust, ".pdf"))
  pdf(paste0(OUT_DIR, "FeaturePlots_single_cluster/", clust, "/Receptors_", clust, ".pdf"), h=4, w=6)
  	for (rec in rec_clust) {
  		print(rec)
  		print(FeaturePlot(ser, features=rec))
  	}
   dev.off()

  print(paste0(OUT_DIR, "FeaturePlots_single_cluster/", clust, "/Ligands_", clust, ".pdf"))
  pdf(paste0(OUT_DIR, "FeaturePlots_single_cluster/", clust, "/Ligands_", clust, ".pdf"), h=4, w=6)
  	for (lig in lig_clust) {
  		print(lig)
  		print(FeaturePlot(ser, features=lig))
  	}
  dev.off()

  
  print(paste0(OUT_DIR, "Incoming_ligands_single_cluster/", clust, "/Incoming_ligand_heatmap_", clust, ".pdf"))
  pdf(paste0(OUT_DIR, "Incoming_ligands_single_cluster/", clust, "/Incoming_ligand_heatmap_", clust, ".pdf"), h=6, w=18)
  	incoming_signaling_heatmap(dom, rec_clust = clust, max_thresh = 2.5)
  dev.off() 

  print(paste0(OUT_DIR, "Incoming_ligands_single_cluster/", clust, "/Incoming_ligands_matrix_", clust, ".csv"))
  write.csv2(dom@cl_signaling_matrices[[clust]], paste0(OUT_DIR, "Incoming_ligands_single_cluster/", clust, "/Incoming_ligands_matrix_", clust, ".csv")) 
}




sink(paste0(OUT_DIR, "C_sessionInfo.txt"))
sessionInfo()
sink()


