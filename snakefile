#To run the file:
#conda activate yuan2

feathers=["hg19-500bp-upstream-10species.mc9nr", "hg19-500bp-upstream-7species.mc9nr", "hg19-tss-centered-10kb-10species.mc9nr", "hg19-tss-centered-10kb-7species.mc9nr", "hg19-tss-centered-5kb-10species.mc9nr", "hg19-tss-centered-5kb-7species.mc9nr"]
CLUSTERS = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30"]

rule all:
	input:
		"C/output/DER/sessionInfo.txt"

# A. Make the input count matrix
rule A:
	input:
		"ARC_neuron_processed.RDS"
	output:
		"A/output/counts.tsv"
	shell:
		"Rscript A/notebooks/A.R $PWD 'ARC_neuron'"

# B. SCENIC
rule B1:
	input:
		"A/output/counts.tsv"
	output:
		"B/output/filtered_scenic.loom"
	shell:
		"python B/notebooks/S1.py $PWD"


rule B2:
	input:
		f_tfs="/scratch/yuan/pyscenicdata/hs_hgnc_curated_tfs.txt",
		f_loom_path_scenic="B/output/filtered_scenic.loom"
	output:
		adjacencies="B/output/adjacencies.csv"	
	shell:
		"~/.conda/envs/yuan2/bin/pyscenic grn {input.f_loom_path_scenic} {input.f_tfs} -o {output.adjacencies} --num_workers 10"


rule B3:
	input:
		adjacencies="B/output/adjacencies.csv",
		f_loom_path_scenic="B/output/filtered_scenic.loom",
		f_db_names=expand("/scratch/yuan/pyscenicdata/{fea}.feather", fea=feathers),
		f_motif_path="/scratch/yuan/pyscenicdata/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"
	output:
		regulon="B/output/regulons.csv"
	shell:
		"~/.conda/envs/yuan2/bin/pyscenic ctx {input.adjacencies} {input.f_db_names} --annotations_fname {input.f_motif_path} --expression_mtx_fname {input.f_loom_path_scenic} --mode 'dask_multiprocessing' --output {output.regulon} --mask_dropouts --num_workers 10"


rule B4:
	input:
		f_loom_path_scenic="B/output/filtered_scenic.loom",
		regulon="B/output/regulons.csv"
	output:
		output_auc_mtx="B/output/auc_mtx.csv"
	shell:
		"~/.conda/envs/yuan2/bin/pyscenic aucell {input.f_loom_path_scenic} {input.regulon} --output {output.output_auc_mtx} --num_workers 10"


rule B5:
	input:
		f_loom_path_scenic="B/output/filtered_scenic.loom",
		regulon="B/output/regulons.csv"
	output:
		f_pyscenic_output="B/output/pyscenic_output.loom"
	shell:
		"~/.conda/envs/yuan2/bin/pyscenic aucell {input.f_loom_path_scenic} {input.regulon} --output {output.f_pyscenic_output} --num_workers 10"


rule B6:
	input:
		"B/output/auc_mtx.csv"
	output:
		"B/output/auc_mtx_binary.csv"
	shell:
		"python B/notebooks/B6.py $PWD"



# C. Reulatory network
### Differentially enriched regulons
rule C1:
	input:
		"ARC_neuron_processed.RDS",
		"ARC_neuron_marker_table_cellType.res.0.7.csv",
		"B/output/auc_mtx.csv",
		"B/output/auc_mtx_binary.csv",
		"B/output/regulons.csv",
		"B/output/adjacencies.csv"
	output:
		"C/output/DER/SCENIC_DER_marker_table_SCT_snn_res.0.7.csv",
		"C/output/DER/SCENIC_DER_Regulon_AUC_heatmap_SCT_snn_res.0.7.pdf",
		"C/output/DER/SCENIC_DER_Regulon_AUC_heatmap_SCT_snn_res.0.7_small.pdf",
		"C/output/DER/SCENIC_DER_cluster_top_regulonAUC_SCT_snn_res.0.7.pdf",
		"C/output/DER/SCENIC_DER_cluster_top_geneExpr_SCT_snn_res.0.7.pdf",
		"C/output/DER/SCENIC_DER_Regulon_AUCbinary_heatmap_SCT_snn_res.0.7.pdf",
		"C/output/DER/sessionInfo.txt"
	shell:
		"Rscript C/notebooks/C1.R $PWD 'ARC_neuron'"

# Target2TF
rule C2:
	input:
		"gw10_processed.RDS",
		"C/output/DER/SCENIC_DER_marker_table_SCT_snn_res.1.3_modi.csv",
		"B/output/auc_mtx.csv",
		"B/output/auc_mtx_binary.csv",
		"B/output/regulons.csv",
		"B/output/adjacencies.csv"
	output:
	     "C/output/target2TF/From_LHX92TF_in_cluster_15.txt"
	shell:
		"Rscript C/notebooks/C2.R"


# TF2target
rule C3:
	input:
		"gw10_processed.RDS",
		"gw10_marker_table_cellType.res.1.3_modi.csv",
		"C/output/DER/SCENIC_DER_marker_table_SCT_snn_res.1.3_modi.csv",
		"B/output/auc_mtx.csv",
		"B/output/auc_mtx_binary.csv",
		"B/output/regulons.csv",
		"B/output/adjacencies.csv"
	output:
	     "C/output/TF2target/TF2targets_Regulon_target_network_for_top_DERs_in15_clr_DEGs_for_regulon_targets15.pdf",
	     "C/output/TF2target/TF2targets_Regulon_target_network_for__clr_DEGs_for_regulon_targets15.pdf",
	     "C/output/TF2target/sessionInfo.txt"
	shell:
		"Rscript C/notebooks/C3.R"





# D. Cell-cell communication
rule D:
	input:
		"gw10_processed.RDS",
		"B/output/auc_mtx.csv"
	output:
		"D/output/Global/Domino.obj.RDS",
		"D/output/Global/signaling_network.pdf",
		"D/output/Global/heatmap_TFscores.pdf",
		"D/output/Global/Correlation_TF_receptor.pdf",
		"D/output/Global/TF_receptor_ligand_network.pdf",
		expand("Network_single_cluster/Network_{clust}.pdf", clust=CLUSTERS),
		expand("FeaturePlots_single_cluster/{clust}/TFs_{clust}.pdf", clust=CLUSTERS),
		expand("FeaturePlots_single_cluster/{clust}/Receptors_{clust}.pdf", clust=CLUSTERS),
		expand("FeaturePlots_single_cluster/{clust}/Ligands_{clust}.pdf", clust=CLUSTERS),
		expand("Incoming_ligands_single_cluster/{clust}/Incoming_ligand_heatmap_{clust}.pdf", clust=CLUSTERS),
		expand("Incoming_ligands_single_cluster/{clust}/Incoming_ligands_matrix_{clust}.pdf", clust=CLUSTERS),
		"D/output/C_sessionInfo.txt"		
	shell:
		"Rscript D/notebooks/D.R"






