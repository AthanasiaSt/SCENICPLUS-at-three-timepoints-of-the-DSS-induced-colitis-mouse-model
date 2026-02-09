# Code for Figure 4 SUPP

#libraries 
library(dplyr)
library(tidyr)
library(tidyverse)
library(pheatmap)
library(viridis)
library(eulerr)


# Download processed data from Zenodo: Processed_datasets_for_scRNA_scATAC_scenicplus 

# set working directory
var<-getwd() 
setwd(paste0(var,'/Processed_datasets_for_scRNA_scATAC_scenicplus/'))


#filtered regulons, can be found in zenodo
regulons_df<-read.table('./SCENICPLUS_Healthy_Inflammation_Regeneration_regulons_results/Direct_Positive_eRegulons_filtered_basedOn_GeneBased_AUC.csv', row.names = 1, sep='\t', header = T)

# keep unique combinations
pos_regulons_filter_df <- regulons_df[!duplicated(regulons_df[c("TF", "Gene")]), ]

#Jaccard similarity function
jaccard <- function(x, y){
  x <- unique(x)
  y <- unique(y)
  intersection <- length(intersect(x, y))
  union <- length(x) + length(y) - intersection
  return (intersection / union)
}


#Create a list of genes for each unique TF
genes_by_tf <- split(pos_regulons_filter_df$Gene, pos_regulons_filter_df$TF)
#Get all combinations of TFs, including both (x, y) and (y, x)
tf_combinations <- expand.grid(names(genes_by_tf), names(genes_by_tf), stringsAsFactors = FALSE)

#Initialize a list to store Jaccard results
jaccard_filtered_results <- list()
for (i in 1:nrow(tf_combinations)) {
  TF1 <- tf_combinations[i, 1]
  TF2 <- tf_combinations[i, 2]
  #Get the respective genes for these TFs
  genes_TF1 <- genes_by_tf[[TF1]]
  genes_TF2 <- genes_by_tf[[TF2]]
  #Calculate Jaccard index for the gene sets of TF1 and TF2
  jaccard_index <- jaccard(genes_TF1, genes_TF2)
  #Store the result, labeled by the TF pair
  jaccard_filtered_results[[paste(TF1, "vs", TF2, sep = "_")]] <- jaccard_index
}

#make a dataframe out of the results
jaccard_filtered_results_df <- data.frame(TF_Comparison = names(jaccard_filtered_results), jaccard_index = unlist(jaccard_filtered_results))
rownames(jaccard_filtered_results_df) <- NULL

#divide the column into TF1 TF2
new_jaccard_results_filtered <- separate(jaccard_filtered_results_df,"TF_Comparison", into = c("TF1", "TF2"), sep = "_vs_")

#Building the heatmap
jaccard_filtered_wide <- new_jaccard_results_filtered %>%
  pivot_wider(names_from = TF2, values_from = jaccard_index)

#Convert the result to a matrix
#Exclude the first column to make it a matrix
#Assign the row names from the TF1 column
jaccard_filtered_matrix <- as.matrix(jaccard_filtered_wide[,-1])
rownames(jaccard_filtered_matrix) <- jaccard_filtered_wide$TF1

# make the colours for the heatmap 
breaks <- c(seq(0, 0.2, length.out = 1000)[-1000],
            seq(0.2, 0.8, length.out = 1000)[-1000],
            seq(0.8, 1, length.out = 1000)[-1000])

custom_colors <- turbo(3000)


pdf('Fig4_c_SUPP_Jaccard_similarity_filtered_regulons.pdf', height = 7, width = 7)

pheatmap(jaccard_filtered_matrix,
         clustering_distance_rows = "euclidean",  # or "correlation"
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",          # or "average", "ward.D2", etc.
         main = "Jaccard Index Heatmap",
         breaks = breaks, color = custom_colors)
dev.off()

# Fig SUPP e)
# Euler plot and Fisher's exact test for TF overlaps 

#filtered regulons for each of the three runs //
# -----------------------------
# Load filtered regulons
# -----------------------------

regulons_all <- read.table(
  "./SCENICPLUS_Healthy_Inflammation_Regeneration_regulons_results/Direct_Positive_eRegulons_filtered_basedOn_GeneBased_AUC.csv",
  row.names = 1, sep = "\t", header = TRUE
)

regulons_kinchen <- read.table(
  "./Revisions_1/scenicplus_revisions_kinchen/scplus_pipeline/Snakemake/output/direct_e_regulon_metadata_filtered_basedOn_max_mean_enrichment_across_categories_GeneBased.csv",
  row.names = 1, sep = "\t", header = TRUE
)

regulons_inhouse <- read.table(
  "./Revisions_1/scenicplus_revisions_D0_D14/scplus_pipeline/Snakemake/output/direct_e_regulon_metadata_filtered_basedOn_max_mean_enrichment_across_categories_GeneBased.csv",
  row.names = 1, sep = "\t", header = TRUE
)


# Load Seurat object
seurat_obj <- readRDS(
  "./scRNA_seurat_integration_Healthy_Inflammation_Regeneration/Seurat5_Harmony_integration_all_samples_Healthy_Inflammation_Regeneration_fibroblasts.rds"
)

# Load motif universe
motifs <- read.table(
  "/home/astavropoulou/SCENICPLUS_runs/scenicplus_D7_D14_D0_fibroblasts/aertslab_motif_colleciton/v10nr_clust_public/snapshots/motifs-v10-nr.mgi-m0.00001-o0.0.tbl",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  comment.char = ""
)

TF_universe <- unique(motifs$gene_name)

# Conservative, biologically defensible background:
TF_background <- intersect(TF_universe, rownames(seurat_obj))
N_background <- length(TF_background)

# Extract unique TFs per run
# -----------------------------
TF_all      <- unique(regulons_all$TF)
TF_kinchen  <- unique(regulons_kinchen$TF)
TF_inhouse  <- unique(regulons_inhouse$TF)

# Fisher overlap test
# -----------------------------
fisher_overlap_test <- function(setA, setB, background) {
  
  setA <- intersect(setA, background)
  setB <- intersect(setB, background)
  
  overlap <- length(intersect(setA, setB))
  onlyA   <- length(setA) - overlap
  onlyB   <- length(setB) - overlap
  neither <- length(background) - (overlap + onlyA + onlyB)
  
  contingency_table <- matrix(
    c(overlap, onlyA,
      onlyB, neither),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      SetA = c("In_SetB", "Not_in_SetB"),
      SetB = c("In_SetA", "Not_in_SetA")
    )
  )
  
  fisher.test(contingency_table, alternative = "greater")
}

test_all_kinchen <- fisher_overlap_test(TF_all, TF_kinchen, TF_background)
test_all_inhouse <- fisher_overlap_test(TF_all, TF_inhouse, TF_background)
test_kinchen_inhouse <- fisher_overlap_test(TF_kinchen, TF_inhouse, TF_background)

# Format p-values for reporting
p_all_kinchen_txt <- paste0("p = ", formatC(test_all_kinchen$p.value, format = "e", digits = 2))
p_all_inhouse_txt <- paste0("p = ", formatC(test_all_inhouse$p.value, format = "e", digits = 2))
p_kinchen_inhouse_txt <- paste0("p = ", formatC(test_kinchen_inhouse$p.value, format = "e", digits = 2))


tf_list <- list(
  Run_all     = TF_all,
  Run_kinchen = TF_kinchen,
  Run_inhouse = TF_inhouse
)

fit <- euler(tf_list)

pdf("Fig4_SUPP_f_Euler_TF_overlap_filtered_SCENICplus_runs.pdf", width = 4, height = 4)
plot(
  fit,
  fills = list(fill = c("#E69F00", "#56B4E9", "#009E73"), alpha = 0.6),
  edges = TRUE,
  quantities = TRUE
)
dev.off()

