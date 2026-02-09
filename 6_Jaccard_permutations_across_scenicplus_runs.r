# Code revisions -- jaccard across different runs of scenicplus 

#libraries 
library(dplyr)
library(tidyr)
library(tidyverse)
library(pheatmap)
library(viridis)
library(ggplot2)

# set working directory
var<-getwd() 
setwd(paste0(var,'/Processed_datasets_for_scRNA_scATAC_scenicplus/'))

#filtered regulons across two Runs
# Test main SCENIC+ with Kinchen separate RUN
regulons_df_run_all<-read.table('./SCENICPLUS_Healthy_Inflammation_Regeneration_regulons_results/Direct_Positive_eRegulons_filtered_basedOn_GeneBased_AUC.csv', row.names = 1, sep='\t', header = T)
regulons_df_run_D7<-read.table('./Revisions_1/scenicplus_revisions_kinchen/scplus_pipeline/Snakemake/output/direct_e_regulon_metadata_filtered_basedOn_max_mean_enrichment_across_categories_GeneBased.csv', row.names = 1, sep='\t', header = T)
#regulons_df_run_D14<-read.table('./Revisions_1/scenicplus_revisions_D0_D14/scplus_pipeline/Snakemake/output/direct_e_regulon_metadata_filtered_basedOn_max_mean_enrichment_across_categories_GeneBased.csv', row.names = 1, sep='\t', header = T)

#seurat object
seurat_obj<-readRDS(file = "./scRNA_seurat_integration_Healthy_Inflammation_Regeneration/Seurat5_Harmony_integration_all_samples_Healthy_Inflammation_Regeneration_fibroblasts.rds")


# check jaccard similarity for each TF regulon across runs 
# keep unique combinations
pairs <- expand.grid(
  TF_run1 = unique(regulons_df_run_all$TF),
  TF_run2 = unique(regulons_df_run_D7$TF),
  stringsAsFactors = FALSE
)

#Jaccard similarity function
jaccard <- function(x, y){
  x <- unique(x)
  y <- unique(y)
  intersection <- length(intersect(x, y))
  union <- length(x) + length(y) - intersection
  return (intersection / union)
}

#Create a list of genes for each unique TF
tmp <- regulons_df_run_all[!duplicated(regulons_df_run_all[c("TF", "Gene")]), ]
genes_by_tf_run_all <- split(tmp$Gene, tmp$TF)

tmp <- regulons_df_run_D7[!duplicated(regulons_df_run_D7[c("TF", "Gene")]), ]
genes_by_tf_run_D7 <- split(tmp$Gene, tmp$TF)

#----------------- Real data
#Initialize a list to store Jaccard results
jaccard_filtered_results <- list()
for (i in 1:nrow(pairs)) {

  TF1 <- pairs[i, 1]#pairs[i, 1]
  TF2 <- pairs[i, 2]#pairs[i, 2]
  #Get the respective genes for these TFs
  genes_TF1 <- genes_by_tf_run_all[[TF1]]
  genes_TF2 <- genes_by_tf_run_D7[[TF2]]
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

#----------------------Random permutations
#-------------- Random Jaccard distributions
set.seed(123)

jaccard_RANDOM <- vector("list", length = 10000)

# first keep only unique relationships with each TF-gene / exclude regions from this analysis 
tmp <- regulons_df_run_D7[!duplicated(regulons_df_run_D7[c("TF", "Gene")]), ]

for (j in 1:10000) {
  
  regulons_df_run_D7_RANDOM <- tmp %>%
    dplyr::mutate(Gene = sample(Gene)) %>%
    dplyr::distinct(TF, Gene)
  
  genes_by_tf_run_D7_RANDOM <- split(
    regulons_df_run_D7_RANDOM$Gene,
    regulons_df_run_D7_RANDOM$TF
  )
  
  perm_res <- data.frame(
    TF1 = pairs[, 1],
    TF2 = pairs[, 2],
    Jaccard = NA_real_,
    perm = j
  )
  
  for (i in seq_len(nrow(pairs))) {
    
    TF1 <- pairs[i, 1]
    TF2 <- pairs[i, 2]
    
    genes_TF1 <- genes_by_tf_run_all[[TF1]]
    genes_TF2 <- genes_by_tf_run_D7_RANDOM[[TF2]]
    
    perm_res$Jaccard[i] <- jaccard(genes_TF1, genes_TF2)
  }
  
  jaccard_RANDOM[[j]] <- perm_res
}

jaccard_RANDOM_df <- dplyr::bind_rows(jaccard_RANDOM)


new_jaccard_results_filtered$comparison <- ifelse(
  new_jaccard_results_filtered$TF1 == new_jaccard_results_filtered$TF2,
  "Same TF",
  "Different TFs"
)

jaccard_RANDOM_df$comparison <- ifelse(
  jaccard_RANDOM_df$TF1 == jaccard_RANDOM_df$TF2,
  "Same TF Random",
  "Different TFs Random"
)

jaccard_RANDOM_df$jaccard_index <- jaccard_RANDOM_df$Jaccard

# compute p values for the same TFs pvalues 
new_jaccard_results_filtered_sameTF <- new_jaccard_results_filtered %>%
  dplyr::filter(TF1 == TF2)

jaccard_RANDOM_sameTF <- jaccard_RANDOM_df %>%
  dplyr::filter(TF1 == TF2)

N_perms<- 10000
# For each TF find pvalues 
pvals_sameTF <- new_jaccard_results_filtered_sameTF %>%
  left_join(jaccard_RANDOM_sameTF, by = "TF1") %>%
  group_by(TF1) %>%
  summarise(
    J_obs = first(jaccard_index.x),
    J_mean_random = mean(jaccard_index.y),
    p_emp = (sum(jaccard_index.y >= J_obs) + 1) / (N_perms + 1),
    .groups = "drop"
  )

write.table(pvals_sameTF, 'Empirical_pvalues_10000_perms_D7_vs_all_runs_scenicplus.csv', sep = '\t')

