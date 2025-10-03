# Code for Figure 4 SUPP

#libraries 
library(dplyr)
library(tidyr)
library(tidyverse)
library(pheatmap)
library(viridis)


# Download processed data from Zenodo: Processed_datasets_for_scRNA_scATAC_scenicplus 

# set working directory
var<-getwd() 
setwd(paste0(var,'/Processed_datasets_for_scRNA_scATAC_scenicplus/'))


#filtered regulons, can be found in zenodo
regulons_df<-read.table('/SCENICPLUS_D0_D7_D14_regulons_results/Direct_Positive_eRegulons_filtered_basedOn_GeneBased_AUC.csv', row.names = 1, sep='\t', header = T)

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

dim(jaccard_filtered_matrix)
