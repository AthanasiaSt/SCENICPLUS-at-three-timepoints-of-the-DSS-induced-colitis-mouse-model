# Code for Main Figure 6

#libraries
library(Seurat)
library("clusterProfiler")
library("enrichplot")
organism = "org.Mm.eg.db"
library(organism, character.only = TRUE)
#BiocManager::install("org.Mm.eg.db")

# Download from zenodo the processed files 
# set working directory
var<-getwd() 
setwd(paste0(var,'/Processed_datasets_for_scRNA_scATAC_scenicplus/'))

# loading the scRNA-seq object with the Harmony integration of D0, D7, D14 fibroblasts
seurat_obj<-readRDS(file = "/scRNA_seurat_integration_D0_D7_D14/Seurat5_Harmony_integration_all_samples_D0_D7_D14_fibroblasts.rds")

#regulons
regulons_df<-read.table('/SCENICPLUS_D0_D7_D14_regulons_results/Direct_Positive_eRegulons_filtered_basedOn_GeneBased_AUC.csv', row.names = 1, sep='\t', header = T)

regulons_df %>% group_by(TF) %>%
  top_n(n = -200, wt = triplet_rank) -> regulons_df_top100

signatures<-list()

# find overlaps of regulons 
signatures$Nfia<- unique(regulons_df[regulons_df$TF == 'Nfia',]$Gene)
signatures$Nfib<- unique(regulons_df[regulons_df$TF == 'Nfib',]$Gene)

#---------------functional enrichment analysis 
ora_analysis_bp <- enrichGO(gene=unique(intersect(signatures$Nfia,signatures$Nfib)) ,  #  unique(top10[top10$TF %in% signatures$Trophocytes,]$Gene)
                            universe = rownames(seurat_obj), 
                            OrgDb = org.Mm.eg.db,  # contains the TAIR/Ensembl id to GO correspondence for A. thaliana
                            keyType = "SYMBOL",
                            ont = "BP",              # either "BP", "CC" or "MF",
                            pAdjustMethod = "BH",
                            qvalueCutoff = 0.05,
                            readable = FALSE, 
                            pool = FALSE)



ora_analysis_bp_simplified <- clusterProfiler::simplify(ora_analysis_bp) 

#plots
pdf('Fig6_b_Ora_Nfia_Nfib_intersection_allGenes.pdf', width = 10, height =10)
dotplot(ora_analysis_bp_simplified, showCategory=12) + theme(text = element_text(size=20),axis.text.y = element_text(size=20),axis.text.x = element_text(size=20)) + scale_size(range = c(5, 15))  # increase min/max dot size
dev.off()

#Fig6. c
#Download from: https://gitlab.com/xavier-lab-computation/public/molecular-cartography-mouse-gut/-/blob/00a1f472546c18a274294e4e4e374b57b19a4428/visium/data/paper_dss_rolls.rds
spatial_data <- readRDS("paper_dss_rolls.rds")  # Load object containing spatial transcriptomics data
DefaultAssay(spatial_data) <- "Spatial"         # Set default assay to "Spatial"
spatial_data <- NormalizeData(spatial_data) %>% ScaleData()  # Normalize and scale the spatial expression data

# Gene-specific range calculation
get_gene_specific_range <- function(seurat_obj, feature, slot = "data") {
  # Extract expression values of the gene from the given assay slot
  all_expr <- GetAssayData(seurat_obj, slot = slot)[feature, ]
  
  # Define min, max, and mid values for color scaling
  min_val <- 0
  max_val <- max(all_expr, na.rm = TRUE)  # Maximum expression across all spots/cells
  mid_val <- median(all_expr[all_expr > 0], na.rm = TRUE)  # Median of nonzero values
  
  # Return values as a list for later use in plotting
  return(list(min = min_val, max = max_val, mid = mid_val, gene = feature))
}

# Create spatial plot for a gene in one sample
create_gene_scaled_plot <- function(seurat_obj, sample_name, feature, 
                                    slot = "data", 
                                    colors = c("grey", "white", "darkred"),
                                    point_size = 0.35, 
                                    alpha = 0.9,
                                    show_background = FALSE) {
  
  # Get min, max, and midpoint values for scaling color
  gene_range <- get_gene_specific_range(seurat_obj, feature, slot)
  
  # Extract spot/cell coordinates and scaling factors for the image
  coords <- seurat_obj@images[[sample_name]]@coordinates
  scale_factor <- seurat_obj@images[[sample_name]]@scale.factors$lowres
  
  # Match cells in the coordinate table with those in the expression matrix
  sample_cells <- intersect(rownames(coords), colnames(seurat_obj))
  
  # Get expression values of the gene for those cells
  expr_data <- GetAssayData(seurat_obj, slot = slot)[feature, sample_cells]
  
  # Build plotting dataframe: x/y coordinates + expression value
  plot_data <- data.frame(
    x = coords[sample_cells, "imagecol"] * scale_factor,
    y = coords[sample_cells, "imagerow"] * scale_factor,
    expression = expr_data
  )
  
  # Base ggplot: scatter plot of spatial spots colored by expression
  p <- ggplot(plot_data, aes(x = x, y = y, color = expression))
  
  # Optionally add background spots in gray for context
  if(show_background) {
    p <- p + 
      geom_point(color = "gray95", size = point_size * 0.8, alpha = 0.3) +
      geom_point(aes(color = expression), size = point_size, alpha = alpha)
  } else {
    p <- p + geom_point(size = point_size, alpha = alpha)
  }
  
  # Apply color scaling with custom min, midpoint, and max
  p <- p +
    scale_color_gradient2(
      low = colors[1], 
      mid = colors[2], 
      high = colors[3],
      midpoint = gene_range$mid,
      limits = c(gene_range$min, gene_range$max),
      name = paste0(feature, "\n(0-", round(gene_range$max, 1), ")"),
      oob = scales::squish  # Squish out-of-bound values into limits
    ) +
    scale_y_reverse() +   # Flip y-axis so it matches image orientation
    theme_void() +        # Remove axes and grid
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      legend.position = "right",
      legend.key.height = unit(0.8, "cm"),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    labs(title = sample_name) +  # Add sample name as title
    coord_fixed()                # Keep aspect ratio fixed
  
  return(p)
}

# Create figure for one gene across multiple samples
create_gene_figure <- function(seurat_obj, gene_name, 
                               colors = c("grey", "white", "darkred"),
                               point_size = 0.35,
                               alpha = 0.9,
                               slot = "data") {
  
  # Get sample names (each corresponds to a tissue section image)
  sample_names <- names(seurat_obj@images)
  
  # Get gene-specific min/max/mid expression values
  gene_range <- get_gene_specific_range(seurat_obj, gene_name, slot)
  
  # Print summary to console for reference
  cat("Gene:", gene_name, "| Range: 0 to", round(gene_range$max, 2), 
      "| Midpoint:", round(gene_range$mid, 2), "\n")
  
  # Generate plots for each sample
  plots <- list()
  for(i in 1:length(sample_names)) {
    plots[[i]] <- create_gene_scaled_plot(
      seurat_obj = seurat_obj,
      sample_name = sample_names[i],
      feature = gene_name,
      colors = colors,
      point_size = point_size,
      alpha = alpha,
      slot = slot,
      show_background = FALSE
    )
  }
  
  # Combine sample plots into one figure (4 per row)
  figure <- wrap_plots(plots, ncol = 4)
  return(figure)
}

# Define list of genes you want to plot
genes_to_plot <- c("Nfia", "Nfib")
# Loop over each gene and save plots
for (g in genes_to_plot) {
  message("Processing gene: ", g)
  
  # Create figure
  gene_fig <- create_gene_figure(
    seurat_obj = spatial_data,   # your Seurat object
    gene_name = g,
    colors = c("grey", "white", "darkred"),
    point_size = 0.35,
    alpha = 0.9,
    slot = "data"
  )
  
  # Save outputs
  ggsave(paste0(g, "_fig6_c_spatial.pdf"), gene_fig, width = 16, height = 4, dpi = 300)
}

