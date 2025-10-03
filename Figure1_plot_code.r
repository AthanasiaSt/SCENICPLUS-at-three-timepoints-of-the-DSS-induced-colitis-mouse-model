# Code for Main Figure 1

#libraries 
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(SCpubr)
library(scater)
library(readxl)
library(RColorBrewer)
library(pals)
library(pheatmap)
library(patchwork)
library(dplyr)

# Download processed data from Zenodo: Processed_datasets_for_scRNA_scATAC_scenicplus 

# set working directory
var<-getwd() 
setwd(paste0(var,'/Processed_datasets_for_scRNA_scATAC_scenicplus/'))

# loading the scRNA-seq object with the Harmony integration of D0, D7, D14 fibroblasts
seurat_obj<-readRDS(file = "/scRNA_seurat_integration_D0_D7_D14/Seurat5_Harmony_integration_all_samples_D0_D7_D14_fibroblasts.rds")

# colours  
colours_conditions = c('D0'='#569B9A','D0_kinchen'="#A1CEC5",'D7_kinchen'='#8e0f62','D14'='#CD5808')
colours_celltypes = c('Trophocytes'="#8E7692",'PDGFRalo'='#a5af37','Telocytes'='#416522')

#-----------------------------Figure 1.b
pdf('Fig1_b_UMAP_celltypes_fibroblasts.pdf', width = 5, height = 5)
SCpubr::do_DimPlot(sample = seurat_obj, pt.size = 0.4, group.by = 'celltypes',colors.use=colours_celltypes)
dev.off()

#-----------------------------Figure 1.c
pdf('Fig1_c_UMAP_samples_fibroblasts.pdf', width = 5, height = 5)
SCpubr::do_DimPlot(sample = seurat_obj, pt.size = 0.4, group.by = 'sample', colors.use=colours_conditions)
dev.off()

#-----------------------------Figure 1.d
genes<-c("Bmp5","Sox6", "Pdgfra",'Col14a1',"Col15a1", 'Fgfr2', "Edil3","C3","Cd81", "Pi16")

pdf('Fig1_dMarkers_celltypes_dotplot.pdf', width = 7, height = 3)

SCpubr::do_DotPlot(sample = seurat_obj, 
                   features = genes, 
                   cluster = F,
                   dot.scale = 8, group.by = 'celltypes', font.size = 10)

dev.off()

#-----------------------------Figure 1.e
df<-as.data.frame(table(seurat_obj$sample,seurat_obj$celltypes)/rowSums(table(seurat_obj$sample,seurat_obj$celltypes)))

pdf('Fig1_e_Percentage_celltypes.pdf', width=2, height = 0.8)
ggplot(df,                  # Stacked barplot using ggplot2
       aes(x = Freq,
           y = Var1,
           fill = Var2)) + 
  labs( 
    y="samples", x='percentage')+
  
  geom_bar(stat = "identity")+   scale_fill_manual(values = colours_celltypes)  + theme(text = element_text(size=3),axis.text  = element_text(size=3)) 

dev.off()

#-----------------------------Figure 1.f
# Heatmaps with up/down regulated gene sets that show significant functional enrichments 

# loading common and condition-specific up/down-regulated enriched functional terms, download Supplementary table 1 in the manuscript.

# upregulated
functional_terms_up<-  read_excel('Supplementary Data 1.xlsx', sheet = 'ORA_upDEGs_top_results')
common_up <- unique(unlist(strsplit(functional_terms_up[functional_terms_up$category == 'common_up_functions',]$geneID,"/")))
d14_up <- unique(unlist(strsplit(functional_terms_up[functional_terms_up$category == 'D14_up_functions',]$geneID,"/")))
d7_up <- unique(unlist(strsplit(functional_terms_up[functional_terms_up$category == 'D7_up_functions',]$geneID,"/")))

#downregulated
d7_down<-  read_excel('Supplementary Data 1.xlsx', sheet = 'ORA_downDEGs_top_results')
d7_down <- unique(unlist(strsplit(d7_down$geneID,"/")))

# extract all the genes corresponding to the functional terms and
#calculate average log2 fold change for each gene, in each cell type and condition against the corresponding control
temp<-c()
temp1<-c()

for (cluster in c('Trophocytes','PDGFRalo','Telocytes')) {
  fold_ch_D14<-FoldChange(seurat_obj, group.by = 'celltypes_timepoint',features =unique(c(d7_up, common_up, d14_up,d7_down)) , ident.1 = paste0('D14_',cluster),ident.2 =  paste0('D0_',cluster))
  temp[[paste0('D14_',cluster)]] <- fold_ch_D14$avg_log2FC
  
  fold_ch_D7<-FoldChange(seurat_obj, group.by = 'celltypes_timepoint',features = unique(c(d7_up, common_up, d14_up,d7_down)), ident.1 =  paste0('D7_kinchen_',cluster),ident.2 = paste0('D0_kinchen_',cluster))
  temp1[[paste0('D7_kinchen_',cluster)]] <- fold_ch_D7$avg_log2FC
  
}

# put the results in a dataframe
df <- cbind(as.data.frame(temp1), as.data.frame(temp)) 
rownames(df) <-unique(c(d7_up, common_up, d14_up,d7_down))
df<-as.data.frame(df)
df <- as.matrix(df)

# Function to create breakpoints for the heatmaps color scale
# based on quantiles of the input values.
quantile_breaks <- function(xs, n = 500) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(df, n = 500)

# Make a single hetmap for each gene_set corresponding to up/down regulated functions
genes_to_plot <- list(d7_up,common_up, d14_up, d7_down)
names_in_titles <- list('D7_up_regulated_genes','Commonly_up_regulated_genes','D14_up_regulated_genes', 'D7_down_regulated_genes')

for (i in 1:4){
  
# multiple levels of annotation on the heatmap / condition X celltype
mat_col_g_c <- data.frame(group = c('D7','D7','D7','D14','D14','D14'), celltype= c('Trophocytes','PDGFRalo','Telocytes','Trophocytes','PDGFRalo','Telocytes'))
rownames(mat_col_g_c) <- colnames(df[genes_to_plot[[i]],])

mat_colors_g_c <- list(group = c('#8e0f62','#CD5808'), celltype = c('#8E7692','#a5af37','#416522'))
names(mat_colors_g_c$group) <- c('D7','D14')
names(mat_colors_g_c$celltype) <- c('Trophocytes','PDGFRalo','Telocytes')

pdf(paste0(names_in_titles[[i]],'_Fig1_f_heatmap.pdf'), width = 3, height = 4)
(
pheatmap(
  mat               = df[genes_to_plot[[i]],],
  cluster_cols      = F,
  cluster_rows      = T,
  color             = coolwarm(length(mat_breaks) - 1),
  breaks            = mat_breaks,
  border_color      = NA,
  show_colnames     = T,
  show_rownames     = T,
  annotation_col    = mat_col_g_c,
  annotation_colors = mat_colors_g_c,
  #annotation_row    = mat_colors_row_g,
  drop_levels       = TRUE,
  fontsize          = 2,
  clustering_distance_rows = "euclidean",
  clustering_method = "ward.D2",
  main              = names_in_titles[[i]]
)
)
dev.off()
}


#----Figure 1.g
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
genes_to_plot <- c("Col5a2", "Igf1", "Cxcl9", "Col18a1",'Marcks','Mmp3','Serpina3n','Slit3')
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
  ggsave(paste0(g, "_fig1_g_spatial.pdf"), gene_fig, width = 16, height = 4, dpi = 300)
}

