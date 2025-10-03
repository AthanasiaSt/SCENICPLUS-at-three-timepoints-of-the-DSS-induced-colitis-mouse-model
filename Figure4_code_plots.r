#code for Main Figure 4 

library(monocle3)
library(Seurat)

# Download processed data from Zenodo: Processed_datasets_for_scRNA_scATAC_scenicplus 

# set working directory
var<-getwd() 
setwd(paste0(var,'/Processed_datasets_for_scRNA_scATAC_scenicplus/'))

# loading the scRNA-seq object with the Harmony integration of D0, D7, D14 fibroblasts
seurat_obj<-readRDS(file = "/scRNA_seurat_integration_D0_D7_D14/Seurat5_Harmony_integration_all_samples_D0_D7_D14_fibroblasts.rds")

# colours  
colours_conditions = c('D0'='#569B9A','D0_kinchen'="#A1CEC5",'D7_kinchen'='#8e0f62','D14'='#CD5808')
colours_celltypes = c('Trophocytes'="#8E7692",'PDGFRalo'='#a5af37','Telocytes'='#416522')

#------------Figure 4.a 
# insert raw data into a monocle object  
counts<- seurat_obj@assays$RNA@layers$counts
expression_matrix <- as(as.matrix(counts), "sparseMatrix")     

metadata<- seurat_obj[[]]

colnames(expression_matrix) <- row.names(metadata)
rownames(expression_matrix) <- rownames(seurat_obj)

gene_annotation <-  as.data.frame(rownames(seurat_obj))
rownames(gene_annotation) <- rownames(seurat_obj)
gene_annotation$gene_short_name <- rownames(seurat_obj)

# Make the CDS object
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = metadata,
                         gene_metadata = gene_annotation)



#monocle preprocessing
monocle_object <- preprocess_cds(cds, num_dim = 40)
monocle_object <- align_cds(monocle_object, alignment_group = "dataset", residual_model_formula_str = '~ percent.mt')
monocle_object <- reduce_dimension(monocle_object)
monocle_object <- cluster_cells(cds = monocle_object)
monocle_object <- learn_graph(monocle_object, use_partition = TRUE)

# a helper function to identify the root principal points: from D0
get_earliest_principal_node <- function(monocle_object, time_bin="D0"){
  cell_ids <- which(colData(monocle_object)[, "sample"] == time_bin)
  
  closest_vertex <-
    monocle_object@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(monocle_object), ])
  root_pr_nodes <-
    igraph::V(principal_graph(monocle_object)[["UMAP"]])$name[as.numeric(names
                                                                         (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
monocle_object <- order_cells(monocle_object, root_pr_nodes=get_earliest_principal_node(monocle_object))

#this is needed for plotting genes on the umaps
rowData( monocle_object ) <- data.frame( 'gene_short_name'=rownames(monocle_object) )

#------------ UMAPs  
pdf('Fig4_a_Monocle_umaps.pdf', width = 6, height = 6)

plot_cells(monocle_object,
           color_cells_by = "pseudotime", cell_size = 1,
           graph_label_size=3,
           label_cell_groups=F,
           show_trajectory_graph = TRUE)

plot_cells(monocle_object,
           color_cells_by = "celltypes", cell_size = 1,
           graph_label_size=3,
           label_cell_groups=F,
           show_trajectory_graph = TRUE) + scale_color_manual(values = colours_celltypes)

plot_cells(monocle_object,
           color_cells_by = "sample", cell_size = 1,
           graph_label_size=3,
           label_cell_groups=F,
           show_trajectory_graph = TRUE)+ scale_color_manual(values = colours_conditions)


dev.off()

#-------------------------------------- Figure 4.j

pdf('Fig4_j_Pi16_grem1_D0_D0kinchen_blended_featureplot.pdf', width = 15, height = 4)
FeaturePlot(
  seurat_obj[,seurat_obj$sample %in% c('D0','D0_kinchen')],
  features = c("Pi16", "Grem1"),
  blend = TRUE, order = TRUE, 
  blend.threshold = 0.3, reduction = 'umap',
  pt.size = 0.2, combine = T, cols = c('lightgrey','darkolivegreen','darkmagenta')
)
dev.off()
