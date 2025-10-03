# Code for Figure5 in the manuscript
#libraries
library(scater)
library(UCell)
library(Seurat)
library(dplyr)
library(circlize)
library(reshape2)
library(ArchR)

# Download processed data from Zenodo: Processed_datasets_for_scRNA_scATAC_scenicplus 

# set working directory
var<-getwd() 
setwd(paste0(var,'/Processed_datasets_for_scRNA_scATAC_scenicplus/'))

# loading the scRNA-seq object with the Harmony integration of D0, D7, D14 fibroblasts
seurat_obj<-readRDS(file = "/scRNA_seurat_integration_D0_D7_D14/Seurat5_Harmony_integration_all_samples_D0_D7_D14_fibroblasts.rds")

# colours  
colours_conditions = c('D0'='#569B9A','D0_kinchen'="#A1CEC5",'D7_kinchen'='#8e0f62','D14'='#CD5808')
colours_celltypes = c('Trophocytes'="#8E7692",'PDGFRalo'='#a5af37','Telocytes'='#416522')


#-----------------------------Figure 5.b
#keep only D7, D14
seurat_obj<- seurat_obj[,seurat_obj$sample %in% c('D7_kinchen','D14')]

# convert seurat object to SingleCellExperiment
df.sce <- as.SingleCellExperiment(seurat_obj)
df.sce$celltypes_timepoint <- factor(df.sce$celltypes_timepoint, levels =c('D7_kinchen_Trophocytes','D7_kinchen_PDGFRalo','D7_kinchen_Telocytes','D14_Trophocytes','D14_PDGFRalo','D14_Telocytes'))


#-------------------plot the mean expression across celltypes 
pdf('fig5_b_TF_expression.pdf',width = 4, height = 4)
plotGroupedHeatmap(df.sce,exprs_values='logcounts', features=c('Nfkb1','Bach1','Smad3','Ets2','Ets1','Nfe2l2','Fosl1','Cebpb','Fosl2','Nfia','Nfib','Egr3','Tcf7l2','Pbx1'),center=T, scale = F,
                   group=c('celltypes_timepoint'), zlim = c(-0.5,0.5), show_rownames=T,display_numbers = F, cluster_cols=F,cluster_rows=T, cutree_rows = 1)# + theme(text = element_text(size=20),axis.text  = element_text(size=40))
#keep the above as a variable to extract the clusters 
dev.off()

#-----------------------------Figure 5.c
# circos plot with TFs and target genes based on scenicplus 

# Load regulon data, find in zenodo
regulons_df<-read.table('/SCENICPLUS_D0_D7_D14_regulons_results/Direct_Positive_eRegulons_filtered_basedOn_GeneBased_AUC.csv', row.names = 1, sep='\t', header = T)


# Define genes and TFs
gene_targets_inf <- c(
  'Mmp3','Ptgs2','Ctsk','Timp1','Ccl2','Saa3','Mmp2','Serpina3n','Serpina3m',
  'Fgf7','Angptl4','Ctsb','Il1r1','Chl1',
  'Lox','Cxcl13','Ccl7','Il6','Mt1','Mt2','Mmp14','Lcn2','Ereg')

tfs_inf <-  rev(c('Cebpb','Nfkb1','Nfe2l2','Fosl1','Ets2','Ets1','Smad3'))

gene_targets_reg <- c(
  'Robo2','Col5a1','Mdk','Marcks','Ror1','Sema3c','Cacna2d1','Adamts17','Adamts19',
  'Col1a1','Ptk2','Fn1','Igf1','Col5a2',
  'Gli2','Ptn','Fbn1','Adamtsl1','Adamtsl3','Slit2','Slit3','Thsd4')

tfs_reg <-  rev(c('Nfia','Nfib','Pbx1','Tcf7l2'))


#set colour groups
gene_colors_inf <- rep("", length(gene_targets_inf))
names(gene_colors_inf) <- gene_targets_inf
gene_colors_inf[names(gene_colors_inf) %in% gene_targets_inf] <- "magenta4"

gene_colors_reg <- rep("", length(gene_targets_reg))
names(gene_colors_reg) <- gene_targets_reg
gene_colors_reg[names(gene_colors_reg) %in% gene_targets_reg] <- "orange"

# for TFs as well
tf_colors_inf <- rep("", length(tfs_inf))
names(tf_colors_inf) <- tfs_inf
tf_colors_inf[names(tf_colors_inf) == 'Cebpb'] <- "darkslateblue"
tf_colors_inf[names(tf_colors_inf) == 'Nfkb1'] <- "darkviolet"
tf_colors_inf[names(tf_colors_inf) == 'Nfe2l2'] <- "darkorchid4"
tf_colors_inf[names(tf_colors_inf) %in% c('Ets1','Ets2','Fosl1','Smad3')] <- "gray"



tf_colors_reg <- rep("", length(tfs_reg))
names(tf_colors_reg) <- tfs_reg
tf_colors_reg[names(tf_colors_reg) == 'Nfia'] <- "darkorange4"
tf_colors_reg[names(tf_colors_reg) == 'Nfib'] <- "darkorange3"
tf_colors_reg[names(tf_colors_reg) == 'Pbx1'] <- "coral"
tf_colors_reg[names(tf_colors_reg) == 'Tcf7l2'] <- "gray"


#---------------------------------------------------#
# Function to plot TF-gene circos diagram
#---------------------------------------------------#
plot_tf_gene_circos <- function(regulons_df, tfs, gene_targets, tf_colors, gene_colors, 
                                file_name = "circos_plot.pdf", title_text = "TF-gene regulation") {
  
  # Define sector colors
  all_colors <- c(tf_colors, gene_colors)
  
  # Define sector order
  sector_order <- c(tfs, gene_targets)
  
  # Subset regulon data
  subset_data <- subset(regulons_df, TF %in% tfs & Gene %in% gene_targets)
  subset_data$signed_importance <- subset_data$importance_TF2G
  
  # Aggregate if TF–gene pair has multiple enhancers
  aggregated <- aggregate(signed_importance ~ Gene + TF, data = subset_data, FUN = mean)
  
  # Convert to circos long format
  links <- aggregated
  colnames(links) <- c("to", "from", "value")  # circos expects from→to
  links <- links[links$value > 0, ]            # keep positive only
  
  # Scale link widths
  max_val <- max(links$value)
  links$value <- links$value / max_val * 5     # max width = 5
  
  # Assign link colors = TF color
  links$col <- tf_colors[links$from]  
  
  #----------------------------- Plot
  pdf(file = file_name, width = 15, height = 15)
  
  circos.clear()
  circos.par(gap.after = c(rep(2, length(tfs) - 1), 10,
                           rep(1, length(gene_targets) - 1), 10))
  
  chordDiagram(
    links,
    order = sector_order,
    grid.col = all_colors,
    col = links$col,
    directional = 1,
    direction.type = "arrows",
    link.arr.type = "big.arrow",
    annotationTrack = "grid",
    link.lwd = links$value,        
    preAllocateTracks = list(track.height = 0.1)
  )
  
  circos.trackPlotRegion(
    track.index = 1,
    panel.fun = function(x, y) {
      sector.name <- get.cell.meta.data("sector.index")
      xlim <- get.cell.meta.data("xlim")
      ylim <- get.cell.meta.data("ylim")
      circos.text(mean(xlim), ylim[1] + 1, sector.name,
                  facing = "clockwise", niceFacing = TRUE,
                  adj = c(0, 0.5), cex = 1.2)
    },
    bg.border = NA
  )
  
  title(title_text)
  
  dev.off()
  circos.clear()
}


# For inflammation set
plot_tf_gene_circos(
  regulons_df, 
  tfs = tfs_inf, 
  gene_targets = gene_targets_inf, 
  tf_colors = tf_colors_inf, 
  gene_colors = gene_colors_inf, 
  file_name = "Fig5c_Circos_Inflammation.pdf",
  title_text = "TF-gene regulation (Inflammation)"
)

# For regeneration set
plot_tf_gene_circos(
  regulons_df, 
  tfs = tfs_reg, 
  gene_targets = gene_targets_reg, 
  tf_colors = tf_colors_reg, 
  gene_colors = gene_colors_reg, 
  file_name = "Fig5c_Circos_Regeneration.pdf",
  title_text = "TF-gene regulation (Regeneration)"
)

#-----------------------------Figure 5.d
# start with ArchR commads
addArchRGenome("mm10")
addArchRThreads(threads = 48) 

#load ArchR project. Find in zenodo
proj_fib<-loadArchRProject('Final_ArchR_object')
proj_fib <- addImputeWeights(proj_fib)

#-- ChromVAr--- deviation scores 
if("Motif" %ni% names(proj_fib@peakAnnotation)){
  proj_fib <- addMotifAnnotations(ArchRProj = proj_fib, motifSet = "cisbp", name = "Motif")
}

#proj_fib <- addBgdPeaks(proj_fib)

#proj_fib <- addDeviationsMatrix(
#  ArchRProj = proj_fib, 
#  peakAnnotation = "Motif",
#  force = TRUE
#)

markerMotifs <- getFeatures(
  proj_fib, 
  select = paste(c('Nfkb1','Cebpb','Fosl1','Nfia','Nfib','Pbx1'), collapse="|"), 
  useMatrix = "MotifMatrix"
)

#Keep the z-scores also know as deviation scores. 
markerMotifs <- grep("z:", markerMotifs, value = TRUE)

proj_fib <- addImputeWeights(proj_fib)

p <- plotGroups(ArchRProj = proj_fib, 
                groupBy = "celltype_timepoint", 
                colorBy = "MotifMatrix", 
                name = markerMotifs,
                imputeWeights = getImputeWeights(proj_fib)
)

pdf('fig5_d_deviations_violin_plots.pdf',width = 5, height = 5)
p
dev.off()

#-----------------------------Figure 5.e
all_atac_genes<-getFeatures(proj_fib)

# load regulons
regulons_df<-read.table('/SCENICPLUS_D0_D7_D14_regulons_results/Direct_Positive_eRegulons_filtered_basedOn_GeneBased_AUC.csv', row.names = 1, sep='\t', header = T)

regulons_df %>% group_by(TF) %>%
  top_n(n = -200, wt = triplet_rank) -> regulons_df_top100

signatures<-list()
# characterize each state
signatures$inflammation<- unique(regulons_df_top100[regulons_df_top100$TF %in% c('Nfkb1','Bach1','Smad3','Ets2','Ets1','Nfe2l2','Fosl1','Cebpb'),]$Gene)
signatures$regeneration<- unique(regulons_df_top100[regulons_df_top100$TF %in% c('Fosl2','Nfia','Nfib','Egr3','Tcf7l2','Pbx1'),]$Gene)

signatures$inflammation<- signatures$inflammation[signatures$inflammation %in% all_atac_genes]
signatures$regeneration<- signatures$regeneration[signatures$regeneration %in% all_atac_genes]

# -------- module score --
proj_fib <- addImputeWeights(proj_fib)


proj_fib <- addModuleScore(proj_fib,
                           useMatrix = "GeneScoreMatrix",
                           name = "Module",
                           features = signatures,
                           imputeWeights = getImputeWeights(proj_fib))


# we change the order , so that order is changed in the plot as well. 
proj_fib$new_celltype_timepoint <- paste0(proj_fib$celltype,'_',proj_fib$Sample)


p <- plotGroups(
  ArchRProj = proj_fib, 
  groupBy = "new_celltype_timepoint", 
  colorBy = "cellColData", 
  name = "Module.inflammation",
  plotAs = "violin",
  alpha = 0.4,
  baseSize = 10,
  addBoxPlot = TRUE,
)

p2 <- plotGroups(
  ArchRProj = proj_fib, 
  groupBy = "new_celltype_timepoint", 
  colorBy = "cellColData", 
  name = "Module.regeneration",
  plotAs = "violin",
  alpha = 0.4,
  baseSize = 10,
  addBoxPlot = TRUE,
)


pdf('Fig5_e_top200_activation_signatures_genescorematrix.pdf', width = 4, height = 4)
p
p2
dev.off()

#-----------------------------Figure 5.f
# umaps in scATAC // 
genes<- c('Nfkb1', 'Ccl2', 'Fosl1', 'Cebpb','Saa3','Serpina3n','Nfia', 'Nfib', 'Pbx1','Marcks','Sema3c', 'Slit3')

proj_fib <- addImputeWeights(proj_fib)

plots <- plotEmbedding(
  ArchRProj = proj_fib, 
  colorBy = "GeneScoreMatrix", 
  name = unique(genes), 
  embedding = "UMAP_unintegrated",
  quantCut = c(0.01, 0.95),
  imputeWeights = getImputeWeights(proj_fib)
)

# Apply theme to every ggplot in the list
plots_clean <- lapply(plots, function(p){
  p + 
    theme_void() + 
    theme(legend.position = "right")
})

# Save them all to a multi-page PDF
pdf("Fig5_f_umaps_scATAC_no_axes.pdf", width = 4, height = 4)
for (p in plots_clean) {
  print(p)
}
dev.off()
