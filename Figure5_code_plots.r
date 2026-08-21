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

# loading the scRNA-seq object with the Harmony integration of D0, D7, Regeneration fibroblasts
seurat_obj<-readRDS(file = "./scRNA_seurat_integration_Healthy_Inflammation_Regeneration/Seurat5_Harmony_integration_all_samples_Healthy_Inflammation_Regeneration_fibroblasts.rds")

# colours  
colours_conditions = c('Healthy'='#569B9A','HealthyK'="#A1CEC5",'InflammationK'='#8e0f62','Regeneration'='#CD5808')
colours_celltypes = c('Trophocytes'="#8E7692",'CD81-stroma'='#a5af37','SEMFs'='#416522')


#-----------------------------Figure 5.b
#keep only D7, Regeneration
seurat_obj<- seurat_obj[,seurat_obj$sample %in% c('InflammationK','Regeneration')]

# convert seurat object to SingleCellExperiment
df.sce <- as.SingleCellExperiment(seurat_obj)
df.sce$celltypes_condition <- factor(df.sce$celltypes_condition, levels =c('InflammationK_Trophocytes','InflammationK_CD81-stroma','InflammationK_SEMFs','Regeneration_Trophocytes','Regeneration_CD81-stroma','Regeneration_SEMFs'))


#-------------------plot the mean expression across celltypes 
pdf('fig5_b_TF_expression.pdf',width =3, height = 4)
plotGroupedHeatmap(df.sce,exprs_values='logcounts', features=c('Nfkb1','Bach1','Smad3','Ets2','Ets1','Nfe2l2','Fosl1','Cebpb','Fosl2','Nfia','Nfib','Egr3','Tcf7l2','Pbx1'),center=T, scale = F,
                   group=c('celltypes_condition'), zlim = c(-0.5,0.5), show_rownames=T,display_numbers = F, cluster_cols=F,cluster_rows=T, cutree_rows = 1)# + theme(text = element_text(size=20),axis.text  = element_text(size=40))
#keep the above as a variable to extract the clusters 
dev.off()

#-----------------------------Figure 5.c
# circos plot with TFs and target genes based on scenicplus 

# Load regulon data, find in zenodo
regulons_df<-read.table('./SCENICPLUS_Healthy_Inflammation_Regeneration_regulons_results/Direct_Positive_eRegulons_filtered_basedOn_GeneBased_AUC.csv', row.names = 1, sep='\t', header = T)


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
proj_fib<-loadArchRProject('./Final_ArchR_object')
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
regulons_df<-read.table('./SCENICPLUS_Healthy_Inflammation_Regeneration_regulons_results/Direct_Positive_eRegulons_filtered_basedOn_GeneBased_AUC.csv', row.names = 1, sep='\t', header = T)

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
                           features = signatures
                           )


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

#-----------------------------Figure 5.g
library(sceasy)
library(reticulate)

sceasy::convertFormat("fibperi_final.h5ad", from = "anndata", to = "seurat",
                      outFile = "fibperi_final.rds")

#import the human dataset 
seurat_obj <- readRDS('fibperi_final.rds')

#exlude samples 
seurat_obj$Remission_status <- as.character(seurat_obj$Remission_status)
seurat_obj$Remission_status[seurat_obj$Remission_status == 'None '] <- 'Healthy'
seurat_obj<-seurat_obj[,seurat_obj$Remission_status != 'Not_avail']
seurat_obj<-seurat_obj[,seurat_obj$Disease %in% c('Healthy','UC')]

#normalize 
seurat_obj <- NormalizeData(seurat_obj)

# compute FB subtype signatures in human clusters to match our celltype annotation based on mouse genes
human <- useEnsembl("ensembl","hsapiens_gene_ensembl", mirror = "useast")

signatures<- list()
signatures$Trophocytes<-unique(getBM(c("hgnc_symbol"), "external_gene_name",values =  c('Grem1','Pi16','Dpt','C3','Pcolce2') , human)$hgnc_symbol)
signatures$SEMFs<- unique(getBM(c("hgnc_symbol"), "external_gene_name",values = c('Pdgfra','Sox6','Bmp5','Bmp7','Wnt5a','F3'), human)$hgnc_symbol)
signatures$CD81_stroma<- unique(getBM(c("hgnc_symbol"), "external_gene_name",values = c('Adamdec1','Tcf21','Fgfr2','Col15a1','Edil3'), human)$hgnc_symbol)

sample.combined <- AddModuleScore_UCell(seurat_obj,  features = signatures)
signature.names <- paste0(names(signatures), "_UCell")

#average UCell per category
df_all <- sample.combined@meta.data[, c('Trophocytes_UCell','CD81_stroma_UCell','SEMFs_UCell')]
df_all$celltypes<-as.character(sample.combined$final_analysis)

df_means <- df_all %>%
  dplyr::group_by(celltypes) %>%
  dplyr::summarise(dplyr::across(dplyr::all_of( c('Trophocytes_UCell','CD81_stroma_UCell','SEMFs_UCell')), ~ mean(.x, na.rm = TRUE), .names = "{.col}_mean"))

df_means<- as.data.frame(df_means)
rownames(df_means) <- df_means$celltypes
df_means$celltypes <- NULL

#Z-score normalization by column (signature)
df_scaled <- scale(df_means)

max_val <- max(abs(df_scaled))  # find the largest absolute z-score
breaks <- seq(-max_val, max_val, length.out = 101)  # 101 breaks for 100 colors

#Plot heatmap
pdf('fig5_g_thomas_et_al_Signatures_MouseSubtypes_categories_detailed.pdf', width = 5, height = 5)
pheatmap(t(df_scaled),
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = breaks,   # <-- this forces white = 0, symmetric scale
         cluster_rows = F,
         cluster_cols = T,
         display_numbers = F,
         fontsize_row = 10,
         fontsize_col = 10,
         border_color = NA)
dev.off()

#inlcude FB subtypes of interest 
seurat_obj<-seurat_obj[,seurat_obj$final_analysis %in% c('ABCA8pos WNT2Bpos FOSlo fibroblast','ABCA8pos WNT2Bpos FOShi fibroblast','C3hi RSPO3pos fibroblast','C3hi CCL19pos fibroblast','SOX6pos POSTNpos NRG1hi NPYpos fibroblast','SOX6pos POSTNpos fibroblast','THY1pos FAPpos PDPNpos fibroblast')]
seurat_obj$Remission_status <- factor(seurat_obj$Remission_status, levels =c('Healthy','Remission','Non_Remission'))
seurat_obj$Treatment <- factor(seurat_obj$Treatment, levels = c('Pre','Post'))
seurat_obj <- seurat_obj[, seurat_obj$Remission_status %in% c('Healthy','Remission','Non_Remission')]

# ---- Detailed fibroblast subtypes ----
celltype_order <- c(
  'C3hi RSPO3pos fibroblast',
  'C3hi CCL19pos fibroblast',
  'ABCA8pos WNT2Bpos FOShi fibroblast',
  'ABCA8pos WNT2Bpos FOSlo fibroblast',
  'SOX6pos POSTNpos NRG1hi NPYpos fibroblast',
  'SOX6pos POSTNpos fibroblast',
  'THY1pos FAPpos PDPNpos fibroblast'
)
seurat_obj$final_analysis <- factor(seurat_obj$final_analysis, levels = celltype_order)
seurat_obj <- seurat_obj[, !is.na(seurat_obj$final_analysis)]

# ---- Simplified status label ----
seurat_obj$status_simple <- ifelse(
  seurat_obj$Remission_status == 'Healthy',
  'Healthy',
  paste(seurat_obj$Remission_status, seurat_obj$Treatment, sep = "_")
)
status_order <- c('Healthy','Non_Remission_Pre','Non_Remission_Post','Remission_Pre','Remission_Post')
seurat_obj$status_simple <- factor(seurat_obj$status_simple, levels = status_order)

colours_status <- c(
  'Healthy'            = '#569B9A',
  'Non_Remission_Pre'  = '#7B2D8B',
  'Non_Remission_Post' = '#B22222',
  'Remission_Pre'      = 'darkorange4',
  'Remission_Post'     = '#CD5808'
)
# ---- Gene lists ----
tfs <- c('NFKB1','BACH1','SMAD3','ETS2','ETS1','NFE2L2','FOSL1','CEBPB','FOSL2','NFIA','NFIB','EGR3','TCF7L2','PBX1')
targets <- c('LCN2','IL6','CXCL9','CXCL13','MMP3','MMP10','TIMP1','SERPINA3','COL1A1','COL5A1','COL5A2','COL18A1','FN1','MDK','IGF1','EREG','ROBO2','ADAMTSL3')
genes_use <- c(tfs, targets)

# ---- Colors (status + celltype only) ----
colours_status <- c(
  'Non_Remission_Pre'  = '#7B2D8B',
  'Non_Remission_Post' = '#B22222',
  'Healthy'            = '#569B9A',
  'Remission_Pre'      = 'darkorange4',
  'Remission_Post'     = '#CD5808'
)

colours_celltype <- c(
  'C3hi RSPO3pos fibroblast'                  = '#8E7692',
  'C3hi CCL19pos fibroblast'                  = '#5E3A66',
  'ABCA8pos WNT2Bpos FOShi fibroblast'        = '#a5af37',
  'ABCA8pos WNT2Bpos FOSlo fibroblast'        = '#6E7A1E',
  'SOX6pos POSTNpos NRG1hi NPYpos fibroblast' = '#416522',
  'SOX6pos POSTNpos fibroblast'               = '#7FA34A',
  'THY1pos FAPpos PDPNpos fibroblast'         = '#4472C4'
)

my_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(256)
my_breaks <- seq(-2, 2, length.out = 257)

# ---- Helper: normalize strings so matching survives Seurat's sanitization ----
normalize <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))

# ---- Helper: build reordered avg_exp + matching annotation, robust to sanitization ----
build_annotated_matrix <- function(seurat_obj, genes_use, group_by_col, grid) {
  avg_exp <- AverageExpression(
    object   = seurat_obj,
    features = genes_use,
    group.by = group_by_col,
    assays   = "RNA",
    slot     = "data"
  )$RNA
  
  raw_names  <- colnames(avg_exp)
  norm_raw   <- normalize(raw_names)
  norm_combo <- normalize(grid$combo)
  
  # match each grid combo to the actual (possibly sanitized) avg_exp column
  match_idx <- match(norm_combo, norm_raw)
  if (any(is.na(match_idx))) {
    stop("Could not match these combos to avg_exp colnames: ",
         paste(grid$combo[is.na(match_idx)], collapse = ", "),
         "\nActual colnames(avg_exp): ", paste(raw_names, collapse = " | "))
  }
  
  avg_exp <- avg_exp[, match_idx]
  colnames(avg_exp) <- grid$combo   # relabel with clean names for the plot
  
  mat_col <- data.frame(celltype = grid$celltype, group = grid$status)
  rownames(mat_col) <- grid$combo
  
  list(avg_exp = avg_exp, mat_col = mat_col)
}

# Ordering 1: celltype-first, status-inner
# ============================================================
grid_ct_first <- expand.grid(celltype = celltype_order, status = status_order, stringsAsFactors = FALSE)
grid_ct_first <- grid_ct_first[order(match(grid_ct_first$celltype, celltype_order),
                                     match(grid_ct_first$status, status_order)), ]
grid_ct_first$combo <- paste(grid_ct_first$celltype, grid_ct_first$status, sep = "_")

seurat_obj$group_combined_celltype_first <- paste(seurat_obj$final_analysis, seurat_obj$status_simple, sep = "_")
grid_ct_first <- grid_ct_first[grid_ct_first$combo %in% unique(seurat_obj$group_combined_celltype_first), ]

res_ct <- build_annotated_matrix(seurat_obj, genes_use, "group_combined_celltype_first", grid_ct_first)

mat_colors_ct <- list(group = colours_status, celltype = colours_celltype)

pdf('Fig5_g_heatmap_celltype_first.pdf', width = 9, height = 7)
pheatmap(
  mat                      = res_ct$avg_exp,
  cluster_cols             = F,
  cluster_rows             = T,
  scale                    = "row",
  color                    = my_colors,
  breaks                   = my_breaks,
  border_color             = NA,
  show_colnames            = T,
  show_rownames            = T,
  annotation_col           = res_ct$mat_col,
  annotation_colors        = mat_colors_ct,
  drop_levels              = TRUE,
  fontsize                 = 6,
  clustering_distance_rows = "euclidean",
  clustering_method        = "ward.D2",
  main                     = "TFs & targets: celltype-first ordering"
)
dev.off()



