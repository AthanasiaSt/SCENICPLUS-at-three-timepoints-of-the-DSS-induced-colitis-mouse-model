# Code for main Figure3

#libraries
library(scater)
library(UCell)
library(Seurat)
library(dplyr)
library(circlize)
library(reshape2)
library(RColorBrewer)
library(ArchR)
library(biomaRt)
library(pheatmap)

# Download processed data from Zenodo: Processed_datasets_for_scRNA_scATAC_scenicplus 

# set working directory
var<-getwd() 
setwd(paste0(var,'/Processed_datasets_for_scRNA_scATAC_scenicplus/'))

# loading the scRNA-seq object with the Harmony integration of D0, D7, D14 fibroblasts
seurat_obj<-readRDS(file = "/scRNA_seurat_integration_D0_D7_D14/Seurat5_Harmony_integration_all_samples_D0_D7_D14_fibroblasts.rds")

# colours  
colours_conditions = c('D0'='#569B9A','D0_kinchen'="#A1CEC5",'D7_kinchen'='#8e0f62','D14'='#CD5808')
colours_celltypes = c('Trophocytes'="#8E7692",'PDGFRalo'='#a5af37','Telocytes'='#416522')


#-----------------------------Figure 3.b
#keep only D0 
seurat_obj<- seurat_obj[,seurat_obj$sample == 'D0']

# convert seurat object to SingleCellExperiment as input to scater
df.sce <- as.SingleCellExperiment(seurat_obj)

order_genes<-c('Ar',"Ebf1" ,"Klf4",'Klf2',"Pbx1",'Egr1','Maf','Tcf21', 'Pitx1', 'Tcf4',"Etv1","Foxf2","Foxf1","Runx1","Runx2")
#-------------------plot the mean expression across celltypes 

pdf('Fig3_b_D0_Tfs_gene_expression.pdf',width = 2, height = 4)

plotGroupedHeatmap(df.sce,exprs_values='logcounts', features=order_genes,center=T, scale = F,
                   group=c('celltypes'), zlim = c(-0.5,0.5), show_rownames=T,display_numbers = F, cluster_cols=F,cluster_rows=F, cutree_rows = 1)# + theme(text = element_text(size=20),axis.text  = element_text(size=40))
#keep the above as a variable to extract the clusters 
dev.off()

#-----------------------------Figure 3.c

#--------------------------- UCell Signatures 
#reload dataset
seurat_obj<-readRDS(file = "/scRNA_seurat_integration_D0_D7_D14/Seurat5_Harmony_integration_all_samples_D0_D7_D14_fibroblasts.rds")

#filtered regulons, can be found in zenodo
regulons_df<-read.table('/SCENICPLUS_D0_D7_D14_regulons_results/Direct_Positive_eRegulons_filtered_basedOn_GeneBased_AUC.csv', row.names = 1, sep='\t', header = T)

#keep the top 100 genes based on triplet ranking per regulon
regulons_df %>% group_by(TF) %>%
  top_n(n = -100, wt = triplet_rank) -> regulons_df_top100

signatures<-list()

# make the signatures per FB subtype, based on the union of the top 100 genes per regulon
signatures$Telocytes<- unique(regulons_df_top100[regulons_df_top100$TF %in% c('Runx1','Runx2','Foxf1','Foxf2','Tcf4','Etv1'),]$Gene)
signatures$PDGFRalo<- unique(regulons_df_top100[regulons_df_top100$TF %in% c('Tcf21'),]$Gene)
signatures$Trophocytes<- unique(regulons_df_top100[regulons_df_top100$TF %in% c('Klf4','Klf2','Ebf1','Egr1','Ar','Pbx1'),]$Gene)

seurat_obj <- AddModuleScore_UCell(seurat_obj,  features = signatures)

pdf('Fig3_c_Homeostasis_TFs_top100_triplet_signatures.pdf', width = 3, height = 3)
VlnPlot(seurat_obj, features = 'Trophocytes_UCell', group.by = 'celltypes', cols =colours_celltypes,  ncol = 1, pt.size = 0, combine = T)
VlnPlot(seurat_obj, features = 'PDGFRalo_UCell', group.by = 'celltypes', cols =colours_celltypes,  ncol = 1, pt.size = 0, combine = T)
VlnPlot(seurat_obj, features = 'Telocytes_UCell', group.by = 'celltypes', cols =colours_celltypes,  ncol = 1, pt.size = 0, combine = T)
dev.off()


pdf('Fig3_c_Homeostasis_split_dataset_TFs_top100_triplet_.pdf', width = 3, height = 3)
VlnPlot(seurat_obj[,seurat_obj$sample %in% c('D0_kinchen','D7_kinchen')], features = 'Trophocytes_UCell',split.by = 'sample', group.by = 'celltypes', cols  = c('D0_kinchen'="#A1CEC5",'D7_kinchen'='#8e0f62'), ncol = 1, pt.size = 0, combine = T)
VlnPlot(seurat_obj[,seurat_obj$sample %in% c('D0','D14')], features = 'Trophocytes_UCell',split.by = 'sample', group.by = 'celltypes', cols  = c('D0'='#569B9A','D14'='#CD5808'), ncol = 1, pt.size = 0, combine = T)

VlnPlot(seurat_obj[,seurat_obj$sample %in% c('D0_kinchen','D7_kinchen')], features = 'PDGFRalo_UCell',split.by = 'sample', group.by = 'celltypes', cols  = c('D0_kinchen'="#A1CEC5",'D7_kinchen'='#8e0f62'), ncol = 1, pt.size = 0, combine = T)
VlnPlot(seurat_obj[,seurat_obj$sample %in% c('D0','D14')], features = 'PDGFRalo_UCell',split.by = 'sample', group.by = 'celltypes', cols  = c('D0'='#569B9A','D14'='#CD5808'), ncol = 1, pt.size = 0, combine = T)

VlnPlot(seurat_obj[,seurat_obj$sample %in% c('D0_kinchen','D7_kinchen')], features = 'Telocytes_UCell',split.by = 'sample', group.by = 'celltypes', cols  = c('D0_kinchen'="#A1CEC5",'D7_kinchen'='#8e0f62'), ncol = 1, pt.size = 0, combine = T)
VlnPlot(seurat_obj[,seurat_obj$sample %in% c('D0','D14')], features = 'Telocytes_UCell',split.by = 'sample', group.by = 'celltypes', cols  = c('D0'='#569B9A','D14'='#CD5808'), ncol = 1, pt.size = 0, combine = T)

dev.off()


#-----------------------------Figure 3.d
# circos plot with TFs and target genes based on scenicplus results

# Load regulon data -- from zenodo
regulons_df <- read.table(
  '/SCENICPLUS_D0_D7_D14_regulons_results/Direct_Positive_eRegulons_filtered_basedOn_GeneBased_AUC.csv',
  row.names = 1, sep = '\t', header = TRUE
)

# Define genes and TFs
gene_targets <- rev(c(
  'Mylk','Tagln','Myl9','Bmp3','Bmp5','Bmp7','Wnt5a','Sox6','Ednrb',
  'Col15a1','Edil3','Igfbp3','Fbln1','Adam23',
  'Atf3','Cd81','Prss23','Pcolce2','Col14a1','Anxa3','Pi16','Igfbp6'
))

tfs <- c('Tcf4','Foxf2','Runx1','Runx2','Foxf1',
         'Tcf21','Klf4','Egr1','Ar','Ebf1')

#set colour groups
group1 <- c('Mylk','Tagln','Myl9','Bmp3','Bmp5','Bmp7','Wnt5a','Sox6','Ednrb')
group2 <- c('Col15a1','Edil3','Igfbp3','Fbln1','Adam23')
group3 <- c('Atf3','Cd81','Prss23','Pcolce2','Col14a1','Anxa3','Pi16','Igfbp6')

gene_colors <- c(
  setNames(rep("#416522", length(group1)), group1),
  setNames(rep("#a5af37", length(group2)), group2),
  setNames(rep("#8E7692", length(group3)), group3)
)
# for TFs as well
group1 <-  c('Tcf4','Foxf2','Runx1','Runx2','Foxf1')
group2 <- c('Tcf21')
group3 <- c('Klf4','Egr1','Ar','Ebf1')

tf_colors <- c(
  setNames(rep("#416522", length(group1)), group1),
  setNames(rep("#a5af37", length(group2)), group2),
  setNames(rep("#8E7692", length(group3)), group3)
)

# Sector colors
all_colors <- c(tf_colors, gene_colors)

# Order
sector_order <- c(tfs, gene_targets)

# Subset and aggregate in case TF is linked with a gene through more than one enhancer
subset_data <- subset(regulons_df, TF %in% tfs & Gene %in% gene_targets)
subset_data$signed_importance <- subset_data$importance_TF2G
aggregated <- aggregate(signed_importance ~ Gene + TF, data = subset_data, FUN = mean)

# Convert to long format for circos
links <- aggregated
colnames(links) <- c("to", "from", "value")  # circos expects from→to
links <- links[links$value > 0, ]  # keep positive only

# Scale link width
max_val <- max(links$value)
links$value <- links$value / max_val * 5  # scale to max width 5

# Assign link colors (use gene color)
links$col <- gene_colors[links$to]

#----------------------------- Plot
pdf(file = "Fig3_d_Circos_homeostasis_TFs_genes.pdf", width = 10, height = 10)

circos.clear()
circos.par(gap.after = c(rep(2, length(tfs)-1), 3, rep(1, length(gene_targets)-1), 3))

chordDiagram(
  links,
  order = sector_order,
  grid.col = all_colors,
  col = links$col,
  directional = 1,
  direction.type = "arrows",
  link.arr.type = "big.arrow",
  annotationTrack = "grid",
  link.lwd = links$value,        # <-- width reflects importance
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

title("TF-gene regulation (link width = importance)")

dev.off()
circos.clear()

# ---------------------------------------- Figure 3.e 
#ArchR -- gene set modules in scATAC plotted in fig3.c

# start with ArchR commads
addArchRGenome("mm10")
addArchRThreads(threads = 48) 

#load ArchR project. Find in zenodo:........
proj_fib<-loadArchRProject('Final_ArchR_object')
all_atac_genes<-getFeatures(proj_fib)

# load regulons found in zenodo
regulons_df<-read.table('/SCENICPLUS_D0_D7_D14_regulons_results/Direct_Positive_eRegulons_filtered_basedOn_GeneBased_AUC.csv', row.names = 1, sep='\t', header = T)

regulons_df %>% group_by(TF) %>%
  top_n(n = -100, wt = triplet_rank) -> regulons_df_top100


signatures<-list()

signatures$Telocytes<- unique(regulons_df_top100[regulons_df_top100$TF %in% c('Runx1','Runx2','Foxf1','Foxf2','Tcf4','Etv1'),]$Gene)
signatures$PDGFRalo<- unique(regulons_df_top100[regulons_df_top100$TF %in% c('Tcf21'),]$Gene)
signatures$Trophocytes<- unique(regulons_df_top100[regulons_df_top100$TF %in% c('Klf4','Klf2','Ebf1','Egr1','Ar','Pbx1'),]$Gene)

#keep genes that are included in ArchR object 
signatures$Telocytes<- signatures$Telocytes[signatures$Telocytes %in% all_atac_genes]
signatures$PDGFRalo<- signatures$PDGFRalo[signatures$PDGFRalo %in% all_atac_genes]
signatures$Trophocytes<- signatures$Trophocytes[signatures$Trophocytes %in% all_atac_genes]

# -------- module score --
proj_fib <- addImputeWeights(proj_fib)


proj_fib <- addModuleScore(proj_fib,
                            useMatrix = "GeneScoreMatrix",
                            name = "Module",
                            features = signatures)


# we change the order , so that order is changed in the plot as well. 
proj_fib$new_celltype_timepoint <- paste0(proj_fib$celltype,'_',proj_fib$Sample)


p <- plotGroups(
  ArchRProj = proj_fib, 
  groupBy = "new_celltype_timepoint", 
  colorBy = "cellColData", 
  name = "Module.Telocytes",
  plotAs = "violin",
  alpha = 0.4,
  baseSize = 10,
  addBoxPlot = TRUE,
)

p2 <- plotGroups(
  ArchRProj = proj_fib, 
  groupBy = "new_celltype_timepoint", 
  colorBy = "cellColData", 
  name = "Module.Trophocytes",
  plotAs = "violin",
  alpha = 0.4,
  baseSize = 10,
  addBoxPlot = TRUE,
)

p3 <- plotGroups(
  ArchRProj = proj_fib, 
  groupBy = "new_celltype_timepoint", 
  colorBy = "cellColData", 
  name = "Module.PDGFRalo",
  plotAs = "violin",
  alpha = 0.4,
  baseSize = 10,
  addBoxPlot = TRUE,
)

pdf('Fig3_e_Homeostasis_top100_homeostasis_signatures_scATAC_genescorematrix.pdf', width = 4, height = 4)
p
p2
p3
dev.off()

#------------------------------- Figure 3.f

# umaps in scATAC // 
genes<- c('Tcf21', 'Edil3', 'Tcf4', 'Bmp5', 'Ebf1', 'Cd81')

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
pdf("Fig3_f_umaps_Genescore_no_axes.pdf", width = 4, height = 4)
for (p in plots_clean) {
  print(p)
}
dev.off()

#------------------------------- Figure 3.g
# start with ArchR commads
addArchRGenome("mm10")
addArchRThreads(threads = 48) 

#load ArchR project. Find in zenodo:........
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
  select = paste(c('Tcf21','Tcf4','Runx1','Runx2','Foxf2','Ebf1','Ar','Klf4','Egr1'), collapse="|"), 
  useMatrix = "MotifMatrix"
)

#Keep the z-scores also know as deviation scores. 
markerMotifs <- grep("z:", markerMotifs, value = TRUE)
markerMotifs <- markerMotifs[markerMotifs %in% c('z:Runx2_713','z:Runx1_712','z:Ar_687','z:Foxf2_308','z:Egr1_190','z:Klf4_143','z:Ebf1_90','z:Tcf4_88','z:Tcf21_79')]

p <- plotGroups(ArchRProj = proj_fib, 
                groupBy = "celltype_timepoint", 
                colorBy = "MotifMatrix", 
                name = markerMotifs,
                imputeWeights = getImputeWeights(proj_fib)
)

pdf('Fig3_g_deviations_violin_plots.pdf',width = 5, height = 5)
p
dev.off()

#------------------------------- Figure 3.i
# load human datasets and make common variables in order to merge  -- objects not included in 

# kinchen 
seurat_obj_kin <- readRDS('seurat5_Kinchen_cellClusters_HC.rds')
seurat_obj_kin$dataset<- c('kinchen')

#Lee
seurat_obj <- readRDS('./public_data/Human_CRC/Lee_NatGen_human_CRC/All_Stromal/seurat_Lee_combined_Stromal_QC.rds')
seurat_obj<- seurat_obj[,seurat_obj$Class == 'Normal'] # keep only normal tissue
seurat_obj_lee<- seurat_obj[,seurat_obj$Cell_subtype %in% c('Myofibroblasts','Stromal 1','Stromal 2','Stromal 3')] # keep only fibroblasts
seurat_obj_lee$celltypes<- seurat_obj_lee$Cell_subtype
seurat_obj_lee$dataset<- c('Lee')

#Pelka
seurat_obj <- readRDS('./public_data/Human_CRC/Pelka_2021_human_CRC_human/analysis/Stromal/seurat_Pelka_Stromal_QC.rds')
seurat_obj<- seurat_obj[,seurat_obj$SPECIMEN_TYPE == 'N']
seurat_obj<- seurat_obj[,seurat_obj$clMidwayPr == 'Fibro']
seurat_obj_pelka<- seurat_obj[,seurat_obj$cl295v11SubFull %in% c('cS21 (Fibro stem cell niche)','cS22 (Fibro stem cell niche)','cS23 (Fibro BMP-producing)','cS24 (Fibro BMP-producing)','cS25 (Fibro CCL8+)')]
seurat_obj_pelka$celltypes<- seurat_obj_pelka$cl295v11SubFull
seurat_obj_pelka$dataset<- c('pelka')

#Qi
seurat_obj <- readRDS('./public_data/Human_CRC/Qi_2022_CRC/Stromal/seurat_Qi_stromal_QC.rds')
seurat_obj<- seurat_obj[,seurat_obj$Tissues == 'N']
seurat_obj_qi<- seurat_obj[,seurat_obj$`Cell Types` %in% c('CD24+ fibroblasts','CD73+ fibroblasts','DES+ myofibroblasts','FGFR2+ fibroblasts','FAP+ fibroblasts','MFAP5+ myofibroblasts','ICAM1- telocytes','ICAM1+ telocytes')]
seurat_obj_qi$celltypes<- seurat_obj_qi$`Cell Types`
seurat_obj_qi$dataset<- c('Qi')

#merge all of them 
seurat_obj <- merge(seurat_obj_lee, y = c(seurat_obj_kin,seurat_obj_pelka, seurat_obj_qi), add.cell.ids = c('Lee','Kinchen','Pelka','Qi'), project = "Human_healthy_colon")

#Do some basic steps of the standard analysis 

# join the layers
seurat_obj <- JoinLayers(seurat_obj)
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj,selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj,verbose = FALSE,vars.to.regress = c("nCount_RNA", "percent.mt"))
seurat_obj <- RunPCA(seurat_obj, npcs = 40)

#make the order 
seurat_obj$celltypes_dataset <- paste0(seurat_obj$celltypes,'_', seurat_obj$dataset)
seurat_obj$celltypes_dataset <- factor(seurat_obj$celltypes_dataset, levels = c('cS21 (Fibro stem cell niche)_pelka','cS22 (Fibro stem cell niche)_pelka','CD73+ fibroblasts_Qi','S3_kinchen','Stromal 3_Lee','cS23 (Fibro BMP-producing)_pelka','cS24 (Fibro BMP-producing)_pelka','S2_kinchen','Stromal 2_Lee','ICAM1+ telocytes_Qi','ICAM1- telocytes_Qi','FGFR2+ fibroblasts_Qi','cS25 (Fibro CCL8+)_pelka','CD24+ fibroblasts_Qi','Myofibroblasts_Lee','DES+ myofibroblasts_Qi','MFAP5+ myofibroblasts_Qi','FAP+ fibroblasts_Qi','S1_kinchen','Stromal 1_Lee','S4_kinchen'))
Idents(seurat_obj) <- seurat_obj$celltypes_dataset

# compute FB subtype signatures in human clusters to match our celltype annotation , based on mouse genes
human <- useEnsembl("ensembl","hsapiens_gene_ensembl", mirror = "useast")

signatures<- list()
signatures$Trophocytes<-unique(getBM(c("hgnc_symbol"), "external_gene_name",values =  c('Grem1','Pi16','Dpt','C3','Pcolce2') , human)$hgnc_symbol)
signatures$Telocytes<- unique(getBM(c("hgnc_symbol"), "external_gene_name",values = c('Pdgfra','Sox6','Bmp5','Bmp7','Wnt5a','F3'), human)$hgnc_symbol)
signatures$PDGFRalo<- unique(getBM(c("hgnc_symbol"), "external_gene_name",values = c('Adamdec1','Tcf21','Fgfr2','Col15a1','Edil3'), human)$hgnc_symbol)

sample.combined <- AddModuleScore_UCell(seurat_obj,  features = signatures)
signature.names <- paste0(names(signatures), "_UCell")

#compute mean enrichment ucell signatures per category and then plot on a heatmap
df_all <- sample.combined@meta.data[, tail(seq_len(ncol(sample.combined@meta.data)), 3)]
df_all$celltypes<-sample.combined$celltypes_dataset

df_means <- df_all %>%
  group_by(celltypes) %>%
  summarise(across(all_of(signature.names), ~ mean(.x, na.rm = TRUE), .names = "{.col}_mean"))

df_means<- as.data.frame(df_means)
rownames(df_means) <- df_means$celltypes
df_means$celltypes <- NULL

#Z-score normalization by column (signature)
df_scaled <- scale(df_means)

#Plot heatmap
pdf('Fig3_i_Signatures_MouseSubtypes_separate_publicData_human_data_categories.pdf', width = 7, height = 3)
pheatmap(t(df_scaled),
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = F,
         fontsize_row = 10,
         fontsize_col = 10,
         border_color = NA)
dev.off()

#plot mean expression of TF genes across categories 
df.sce <- as.SingleCellExperiment(seurat_obj)
pdf('Fig3_i_TF_mean_expression_categories.pdf',width = 5, height = 5)

plotGroupedHeatmap(df.sce,exprs_values='logcounts', features=c('EBF1','KLF2','KLF4','PBX1','TCF4','FOXF1','FOXF2','RUNX1','RUNX2','PITX1','TCF21'),center=T, scale = F,
                   group=c('celltypes_dataset'), zlim = c(-0.5,0.5), show_rownames=T,display_numbers = F, cluster_cols=F,cluster_rows=F, cutree_rows = 1)# + theme(text = element_text(size=20),axis.text  = element_text(size=40))

dev.off()

# compute regulon signatures in human clusters
#regulons , find in zenodo
regulons_df<-read.table('/SCENICPLUS_D0_D7_D14_regulons_results/Direct_Positive_eRegulons_filtered_basedOn_GeneBased_AUC.csv', row.names = 1, sep='\t', header = T)

regulons_df %>% group_by(TF) %>%
  top_n(n = -100, wt = triplet_rank) -> regulons_df_top100


# Separately for each regulon
# initialize empty list
signatures_regs <- list()

# TFs of interest
lst <- c('Ebf1','Klf2','Klf4','Pbx1','Tcf4','Foxf1','Foxf2','Runx1','Runx2','Pitx1','Tcf21')

for (TF in lst) {
  # extract target genes for this TF
  targets <- regulons_df_top100[regulons_df_top100$TF == TF, ]$Gene
  
  # query biomart for HGNC symbols
  gene_symbols <- unique(getBM(
    attributes = c("hgnc_symbol"),
    filters = "external_gene_name",
    values = unique(targets),
    mart = human
  )$hgnc_symbol)
  
  # store in the list with TF name as element name
  signatures_regs[[TF]] <- gene_symbols
}

sample.combined <- AddModuleScore_UCell(seurat_obj,  features = signatures_regs)
signature.names <- paste0(names(signatures_regs), "_UCell")

#take the 11 last columns with regulon signatures
df_all <- sample.combined@meta.data[, tail(seq_len(ncol(sample.combined@meta.data)), 11)]
df_all$celltypes<-sample.combined$celltypes_dataset

df_means <- df_all %>%
  group_by(celltypes) %>%
  summarise(across(all_of(signature.names), ~ mean(.x, na.rm = TRUE), .names = "{.col}_mean"))

df_means<- as.data.frame(df_means)
rownames(df_means) <- df_means$celltypes
df_means$celltypes <- NULL

#Z-score normalization by column (signature)
df_scaled <- scale(df_means)

#Plot heatmap
pdf('Fig3_i_Signatures_regulons_separate_Human_normal_publicData_categories.pdf', width = 7, height = 7)

pheatmap(t(df_scaled),
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = F,
         fontsize_row = 10,
         fontsize_col = 10,
         border_color = NA)
dev.off()

