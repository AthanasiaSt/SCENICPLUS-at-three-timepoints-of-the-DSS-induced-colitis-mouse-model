# Code for Figure 5 SUPP

#libraries 
library(Seurat)
library(scater)
library(ArchR)

setwd('/yout pathway/')


# figure 5 SUPP a 
# Heatmap with mean expression of TF genes across other mouse datases of healthy colon

# --------------------- Load other D0 Mouse datasets -- not included in Zenodo processed data
#Fazilaty et al
df <- readRDS(file = "/home/astavropoulou/D0_D7_D14_PAPER_FINAL_FILES_RUNS/Fleming_PC/D0_Fazilaty.rds")
df$celltypes <- factor(df$celltypes, levels = c('Trophocytes','PDGFRalo','Telocytes'))
df.sce_faz <- as.SingleCellExperiment(df)

#Ho et al
df <- readRDS(file = "/home/astavropoulou/D0_D7_D14_PAPER_FINAL_FILES_RUNS/Fleming_PC/D0_H0_et_al.rds")
df$celltypes <- factor(df$celltypes, levels = c('Trophocytes','PDGFRalo','Telocytes'))
df.sce_ho <- as.SingleCellExperiment(df)

#Ho et al
df <- readRDS(file = "/home/astavropoulou/D0_D7_D14_PAPER_FINAL_FILES_RUNS/Fleming_PC/H20_chr_DSS.rds")
df$celltypes <- factor(df$celltypes, levels = c('Trophocytes','PDGFRalo','Telocytes'))
df.sce_Jasso <- as.SingleCellExperiment(df)


#Load our integrated datase for D0 kinchen and in house D0 
seurat_obj<-readRDS(file = "/home/astavropoulou/seurat5_Harmony_integration_all_samples_treated_together_D7_D14_Subclustering.rds")
df.sce_in_house <- as.SingleCellExperiment(seurat_obj[,seurat_obj$orig.ident %in% c('D0')])
df.sce_in_house$celltypes<- df.sce_in_house$celltypes_origIdent_older
df.sce_kinchen <- as.SingleCellExperiment(seurat_obj[,seurat_obj$orig.ident %in% c('D0_kinchen')])
df.sce_kinchen$celltypes<- df.sce_kinchen$celltypes_origIdent_older

#-------------------plot the mean of groups across cells 

# Define datasets and labels
datasets <- list(
  Kinchen    = df.sce_kinchen,
  InHouse    = df.sce_in_house,
  Jasso      = df.sce_Jasso,
  Ho         = df.sce_ho,
  Fazilaty   = df.sce_faz
)

order_genes <- c('Ar',"Ebf1" ,"Klf4",'Klf2',"Pbx1",'Egr1','Maf','Tcf21', 
                 'Pitx1', 'Tcf4',"Etv1","Foxf2","Foxf1","Runx1","Runx2")

# Open one PDF
pdf("Fig5_SUPP_a_TFs_mean_expression_D0_datasets.pdf", width = 2, height = 4)

# Loop and plot
for (name in names(datasets)) {
  df.sce <- datasets[[name]]
  
  p <- plotGroupedHeatmap(
    df.sce,
    exprs_values = "logcounts",
    features = order_genes,
    center = TRUE,
    scale = FALSE,
    group = "celltypes",
    zlim = c(-0.5, 0.5),
    show_rownames = TRUE,
    display_numbers = FALSE,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    cutree_rows = 1,
    main = paste("TFs in", name)
  )
  
  print(p)  # <-- ensure it renders into the PDF
}

# Close device
dev.off()

#Figure 5 SUPP b. 

#Load our integrated dataset -- find in zenodo processed files 
seurat_obj<-readRDS(file = "/scRNA_seurat_integration_D0_D7_D14/Seurat5_Harmony_integration_all_samples_D0_D7_D14_fibroblasts.rds")
df.sce <- as.SingleCellExperiment(seurat_obj[,seurat_obj$orig.ident %in% c('D14','D7_kinchen')])
df.sce$celltypes<- df.sce$celltypes_origIdent_older


pdf("Fig5_SUPP_b_TFs_mean_expression_D7_D14_datasets.pdf", width = 4, height = 5)

plotGroupedHeatmap(
  df.sce,
  exprs_values = "logcounts",
  features = order_genes,
  center = TRUE,
  scale = FALSE,
  group = "celltype_stim_2",
  zlim = c(-0.5, 0.5),
  show_rownames = TRUE,
  display_numbers = FALSE,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  cutree_rows = 1,
  main = ('TF mean expression in D7 and D14')
)
dev.off()

#Figure 5 SUPP c.

#------------------------------- Figure 3.f

# set working directory
var<-getwd() 
setwd(paste0(var,'/Processed_datasets_for_scRNA_scATAC_scenicplus/'))

# start with ArchR commads
addArchRGenome("mm10")
addArchRThreads(threads = 48) 

#load ArchR project. Find in zenodo:........
proj_fib<-loadArchRProject('Final_ArchR_object')

# umaps in scATAC // 
genes<- c('Egr1', 'C3', 'Runx1', 'Pdgfra', 'Ar', 'Col14a1','Runx2','Sox6','Klf4','Pcolce2','Foxf1','Wnt5a','Col15a1','Fgfr2','Foxf2','Mylk')

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
pdf("Fig5_supp_c_umaps_Genescore_no_axes.pdf", width = 4, height = 4)
for (p in plots_clean) {
  print(p)
}
dev.off()
