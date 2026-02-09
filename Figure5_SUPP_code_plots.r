# Code for Figure 5 SUPP

#libraries 
library(Seurat)
library(scater)
library(ArchR)

setwd('/yout pathway/')


# figure 5 SUPP a 
# Heatmap with mean expression of TF genes across other mouse datases of healthy colon
colours_celltypes = c('Trophocytes'="#8E7692",'CD81-stroma'='#a5af37','SEMFs'='#416522')

# --------------------- Load other D0 Mouse datasets -- not included in Zenodo processed data
#Fazilaty et al
df <- readRDS(file = "/home/astavropoulou/D0_D7_D14_PAPER_FINAL_FILES_RUNS/Fleming_PC/D0_Fazilaty.rds")

df$celltypes<- as.character(df$celltypes)
df[[]] %>%
  mutate(celltypes = case_when(
    endsWith(celltypes, "rophocytes") ~ "Trophocytes",
    endsWith(celltypes, "DGFRalo") ~ "CD81-stroma",
    endsWith(celltypes, "elocytes") ~ "SEMFs"
  )) -> df[[]]

df$celltypes <- factor(df$celltypes, levels = c('Trophocytes','CD81-stroma','SEMFs'))
df$sample <- 'Fazilaty'
df.sce_faz <- as.SingleCellExperiment(df)

#Ho et al
df <- readRDS(file = "/home/astavropoulou/D0_D7_D14_PAPER_FINAL_FILES_RUNS/Fleming_PC/D0_H0_et_al.rds")
df$celltypes <- factor(df$celltypes, levels = c('Trophocytes','CD81-stroma','SEMFs'))
df$sample <- 'Ho'
df.sce_ho <- as.SingleCellExperiment(df)

#Ho et al
df <- readRDS(file = "/home/astavropoulou/D0_D7_D14_PAPER_FINAL_FILES_RUNS/Fleming_PC/H20_chr_DSS.rds")
df$celltypes <- factor(df$celltypes, levels = c('Trophocytes','CD81-stroma','SEMFs'))
df$sample <- 'Jasso'
df.sce_Jasso <- as.SingleCellExperiment(df)


#Load our integrated datase for D0 kinchen and in house D0 
seurat_obj<-readRDS(file = "./scRNA_seurat_integration_Healthy_Inflammation_Regeneration/Seurat5_Harmony_integration_all_samples_Healthy_Inflammation_Regeneration_fibroblasts.rds")
seurat_obj$celltypes <- factor(seurat_obj$celltypes, levels = c('Trophocytes','CD81-stroma','SEMFs'))
df.sce_in_house <- as.SingleCellExperiment(seurat_obj[,seurat_obj$sample %in% c('Healthy')])
df.sce_kinchen <- as.SingleCellExperiment(seurat_obj[,seurat_obj$sample %in% c('HealthyK')])

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
dev.off()
# Open one PDF
pdf("Percentage_celltypes_D0_datasets.pdf", width = 4, height = 1.5)

# Loop and plot
for (name in names(datasets)) {
  df.sce <- datasets[[name]]
  df<-as.data.frame(table(df.sce$sample,df.sce$celltypes)/rowSums(table(df.sce$sample,df.sce$celltypes)))
  
  p <- ggplot(df,                  # Stacked barplot using ggplot2
              aes(x = Freq,
                  y = Var1,
                  fill = Var2)) + 
    labs( 
      y="samples", x='percentage')+
    
    geom_bar(stat = "identity")+  scale_fill_manual(values = colours_celltypes)  + theme(text = element_text(size=3),axis.text  = element_text(size=3)) 
  
  print(p)  # <-- ensure it renders into the PDF
}
# Close device
dev.off()

# Combine all datasets first
df_all <- do.call(rbind, lapply(names(datasets), function(name) {
  df.sce <- datasets[[name]]
  tab <- table(df.sce$sample, df.sce$celltypes)
  prop <- prop.table(tab, margin = 1)  # percentages per sample
  
  df <- as.data.frame(prop)
  colnames(df) <- c("sample", "celltype", "percentage")
  df$dataset <- name
  df
}))

pdf("Percentage_celltypes_D0_datasets.pdf", width = 4, height = 4)
ggplot(df_all,
            aes(x = percentage,
                y = sample,
                fill = celltype)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colours_celltypes) +
  labs(x = "Percentage", y = "Sample")
dev.off()

#Figure 5 SUPP b. 

#Load our integrated dataset -- find in zenodo processed files 
seurat_obj<-readRDS(file = "/scRNA_seurat_integration_D0_D7_D14/Seurat5_Harmony_integration_all_samples_D0_D7_D14_fibroblasts.rds")
df.sce <- as.SingleCellExperiment(seurat_obj[,seurat_obj$orig.ident %in% c('Regeneration','InflammationK')])
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
