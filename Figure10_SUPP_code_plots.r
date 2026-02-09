# Code for Figure 9 SUPP

#libraries 
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(SCpubr)
library(UCell)
library(ggplot2)
library(ggpubr)
library(FSA)      # for dunnTest
library(dplyr)
library(purrr)

# set working directory 
setwd('/home/astavropoulou/D0_D7_D14_PAPER_FINAL_FILES_RUNS/analysis_and_plots/')

# loading the scRNA-seq object with the Harmony integration of Liu et al dataset// not included in processed files
seurat_obj<-readRDS(file = "./11_Nfia_in_Liu_dataset/2_Fibroblasts/Seurat_fibroblasts_DSS_Liu.rds")

# colours  
colours_celltypes = c('Trophocytes'="#8E7692",'PDGFRalo'='#a5af37','Telocytes'='#416522')

#------------------ Fig9.a supp
pdf('Fig9_a_Supp_UMAPs_fibroblasts.pdf', width = 5, height = 5)
SCpubr::do_DimPlot(sample = seurat_obj, pt.size = 0.4, group.by = 'celltypes_2', colors.use= colours_celltypes)
dev.off()

#------------------ Fig9.b supp
pdf('Fig9_b_Supp_UMAPs_fibroblasts_days.pdf', width = 5, height = 5)
SCpubr::do_DimPlot(sample = seurat_obj, pt.size = 0.4, group.by = 'day')
dev.off()

#------------------ Fig9.c supp
genes<-c("Pi16", "Cd81","C3","Edil3","Fgfr2", "Col15a1","Pdgfra", "F3", "Bmp5")

pdf('Fig9_c_Supp_Markers_celltypes_dotplot.pdf', width = 6, height = 4)

SCpubr::do_DotPlot(sample = seurat_obj, 
                   features = genes, 
                   cluster = F,
                   dot.scale = 10, group.by = 'celltypes_2', font.size = 10)

dev.off()


#------------------ Fig9.d supp
#-----------------------------Figure 5.b
# convert seurat object to SingleCellExperiment
df.sce <- as.SingleCellExperiment(seurat_obj)
df.sce$celltype_stim <- factor(df.sce$celltype_stim, levels =c('0_Trophocytes','0_PDGFRalo','0_Telocytes','6_Trophocytes','6_PDGFRalo','6_Telocytes','7_Trophocytes','7_PDGFRalo','7_Telocytes','9_Trophocytes','9_PDGFRalo','9_Telocytes','11_Trophocytes','11_PDGFRalo','11_Telocytes','12_Trophocytes','12_PDGFRalo','12_Telocytes','16_Trophocytes','16_PDGFRalo','16_Telocytes','18_Trophocytes','18_PDGFRalo','18_Telocytes','19_Trophocytes','19_PDGFRalo','19_Telocytes','23_Trophocytes','23_PDGFRalo','23_Telocytes','42_Trophocytes','42_PDGFRalo','42_Telocytes'))

#-------------------plot the mean expression across celltypes 
pdf('Fig9_d_Supp_TF_expression.pdf',width = 6, height = 3)
plotGroupedHeatmap(df.sce,exprs_values='logcounts', features=c('Nfkb1','Bach1','Smad3','Ets2','Ets1','Nfe2l2','Fosl1','Cebpb','Fosl2','Nfia','Nfib','Egr3','Tcf7l2','Pbx1'),center=T, scale = F,
                   group=c('celltype_stim'), zlim = c(-0.5,0.5), show_rownames=T,display_numbers = F, cluster_cols=F,cluster_rows=T, cutree_rows = 1)# + theme(text = element_text(size=20),axis.text  = element_text(size=40))
dev.off()

#------------------ Fig9.e supp

#regulons
regulons_df<-read.table('/home/astavropoulou/D0_D7_D14_PAPER_FINAL_FILES_RUNS/Processed_datasets_for_scRNA_scATAC_scenicplus/SCENICPLUS_D0_D7_D14_regulons_results/Direct_Positive_eRegulons_filtered_basedOn_GeneBased_AUC.csv', row.names = 1, sep='\t', header = T)

regulons_df %>% group_by(TF) %>%
  top_n(n = -200, wt = triplet_rank) -> regulons_df_top100

signatures<-list()

# characterize each state
signatures$inflammation<- unique(regulons_df_top100[regulons_df_top100$TF %in% c('Nfkb1','Bach1','Smad3','Ets2','Ets1','Nfe2l2','Fosl1','Cebpb'),]$Gene)
signatures$regeneration<- unique(regulons_df_top100[regulons_df_top100$TF %in% c('Fosl2','Nfia','Nfib','Egr3','Tcf7l2','Pbx1'),]$Gene)
signatures$NFIA<- unique(regulons_df_top100[regulons_df_top100$TF %in% c('Nfia'),]$Gene)
signatures$NFIB<- unique(regulons_df_top100[regulons_df_top100$TF %in% c('Nfib'),]$Gene)

seurat_obj <- AddModuleScore_UCell(seurat_obj,  features = signatures, maxRank = 3000)
signature.names <- paste0(names(signatures), "_UCell")

#Figure 9 SUPP e.
seurat_obj$celltype_stim <- factor(seurat_obj$celltype_stim, levels =c('0_Trophocytes','0_PDGFRalo','0_Telocytes','6_Trophocytes','6_PDGFRalo','6_Telocytes','7_Trophocytes','7_PDGFRalo','7_Telocytes','9_Trophocytes','9_PDGFRalo','9_Telocytes','11_Trophocytes','11_PDGFRalo','11_Telocytes','12_Trophocytes','12_PDGFRalo','12_Telocytes','16_Trophocytes','16_PDGFRalo','16_Telocytes','18_Trophocytes','18_PDGFRalo','18_Telocytes','19_Trophocytes','19_PDGFRalo','19_Telocytes','23_Trophocytes','23_PDGFRalo','23_Telocytes','42_Trophocytes','42_PDGFRalo','42_Telocytes'))

pdf('Fig9_e_SUPP_Activation_regulons_top200_triplet_signature.pdf', width = 12, height =6)
VlnPlot(seurat_obj, features = 'inflammation_UCell', group.by = 'celltype_stim',  ncol = 1, pt.size = 0.2, alpha=0.3, combine = T)
VlnPlot(seurat_obj, features = 'regeneration_UCell', group.by = 'celltype_stim',  ncol = 1,  pt.size = 0.2, alpha=0.3,  combine = T)
VlnPlot(seurat_obj, features = 'NFIA_UCell', group.by = 'celltype_stim',  ncol = 1,  pt.size = 0.2, alpha=0.3,  combine = T)
VlnPlot(seurat_obj, features = 'NFIB_UCell', group.by = 'celltype_stim',  ncol = 1,  pt.size = 0.2, alpha=0.3,  combine = T)
dev.off()

# different ordering
timepoints <- c(0, 6, 7, 9, 11, 12, 16, 18, 19, 23, 42)
celltypes  <- c("Trophocytes", "PDGFRalo", "Telocytes")

levels_new <- as.vector(outer(timepoints, celltypes, paste, sep = "_"))

seurat_obj$celltype_stim <- factor(
  seurat_obj$celltype_stim,
  levels = levels_new
)

# alternative 
#Figure 9 SUPP e.
pdf('Fig9_e_SUPP_Activation_regulons_top200_triplet_signature_2.pdf', width = 12, height =6)
VlnPlot(seurat_obj, features = 'inflammation_UCell', group.by = 'celltype_stim',  ncol = 1, pt.size = 0.2, alpha=0.3, combine = T)
VlnPlot(seurat_obj, features = 'regeneration_UCell', group.by = 'celltype_stim',  ncol = 1,  pt.size = 0.2, alpha=0.3,  combine = T)
VlnPlot(seurat_obj, features = 'NFIA_UCell', group.by = 'celltype_stim',  ncol = 1,  pt.size = 0.2, alpha=0.3,  combine = T)
VlnPlot(seurat_obj, features = 'NFIB_UCell', group.by = 'celltype_stim',  ncol = 1,  pt.size = 0.2, alpha=0.3,  combine = T)
dev.off()
