# Code for Figure2 SUPP

#libraries 
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(SCpubr)

# set working directory 
setwd('/your home/')

# loading the scRNA-seq object with the Harmony integration of D0, D7, D14 all stromal cells // not included in processed files
seurat_obj<-readRDS(file = "seurat5_Harmony_integration_all_samples_treated_together_D7_D14.rds")

unique(seurat_obj$celltypes_origIdent)

# colours  
colours_conditions = c('D0'='#569B9A','D0_kinchen'="#A1CEC5",'D7_kinchen'='#8e0f62','D14'='#CD5808')
colours_celltypes = c('Trophocytes'="#8E7692",'PDGFRalo'='#a5af37','Telocytes'='#416522','SMCs'='gold','BECs'='firebrick','LECs'='mediumorchid4','Pericytes'='darkorange3','Proliferating'='turquoise4')

#------------------ Fig2.a supp
pdf('Fig2_a_SUPP_Quality_control_metrics.pdf', width =9, height = 7)
VlnPlot(seurat_obj, alpha = 0.1, group.by = 'sample', features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp","log10GenesPerUMI",'percent.Largest.Gene'), cols = colours_conditions, ncol = 3)
dev.off()

#------------------ Fig2.b supp
pdf('Fig2_b_Supp_UMAPs_allCells_sample.pdf', width = 5, height = 5)
SCpubr::do_DimPlot(sample = seurat_obj, pt.size = 0.4, group.by = 'sample', colors.use= colours_conditions)
dev.off()

#------------------ Fig2.c supp
pdf('Fig2_c_Supp_UMAPs_allCells.pdf', width = 5, height = 5)
SCpubr::do_DimPlot(sample = seurat_obj, pt.size = 0.4, group.by = 'celltypes',colors.use=colours_celltypes)
dev.off()

#------------------ Fig2.d supp
genes<-c("Pi16", "Cd81","C3","Edil3","Fgfr2", "Col15a1","Pdgfra", "Sox6", "Bmp5",'Myh11','Acta2','Tagln','Pecam1','Flt1','Eng','Lyve1','Prox1','Pdpn','Rgs5','Mcam','Pdgfrb','Stmn1','Pcna','Rrm2')

pdf('Fig2_d_Supp_Markers_celltypes_dotplot.pdf', width = 7, height = 5)

SCpubr::do_DotPlot(sample = seurat_obj, 
                   features = genes, 
                   cluster = F,
                   dot.scale = 10, group.by = 'celltypes', font.size = 10)

dev.off()
