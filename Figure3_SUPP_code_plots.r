# Code for Figure3 SUPP

#libraries
library(ArchR)

#set working directory
setwd('.../ArchR_analysis/')

# start with ArchR commads
addArchRGenome("mm10")
addArchRThreads(threads = 48) 

#load ArchR project. Not included in Zenodo processed files / 
proj_fib<-loadArchRProject('3_Save-Allcells_initial_analysis_3')


#----------- Make plots with specific colours ------------   
colours_conditions = c('cnt'='#569B9A','D7'='#8e0f62','D14'='#CD5808')
colours_celltypes = c('Trophocytes'="#8E7692",'PDGFRalo'='#a5af37','Telocytes'='#416522','SMCs'='gold','BECs'='firebrick','LECs'='mediumorchid4','Pericytes'='darkorange3','Glial'='slateblue3')

p1<-plotEmbedding(ArchRProj = proj_fib, baseSize = 0, size = 0.7, labelSize=0, colorBy = "cellColData", name = "Sample", embedding = "UMAP_unintegrated", pal = colours_conditions) +theme(axis.line =  element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(), panel.background = element_blank()) 
p2<-plotEmbedding(ArchRProj = proj_fib, baseSize = 0, size = 1, labelSize=0, colorBy = "cellColData", name = "all_celltypes", embedding = "UMAP_unintegrated", pal=colours_celltypes) + theme(axis.line =  element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(), panel.background = element_blank())

pdf('Fig3_a_SUPP_All_cells_unintegrated_UMAP_sample.pdf', width = 5, height = 5)
p1
dev.off()

pdf('Fig3_b_SUPP_All_cells_unintegrated_UMAP_celltypes.pdf', width = 5, height = 5)
p2
dev.off()

#----------------------------Figure SUPP 3

# -------- module score 
proj_fib <- addImputeWeights(proj_fib)

features <- list(
  Pericytes = c("Rgs5", "Notch3", "Pdgfrb"),
  SMCs = c("Tagln", "Hhip", "Cnn1"),
  BECs = c('Pecam1','Kdr','Plvap'),
  LECs = c('Lyve1','Prox1','Kdr'),
  Glial = c('S100b','Plp1','Gfap'),
  Trophocytes = c('Pi16','C3','Cd81'),
  Telocytes = c('Pdgfra','Sox6','Bmp5','F3'),
  PDGFRalo = c('Edil3','Fgfr2','Col15a1')
)

proj_fib <- addModuleScore(proj_fib,
                            useMatrix = "GeneScoreMatrix",
                            name = "Module",
                            features = features)

p1 <- plotEmbedding(proj_fib,
                    embedding = "UMAP_unintegrated",
                    colorBy = "cellColData",
                    name="Module.Pericytes", size=2,
                    imputeWeights = getImputeWeights(proj_fib),plotAs='points')
p2 <- plotEmbedding(proj_fib,
                    embedding = "UMAP_unintegrated",
                    colorBy = "cellColData",
                    name="Module.BECs",size=2,
                    imputeWeights = getImputeWeights(proj_fib),plotAs='points')
p3 <- plotEmbedding(proj_fib,
                    embedding = "UMAP_unintegrated",
                    colorBy = "cellColData",
                    name="Module.LECs",
                    imputeWeights = getImputeWeights(proj_fib),plotAs='points', size=2)
p4 <- plotEmbedding(proj_fib,
                    embedding = "UMAP_unintegrated",
                    colorBy = "cellColData",
                    name="Module.Trophocytes",
                    imputeWeights = getImputeWeights(proj_fib),plotAs='points', size=2)
p5 <- plotEmbedding(proj_fib,
                    embedding = "UMAP_unintegrated",
                    colorBy = "cellColData",
                    name="Module.Telocytes",
                    imputeWeights = getImputeWeights(proj_fib),plotAs='points', size=2)
p6 <- plotEmbedding(proj_fib,
                    embedding = "UMAP_unintegrated",
                    colorBy = "cellColData",
                    name="Module.PDGFRalo",
                    imputeWeights = getImputeWeights(proj_fib),plotAs='points', size=2)
p7 <- plotEmbedding(proj_fib,
                    embedding = "UMAP_unintegrated",
                    colorBy = "cellColData",
                    name="Module.SMCs",
                    imputeWeights = getImputeWeights(proj_fib),plotAs='points', size=2)
p8 <- plotEmbedding(proj_fib,
                    embedding = "UMAP_unintegrated",
                    colorBy = "cellColData",
                    name="Module.Glial",
                    imputeWeights = getImputeWeights(proj_fib),plotAs='points', size=2)



pdf('Fig3_c_SUPP_UMAPs_signatures_all_celltypes.pdf', width = 5, height = 5)
p1
p2
p3
p4
p5
p6
p7
p8
dev.off()

