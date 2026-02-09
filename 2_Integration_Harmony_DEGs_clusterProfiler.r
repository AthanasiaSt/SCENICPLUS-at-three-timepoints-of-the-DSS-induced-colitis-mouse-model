# 2. Code for dataset integration
# integration -- Harmony
# cell cluster annotation -- seurat 
# DEGs -- MAST
# Over-representation analysis -- ClusterProfiler
# Gene set enrichment analysis -- fgsea

#libraries
library(dplyr)
library(Seurat)
library(ggplot2)
library(Matrix)
library(harmony)
library(writexl)
library(readxl)
library(RColorBrewer)
library("clusterProfiler")
library("enrichplot")
organism = "org.Mm.eg.db"
library(organism, character.only = TRUE)
library(fgsea)
library(data.table)

setwd('/home/astavropoulou/')

set.seed(38)

# After perfomring QC and initial analysis based on the pipeline in 1. 
# continue by merging the two datasets (four samples) and addressing batch effects 

#------------------ Merging the four samples ----------------------------
#------------------------------kinchen
healthyK <- readRDS('Initial_QC_analysis/healthyK/seurat_normalized_healthyK_colon.rds')
inflammationK<- readRDS('Initial_QC_analysis/inflammationK/seurat_normalized_inflammationK_colon.rds')

# in-house 
healthy <- readRDS('Initial_QC_analysis/healthy_colon/seurat_normalized_healthy_colon.rds')
regeneration<- readRDS('Initial_QC_analysis/regeneration_colon/seurat_normalized_regeneration_colon.rds')

#find common genes between the datasets
common.features_1 <- intersect(rownames(healthyK), rownames(inflammationK))
common.features_2 <- intersect(rownames(healthy), rownames(regeneration))
common.features <- intersect(common.features_1, common.features_2) 

#--------- merge each sample based on common features
seurat_obj <- merge(healthy[common.features, ], y = c(regeneration[common.features, ],healthyK[common.features, ],inflammationK[common.features, ]), add.cell.ids = c('healthy','regeneration','healthyK','inflammationK'), project = "healthy_inflammation_regeneration")


# ---------- filter out genes expressed in less that 1% of cells 
# join the layers
seurat_obj <- JoinLayers(seurat_obj) # this is needed for seurat 5 , each sample is a different layer  

counts <- seurat_obj[['RNA']]$counts
genes.percent.expression <- rowMeans(counts>0 )*100   
genes.use <- names(genes.percent.expression[genes.percent.expression>1])

# exclude lowly expressed features
seurat_obj<-seurat_obj[genes.use,]

#------------------ make dataset label 
seurat_obj[[]] %>%
  mutate(group = case_when(
    endsWith(orig.ident, "healthy") ~ "healthy_regeneration",
    endsWith(orig.ident, "regeneration") ~ "healthy_regeneration",
    endsWith(orig.ident, "healthyK") ~ "healthyK_inflammationK",
    endsWith(orig.ident, "inflammationK") ~ "healthyK_inflammationK",
  )) -> seurat_obj@meta.data


#--------------------------------- processing the samples all together 
# join the layers
seurat_obj <- JoinLayers(seurat_obj) # do this at this point if skipped before 

#-------------- inital analysis -- no integration 
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj,selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj,verbose = FALSE,vars.to.regress = c("nCount_RNA", "percent.mt"))
seurat_obj <- RunPCA(seurat_obj, npcs = 40)

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:40, reduction = "pca")
seurat_obj <- FindClusters(seurat_obj, resolution = 0.1, cluster.name = "unintegrated_clusters")

seurat_obj <- RunUMAP(seurat_obj, dims = 1:40, reduction = "pca", reduction.name = "umap.unintegrated")

# inspect resulting umap without integration
# cells are mostly grouped by dataset of origin
pdf('Unintegrated_treat_allSamples_together.pdf', width = 25, height = 7)
DimPlot(seurat_obj, reduction = "umap.unintegrated", group.by = c("orig.ident",'group','unintegrated_clusters'))
dev.off()


#------------------ integration at the dataset level (batch)
seurat_obj <- seurat_obj %>% 
  RunHarmony(group.by.vars = "group", theta = 2, sigma = 0.1, assay.use = 'RNA', max_iter=20,  plot_convergence = TRUE, reduction.save = 'Theta2_sigma0.1_commonVargenes_onGroup')

seurat_obj <- seurat_obj %>% 
  RunUMAP(reduction = "Theta2_sigma0.1_commonVargenes_onGroup", dims = 1:40) %>% 
  FindNeighbors(reduction = "Theta2_sigma0.1_commonVargenes_onGroup", dims = 1:40) %>% 
  FindClusters(resolution = 0.1) %>%  # change for different number of clusters based on the data inspection
  identity()


#inspecting integration results
pdf('Harmony_integration.pdf', width = 6, height = 6)
DimPlot(seurat_obj, reduction = "umap", label = T, group.by = c("RNA_snn_res.0.1"))
DimPlot(seurat_obj, reduction = "umap", group.by = c('orig.ident','group'), label = F)
dev.off()

#------------------------------------- marker genes in order to annotate celltypes
my_data.markers <- FindAllMarkers(seurat_obj, group.by = 'seurat_clusters', only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
my_data.markers <- my_data.markers %>% filter (p_val_adj <= 0.05) 

my_data.markers %>%
  group_by(cluster) %>%
  top_n(n = 8, wt = avg_log2FC) -> top10

write_xlsx(my_data.markers,'MarkerGenes_onGroup_Celltypes_padj_log2fc_0_25.xlsx')

#-------------------change the names of the clusters after inspecting marker genes 
new.cluster.ids <- c("CD81_stroma","Trophocytes",'SEMFs',"CD81_stroma",'CD81_stroma','CD81_stroma','Trophocytes', "SMCs",'LECs','BECs','Proliferating','Pericytes')
names(new.cluster.ids) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)
seurat_obj$celltypes<- Idents(seurat_obj)

#------------------------------------------------------------------
#---------------change the order of samples
seurat_obj$celltypes <- factor(seurat_obj$celltypes, levels = c('Trophocytes','CD81_stroma','SEMFs'))
seurat_obj$orig.ident <- factor(seurat_obj$orig.ident, levels = c('healthy','healthyK','inflammationK','regeneration'))
seurat_obj$celltype_stim_2<- paste0(seurat_obj$orig.ident,'_',seurat_obj$celltypes)

#save final object
saveRDS(seurat_obj, "seurat_harmony_all_samples_fibroblasts.rds")

# -------------------- Differential expression analysis
#----------------------- for inflammationK against healthyK and regeneration against healthy and only on Fibroblasts
# all fibroblasts 
temp1<-c()

sample='regeneration'
sample1='healthy'

# find the group with the smallest number of cells
min_number<-min(dim(seurat_obj[,seurat_obj$orig.ident == sample])[2],dim(seurat_obj[,seurat_obj$orig.ident == sample1])[2])

# use MAST with a specific number of cells per category
temp <- FindMarkers(seurat_obj, test.use ='MAST', random.seed = 2 ,max.cells.per.ident=min_number, assay = 'RNA', group.by = 'orig.ident', ident.1 = sample, ident.2 = sample1,logfc.threshold = 0.25, min.pct = 0.01, min.cells.feature = 10, min.cells.group = 10, verbose = FALSE, latent.vars = c("nCount_RNA",'percent.mt'))
temp<-temp %>% filter (p_val_adj <= 0.05) 

if (length(rownames(temp)) > 0){
  temp$genes <- rownames(temp)
  temp$celltype <- c(paste0(sample,'_vs_',sample1))
  rownames(temp) <- NULL
  temp1<-rbind(temp1, temp)
}
write_xlsx(temp1,'DEGs_MAST_Harmony_origIdent_regeneration_VS_healthy_pAdj_filter_BalancedNumbers_AllFibroblasts.xlsx')

# for each fibroblast subtype separately
temp1<-c()
for (cluster in unique(seurat_obj$celltypes)){ # for each cell cluster
  
  print(cluster)
  sample=paste0('regeneration_',cluster)
  sample1=paste0('healthy_',cluster)

  # find the group with the smallest number of cells
  min_number<-min(dim(seurat_obj[,seurat_obj$celltype_stim_2 == sample])[2],dim(seurat_obj[,seurat_obj$celltype_stim_2 == sample1])[2])

  # use MAST with a specific number of cells per category
  temp <- FindMarkers(seurat_obj, test.use ='MAST', random.seed = 2 ,max.cells.per.ident=min_number, assay = 'RNA', group.by = 'celltype_stim_2', ident.1 = sample, ident.2 = sample1,logfc.threshold = 0.25, min.pct = 0.01, min.cells.feature = 10, min.cells.group = 10, verbose = FALSE, latent.vars = c("nCount_RNA",'percent.mt'))
  temp<-temp %>% filter (p_val_adj <= 0.05) 

  if (length(rownames(temp)) > 0){
    temp$genes <- rownames(temp)
    temp$celltype <- c(paste0(sample,'_vs_',sample1))
    rownames(temp) <- NULL
    temp1<-rbind(temp1, temp)
  }
}

write_xlsx(temp1,'DEGs_MAST_Harmony_origIdent_regeneration_VS_healthy_pAdj_filter_BalancedNumbers.xlsx')

# ------------------------ Over-representation analysis -------------------------------  
# ------------------------ Cluster profiler on the DEGs 

# load files with DEGs
DEGs_regeneration<-read_excel('DEGs_origIdent/DEGs_MAST/origIdent_integration_BALANCED/DEGs_MAST_Harmony_origIdent_regeneration_VS_healthy_pAdj_filter_BalancedNumbers_AllFibroblasts.xlsx') 
DEGs_inflammation<-read_excel('DEGs_origIdent/DEGs_MAST/origIdent_integration_BALANCED/DEGs_MAST_Harmony_origIdent_inflammation_VS_healthy_pAdj_filter_BalancedNumbers_AllFibroblasts.xlsx') 

# filter based on Log2fc / Do the same for downregulation as well
regeneration<-DEGs_regeneration[DEGs_regeneration$avg_log2FC >= 1,]
inflammation<-DEGs_inflammation[DEGs_inflammation$avg_log2FC >= 1,]

#divide them in common and condition specific lists
common<-intersect(regeneration$genes,inflammation$genes)
regeneration_specif<-setdiff(regeneration$genes,inflammation$genes)
inflammation_specif<-setdiff(inflammation$genes,regeneration$genes)

write.table(common, 'common_downDEGs_regeneration_inflammation_ALLfibs_log2fc_1.txt', quote = FALSE,row.names = FALSE)
write.table(regeneration_specif, 'regeneration_specific_downDEGs_ALLfibs_log2fc_1.txt', quote = FALSE,row.names = FALSE)
write.table(inflammation_specif, 'inflammation_specific_downDEGs_ALLfibs_log2fc_1.txt', quote = FALSE,row.names = FALSE)


#---------------functional enrichment analysis 
ora_analysis_bp <- enrichGO(gene=unique(regeneration_specif) , # put the other lists as well 
                            universe = rownames(seurat_obj), # set as universe the genes expressed in your dataset
                            OrgDb = org.Mm.eg.db,
                            keyType = "SYMBOL",
                            ont = "BP",              # either "BP", "CC" or "MF",
                            pAdjustMethod = "BH",
                            qvalueCutoff = 0.05,
                            readable = FALSE, 
                            pool = FALSE)



#-----------------export
ora_analysis_bp_simplified <- clusterProfiler::simplify(ora_analysis_bp) 
write.csv(x = as.data.frame(ora_analysis_bp_simplified@result), 
          file = "Ora_upDEGs_regeneration_specific_allFibs_avgLog2fc_1.csv")


#---------------GSEA 
pathways <- gmtPathways("/home/astavropoulou/D0_D7_D14_PAPER_FINAL_FILES_RUNS/analysis_and_plots/msidb_gene_sets/m5.go.bp.v2025.1.Mm.symbols.gmt")

#compute log2FC between all genes across two conditions
markers <- FindMarkers(
  seurat_obj,
  ident.1 = "inflammationK",
  ident.2 = "healthyK",
  group.by = "orig.ident",
  logfc.threshold = -Inf,   # no logFC filter
  min.pct = 0,              # include all genes
  min.diff.pct = -Inf,      # include all genes
  test.use = "wilcox"       
)

# we want the log2 fold change 
original_gene_list <- markers$avg_log2FC
names(original_gene_list) <- rownames(markers)
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order
gene_list_ordered <- sort(gene_list, decreasing = TRUE)

fgseaRes <- fgsea(pathways = pathways, 
                  stats    = gene_list_ordered,
                  minSize  = 15,
                  maxSize  = 500)

fwrite(fgseaRes, file="fgseaRes_inflammation_BP_GO.txt", sep="\t", sep2=c("", " ", ""))

