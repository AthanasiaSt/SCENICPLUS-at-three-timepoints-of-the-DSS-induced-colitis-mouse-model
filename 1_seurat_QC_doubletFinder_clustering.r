# 1. Code with the basic processing of scRNA analysis with 
# - quality control -- seurat
# - doublet extraction -- DoubletFinder
# - initial clustering and basic preprocessing - seurat 

# ----------------------- Loading packages 
library(dplyr)
library(Seurat)
library(ggplot2)
library(readxl)
library(writexl)
library(Matrix)
library(DoubletFinder)

setwd('/home/astavropoulou/')

set.seed(38)

#---------- Run this separately for each sample healthy, healthyK, inflammationK, regeneration // 
# 1. Quality control
#---------------Load the dataset -- in house 
healthy<-Read10X(data.dir = "./Cellranger_Results_various_datasets/control_Colon_scRNA/VK_C_count/outs/filtered_feature_bc_matrix/")
#regeneration<-Read10X(data.dir = "./Cellranger_Results_various_datasets/regeneration/filtered_feature_bc_matrix/")

#----------------Loading kinchen et al dataset
#healthyK<-Read10X(data.dir = "./Cellranger_Results_various_datasets/Mouse_HC_aggr/outs/count/filtered_feature_bc_matrix/")
#inflammationK<-Read10X(data.dir = "./Cellranger_Results_various_datasets/Mouse_DSS_aggr/outs/count/filtered_feature_bc_matrix/")

# --------------------------------------------------------------------------------------
# ------------------------only for Kinchen et al dataset
samples <- read.csv(file.path("./Cellranger_Results_various_datasets/Mouse_DSS_aggr/outs/", "aggregation.csv"), stringsAsFactors=F)
cells <- new("seurat", raw.data=healthyK)

# -------- include the information of each mouse  
cellcodes <- as.data.frame(cells@raw.data@Dimnames[[2]])
colnames(cellcodes) <- "barcodes"
rownames(cellcodes) <- cellcodes$barcodes

cellcodes$libcodes <- as.factor(gsub(pattern=".+-", replacement="", cellcodes$barcodes))
cellcodes$samples <- as.vector(samples$sample_id[cellcodes$libcodes])

sampleidentity <- cellcodes["samples"]

healthyK <- CreateSeuratObject(counts = healthyK, project = "healthy_kinchen", min.cells = 10, min.features = 200)
healthyK<-healthyK[,cellcodes$barcodes]
healthyK$individuals<- cellcodes$samples
# --------------------------------------------------------------------------------------
# Initialize the Seurat object with the raw (non-normalized data). Change the min.cells to filter out more genes 
seurat_initial <- CreateSeuratObject(counts = healthy, project = "healthy", min.cells = 10, min.features = 200)

#---------------manually filter the genes based on percentage of cells
counts <- GetAssayData(seurat_initial, layer="counts", assay="RNA")   
genes.percent.expression <- rowMeans(counts>0 )*100   
genes.filter <- names(genes.percent.expression[genes.percent.expression>1]) # exclude genes that are expressed in less than 1% of cells
counts.sub <- counts[genes.filter,]
new_seurat_object <- CreateSeuratObject(counts=counts.sub) # new object with filtered features/genes

# ---------------- Compute metrics to use for the QC --------------------------
#compute the mitochondrial percentage of each cell
new_seurat_object[["percent.mt"]] <- PercentageFeatureSet(new_seurat_object, pattern = "^mt-")

#compute the percentage of ribosomal genes
new_seurat_object[["percent.rp"]] <- PercentageFeatureSet(new_seurat_object,pattern="^Rpl|^Rps") 

#complexity
new_seurat_object$log10GenesPerUMI <- log10(new_seurat_object$nFeature_RNA) / log10(new_seurat_object$nCount_RNA)

#calculate largest gene percentage of reads 
dataset<-new_seurat_object
dataset[rownames(dataset) != "Malat1",] -> dataset # exclude Malat1
apply(
  dataset@assays$RNA@layers$counts,
  2,
  max
) -> dataset$largest_count

apply(
  dataset@assays$RNA@layers$counts,
  2,
  which.max
) -> dataset$largest_index

rownames(dataset)[dataset$largest_index] -> dataset$largest_gene

100 * dataset$largest_count / dataset$nCount_RNA -> dataset$percent.Largest.Gene

dataset$largest_gene -> new_seurat_object$largest_gene
dataset$percent.Largest.Gene -> new_seurat_object$percent.Largest.Gene


#--------------Choose thresholds based on upper/bottom quantiles for each metric -----------------
#new_seurat_object$individuals<- healthyK$individuals #only for kinchen et al
data<-new_seurat_object # new object in order to test various quality thresholds

#set thresholds
mt<-quantile(data$percent.mt, .95)

largest<-quantile(data$percent.Largest.Gene, .95)

up_rp<-quantile(data$percent.rp, .95)
down_rp<-quantile(data$percent.rp, .05)

complexity<-quantile(data$log10GenesPerUMI, .05)

up_genes<-quantile(data$nFeature_RNA, .95)
down_genes<-quantile(data$nFeature_RNA, .05)

up_counts<-quantile(data$nCount_RNA, .93)
down_counts<-quantile(data$nCount_RNA, .07)

# Visualize QC metrics before and after filtering
pdf('testing_QC_metrics_healthy.pdf', width = 20, height = 15)

#plots
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp","percent.Largest.Gene","log10GenesPerUMI"), ncol = 3)

#scatter plots  
plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# make scatterplots with colours 
# initial dataset 
data[[]] %>%
  arrange(percent.mt) %>%
  ggplot(aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) + 
  geom_point() + 
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
  ggtitle("Plotting QC metrics") +
  
  
  data[[]] %>%
  arrange(percent.mt) %>%
  ggplot(aes(nCount_RNA,nFeature_RNA,colour=log10GenesPerUMI)) + 
  geom_point() + 
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
  ggtitle("Plotting QC metrics")


#filter dataset based on thresholds
data <- subset(data, subset = nFeature_RNA >= down_genes  &  nCount_RNA >= down_counts &  nCount_RNA <= up_counts & percent.mt <= mt & percent.rp >= down_rp & percent.Largest.Gene<= largest & log10GenesPerUMI >= complexity)

#------------------after filtering
# Visualize QC metrics as a violin plot
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp","log10GenesPerUMI",'percent.Largest.Gene'), ncol = 3)

data[[]] %>%
  arrange(percent.mt) %>%
  ggplot(aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) + 
  geom_point(size=0.7) + 
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
  ggtitle("Plotting QC metrics after thresholds for Tumor") +
  scale_x_log10() + scale_y_log10() 


#scatter plots
plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(plot1 + plot2)

dev.off()
# rerun the previous step of changing thresholds and filtering while doing manual inspection 

#---------------------- Starting the initial analysis to remove doublets - doublet finder ---------
# use the corrected dataset 
name_sample='healthy.pdf'

#preprocessing 
norm='default'
if (norm == 'default'){
  value_in_doubletFinder<-FALSE
  my_data <- NormalizeData(data,assay = 'RNA', normalization.method = "LogNormalize", scale.factor = 10000)
  my_data <- FindVariableFeatures(my_data, selection.method = "vst", nfeatures = 2000)
  my_data <- ScaleData(my_data,  verbose = F, vars.to.regress = c("nFeature_RNA", "percent.mt")) # regress out the effect of these variables , test also this vars.to.regress = c("nFeature_RNA", "percent.mt") maybe not in this stage
}else{
  my_data <- SCTransform(data,assay = 'RNA', verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mt")) # use vars.to.regress = "percent.mt" ? if you want 
  value_in_doubletFinder<-TRUE
}

#-----------------PCA-----------------Find the best number of PCS 
dim<-40
my_data <- RunPCA(my_data, verbose = F, npcs = dim, features = VariableFeatures(object = my_data))

#plots
ElbowPlot(my_data, ndims = dim)

#-------------------------UMAP 
my_data <- RunUMAP(my_data, dims = 1:dim)
my_data <- FindNeighbors(my_data, dims = 1:dim)
my_data <- FindClusters(my_data, resolution = 0.3, algorithm = 1) #change resolution to change the number of clusters in the UMAP, change the number to use the leiden algorithm

#find marker genes to check if all clusters are meaningful
my_data.markers <- FindAllMarkers(my_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
my_data.markers <- my_data.markers %>% filter (p_val_adj <= 0.05) 

my_data.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top10

#---------------plots for metrics and markers across clusters 
DimPlot(my_data, reduction = "umap", pt.size = 1, label=TRUE, label.size = 10)
VlnPlot(my_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp","percent.Largest.Gene","log10GenesPerUMI"), ncol = 3)
FeaturePlot(my_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp","percent.Largest.Gene","log10GenesPerUMI"), ncol = 3)

#plot the markers
DoHeatmap(my_data, features = top10$gene) + NoLegend()
DotPlot(my_data, features = unique(top10$gene), cols = c("blue", "red","pink",'yellow','black','purple'), dot.scale = 8,cluster.idents=F) + RotatedAxis() + FontSize(10)  + ggtitle('Gene expression across tumor clusters')

#--------------filter out clusters with bad quality metrics and no meaningful marker genes before running doublet finder //
# in case you discard a cluster repeat the pre-processing steps of seurat analysis until cell clusters are meaningful 
#my_data<-my_data[,my_data$seurat_clusters != 11]


# ---------------- Doublet finder -- based on the corresponding tutorial ------------------------------------------------------------------------------

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_pbmc <- paramSweep(my_data, PCs = 1:dim, sct = value_in_doubletFinder)
sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)

pK <- bcmvn_pbmc %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- my_data@meta.data$seurat_clusters #put the different clusters you identified 
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.076*nrow(my_data@meta.data))  ## Assuming 7.6% doublet formation rate - tailor for your dataset --> depending on the number of initial cells
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


# run doubletFinder 
my_data <- doubletFinder(my_data, 
                         PCs = 1:dim, 
                         pN = 0.25, 
                         pK = pK, 
                         nExp = nExp_poi.adj,
                         reuse.pANN = FALSE, sct = value_in_doubletFinder)

#---------plotting 
pdf(paste0('singlets_doublets_', name_sample))

# visualize doublets
DimPlot(my_data, reduction = 'umap', group.by = rev(names(my_data[[]]))[1])
VlnPlot(my_data, features = "nFeature_RNA", group.by = rev(names(my_data[[]]))[1], pt.size = 0.1)

# number of singlets and doublets
barplot(table(rev(my_data[[]])[1]))

dev.off()

#keep the singlets based on DoubletFinder
my_data = my_data[, my_data@meta.data[, rev(names(my_data[[]]))[1]] == "Singlet"]

#reload the raw initial dataset 
healthy<-Read10X(data.dir = "./Cellranger_Results_various_datasets/control_Colon_scRNA/VK_C_count/outs/filtered_feature_bc_matrix/")
# Initialize the Seurat object with the raw (non-normalized data). Change the min.cells to filter out more genes 
seurat_initial <- CreateSeuratObject(counts = healthy, project = "Healthy", min.cells = 10, min.features = 200)

#keep cells/genes after quality control and doublet finder
seurat_initial<-seurat_initial[rownames(my_data),colnames(my_data)]
#seurat_initial$individuals<- my_data$individuals # only for Kinchen
#save
saveRDS(seurat_initial, "seurat_healthy_colon_qc_doublet.rds")
#--------------------------------------------------------------------------------------------

# After quality control and doublet extraction, redo basic steps of seurat analysis
# take output from the previous section
# Do this for each sample separately 

seurat_initial<- readRDS('./healthy_colon/initial_QC/seurat_healthy_colon_qc_doublet.rds')

#plot the final quality metrics 
#compute the mitochondrial percentage of each cell  --> should not be too high
seurat_initial[["percent.mt"]] <- PercentageFeatureSet(seurat_initial, pattern = "^mt-")

#compute the percentage of ribosomal genes --> should not be too low nor too high
seurat_initial[["percent.rp"]] <- PercentageFeatureSet(seurat_initial,pattern="^Rpl|^Rps|^Rps") 

#complexity
seurat_initial$log10GenesPerUMI <- log10(seurat_initial$nFeature_RNA) / log10(seurat_initial$nCount_RNA)

#calculate largest gene percentage of reads 
dataset<-seurat_initial
dataset[rownames(dataset) != "Malat1",] -> dataset
apply(
  dataset@assays$RNA@layers$counts,
  2,
  max
) -> dataset$largest_count

apply(
  dataset@assays$RNA@layers$counts,
  2,
  which.max
) -> dataset$largest_index

rownames(dataset)[dataset$largest_index] -> dataset$largest_gene

100 * dataset$largest_count / dataset$nCount_RNA -> dataset$percent.Largest.Gene

dataset$largest_gene -> seurat_initial$largest_gene
dataset$percent.Largest.Gene -> seurat_initial$percent.Largest.Gene

pdf('2_final_cells_metrics_healthy.pdf',width = 10, height = 10)
# Visualize QC metrics as a violin plot
VlnPlot(seurat_initial, group.by = 'orig.ident', features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp","log10GenesPerUMI",'percent.Largest.Gene'), ncol = 3)
VlnPlot(seurat_initial, group.by = 'individuals', features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp","log10GenesPerUMI",'percent.Largest.Gene'), ncol = 3)

seurat_initial[[]] %>%
  arrange(percent.mt) %>%
  ggplot(aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) + 
  geom_point(size=0.7) + 
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
  ggtitle("Plotting QC metrics after thresholds for Tumor") +
  scale_x_log10() + scale_y_log10() 


#scatter plots
plot1 <- FeatureScatter(seurat_initial, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_initial, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(plot1 + plot2)

dev.off()

#--------------------------------------------------------------------------------------------
#------------------------------ Clustering/cell type annotation/Markers ---------------------------------------------------------------

#------------------------------ preprocessing 
norm='default'
if (norm == 'default'){
  my_data <- NormalizeData(seurat_initial,assay = 'RNA', normalization.method = "LogNormalize", scale.factor = 10000)
  my_data <- FindVariableFeatures(my_data, selection.method = "vst", nfeatures = 2000) # you can try here higher numbers of most variable genes 
  my_data <- ScaleData(my_data,  verbose = F, vars.to.regress = c("nCount_RNA", "percent.mt")) # regress out the effect of these variables , test also this vars.to.regress = c("nFeature_RNA", "percent.mt") maybe not in this stage
}else{
  my_data <- SCTransform(data,assay = 'RNA',vst.flavor = "v2", verbose = FALSE,  vars.to.regress = c("nCount_RNA", "percent.mt")) # use vars.to.regress = "percent.mt" ? if you want 
}

#-------------------------PCA -- Find the best number of PCS
dim<-40
my_data <- RunPCA(my_data, verbose = F, npcs = dim, features = VariableFeatures(object = my_data))

#----------------------UMAP 
my_data <- RunUMAP(my_data,reduction = "pca", dims = 1:dim)
my_data <- FindNeighbors(my_data,reduction = "pca", dims = 1:dim)
my_data <- FindClusters(my_data, resolution = 0.3, algorithm = 1) # change the number for a different algorithm

#------------------- Markers
#find marker genes to check if you have to filter out a cluster 
my_data.markers <- FindAllMarkers(my_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
my_data.markers<-my_data.markers[my_data.markers$p_val_adj <= 0.05,] 

my_data.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

#write_xlsx(my_data.markers,'Allmarkers_Padjvalue_0_05_D7_Kinchen_colon.xlsx')

pdf('3_healthy_colon_qc_doublet_final_clusters.pdf', width = 20, height = 10)

# plot top markers
DoHeatmap(my_data, features = top10$gene) + NoLegend()
DotPlot(my_data, features = unique(top10$gene), cols = c("blue", "red","pink",'yellow','black','purple'), dot.scale = 8,cluster.idents=F) + RotatedAxis() + FontSize(10)  + ggtitle('Gene expression across tumor clusters')

DimPlot(my_data, reduction = "umap", pt.size = 1, label = TRUE, label.size = 10, repel = TRUE)
DimPlot(my_data, split.by = 'individuals', reduction = "umap", pt.size = 1, label = TRUE, label.size = 10, repel = TRUE)

DimPlot(my_data, reduction = "pca", pt.size = 1, label = TRUE, label.size = 10, repel = TRUE)

VlnPlot(my_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp","percent.Largest.Gene","log10GenesPerUMI"), ncol = 3)

#--------remove some clusters if it is needed and rerun the analysis after quality control 
#my_data<-my_data[,(my_data$seurat_clusters != 8) ]

#save final object
saveRDS(my_data, "seurat_QC_doublet_normalized_clusters_healthy_colon.rds")
