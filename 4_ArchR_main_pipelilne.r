# 4. Code with main ArchR analysis following quality control 
# - basic tutorial: https://www.archrproject.com/bookdown/index.html

#----------------- Libraries
library(ArchR)
ref_genome="BSgenome.Mmusculus.UCSC.mm10"
BSgenome::getBSgenome(genome = ref_genome)       # load reference genome
library(BSgenome.Mmusculus.UCSC.mm10)            # mm10 genome package

setwd('/home/astavropoulou/D0_D7_D14_PAPER_FINAL_FILES_RUNS/ArchR_new/')
set.seed(38) # set seed for reproducibility

#----------------- Load project from QC step
projHeme1<-loadArchRProject('1_Archr_Cnt_Day7_Day14_qc')

#----------------------------------------- Initial analysis 
# ---------------------------- dimensionality reduction without batch correction
projHeme1 <- addIterativeLSI(
  ArchRProj = projHeme1,
  useMatrix = "TileMatrix",    # use genome-wide bin matrix
  name = "IterativeLSI", 
  iterations = 6, 
  clusterParams = list(        # clustering parameters (similar to Seurat)
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 15000,         # number of features to use
  dimsToUse = 1:30,            # number of LSI dimensions
  force=T
)

# Add clusters (unsupervised)
projHeme1 <- addClusters(
  input = projHeme1,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters_unintegrated",
  resolution = 0.5,
  force = TRUE
)

# Add UMAP embedding
projHeme1 <- addUMAP(
  ArchRProj = projHeme1, 
  reducedDims = "IterativeLSI", 
  name = "UMAP_unintegrated", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine", 
  force = T
)

# Plot UMAPs colored by sample and cluster
getCellColData(ArchRProj = projHeme1, select = NULL, drop = FALSE)

p1 <- plotEmbedding(ArchRProj = projHeme1, colorBy = "cellColData", name = "Sample", embedding = "UMAP_unintegrated")
p2 <- plotEmbedding(ArchRProj = projHeme1, colorBy = "cellColData", name = "Clusters_unintegrated", embedding = "UMAP_unintegrated")

pdf('1_initial_clustering_noBatch_correction.pdf', width = 5, height = 5)
p1
p2
dev.off()

#------------------------------------- Batch effect correction (Harmony)
projHeme1 <- addHarmony(
  ArchRProj = projHeme1,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample",
  force = T
)

# Clustering with Harmony correction
projHeme1 <- addClusters(
  input = projHeme1,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "Clusters_Harmony",
  resolution = 0.1,
  force = TRUE
)

# UMAP with Harmony correction
projHeme1 <- addUMAP(
  ArchRProj = projHeme1, 
  reducedDims = "Harmony", 
  name = "UMAPHarmony", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)

p1 <- plotEmbedding(ArchRProj = projHeme1, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
p2 <- plotEmbedding(ArchRProj = projHeme1, colorBy = "cellColData", name = "Clusters_Harmony", embedding = "UMAPHarmony")

pdf('1_initial_clustering_Harmony_correction.pdf', width = 5, height = 5)
p1
p2
dev.off()

#---------------------------------------- Identify marker genes
# Marker features from GeneScoreMatrix (accessibility near genes)
markersGS <- getMarkerFeatures(
  ArchRProj = projHeme1, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters_unintegrated",
  bias = c("TSSEnrichment", "log10(nFrags)"), # control for depth / quality
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1")

write.table(markerList,"Markers_unintegrated_clusters_allCells_log2fc1_padj_0_01_celltypes.csv",quote=F,sep="\t",row.names=F,col.names=T)

#---------------------------------------- Keep only mesenchymal cells
# Filter out specific clusters (retain mesenchymal)
idxPass <- which(!(projHeme1$Clusters_unintegrated %in% c('C1','C2','C3','C4','C10')))
cellsPass <- projHeme1$cellNames[idxPass]
projHeme1<-projHeme1[cellsPass, ]

#---------------------------------------- Imputation weights for smoothing
projHeme1 <- addImputeWeights(projHeme1)

#-------------------------- Make pseudoreplicates
projHeme1$celltype.stim<- paste0(projHeme1$Sample,'_',projHeme1$predictedGroup_Un)

# Export metadata (needed for downstream tools like SCENIC+)
write.table(getCellColData(projHeme1), file = 'Cnt_Day7_Day14_ATAC_analysis_mesenchymal.csv', sep = '\t', row.names = T)

# Restrict peak calling to fibroblast subtypes with enough cells
table(projHeme1$Sample, projHeme1$Clusters_unintegrated)

projHeme2 <- addGroupCoverages(
  ArchRProj = projHeme1, 
  groupBy = "Clusters_unintegrated", 
  minRep = 2, maxRep = 2, 
  minCells = 300, maxCells = 1000, 
  sampleRatio = 0.8
)

#-------------------------- Peak calling
pathToMacs2 <- findMacs2() # locate MACS2
projHeme2 <- addReproduciblePeakSet(
  ArchRProj = projHeme2, 
  groupBy = "celltype.stim", 
  pathToMacs2 = pathToMacs2
)

# Export peak list
write.table(getPeakSet(projHeme2), file = 'All_peaks_called_mesenchymal_clusters_unintegrated.csv', sep = '\t', row.names = F)

# Add peak accessibility matrix
projHeme2 <- addPeakMatrix(projHeme2)

# ------------------- Find marker peaks
markersPeaks <- getMarkerFeatures(
  ArchRProj = projHeme2, 
  useMatrix = "PeakMatrix", 
  groupBy = "celltype.stim",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
write.table(markerList,"Markers_fromPeaks_log2fc1_padj_0_01_Clusters_unintegrated_mesenchymal.csv",quote=F,sep="\t",row.names=F,col.names=T)

#-------------------------- Motif analysis
projHeme2 <- addMotifAnnotations(ArchRProj = projHeme2, motifSet = "cisbp", name = "Motif")

# Motif enrichment in differential peaks
motifsUp <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = projHeme2,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.01 & Log2FC >= 1"
)

saveArchRProject(ArchRProj = projHeme2, outputDirectory = "7_pseudoreplicates_peaks_celltype_stim_motifs", load = TRUE)

#----------------------------- ChromVAR deviations
if("Motif" %ni% names(projHeme2@peakAnnotation)){
  projHeme2 <- addMotifAnnotations(ArchRProj = projHeme2, motifSet = "cisbp", name = "Motif")
}
projHeme2 <- addBgdPeaks(projHeme2)

projHeme2 <- addDeviationsMatrix(
  ArchRProj = projHeme2, 
  peakAnnotation = "Motif",
  force = TRUE
)

plotVarDev <- getVarDeviations(projHeme2, name = "MotifMatrix", plot = F)

projHeme2$celltype.stim<-as.character(projHeme2$celltype.stim)

# Filter motifs (custom exclusion list here)
markerMotifs <- getFeatures(projHeme2, select = paste(unlist(features), collapse="|"), useMatrix = "MotifMatrix")
markerMotifs <- grep("z:", markerMotifs, value = TRUE)
markerMotifs <- markerMotifs[markerMotifs %ni% c('z:Ppara_670','z:Rara_685','z:Arx_500','z:Rarb_663','z:Ppard_656','z:Rarg_654','z:Pparg_653','z:Smarcc2_646','z:Arx_500','z:Barhl2_492','z:Barx2_488','z:Barhl1_457','z:Barx1_432','z:Arnt2_26','z:Arid5a_9',"z:Srebf1_32",'z:Arid5b_8','z:Arid3a_7','z:Arid3b_6','z:Arid3c_10','z:Arntl_89','z:Arnt_25','z:Arid2_11')]

# Extract motif matrix
MotifMatrix<-getMatrixFromProject(
  ArchRProj = projHeme2,
  useMatrix = "MotifMatrix",
  useSeqnames = 'z',
  verbose = TRUE,
  binarize = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)
