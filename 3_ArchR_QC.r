# 3. Code with initial preprocessing and quality control of scATAC data 
# - quality control -- ArchR ,following basic tutorial: https://www.archrproject.com/bookdown/index.html

#----------------- Libraries
library(ArchR)

setwd('/home/astavropoulou/D0_D7_D14_PAPER_FINAL_FILES_RUNS/ArchR_new/')
set.seed(38) # set seed for reproducibility

#----------------- Initialize ArchR
addArchRGenome("mm10")              # set genome to mouse mm10
addArchRThreads(threads = 48)       # use multiple threads (speed-up)


# Valid barcodes from CellRanger outputs (filtered barcodes only)
valid_cells_cnt<-as.list(read.csv('/home/astavropoulou/Cellranger_Results_various_datasets/Control_colon_scATAC/outs/filtered_peak_bc_matrix/barcodes.tsv', header = F,sep = '\n'))
names(valid_cells_cnt) <- 'cnt'

valid_cells_D7<-as.list(read.csv('/home/astavropoulou/Cellranger_Results_various_datasets/scATAC_D7_DSS/outs/filtered_peak_bc_matrix/barcodes.tsv', header = F,sep = '\n'))
names(valid_cells_D7) <- 'D7'

valid_cells_D14<-as.list(read.csv('/home/astavropoulou/Cellranger_Results_various_datasets/scATAC_D14_DSS/outs/filtered_peak_bc_matrix/barcodes.tsv', header = F,sep = '\n'))
names(valid_cells_D14) <- 'D14'


# Create Arrow file for Control
ArrowFiles_cnt <- createArrowFiles(
  inputFiles = c('/home/astavropoulou/Cellranger_Results_various_datasets/Control_colon_scATAC/outs/fragments.tsv.gz'),
  sampleNames = c('cnt'),
  outputNames = c('cnt'),
  minTSS = 8,      # minimum TSS enrichment score
  minFrags = 5500, # minimum fragments per cell
  addTileMat = TRUE,      # create TileMatrix
  validBarcodes = valid_cells_cnt, # use filtered barcodes
  addGeneScoreMat = TRUE  # create GeneScoreMatrix
)

# Create Arrow file for inflammation
ArrowFiles_D7 <- createArrowFiles(
  inputFiles = c('/home/astavropoulou/Cellranger_Results_various_datasets/scATAC_D7_DSS/outs/fragments.tsv.gz'),
  sampleNames = c('D7'),
  outputNames = c('D7'),
  minTSS = 6,
  minFrags = 5500,
  addTileMat = TRUE,
  validBarcodes = valid_cells_D7,
  addGeneScoreMat = TRUE
)

# Create Arrow file for regeneration
ArrowFiles_D14 <- createArrowFiles(
  inputFiles = c('/home/astavropoulou/Cellranger_Results_various_datasets/scATAC_D14_DSS/outs/fragments.tsv.gz'),
  sampleNames = c('D14'),
  outputNames = c('D14'),
  minTSS = 8,
  minFrags = 5500,
  addTileMat = TRUE,
  validBarcodes = valid_cells_D14,
  addGeneScoreMat = TRUE
)


###---------------------------------------------------------------------------

#----------------- Doublet detection
# Compute doublet scores for each dataset
doubScores <- addDoubletScores(
  input = c('./cnt.arrow','./D7.arrow','./D14.arrow'),
  k = 10,          # number of neighbors for pseudo-doublet projection
  knnMethod = "UMAP", # embedding for NN search
  LSIMethod = 1
)

#----------------- Merge Arrow files into one ArchRProject
cnt_D7_D14_AOM_DSS_archR <- ArchRProject(
  ArrowFiles =  c('./cnt.arrow','./D7.arrow','./D14.arrow'), 
  outputDirectory = "D0_D7_D14_scATAC_analysis",
  copyArrows = T # keeps original Arrow files safe
)


#----------------- Plot QC metrics
# Ridge plots of TSS enrichment per sample
p1 <- plotGroups(
  ArchRProj = cnt_D7_D14_AOM_DSS_archR, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "ridges"
)

# Save QC plot
png('TSSEnrichment.png', width = 5, height = 5, res=1200, unit='in')
p1
dev.off()

#---------- plot QC metrics for each sample
idxSample <- BiocGenerics::which(cnt_D7_D14_AOM_DSS_archR$Sample %in% "D7")
cellsSample <- cnt_D7_D14_AOM_DSS_archR$cellNames[idxSample]

df <- getCellColData(cnt_D7_D14_AOM_DSS_archR[cellsSample, ], select = c("log10(nFrags)", "TSSEnrichment"))

p <- ggPoint(
  x = df[,1], 
  y = df[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
)

png('Quality_metrics_D7.png', width = 8, height = 8, res=1200, unit='in')
p
dev.off()


#---------------------------------------------------------------------------------------
# Filtering cells based on thresholds chosen from outliers
# manually chose the samples with most extreme outliers per variable 
# We used the same thresholds across the three samples 

# - Use D14 distribution for doublet enrichment
# - Use D7 distribution for TSS enrichment and fragment counts

# D14 cells
idxSample <- BiocGenerics::which(cnt_D7_D14_AOM_DSS_archR$Sample %in% "D14")
cellsSample <- cnt_D7_D14_AOM_DSS_archR$cellNames[idxSample]

# Filter doublet enrichment (keep ≤ 95th percentile of D14)
quantile(cnt_D7_D14_AOM_DSS_archR[cellsSample, ]$DoubletEnrichment, probs=0.95)
idxPass <- which(cnt_D7_D14_AOM_DSS_archR$DoubletEnrichment <= quantile(cnt_D7_D14_AOM_DSS_archR[cellsSample, ]$DoubletEnrichment, probs=0.95))
cellsPass <- cnt_D7_D14_AOM_DSS_archR$cellNames[idxPass]
cnt_D7_D14_AOM_DSS_archR<-cnt_D7_D14_AOM_DSS_archR[cellsPass, ]

# Get D7 cells
idxSample <- BiocGenerics::which(cnt_D7_D14_AOM_DSS_archR$Sample %in% "D7")
cellsSample <- cnt_D7_D14_AOM_DSS_archR$cellNames[idxSample]

# Filter TSS enrichment (keep ≥ 5th percentile of D7)
idxPass <- which(cnt_D7_D14_AOM_DSS_archR$TSSEnrichment >= quantile(cnt_D7_D14_AOM_DSS_archR[cellsSample, ]$TSSEnrichment, probs=0.05))
cellsPass <- cnt_D7_D14_AOM_DSS_archR$cellNames[idxPass]
cnt_D7_D14_AOM_DSS_archR<-cnt_D7_D14_AOM_DSS_archR[cellsPass, ]

# Filter nFrags (keep ≥ 5th percentile of D7)
idxPass <- which(cnt_D7_D14_AOM_DSS_archR$nFrags >= quantile(cnt_D7_D14_AOM_DSS_archR[cellsSample, ]$nFrags, probs=0.05))
cellsPass <- cnt_D7_D14_AOM_DSS_archR$cellNames[idxPass]
cnt_D7_D14_AOM_DSS_archR<-cnt_D7_D14_AOM_DSS_archR[cellsPass, ]

# Remove doublets based on computed scores
cnt_D7_D14_AOM_DSS_archR <- filterDoublets(cnt_D7_D14_AOM_DSS_archR, filterRatio = 1)


#----------------- Final QC visualization
# Example: plot for D7
idxSample <- BiocGenerics::which(cnt_D7_D14_AOM_DSS_archR$Sample %in% "D7")
cellsSample <- cnt_D7_D14_AOM_DSS_archR$cellNames[idxSample]

# Scatter of log10(nFrags) vs TSS enrichment
df <- getCellColData(cnt_D7_D14_AOM_DSS_archR[cellsSample, ], select = c("log10(nFrags)", "TSSEnrichment"))

p <- ggPoint(
  x = df[,1], 
  y = df[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
)

# Save plot
png('Quality_metrics_D7.png', width = 8, height = 8, res=1200, unit='in')
p
dev.off()

#----------------- Save final filtered project
cnt_D7_D14_AOM_DSS_archR <-saveArchRProject(
  ArchRProj = cnt_D7_D14_AOM_DSS_archR, 
  outputDirectory = "1_Archr_Cnt_Day7_Day14_qc", 
  load = T
)

