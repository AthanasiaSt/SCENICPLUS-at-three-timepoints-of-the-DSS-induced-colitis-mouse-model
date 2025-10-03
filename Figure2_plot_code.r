# Code for main Figure 2

#libraries
library(ArchR)
library(pheatmap)
library(RColorBrewer)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggraph)
library(grid)
library(igraph)

# Download processed data from Zenodo: Processed_datasets_for_scRNA_scATAC_scenicplus 

# set working directory
var<-getwd() 
setwd(paste0(var,'/Processed_datasets_for_scRNA_scATAC_scenicplus/'))

# start with ArchR commads
addArchRGenome("mm10")
addArchRThreads(threads = 48) 

#load ArchR project
proj_fib<-loadArchRProject('Final_ArchR_object')

# colours 
colours_conditions = c('cnt'='#569B9A','D7'='#8e0f62','D14'='#CD5808')
colours_celltypes = c('Trophocytes'="#8E7692",'PDGFRalo'='#a5af37','Telocytes'='#416522')

#figure2.a,d,e
p1 <- plotEmbedding(ArchRProj = proj_fib,baseSize = 0, size = 0.7, labelSize=0, colorBy = "cellColData", name = "Clusters_unintegrated", embedding = "UMAP_unintegrated")+theme(axis.line =  element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(), panel.background = element_blank()) 
p2<-plotEmbedding(ArchRProj = proj_fib, baseSize = 0, size = 1, labelSize=0, colorBy = "cellColData", name = "celltype", embedding = "UMAP_unintegrated", pal=colours_celltypes) + theme(axis.line =  element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(), panel.background = element_blank())
p3<-plotEmbedding(ArchRProj = proj_fib, baseSize = 0, size = 0.7, labelSize=0, colorBy = "cellColData", name = "Sample", embedding = "UMAP_unintegrated", pal = colours_conditions) +theme(axis.line =  element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(), panel.background = element_blank()) 

pdf('Fig2_a_d_e_ArchR_UMAPs.pdf', width = 5, height = 5)
p1
p2
p3
dev.off()

#figure2.b
proj_fib <- addImputeWeights(proj_fib)

features <- list(
  Trophocytes = c('Pi16','C3','Cd81','Ebf1','Ar', 'Klf4', 'Egr1'),
  PDGFRalo = c('Tcf21', 'Edil3', 'Fgfr2', 'Col15a1'),
  Telocytes = c('Pdgfra', 'Sox6', 'Bmp5', 'Tcf4', 'Foxf2', 'Runx1', 'Runx2','F3')
)

#add module scores based on gene scores 
proj_fib <- addModuleScore(proj_fib,
                            useMatrix = "GeneScoreMatrix",
                            name = "Module",
                            features = features)


p1 <- plotEmbedding(proj_fib,
                    embedding = "UMAP_unintegrated",
                    colorBy = "cellColData",
                    name="Module.Trophocytes",
                    imputeWeights = getImputeWeights(proj_fib),plotAs='points', size=2)

p2 <- plotEmbedding(proj_fib,
                    embedding = "UMAP_unintegrated",
                    colorBy = "cellColData",
                    name="Module.PDGFRalo",
                    imputeWeights = getImputeWeights(proj_fib),plotAs='points', size=2)

p3 <- plotEmbedding(proj_fib,
                    embedding = "UMAP_unintegrated",
                    colorBy = "cellColData",
                    name="Module.Telocytes",
                    imputeWeights = getImputeWeights(proj_fib),plotAs='points', size=2)

pdf('Fig2_b_ArchR_modules_FB_subtypes.pdf', width = 5, height = 5)
p1
p2
p3
dev.off()

#figure2.c
# common peaks across clusters of 2.a

# load peaks from PeakCalls output in ArchR folder corresponding to initial clustering
# List only .rds files starting with "C"
peak_files <- list.files(
  "~/Final_ArchR_object/PeakCalls/",
  pattern = "^C.*\\.rds$",     # ^C means "starts with C"
  full.names = TRUE
)

# Clean filenames for naming: remove suffix "-reproduciblePeaks.gr.rds"
names(peak_files) <- gsub("-reproduciblePeaks\\.gr\\.rds$", "", basename(peak_files))


# Extract peak sets as GRanges
peak_list <- lapply(peak_files, readRDS)  # Each is a GRanges


#Compute shared peak ratio matrix across cell clusters
celltypes <- names(peak_list)
n <- length(celltypes)
ratio_matrix <- matrix(NA, nrow = n, ncol = n, dimnames = list(celltypes, celltypes))

for (i in seq_len(n)) {
  for (j in seq_len(n)) {
    peaks_i <- peak_list[[i]]
    peaks_j <- peak_list[[j]]
    
    shared <- sum(countOverlaps(peaks_i, peaks_j) > 0)
    denom <- min(length(peaks_i), length(peaks_j))
    
    ratio_matrix[i, j] <- shared / denom
  }
}

# heatmap with peak overlap across clusters
my_colors <- colorRampPalette(brewer.pal(9, "Reds"))(1000)  
pdf('Fig2_c_overlap_of_ATAC_peaks.pdf', width = 4, height = 4)
pheatmap(ratio_matrix,
         color=my_colors,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         display_numbers = F,
         main = "Shared Peak Ratio Between Cell Types")

dev.off()

#figure2.f
# compute again marker peaks per category
markersPeaks <- getMarkerFeatures(
  ArchRProj = proj_fib, 
  useMatrix = "PeakMatrix", 
  groupBy = "celltype_timepoint",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)


#make the heatmap with specific order of categories
heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks[,c('cnt_Trophocytes','D7_Trophocytes','D14_Trophocytes','cnt_PDGFRalo','D7_PDGFRalo','D14_PDGFRalo','cnt_Telocytes','D7_Telocytes','D14_Telocytes')], 
  cutOff = "FDR <= 0.01 & Log2FC >= 1",
  transpose = T,
  nLabel = 4,
  nPrint = 3,
  clusterCols = FALSE)


pdf('Fig2_f_Marker_Peaks_heatmap.pdf', width = 5.2, height = 4)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()

#not needed
#proj_fib <- addMotifAnnotations(ArchRProj = proj_fib, motifSet = "cisbp", name = "Motif",force = TRUE)

#Compute motif enrichments in the marker peaks 
motifsUp <- peakAnnoEnrichment(
  seMarker = markersPeaks[,c('cnt_Trophocytes','D7_Trophocytes','D14_Trophocytes','cnt_PDGFRalo','D7_PDGFRalo','D14_PDGFRalo','cnt_Telocytes','D7_Telocytes','D14_Telocytes')],
  ArchRProj = proj_fib,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.01 & Log2FC >= 1"
)

#make the heatmap
heatmapEM <- plotEnrichHeatmap(motifsUp,clusterCols=FALSE,returnMatrix=T, n=15, transpose = T) 
heatmapEM_scaled<-t(scale(t(heatmapEM), center = T, scale=T)) # scale per row (category) and then transpose again

pdf('Fig2_f_Upmotifs_in_markerPeaks.pdf', width = 15, height = 5)
ComplexHeatmap::Heatmap(heatmapEM_scaled, cluster_rows =F, cluster_columns =F)
dev.off()


#figure 2.h

# network of TF regulons regulating other TFs 

# load filtered regulons -- from zenodo processed files 
regulons_df<-read.table('/SCENICPLUS_D0_D7_D14_regulons_results/Direct_Positive_eRegulons_filtered_basedOn_GeneBased_AUC.csv', row.names = 1, sep='\t', header = T)

# keep gene targets that are themselves TFs
regulons_df<-regulons_df[regulons_df$Gene %in% unique(regulons_df$TF),]

#keep unique combinations of the two
df<-unique(regulons_df[,c('TF','Gene')])

#make a directed graph
g1 <- graph_from_data_frame(d = df, vertices = unique(c(df$TF,df$Gene)), directed = T)

# compute a clustering for node colors
g1_undir <- as.undirected(g1, mode = "collapse")  # merge parallel edges
V(g1)$clu <- as.character(membership(cluster_louvain(g1_undir, resolution = 1)))

# compute degree as node size
V(g1)$size <- degree(g1)

got_palette <- c(
  '2' = '#8e0f62', '1' ='#CD5808', '3'='#416522')

# Network plot 
pdf("Fig2_h_Directed_Network.pdf", width = 9, height =9)

ggraph(g1, layout = "stress") +
  # directed edges with arrows
  geom_edge_link(
    edge_colour = "grey66",
    arrow = arrow(length = unit(3, "mm")),
    end_cap = circle(3, "mm")
  ) +
  # nodes colored by cluster and sized by degree
  geom_node_point(aes(fill = clu, size = size), shape = 21) +
  # node labels
  geom_node_text(
    aes(label = name),
    hjust = 0.2,
    vjust = -0.8,
    size = 4,
    family = "sans"
  ) +
  # node fill colors
  scale_fill_manual(values = got_palette) +
  # theme
  theme_graph() +
  theme(legend.position = "none")

dev.off()

