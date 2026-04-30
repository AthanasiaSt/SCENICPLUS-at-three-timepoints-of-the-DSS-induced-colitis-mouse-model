# Code for Figure11 SUPP in the manuscript

#libraries 
library(Seurat)
library(UCell)
library("clusterProfiler")
library("enrichplot")
organism = "org.Mm.eg.db"
library(organism, character.only = TRUE)

# Download from zenodo the processed files 
# set working directory
var<-getwd() 
setwd(paste0(var,'/Processed_datasets_for_scRNA_scATAC_scenicplus/'))

# loading the scRNA-seq object with the Harmony integration of D0, D7, D14 fibroblasts
seurat_obj<-readRDS(file = "./scRNA_seurat_integration_Healthy_Inflammation_Regeneration/Seurat5_Harmony_integration_all_samples_Healthy_Inflammation_Regeneration_fibroblasts.rds")

#regulons
regulons_df<-read.table('./SCENICPLUS_Healthy_Inflammation_Regeneration_regulons_results/Direct_Positive_eRegulons_filtered_basedOn_GeneBased_AUC.csv', row.names = 1, sep='\t', header = T)

regulons_df %>% group_by(TF) %>%
  top_n(n = -200, wt = triplet_rank) -> regulons_df_top100

signatures<-list()

# characterize each state
signatures$inflammation<- unique(regulons_df_top100[regulons_df_top100$TF %in% c('Nfkb1','Bach1','Smad3','Ets2','Ets1','Nfe2l2','Fosl1','Cebpb'),]$Gene)
signatures$regeneration<- unique(regulons_df_top100[regulons_df_top100$TF %in% c('Fosl2','Nfia','Nfib','Egr3','Tcf7l2','Pbx1'),]$Gene)

seurat_obj <- AddModuleScore_UCell(seurat_obj,  features = signatures, maxRank = 3000)
signature.names <- paste0(names(signatures), "_UCell")

#figure 11 SUPP a.
pdf('fig9_supp_a_Activation_regulons_top200_triplet_signature.pdf', width = 6, height = 4)
VlnPlot(seurat_obj, features = 'inflammation_UCell', group.by = 'celltypes_condition', cols = alpha(c('red2','blue3','forestgreen','red2','blue3','forestgreen','magenta4','orange','yellow','steelblue2','palevioletred','peachpuff2'), 0.6),  ncol = 1, pt.size = 0, combine = T)
VlnPlot(seurat_obj, features = 'regeneration_UCell', group.by = 'celltypes_condition', cols = alpha(c('red2','blue3','forestgreen','red2','blue3','forestgreen','magenta4','orange','yellow','steelblue2','palevioletred','peachpuff2'), 0.6),  ncol = 1, pt.size = 0, combine = T)
dev.off()

#---------------functional enrichment analysis 
ora_analysis_bp_inf <- enrichGO(gene=unique(signatures$inflammation) ,
                            universe = rownames(seurat_obj), 
                            OrgDb = org.Mm.eg.db, 
                            keyType = "SYMBOL",
                            ont = "BP",              # either "BP", "CC" or "MF",
                            pAdjustMethod = "BH",
                            qvalueCutoff = 0.05,
                            readable = FALSE, 
                            pool = FALSE)
ora_analysis_bp_simplified_inf <- clusterProfiler::simplify(ora_analysis_bp_inf) 

ora_analysis_bp_reg <- enrichGO(gene=unique(signatures$regeneration) ,  #  unique(top10[top10$TF %in% signatures$Trophocytes,]$Gene)
                            universe = rownames(seurat_obj), 
                            OrgDb = org.Mm.eg.db,  # contains the TAIR/Ensembl id to GO correspondence for A. thaliana
                            keyType = "SYMBOL",
                            ont = "BP",              # either "BP", "CC" or "MF",
                            pAdjustMethod = "BH",
                            qvalueCutoff = 0.05,
                            readable = FALSE, 
                            pool = FALSE)

#-----------------export
ora_analysis_bp_simplified_reg <- clusterProfiler::simplify(ora_analysis_bp_reg) 


#figure 11 SUPP b. 
pdf('FigS9_b_Functional_enrichments_activation_regulons_top200Genes.pdf', width = 10, height =10)
dotplot(ora_analysis_bp_simplified_inf, showCategory=15) + theme(text = element_text(size=20),axis.text.y = element_text(size=20),axis.text.x = element_text(size=20)) + scale_size(range = c(5, 15))  # increase min/max dot size
dotplot(ora_analysis_bp_simplified_reg, showCategory=15) + theme(text = element_text(size=20),axis.text.y = element_text(size=20),axis.text.x = element_text(size=20)) + scale_size(range = c(5, 15))  # increase min/max dot size
dev.off()

# Fig 11 SUPP c. 

# loading the scRNA-seq object with the Harmony integration of Healthy, InflammationK, Regeneration fibroblasts
seurat_obj<-readRDS(file = "./scRNA_seurat_integration_Healthy_Inflammation_Regeneration/Seurat5_Harmony_integration_all_samples_Healthy_Inflammation_Regeneration_fibroblasts.rds")

# colours  
colours_conditions = c('Healthy'='#569B9A','HealthyK'="#A1CEC5",'InflammationK'='#8e0f62','Regeneration'='#CD5808')
colours_celltypes = c('Trophocytes'="#8E7692",'CD81-stroma'='#a5af37','SEMFs'='#416522')

# ── 1. Define IAF marker gene sets ────────────────────────────────────────────
signatures <- list(
  IAF1   = c("Egr1", "Fos", "Fosb", "Nr4a1",'Sfrp2'),
  IAF2_3 = c("Igfbp5", "Grem1", "Col18a1",'Il11'),
  IAF3   = c("Col18a1", "Mmp2",'Sdc1'),
  IAF4   = c("Mmp3", "Mmp10", "Il11", 'Cd14', 'Serpina3g', "Plau"),
  IAF5 = c('C3','Thbs1','Egr1')
)

sample.combined <- AddModuleScore_UCell(seurat_obj,  features = signatures)
signature.names <- paste0(names(signatures), "_UCell")


pdf('Cadinu_signatures_IAFs_vlplot.pdf', width = 8, height = 15)
VlnPlot(sample.combined, features = signature.names, group.by = 'celltypes_condition',  ncol = 2, pt.size = 0, combine = T)
dev.off()

# Smillie et al.
# covert Human symbols to mouse symbols
human <- useEnsembl("ensembl", "hsapiens_gene_ensembl", mirror = "useast")

genes <- unique(getBM(
  c("hgnc_symbol"),
  attributes = "mmusculus_homolog_associated_gene_name",
  values     = c('PLAU','CHI3L1','MMP3','IL1R1','IL13RA2','TNFSF11',
                 'MMP10','OSMR','IL11','STRA6','FAP','WNT2','TWIST1','IL24'),
  mart       = human
)$mmusculus_homolog_associated_gene_name)

# Remove empty strings (human genes with no mouse ortholog)
all_genes <- genes[genes != ""]
all_genes <- intersect(rownames(seurat_obj), all_genes)
signatures <- list('Smillie'=all_genes)

sample.combined <- AddModuleScore_UCell(seurat_obj,  features = signatures)
signature.names <- paste0(names(signatures), "_UCell")


pdf('Smillie_signatures_IAFs_vlplot.pdf', width = 7, height = 5)
VlnPlot(sample.combined, features = signature.names, group.by = 'celltypes_condition',  ncol = 1, pt.size = 0, combine = T)
dev.off()

#Fig 11 d.
# Log2FC comparing accessibility as genescores across groups // Nfia , Nfib // 

colours_celltypes = c('Trophocytes'="#9E7692",'PDGFRalo'='#a5af37','Telocytes'='#416522')

#marker genes based on gene scores -- output from ArchR can be dowloaded as supplementary table from the manuscript (supplementary table 4)
markers_atac<-read.table('Markers_celltype_Stim_padj_0_05_genescorematrix.csv', header = T)


genes<- c('Nfia','Nfib')

pdf("Fig9_SUPP_c_genes_barplots_log2fc_genescores.pdf", width = 7, height = 5)

for (gene in genes) {
  
  tmp<-markers_atac[markers_atac$name == gene,]
  
  # --- prepare dataframe ---
  df2 <- tmp %>%
    mutate(
      # split "cnt_PDGFRalo" etc. into two parts
      timepoint = str_extract(group_name, "Healthy|Inflammation|Regeneration"),
      celltype  = str_remove(group_name, "Healthy_|Inflammation_|Regeneration_")
    )
  
  # order timepoints
  df2$timepoint <- factor(df2$timepoint, levels = c("Healthy", "Inflammation", "Regeneration"))
  
  # order celltypes if you want a consistent order
  df2$celltype <- factor(df2$celltype, levels = c('Trophocytes',"CD81- stroma", "SEMFs"))
  df2$FDR <- as.numeric(df2$FDR)
  df2 <- df2 %>%
    mutate(
      # set offset as a fraction of the maximum absolute Log2FC across all bars
      offset = 0.05 * max(abs(Log2FC), na.rm = TRUE),
      y_pos  = ifelse(Log2FC >= 0, Log2FC + offset, Log2FC - offset),
      sig = case_when(
        FDR < 0.001 ~ "***",
        FDR < 0.01  ~ "**",
        FDR < 0.05  ~ "*",
        TRUE        ~ ""
      )
    )
  # --- plot ---
  p <- ggplot(df2, aes(x = timepoint, y = Log2FC, fill = celltype)) +
    geom_col(position = position_dodge2(width = 0.9, preserve = "single"), width = 0.7)+
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_text(aes(y = y_pos, label = sig),
              position = position_dodge2(width = 0.9, preserve = "single"),
              size = 6)+
    scale_fill_manual(values = c('Trophocytes'="#9E7692",'CD81- stroma'='#a5af37','SEMFs'='#416522')) +
    theme_minimal(base_size = 16) +
    theme(
      axis.text.x = element_text(size = 14),
      legend.title = element_blank()
    ) +
    labs(
      x = "Timepoint",
      y = paste0("Log2 Fold Change (", gene, ")"),
      title = gene
    )
  
  # print to PDF (each plot is a page)
  print(p)
}

# close PDF
dev.off()

