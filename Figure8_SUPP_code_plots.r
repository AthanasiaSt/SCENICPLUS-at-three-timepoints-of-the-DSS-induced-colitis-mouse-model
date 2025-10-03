# Code for Figure8 SUPP in the manuscript

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
seurat_obj<-readRDS(file = "/scRNA_seurat_integration_D0_D7_D14/Seurat5_Harmony_integration_all_samples_D0_D7_D14_fibroblasts.rds")

#regulons
regulons_df<-read.table('/SCENICPLUS_D0_D7_D14_regulons_results/Direct_Positive_eRegulons_filtered_basedOn_GeneBased_AUC.csv', row.names = 1, sep='\t', header = T)

regulons_df %>% group_by(TF) %>%
  top_n(n = -200, wt = triplet_rank) -> regulons_df_top100

signatures<-list()

# characterize each state
signatures$inflammation<- unique(regulons_df_top100[regulons_df_top100$TF %in% c('Nfkb1','Bach1','Smad3','Ets2','Ets1','Nfe2l2','Fosl1','Cebpb'),]$Gene)
signatures$regeneration<- unique(regulons_df_top100[regulons_df_top100$TF %in% c('Fosl2','Nfia','Nfib','Egr3','Tcf7l2','Pbx1'),]$Gene)

seurat_obj <- AddModuleScore_UCell(seurat_obj,  features = signatures, maxRank = 3000)
signature.names <- paste0(names(signatures), "_UCell")

#Figure 8 SUPP a.
pdf('Activation_regulons_top200_triplet_signature_V2.pdf', width = 6, height = 4)
VlnPlot(seurat_obj, features = 'inflammation_UCell', group.by = 'celltype_stim_2', cols = alpha(c('red2','blue3','forestgreen','red2','blue3','forestgreen','magenta4','orange','yellow','steelblue2','palevioletred','peachpuff2'), 0.6),  ncol = 1, pt.size = 0, combine = T)
VlnPlot(seurat_obj, features = 'regeneration_UCell', group.by = 'celltype_stim_2', cols = alpha(c('red2','blue3','forestgreen','red2','blue3','forestgreen','magenta4','orange','yellow','steelblue2','palevioletred','peachpuff2'), 0.6),  ncol = 1, pt.size = 0, combine = T)
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


#Figure 8 SUPP b. 
pdf('FigS8_b_Functional_enrichments_activation_regulons_top200Genes.pdf', width = 10, height =10)
dotplot(ora_analysis_bp_simplified_inf, showCategory=15) + theme(text = element_text(size=20),axis.text.y = element_text(size=20),axis.text.x = element_text(size=20)) + scale_size(range = c(5, 15))  # increase min/max dot size
dotplot(ora_analysis_bp_simplified_reg, showCategory=15) + theme(text = element_text(size=20),axis.text.y = element_text(size=20),axis.text.x = element_text(size=20)) + scale_size(range = c(5, 15))  # increase min/max dot size
dev.off()

# Fig 8 SUPP c. 

# Log2FC comparing accessibility as genescores across groups // Nfia , Nfib // 

colours_celltypes = c('Trophocytes'="#8E7692",'PDGFRalo'='#a5af37','Telocytes'='#416522')

#marker genes based on gene scores -- output from ArchR can be dowloaded as supplementary table from the manuscript
markers_atac<-read.table('Markers_celltype_Stim_padj_0_05_genescorematrix.csv', header = T)

genes<- c('Nfia','Nfib')

pdf("Fig8_SUPP_c_genes_barplots_log2fc_genescores.pdf", width = 7, height = 5)

for (gene in genes) {
  
  tmp<-markers_atac[markers_atac$name == gene,]
  
  # --- prepare dataframe ---
  df2 <- tmp %>%
    mutate(
      # split "cnt_PDGFRalo" etc. into two parts
      timepoint = str_extract(group_name, "cnt|D7|D14"),
      celltype  = str_remove(group_name, "cnt_|D7_|D14_")
    )
  
  # order timepoints
  df2$timepoint <- factor(df2$timepoint, levels = c("cnt", "D7", "D14"))
  
  # order celltypes if you want a consistent order
  df2$celltype <- factor(df2$celltype, levels = c('Trophocytes',"PDGFRalo", "Telocytes"))
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
    geom_col(position = position_dodge2(width = 0.8, preserve = "single"), width = 0.7)+
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_text(aes(y = y_pos, label = sig),
              position = position_dodge2(width = 0.8, preserve = "single"),
              size = 6)+
    scale_fill_manual(values = c('Trophocytes'="#8E7692",'PDGFRalo'='#a5af37','Telocytes'='#416522')) +
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

