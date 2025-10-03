#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 15 17:08:18 2025

@author: Athanasia Stavropoulou
"""

# code for Figure3.a

# libraries 
import pandas as pd
import numpy as np
import seaborn as sns


# find enrichment of D0 marker gene sets in the regulons
#marker genes D0-specific , see supplementary data 2 in the manuscript
markers=pd.read_excel('Supplementary Data 2.xlsx',  sheet_name='D0_MarkerGenes')

#filtered positive regulons // uploaded in zenodo processed files 
regulons_df=pd.read_csv('/SCENICPLUS_D0_D7_D14_regulons_results/Direct_Positive_eRegulons_filtered_basedOn_GeneBased_AUC.csv', sep='\t')

# ---------------------------------filtering marker genes based on log2FC and percentage of cells 
# Filter markers based on effect size and percentage threshold
markers = markers[markers['avg_log2FC'] >= 0.58]
markers = markers[markers['pct.1'] >= 0.3]

# Create explicit copies of gene names for downstream use
markers['Gene'] = markers['genes'].copy()
markers['TF'] = markers['genes'].copy()


# -------------------------------------------------------------------------
# Find enrichments of marker genes in regulons
# -------------------------------------------------------------------------
df = pd.DataFrame()
for category in markers['celltype'].unique():
    # Get marker genes for this cell type
    genes = markers[markers['celltype'] == category]['Gene'].unique()
    
    # Observed: regulons containing those marker genes
    observed_regulon = regulons_df[regulons_df['Gene'].isin(genes)]
    observed_regulon = observed_regulon.groupby(['TF'])['Gene'].nunique()
    
    # Expected: all genes in regulons for each TF
    expected_regulon = regulons_df.groupby(['TF'])['Gene'].nunique()
    expected_regulon = expected_regulon[observed_regulon.index]
    
    # Compute enrichment as observed / expected
    observed_regulon_enrich = observed_regulon / expected_regulon
    
    # Baseline probability of a DEG (normalized by total features)
    prob_of_deg = len(genes) / 12263  # total number of genes in dataset
    
    # Final enrichment score normalized by baseline probability
    enrich = pd.DataFrame(observed_regulon_enrich / prob_of_deg)
    enrich['category'] = category
    enrich['TF'] = enrich.index
    
    # Reset index, sort by enrichment
    enrich = enrich.reset_index(drop=True).sort_values(['Gene'], ascending=False)    
    
    # Add results for this cell type
    df = pd.concat([df, enrich], axis=0)

# -------------------------------------------------------------------------
# Select top enrichments across regulons for each category
# -------------------------------------------------------------------------
df_temp = df.groupby('category').head(10).reset_index(drop=True)  # top 10 TFs per category
TFs = df_temp['TF'].unique()

# Restrict to these top TFs
df = df[df['TF'].isin(TFs)]

# Pivot to create matrix (cell type x TF)
matrix = df.pivot_table(index='category', columns='TF', values='Gene')
matrix[np.isnan(matrix)] = 0  # replace NaNs with 0


# -------------------------------------------------------------------------
# Reorder rows and columns manually 
# -------------------------------------------------------------------------
d = pd.DataFrame(matrix).reindex(
    index=['D0_Trophocytes','D0_PDGFRalo','D0_Telocytes'], 
    columns=['Ar', 'Ebf1', 'Klf4', 'Klf2','Hoxb4', 'Ebf3', 'Egr1', 'Pbx1',
             'Nfix','Tcf7l2', 'Egr3', 'Maf', 'Tcf21', 'Frem1', 'Pitx1',
             'Stat5b', 'Foxf2', 'Tcf4', 'Runx1', 'Runx2', 'Batf','Foxf1','Etv1']
)


# -------------------------------------------------------------------------
# Plot clustered heatmap
# -------------------------------------------------------------------------
sns.set_theme(style="darkgrid", font_scale=0.8)
g = sns.clustermap(
    d, mask=d==0, z_score=1, cmap='coolwarm',
    col_cluster=False, row_cluster=False, standard_scale=None,
    method='ward', metric='euclidean',
    annot=False, figsize=(6,3), annot_kws={"size": 2}
)

# Save figure
g.savefig('/Fig3_a_heatmap_enrichments_D0_FB_subtype_Tf_regulons.pdf')

