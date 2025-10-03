#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 17 15:02:42 2025

@author: Athanasia Stavropoulou
"""

# code for Figure5.a

#----------------------------- Libraries
import pandas as pd
import numpy as np
import seaborn as sns

#-------------------------------------- Load DEGs for D7 and D14 --- see supplementary data 5 in the manuscript
D7_DEGs = pd.read_excel('Supplementary Data 5.xlsx', sheet_name='D7k_vs_D0k_upDEGs_subtypeLevel')
D14_DEGs = pd.read_excel('Supplementary Data 5.xlsx', sheet_name='D14_vs_D0_upDEGs_subtypeLevel')


#-------------------------------------- Filter DEGs for strong upregulation (log2FC >= 1)
# Only positive DEGs because regulons describe only positive relationships, based on our filtering
D7_DEGs = D7_DEGs[D7_DEGs['avg_log2FC'] >= 1]
D7_DEGs['Gene'] = D7_DEGs['genes'].copy()
D7_DEGs['TF'] = D7_DEGs['genes'].copy()

D14_DEGs = D14_DEGs[D14_DEGs['avg_log2FC'] >= 1]
D14_DEGs['Gene'] = D14_DEGs['genes'].copy()
D14_DEGs['TF'] = D14_DEGs['genes'].copy()

#-------------------------------------- Load filtered positive regulons, find in zenodo
regulons_df = pd.read_csv(
    '/SCENICPLUS_D0_D7_D14_regulons_results/Direct_Positive_eRegulons_filtered_basedOn_GeneBased_AUC.csv',
    sep='\t'
)

#-------------------------------------- Function to compute regulon enrichment
def compute_enrichment(DEGs_subset):
    """
    Compute enrichment of DEGs in each regulon per cell type.
    Returns a dataframe with TF, category, and enrichment value ('Gene').
    """
    df = pd.DataFrame()
    for category in DEGs_subset['celltype'].unique():
        # Get all genes in this category
        genes = DEGs_subset[DEGs_subset['celltype'] == category]['Gene'].unique()
        
        # Observed number of regulon genes in this category
        observed_regulon = regulons_df[regulons_df['Gene'].isin(genes)]
        observed_regulon = observed_regulon.groupby('TF')['Gene'].nunique()
        
        # Expected number of genes per TF
        expected_regulon = regulons_df.groupby('TF')['Gene'].nunique()
        expected_regulon = expected_regulon[observed_regulon.index]
        
        # Compute enrichment (observed / expected / probability of DEG)
        observed_regulon_enrich = observed_regulon / expected_regulon
        prob_of_deg = len(genes) / 12263  # total number of genes in dataset
        enrich = pd.DataFrame(observed_regulon_enrich / prob_of_deg, columns=['Gene'])
        
        # Add category and TF columns
        enrich['category'] = category
        enrich['TF'] = enrich.index
        enrich = enrich.reset_index(drop=True).sort_values(['Gene'], ascending=False)
        
        # Combine into final dataframe
        df = pd.concat([df, enrich], axis=0)
    return df

#-------------------------------------- Compute enrichment for D7 and D14
df_d7 = compute_enrichment(D7_DEGs)
df_d14 = compute_enrichment(D14_DEGs)

# Get all unique TFs from D7 and D14 enrichments
TFs_d7 = df_d7['TF'].unique()
TFs_d14 = df_d14['TF'].unique()

#-------------------------------------- Filter so as TFs are upregulated themselves (log2FC >= 0.58) in at least one category
D7_DEGs = pd.read_excel('Supplementary Data 5.xlsx', sheet_name='D7k_vs_D0k_upDEGs_subtypeLevel')
D14_DEGs = pd.read_excel('Supplementary Data 5.xlsx', sheet_name='D14_vs_D0_upDEGs_subtypeLevel')

D7_DEGs = D7_DEGs[D7_DEGs['avg_log2FC'] >= 0.58]
D7_DEGs['Gene'] = D7_DEGs['genes'].copy()
D7_DEGs['TF'] = D7_DEGs['genes'].copy()
D7_DEGs = D7_DEGs[D7_DEGs['Gene'].isin(set(list(TFs_d14) + list(TFs_d7)))]

D14_DEGs = D14_DEGs[D14_DEGs['avg_log2FC'] >= 0.58]
D14_DEGs['Gene'] = D14_DEGs['genes'].copy()
D14_DEGs['TF'] = D14_DEGs['genes'].copy()
D14_DEGs = D14_DEGs[D14_DEGs['Gene'].isin(set(list(TFs_d14) + list(TFs_d7)))]

#-------------------------------------- Combine TFs that are DEGs in at least one comparison
tfs = set(list(D14_DEGs['genes'].unique()) + list(D7_DEGs['genes'].unique()))

# Filter enrichment tables to only include these TFs
df_d7_temp = df_d7[df_d7['TF'].isin(tfs)]
df_d14_temp = df_d14[df_d14['TF'].isin(tfs)]

# Keep the top 5 enriched TFs per category
df_d7_temp = df_d7_temp.groupby('category').head(5).reset_index(drop=True)
df_d14_temp = df_d14_temp.groupby('category').head(5).reset_index(drop=True)

# Get TFs that are both upregulated DEGs and top5 enriched
TFs_d7 = df_d7_temp['TF'].unique()
TFs_d14 = df_d14_temp['TF'].unique()

# Filter enrichment tables again for these TFs
df_d7_temp = df_d7[df_d7['TF'].isin(set(list(TFs_d14) + list(TFs_d7)))]
df_d14_temp = df_d14[df_d14['TF'].isin(set(list(TFs_d14) + list(TFs_d7)))]

# Combine D7 and D14 enrichment tables
df = pd.concat([df_d7_temp, df_d14_temp], axis=0)

#-------------------------------------- Pivot to matrix for heatmap
matrix = df.pivot_table(index='category', columns='TF', values='Gene')
matrix[np.isnan(matrix)] = 0  # replace NaNs with 0

# Reorder rows and columns based on original category and TF order
d = pd.DataFrame(matrix).reindex(index=df['category'].unique(), columns=df['TF'].unique())

#-------------------------------------- Plot heatmap
sns.set_theme(style="darkgrid", font_scale=0.8)

g = sns.clustermap(
    d,
    mask=d==0,               # Hide zero enrcihments
    vmax=1,                  # Max value for color 
    center=-0.3,             # Center the color scale (so colors diverge around -0.3)
    z_score=1,               # Normalize each column (TF) by z-score
    cmap='coolwarm',       
    col_cluster=True,        
    row_cluster=False,      
    standard_scale=None,     
    method='ward',          
    metric='euclidean',      
    annot=False,            
    figsize=(7, 2),       
    annot_kws={"size": 2}    
)

# Save the heatmap as PDF
g.savefig('/Fig5_a_DEGs_enrichments.pdf')
