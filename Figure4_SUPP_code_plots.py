#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 19 17:03:35 2025

@author: Athanasia Stavropoulou
"""

# libraries 
import pandas as pd
import mudata
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

# -------------------- Figure SUPP 4 

# -------------------- loading the scplusmdata.h5mu SCENIC+ main output // download from zenodo 
#--------------------- and the filtered direct regulons that we used in our analysis 
scplus_mdata = mudata.read("./SCENICPLUS_D0_D7_D14_regulons_results/scplusmdata.h5mu")
filt_regulons=pd.read_csv('./SCENICPLUS_D0_D7_D14_regulons_results/Direct_Positive_eRegulons_filtered_basedOn_GeneBased_AUC.csv',index_col=0 , sep='\t')
filt_regulons=filt_regulons['Gene_signature_name'].unique()

# ----------------------- keep the filtered regulons
scplus_mdata_tmp=scplus_mdata.copy()
scplus_mdata_tmp.uns['direct_e_regulon_metadata']=scplus_mdata_tmp.uns['direct_e_regulon_metadata'][scplus_mdata_tmp.uns['direct_e_regulon_metadata']['Gene_signature_name'].isin(filt_regulons)]

# make variable with celltype and timepoint from indices of metacells 
scplus_mdata_tmp.obs['Celltype_Timepoint']=scplus_mdata_tmp.obs.index.str.replace(r'.(\d*$)','').copy()

# ---------------------------- Gene or Region based enrichments
eRegulon_gene_AUC = scplus_mdata_tmp["direct_gene_based_AUC"]
eRegulon_region_AUC = scplus_mdata_tmp["direct_region_based_AUC"]

# ------------------------------ Keep filtered regulons 
#--------------- Gene based
regs_to_keep = scplus_mdata_tmp.uns['direct_e_regulon_metadata']['Gene_signature_name'].unique()
eRegulon_gene_AUC = eRegulon_gene_AUC[:, eRegulon_gene_AUC.var_names.isin(regs_to_keep)].copy()

eRegulon_gene_AUC.obs = scplus_mdata_tmp.obs.loc[eRegulon_gene_AUC.obs_names]
eRegulon_gene_AUC.obs['Celltype_Timepoint']=eRegulon_gene_AUC.obs.index.str.replace(r'.(\d*$)','').copy()


#--------------- Region based
regs_to_keep = scplus_mdata_tmp.uns['direct_e_regulon_metadata']['Region_signature_name'].unique()
eRegulon_region_AUC = eRegulon_region_AUC[:, eRegulon_region_AUC.var_names.isin(regs_to_keep)].copy()

eRegulon_region_AUC.obs = scplus_mdata_tmp.obs.loc[eRegulon_region_AUC.obs_names]
eRegulon_region_AUC.obs['Celltype_Timepoint']=eRegulon_region_AUC.obs.index.str.replace(r'.(\d*$)','').copy()


# set defaults for all plots
sc.settings.set_figure_params(
    dpi=300,          # resolution
    dpi_save=300,     # resolution when saving
    figsize=(5, 5),   # default figure size
    facecolor='white' # background
)

sc.tl.pca(eRegulon_gene_AUC)
sc.pp.neighbors(eRegulon_gene_AUC, use_rep="X_pca")
sc.tl.umap(eRegulon_gene_AUC)

sc.pl.umap(
    eRegulon_gene_AUC,
    color="Celltype_Timepoint",
    size=30,              # dot size (default ~20)
    #frameon=False,        # remove axis frame
    save="_Fig4_SUPP_a_gene_based.pdf"
)

sc.tl.pca(eRegulon_region_AUC)
sc.pp.neighbors(eRegulon_region_AUC, use_rep="X_pca")
sc.tl.umap(eRegulon_region_AUC)

sc.pl.umap(
    eRegulon_region_AUC,
    color="Celltype_Timepoint",
    size=30,              # dot size (default ~20)
    #frameon=False,        # remove axis frame
    save="_Fig4_SUPP_a_region_based.pdf"
)


#------------------------------ Fig4.b SUPP 
# --------------------------plot numbers of genes and enhancers per regulon after filtering
regulons_df=pd.read_csv('./SCENICPLUS_D0_D7_D14_regulons_results/Direct_Positive_eRegulons_filtered_basedOn_GeneBased_AUC.csv',index_col=0 , sep='\t')

# number of genes and enhancers -- keep combinations once
genes=regulons_df.groupby(['TF','Gene']).size().reset_index().rename(columns={0:'count'}).groupby('TF').size()
regions=regulons_df.groupby(['TF','Region']).size().reset_index().rename(columns={0:'count'}).groupby('TF').size()

df=pd.concat([genes,regions],  axis=1)
df.columns = ['nGenes', 'nRegions']
df=df.sort_values('nGenes', ascending=False)

df_long = df.reset_index().melt(id_vars='TF', var_name='Type', value_name='Count')

# barplots
sns.set_theme(style="whitegrid")

plt.figure(figsize=(9, 6))
sns.barplot(data=df_long, x='TF', y='Count', hue='Type')
plt.xticks(rotation=45, ha='right')
plt.title('Number of Genes and Enhancers per TF-regulon')
plt.tight_layout()

# Save as PDF
plt.savefig("Fig4_b_SUPP_Number_Enhancers_positive_regulons.pdf", dpi=500)  