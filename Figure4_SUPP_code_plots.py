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
import mudata
from scenicplus.plotting.dotplot import heatmap_dotplot
from plotnine import theme, element_text
import mudata
from scenicplus.scenicplus_class import mudata_to_scenicplus

# -------------------- Figure SUPP 4 

# -------------------- loading the scplusmdata.h5mu SCENIC+ main output // download from zenodo 
#--------------------- and the filtered direct regulons that we used in our analysis 
scplus_mdata = mudata.read("./SCENICPLUS_Healthy_Inflammation_Regeneration_regulons_results/scplusmdata.h5mu")
filt_regulons=pd.read_csv('./SCENICPLUS_Healthy_Inflammation_Regeneration_regulons_results/Direct_Positive_eRegulons_filtered_basedOn_GeneBased_AUC.csv',index_col=0 , sep='\t')
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

#change names 
scplus_mdata_tmp.obs['Celltype_Timepoint'] = (
    scplus_mdata_tmp.obs['Celltype_Timepoint']
    .astype(str)
    .str.replace('D0', 'Healthy', regex=False)
    .str.replace('D7', 'Inflammation', regex=False)
    .str.replace('D14', 'Regeneration', regex=False)
    .str.replace('Telocytes', 'SEMFs', regex=False)
    .str.replace('PDGFRalo', 'CD81-stroma', regex=False)
    .astype('category')
)


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
regulons_df=pd.read_csv('./SCENICPLUS_Healthy_Inflammation_Regeneration_regulons_results/Direct_Positive_eRegulons_filtered_basedOn_GeneBased_AUC.csv',index_col=0 , sep='\t')

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
plt.savefig("./Fig4_b_SUPP_Number_Enhancers_positive_regulons.pdf", dpi=500)  

#------------------------------ Fig4.e SUPP 
#Making SCENIC+ heatmaps for the two separate runs for the filtered regulons // not provided

# Kinchen / Healthy(K)-Inflammation(K)
scplus_mdata = mudata.read("./Revisions_1/scenicplus_revisions_kinchen/scplus_pipeline/Snakemake/scplusmdata.h5mu")
filt_regulons=pd.read_csv('./Revisions_1/scenicplus_revisions_kinchen/scplus_pipeline/Snakemake/output/direct_e_regulon_metadata_filtered_basedOn_max_mean_enrichment_across_categories_GeneBased.csv',index_col=0 , sep='\t')
filt_regulons=filt_regulons['Gene_signature_name'].unique()

len(filt_regulons) # this should be 31 
# ----------------------- keep the filtered regulons
scplus_mdata_tmp=scplus_mdata.copy()
scplus_mdata_tmp.uns['direct_e_regulon_metadata']=scplus_mdata_tmp.uns['direct_e_regulon_metadata'][scplus_mdata_tmp.uns['direct_e_regulon_metadata']['Gene_signature_name'].isin(filt_regulons)]

# make variable with celltype and timepoint from indices of metacells 
scplus_mdata_tmp.obs['Celltype_Timepoint']=scplus_mdata_tmp.obs.index.str.replace(r'.(\d*$)','').copy()

#change names 
scplus_mdata_tmp.obs['Celltype_Timepoint'] = (
    scplus_mdata_tmp.obs['Celltype_Timepoint']
    .astype(str)
    .str.replace('D0', 'Healthy', regex=False)
    .str.replace('D7', 'Inflammation', regex=False)
    .astype('category')
)

#----- scenicplus plot
# ------------------------- heatmap 
p=heatmap_dotplot(
    scplus_mudata = scplus_mdata_tmp,
    color_modality = "direct_gene_based_AUC",
    size_modality = "direct_region_based_AUC",
    group_variable = "Celltype_Timepoint",
    eRegulon_metadata_key = "direct_e_regulon_metadata",
    color_feature_key = "Gene_signature_name",
    size_feature_key = "Region_signature_name",
    feature_name_key = "eRegulon_name",
    sort_data_by = "direct_gene_based_AUC",
    orientation = "vertical",
    scale_size_matrix=True,
    scale_color_matrix=True,
    group_variable_order = ['Healthy_Trophocytes','Healthy_CD81_stroma','Healthy_SEMFs','Inflammation_Trophocytes','Inflammation_CD81_stroma','Inflammation_SEMFs'],
    figsize = (6, 10), 
    save=None,   # don’t save yet, adjust first
)

p = p + theme(axis_text_x=element_text(rotation=45, ha='right'))
p.save("./Fig4_supp_e_Heatmap_scenicplus_results_Kinchen.pdf")

# In-house / Healthy-Regeneration
# -------------------- loading the scplusmdata.h5mu SCENIC+ main output 
#--------------------- and the filtered direct regulons that we used in our analysis 
scplus_mdata = mudata.read("./Revisions_1/scenicplus_revisions_D0_D14/scplus_pipeline/Snakemake/scplusmdata.h5mu")
filt_regulons=pd.read_csv('./Revisions_1/scenicplus_revisions_D0_D14/scplus_pipeline/Snakemake/output/direct_e_regulon_metadata_filtered_basedOn_max_mean_enrichment_across_categories_GeneBased.csv',index_col=0 , sep='\t')
filt_regulons=filt_regulons['Gene_signature_name'].unique()

len(filt_regulons) # this should be 52 

# ----------------------- keep the filtered regulons
scplus_mdata_tmp=scplus_mdata.copy()
scplus_mdata_tmp.uns['direct_e_regulon_metadata']=scplus_mdata_tmp.uns['direct_e_regulon_metadata'][scplus_mdata_tmp.uns['direct_e_regulon_metadata']['Gene_signature_name'].isin(filt_regulons)]

# make variable with celltype and timepoint from indices of metacells 
scplus_mdata_tmp.obs['Celltype_Timepoint']=scplus_mdata_tmp.obs.index.str.replace(r'.(\d*$)','').copy()

#change names 
scplus_mdata_tmp.obs['Celltype_Timepoint'] = (
    scplus_mdata_tmp.obs['Celltype_Timepoint']
    .astype(str)
    .str.replace('D0', 'Healthy', regex=False)
    .str.replace('D14', 'Regeneration', regex=False)
    .astype('category')
)

#----- scenicplus plot
# ------------------------- heatmap 
p=heatmap_dotplot(
    scplus_mudata = scplus_mdata_tmp,
    color_modality = "direct_gene_based_AUC",
    size_modality = "direct_region_based_AUC",
    group_variable = "Celltype_Timepoint",
    eRegulon_metadata_key = "direct_e_regulon_metadata",
    color_feature_key = "Gene_signature_name",
    size_feature_key = "Region_signature_name",
    feature_name_key = "eRegulon_name",
    sort_data_by = "direct_gene_based_AUC",
    orientation = "vertical",
    scale_size_matrix=True,
    scale_color_matrix=True,
    group_variable_order = ['Healthy_Trophocytes','Healthy_CD81_stroma','Healthy_SEMFs','Regeneration_Trophocytes','Regeneration_CD81_stroma','Regeneration_SEMFs'],
    figsize = (6, 10), 
    save=None,   # don’t save yet, adjust first
)

p = p + theme(axis_text_x=element_text(rotation=45, ha='right'))
p.save("./Fig4_supp_e_Heatmap_scenicplus_results_In_house.pdf")

# Fig 4 supp g)
#Load one SCENIC+ run and return the max of average gene-based AUC scores per Positive-positive regulon.


def load_and_compute_max_gene_auc(scplus_mdata_path, ctx_h5_path, dem_h5_path):
  # ---- Load mudata
  scplus_mdata = mudata.read(scplus_mdata_path)
  # ---- Filter regulon metadata
  meta = scplus_mdata.uns['direct_e_regulon_metadata']
  meta = meta[meta['is_extended'] == False]
  meta = meta[meta['regulation'] == 1]
  meta = meta[meta['rho_R2G'] > 0]
  scplus_mdata.uns['direct_e_regulon_metadata'] = meta
  # ---- Convert to scenicplus object
  scplus_obj = mudata_to_scenicplus(mdata=scplus_mdata, path_to_cistarget_h5=ctx_h5_path, path_to_dem_h5=dem_h5_path)
  # ---- Get regulon / region names
  pos_regulons_names = meta['Gene_signature_name'].unique()
  pos_regions_names = meta['Region_signature_name'].unique()
  # ---- Subset AUC matrices
  scplus_obj.uns['eRegulon_AUC']['Gene_based'] =  scplus_obj.uns['eRegulon_AUC']['Gene_based'][pos_regulons_names]
  scplus_obj.uns['eRegulon_AUC']['Region_based'] = scplus_obj.uns['eRegulon_AUC']['Region_based'][pos_regions_names]
  # ---- Add celltype_timepoint
  gene_auc = scplus_obj.uns['eRegulon_AUC']['Gene_based'].copy()
  gene_auc['celltype_timepoint'] = gene_auc.index.str.replace(r'.(\d*$)', '', regex=True)
  # ---- Mean per category
  matrix_gene_based = (gene_auc.groupby(['celltype_timepoint']).mean()[pos_regulons_names])
  # ---- Max of average AUC per regulon
  max_gene_auc = matrix_gene_based.max()
  return max_gene_auc

D0_D7 = load_and_compute_max_gene_auc(
    scplus_mdata_path="/home/astavropoulou/D0_D7_D14_PAPER_FINAL_FILES_RUNS/Processed_datasets_for_scRNA_scATAC_scenicplus/Revisions_1/scenicplus_revisions_kinchen/scplus_pipeline/Snakemake/scplusmdata.h5mu",
    ctx_h5_path="/home/astavropoulou/D0_D7_D14_PAPER_FINAL_FILES_RUNS/Processed_datasets_for_scRNA_scATAC_scenicplus/Revisions_1/scenicplus_revisions_kinchen/scplus_pipeline/Snakemake/ctx_results.hdf5",
    dem_h5_path="/home/astavropoulou/D0_D7_D14_PAPER_FINAL_FILES_RUNS/Processed_datasets_for_scRNA_scATAC_scenicplus/Revisions_1/scenicplus_revisions_kinchen/scplus_pipeline/Snakemake/dem_results.hdf5"
)

D0_D14 = load_and_compute_max_gene_auc(
    scplus_mdata_path="/home/astavropoulou/D0_D7_D14_PAPER_FINAL_FILES_RUNS/Processed_datasets_for_scRNA_scATAC_scenicplus/Revisions_1/scenicplus_revisions_D0_D14/scplus_pipeline/Snakemake/scplusmdata.h5mu",
    ctx_h5_path="/home/astavropoulou/D0_D7_D14_PAPER_FINAL_FILES_RUNS/Processed_datasets_for_scRNA_scATAC_scenicplus/Revisions_1/scenicplus_revisions_D0_D14/scplus_pipeline/Snakemake/ctx_results.hdf5",
    dem_h5_path="/home/astavropoulou/D0_D7_D14_PAPER_FINAL_FILES_RUNS/Processed_datasets_for_scRNA_scATAC_scenicplus/Revisions_1/scenicplus_revisions_D0_D14/scplus_pipeline/Snakemake/dem_results.hdf5"
)

D0_D7_D14 = load_and_compute_max_gene_auc(
    scplus_mdata_path="/home/astavropoulou/D0_D7_D14_PAPER_FINAL_FILES_RUNS/scenicplus_output_from_my_PC_the_run_in_server/scplusmdata.h5mu",
    ctx_h5_path="/home/astavropoulou/D0_D7_D14_PAPER_FINAL_FILES_RUNS/scenicplus_output_from_my_PC_the_run_in_server/ctx_results.hdf5",
    dem_h5_path="/home/astavropoulou/D0_D7_D14_PAPER_FINAL_FILES_RUNS/scenicplus_output_from_my_PC_the_run_in_server/dem_results.hdf5"
)

# Plot
plt.figure(figsize=(6,4))

sns.kdeplot(D0_D7, clip=(0, None), linewidth=2, label="D0_D7")
sns.kdeplot(D0_D14, clip=(0, None), linewidth=2, label="D0_D14")
sns.kdeplot(D0_D7_D14, clip=(0, None), linewidth=2, label="D0_D7_D14")

plt.xlabel("Max of average gene-based AUC scores across categories")
plt.ylabel("Density")
plt.legend()
plt.tight_layout()

plt.savefig(
    "Fig_4_supp_g_Max_geneBased_AUC_score_across_categories_KDE_3runs.pdf",
    dpi=1000
)
plt.show()


#compare distributions with Kolmogorov-Smirnov test
from scipy.stats import ks_2samp

ks_12 = ks_2samp(D0_D7, D0_D14)
ks_13 = ks_2samp(D0_D7, D0_D7_D14)
ks_23 = ks_2samp(D0_D14, D0_D7_D14)

print("Run1 vs Run2:", ks_12)
print("Run1 vs Run3:", ks_13)
print("Run2 vs Run3:", ks_23)

