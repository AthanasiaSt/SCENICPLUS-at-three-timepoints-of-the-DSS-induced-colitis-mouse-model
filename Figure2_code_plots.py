#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 12 13:51:06 2025

@author: Athanasia Stavropoulou
"""
# code for main Figure2.g

# libraries 
import pandas as pd
import mudata
from scenicplus.plotting.dotplot import heatmap_dotplot
from plotnine import theme, element_text

# download processed files form zenodo: Processed_datasets_for_scRNA_scATAC_scenicplus
# -------------------- loading the scplusmdata.h5mu SCENIC+ main output 
#--------------------- and the filtered direct regulons that we used in our analysis 
scplus_mdata = mudata.read("./SCENICPLUS_Healthy_Inflammation_Regeneration_regulons_results/scplusmdata.h5mu")
filt_regulons=pd.read_csv('./SCENICPLUS_Healthy_Inflammation_Regeneration_regulons_results/Direct_Positive_eRegulons_filtered_basedOn_GeneBased_AUC.csv',index_col=0 , sep='\t')
filt_regulons=filt_regulons['Gene_signature_name'].unique()

len(filt_regulons) # this should be 43 

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
    .str.replace('D14', 'Regeneration', regex=False)
    .str.replace('Telocytes', 'SEMFs', regex=False)
    .str.replace('PDGFRalo', 'CD81-stroma', regex=False)
    .astype('category')
)

#figure 2.g ----- scenicplus plot
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
    group_variable_order = ['Healthy_Trophocytes','Healthy_CD81-stroma','Healthy_SEMFs','Inflammation_Trophocytes','Inflammation_CD81-stroma','Inflammation_SEMFs','Regeneration_Trophocytes','Regeneration_CD81-stroma','Regeneration_SEMFs'],
    figsize = (6, 10), 
    save=None,   # don’t save yet, adjust first
)

p = p + theme(axis_text_x=element_text(rotation=45, ha='right'))
p.save("./Fig2_g_Heatmap_scenicplus_results.pdf")

