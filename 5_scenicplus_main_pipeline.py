#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 13:08:19 2024

@author: Athanasia Stavropoulou
"""

#5. Main pipeline for using SCENIC+ and filtering regulons 
# see tutorials: https://scenicplus.readthedocs.io/en/latest/tutorials.html

#------------------- Libraries and setup -------------------
import warnings
import sys
import os
import scenicplus
import pycisTopic
import scanpy as sc
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import mudata
import scanpy as sc
from matplotlib import pyplot as plt
from scenicplus.RSS import (regulon_specificity_scores, plot_rss)
from scenicplus.plotting.dotplot import heatmap_dotplot
from scenicplus.scenicplus_class import mudata_to_scenicplus
import numpy as np

# Working and temporary directories
work_dir = '/home/astavropoulou/D0_D7_D14_PAPER_FINAL_FILES_RUNS/new_scenicplus/'
tmp_dir = '/home/astavropoulou/D0_D7_D14_PAPER_FINAL_FILES_RUNS/new_scenicplus/temp/'

if not os.path.exists(work_dir):
    os.makedirs(work_dir)

#-------------------------------------------------- Processing scRNA ---------------------------------------------
# Load scRNA-seq data that was exported from Seurat to AnnData format
adata = sc.read_h5ad("./scRNA/seurat5_Harmony_integration_all_samples_treated_together_D7_D14_Subclustering.h5ad")

# keep all data or specific samples
#adata=adata[adata.obs['orig.ident'].isin(['D0_kinchen','D7_kinchen'])]
#adata=adata[adata.obs['orig.ident'].isin(['D0','D14'])]

# Peek at raw counts (sparse matrix) to verify data integrity
adata.X[20:30,20:30].toarray()

# ---------------- Quality Control ----------------
# Basic filtering by genes and cells
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# QC metrics: add mitochondrial percentage, total counts, number of genes, etc.
adata.var['mt'] = adata.var_names.str.startswith('mt-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Plot QC metrics and overlay thresholds (variables mito_filter and n_counts_filter should be defined earlier!)
fig, axs = plt.subplots(ncols = 2, figsize = (8,4))
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', ax = axs[0], show=False)
sc.pl.scatter(adata, x='total_counts', y='percent.mt', ax = axs[0], show=False)  # possible redundancy: 'pct_counts_mt' vs 'percent.mt'
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', ax = axs[1], show = False)
axs[0].hlines(y = mito_filter, xmin = 0, xmax = max(adata.obs['total_counts']), color = 'red', ls = 'dashed')
axs[1].hlines(y = n_counts_filter, xmin = 0, xmax = max(adata.obs['total_counts']), color = 'red', ls = 'dashed')
fig.tight_layout()
fig.savefig('/home/astavropoulou/D0_D7_D14_PAPER_FINAL_FILES_RUNS/new_scenicplus/QC_metrics.png', dpi=1200)

# ---------------- Cleanup ----------------
# Remove precomputed reductions/layers from Seurat export that can cause conflicts
del adata.obsm['PCA']
del adata.obsm['UMAP']
del adata.obsm['THETA2_SIGMA0.1_COMMONVARGENES_ONGROUP']
del adata.obsm['THETA2_SIGMA0.1_COMMONVARGENES_ONORIGIDENT']
del adata.obsm['THETA2_SIGMA0.1_COMMONVARGENES_ONORIG_IDENT']
del adata.obsm['UMAP.UNINTEGRATED']
del adata.layers['logcounts']

# ---------------- Normalization ----------------
# Keep raw counts in adata.raw for SCENIC+ later
adata.raw = adata.copy()

# Normalize counts and log-transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Identify highly variable genes (HVGs)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# Save full dataset (no HVG filtering yet)
adata.write(os.path.join(work_dir, 'scRNA/adata_noVariableGenes.h5ad'), compression='gzip')

# ---------------- Exploration ----------------
# adata = adata[:, adata.var.highly_variable]  # subset to HVGs
# sc.pp.scale(adata, max_value=10)
# sc.pp.pca(adata)
# sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50, save='') 
# 
# # Harmony integration on 'group' metadata
# sc.external.pp.harmony_integrate(adata, 'group')
# adata.obsm['X_pca'] = adata.obsm['X_pca_harmony']  # replace PCA with Harmony output
# 
# # Neighbors, UMAP, clustering
# sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)
# sc.tl.umap(adata)
# sc.tl.leiden(adata, resolution=0.5)
# 
# # Plot UMAP by donor/sample
# donor_umap = sc.pl.umap(adata, color=['orig.ident'],
#     show=False, palette=sns.color_palette("husl", 24),
#     legend_fontsize=6, frameon=True, title='Donor')
# fig = donor_umap.get_figure()
# fig.set_size_inches(5, 5)
# fig.savefig(str(sc.settings.figdir) + '/umap_lgd_harmony_sample', dpi=400, bbox_inches='tight')
# 
# # Plot UMAP by celltypes
# leiden_umap = sc.pl.umap(adata, color=['celltype_stim_2'],
#     show=False, palette=sns.color_palette("husl", 24),
#     legend_fontsize=6, frameon=True, title='celltypes')
# fig = leiden_umap.get_figure()
# fig.set_size_inches(5, 5)
# fig.savefig(str(sc.settings.figdir) + '/umap_lgd_harmony_leiden', dpi=400, bbox_inches='tight')

# ---------------- Merge modalities ----------------
# Reload saved AnnData with full genes (no HVG filtering)
adata = sc.read_h5ad("/home/astavropoulou/scenicplGEX:annotation_scenicPlusus_D7_D14_D0_fibroblasts/scRNA/adata_noVariableGenes.h5ad")

# Sanity check again on raw counts
adata.raw.X[20:30,20:30].toarray()
adata.X[20:30,20:30].toarray()

# Adjust metadata for SCENIC+ integration
# merge D0 and D0_kinchen --> D0
adata.obs['GEX:annotation_scenicPlus'] = adata.obs['celltype_stim_2'].copy()

adata.obs['GEX:annotation_scenicPlus']=adata.obs['GEX:annotation_scenicPlus'].str.replace(r'D7_kinchen_Trophocytes','D7_Trophocytes').copy()
adata.obs['GEX:annotation_scenicPlus']=adata.obs['GEX:annotation_scenicPlus'].str.replace(r'D7_kinchen_PDGFRalo','D7_CD81_stroma').copy()
adata.obs['GEX:annotation_scenicPlus']=adata.obs['GEX:annotation_scenicPlus'].str.replace(r'D7_kinchen_Telocytes','D7_SEMFs').copy()
adata.obs['GEX:annotation_scenicPlus']=adata.obs['GEX:annotation_scenicPlus'].str.replace(r'D0_kinchen_Telocytes','D0_SEMFs').copy()
adata.obs['GEX:annotation_scenicPlus']=adata.obs['GEX:annotation_scenicPlus'].str.replace(r'D0_kinchen_PDGFRalo','D0_CD81_stroma').copy()
adata.obs['GEX:annotation_scenicPlus']=adata.obs['GEX:annotation_scenicPlus'].str.replace(r'D0_kinchen_Trophocytes','D0_Trophocytes').copy()

adata.obs['GEX:annotation_scenicPlus']=adata.obs['GEX:annotation_scenicPlus'].str.replace(r'D0_Trophocytes','D0_Trophocytes').copy()
adata.obs['GEX:annotation_scenicPlus']=adata.obs['GEX:annotation_scenicPlus'].str.replace(r'D0_PDGFRalo','D0_CD81_stroma').copy()
adata.obs['GEX:annotation_scenicPlus']=adata.obs['GEX:annotation_scenicPlus'].str.replace(r'D0_Telocytes','D0_SEMFs').copy()
adata.obs['GEX:annotation_scenicPlus']=adata.obs['GEX:annotation_scenicPlus'].str.replace(r'D14_Telocytes','D14_SEMFs').copy()
adata.obs['GEX:annotation_scenicPlus']=adata.obs['GEX:annotation_scenicPlus'].str.replace(r'D14_PDGFRalo','D14_CD81_stroma').copy()
adata.obs['GEX:annotation_scenicPlus']=adata.obs['GEX:annotation_scenicPlus'].str.replace(r'D14_Trophocytes','D14_Trophocytes').copy()

adata.write(os.path.join(work_dir, 'scRNA/adata_to_scenicplus.h5ad'), compression='gzip')


#-------------------------------------------------------------------------------------

#---------------------- scATAC --------------------------------
# Load the metadata for scATAC cells.
# This was exported from ArchR after cell annotation (mesenchymal-only subset).
celldata = pd.read_csv(
    '/home/astavropoulou/scenicplus_D7_D14_D0_fibroblasts/scATAC/Cnt_Day7_Day14_ATAC_analysis_mesenchymal.csv',
    sep=','
)

# Clean up cell names by removing ArchR prefixes (cnt#, D7#, D14#) from the cell IDs
celldata['Cells'] = celldata['Cells'].str.replace(r'cnt#','')
celldata['Cells'] = celldata['Cells'].str.replace(r'D7#','')
celldata['Cells'] = celldata['Cells'].str.replace(r'D14#','')

# Create/ensure useful columns for downstream use
celldata['Sample'] = celldata['Sample'].astype(str)    # make sure 'Sample' is string
celldata['barcode'] = celldata['Cells'].astype(str)    # store cell IDs as barcodes

# Add a SCENIC+ annotation column based on ArchR integration metadata
celldata['GEX:annotation_scenicPlus'] = celldata['celltype.stim'].copy  

# Quick check: number of cells per sample × annotation group
celldata.groupby(['Sample','GEX:annotation_scenicPlus']).size()

# ------------------------ Starting the procedure based on SCENIC+ tutorial ------------------------

# ---------------------------- scATAC fragments locations  
# Dictionary linking sample IDs to their corresponding fragment files (from Cell Ranger output)
fragments_dict = {
    'cnt': '/home/astavropoulou/cellranger_results/Control_colon_scATAC/outs/fragments.tsv.gz',
    'D14':'/home/astavropoulou/cellranger_results/scATAC_D14_DSS/outs/fragments.tsv.gz',
    'D7':'/home/astavropoulou/cellranger_results/scATAC_D7_DSS/outs/fragments.tsv.gz'
}

# ---------------------------- chromosome sizes 
# Load mm10 genome chromosome sizes from UCSC, required for defining peak regions
chromsizes = pd.read_table(
    "http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes",
    header = None,
    names = ["Chromosome", "End"]
)
chromsizes.insert(1, "Start", 0)  # add start column (always 0)
chromsizes.head()

# ------------------------------------ Create pseudobulk tracks (per condition and annotation)
# Export pseudobulk BED and bigWig files, grouping by SCENIC+ annotation and sample
from pycisTopic.pseudobulk_peak_calling import export_pseudobulk

# Create directories for pseudobulk outputs
os.makedirs(os.path.join(work_dir, "consensus_peak_calling"), exist_ok = True)
os.makedirs(os.path.join(work_dir, "consensus_peak_calling/pseudobulk_bed_files"), exist_ok = True)
os.makedirs(os.path.join(work_dir, "consensus_peak_calling/pseudobulk_bw_files"), exist_ok = True)

# Run pseudobulk export (generates pooled signal per annotation x sample)
bw_paths, bed_paths = export_pseudobulk(
    input_data = celldata,
    variable = "GEX:annotation_scenicPlus",  # annotation column to group cells
    sample_id_col = "Sample",                # column containing sample IDs
    chromsizes = chromsizes,
    bed_path = os.path.join(work_dir, "consensus_peak_calling/pseudobulk_bed_files"),
    bigwig_path = os.path.join(work_dir, "consensus_peak_calling/pseudobulk_bw_files"),
    path_to_fragments = fragments_dict,
    n_cpu = 10,
    normalize_bigwig = True,
    temp_dir = "/tmp",
    split_pattern = "-"                      # pattern to split cell IDs if needed
)

# Save pseudobulk bigWig and BED file paths for later reference
with open(os.path.join(work_dir, "consensus_peak_calling/bw_paths.tsv"), "wt") as f:
    for v in bw_paths:
        _ = f.write(f"{v}\t{bw_paths[v]}\n")
        
with open(os.path.join(work_dir, "consensus_peak_calling/bed_paths.tsv"), "wt") as f:
    for v in bed_paths:
        _ = f.write(f"{v}\t{bed_paths[v]}\n")

# ---------------------- Peak calling ----------------------
# Use MACS2 on pseudobulk BED files to identify narrow peaks
from pycisTopic.pseudobulk_peak_calling import peak_calling
macs_path = "macs2"  # path to MACS2 binary

os.makedirs(os.path.join(work_dir, "consensus_peak_calling/MACS"), exist_ok = True)

narrow_peak_dict = peak_calling(
    macs_path = macs_path,
    bed_paths = bed_paths,
    outdir = os.path.join(os.path.join(work_dir, "consensus_peak_calling/MACS")),
    genome_size = 'mm',          # genome size for mouse
    n_cpu = 10,
    input_format = 'BEDPE',
    shift = 73,                  # ATAC-seq shift (recommended)
    ext_size = 146,              # fragment length
    keep_dup = 'all',
    q_value = 0.05,              # MACS2 FDR threshold
    _temp_dir = "/tmp"
    #ignore_reinit_error=True
)

# ---------------------- Consensus peaks ----------------------
from pycisTopic.iterative_peak_calling import get_consensus_peaks

# Parameters
peak_half_width = 250
path_to_blacklist = "/home/astavropoulou/pycisTopic/blacklist/mm10-blacklist.v2.bed"

# Merge peaks across pseudobulks into a consensus set
consensus_peaks = get_consensus_peaks(
    narrow_peaks_dict = narrow_peak_dict,
    peak_half_width = peak_half_width,
    chromsizes = chromsizes,
    path_to_blacklist = path_to_blacklist
)

# Save consensus peaks to BED
consensus_peaks.to_bed(
    path = os.path.join(work_dir, "consensus_peak_calling/consensus_regions.bed"),
    keep = True,
    compression = 'infer',
    chain = False
)

# --------------------- QC --------------------------
# In the tutorial, they filter barcodes based on QC metrics. 
# Since QC was already performed earlier, here we only select barcodes that passed filters from celldata.

sample_id_to_barcodes_passing_filters = {'cnt':(), 'D14':(), 'D7':()}

# Assign barcodes from celldata by sample
sample_id_to_barcodes_passing_filters['cnt'] = celldata[celldata['Sample'] == 'cnt']['barcode'].to_list()
sample_id_to_barcodes_passing_filters['D14'] = celldata[celldata['Sample'] == 'D14']['barcode'].to_list()
sample_id_to_barcodes_passing_filters['D7'] = celldata[celldata['Sample'] == 'D7']['barcode'].to_list()

# Save passing barcodes
#import pickle
#pickle.dump(sample_id_to_barcodes_passing_filters, '/home/astavropoulou/scenicplus_try3_onlyfibroblasts/outs/qc/bc_passing_filters.pkl', 'wb')

# Save final celldata object with cleaned barcodes and annotations
celldata.to_csv('celldata_ATAC_final.csv', sep='\t')

# do the QC metric computation by using pycistopic
# BASH 
#pycistopic qc -f /home/astavropoulou/Cellranger_Results_various_datasets/Control_colon_scATAC/outs/fragments.tsv.gz -r /home/astavropoulou/D0_D7_D14_PAPER_FINAL_FILES_RUNS/Processed_datasets_for_scRNA_scATAC_scenicplus/Revisions_1/scenicplus_revisions_D0_D14/consensus_peak_calling/consensus_regions.bed -t /home/astavropoulou/D0_D7_D14_PAPER_FINAL_FILES_RUNS/Processed_datasets_for_scRNA_scATAC_scenicplus/Revisions_1/scenicplus_revisions_D0_D14/outs/qc/tss.bed -o /home/astavropoulou/D0_D7_D14_PAPER_FINAL_FILES_RUNS/Processed_datasets_for_scRNA_scATAC_scenicplus/Revisions_1/scenicplus_revisions_D0_D14/outs/qc/cnt

# Define paths
path_to_regions = os.path.join(work_dir, "consensus_peak_calling/consensus_regions.bed")  # consensus peak regions
path_to_blacklist = "/home/astavropoulou/pycisTopic/blacklist/mm10-blacklist.v2.bed"     # blacklist regions
pycistopic_qc_output_dir = "outs/qc"                                                     # QC output directory

# ---------------------- Create cisTopic objects per sample ----------------------
from pycisTopic.cistopic_class import create_cistopic_object_from_fragments
from pycisTopic.cistopic_class import *
import polars as pl

cistopic_obj_list = []
for sample_id in fragments_dict:
    # Load precomputed fragment QC stats for each sample (from parquet file)
    sample_metrics = pl.read_parquet(
        os.path.join(pycistopic_qc_output_dir, f'{sample_id}.fragments_stats_per_cb.parquet')
    ).to_pandas().set_index("CB").loc[sample_id_to_barcodes_passing_filters[sample_id]]
    
    # Create cisTopic object for this sample (fragments + peaks + blacklist + QC metrics)
    cistopic_obj = create_cistopic_object_from_fragments(
        path_to_fragments = fragments_dict[sample_id],
        path_to_regions = path_to_regions,
        path_to_blacklist = path_to_blacklist,
        metrics = sample_metrics,
        valid_bc = sample_id_to_barcodes_passing_filters[sample_id],
        n_cpu = 1,
        project = sample_id,
        split_pattern = '-'  # split barcode if needed
    )
    cistopic_obj_list.append(cistopic_obj)
    
# Merge all sample-level cisTopic objects into a single object
cistopic_obj = merge(cistopic_obj_list)

# Save merged cisTopic object
import pickle
pickle.dump(
    cistopic_obj,
    open(os.path.join(work_dir, "cistopic_obj.pkl"), "wb")
)

# ---------------------- Add metadata from celldata ----------------------
# Adjust index in celldata so it matches cisTopic barcodes (Cell-Sample___Sample format)
celldata['new_index'] = celldata['Cells'] + '-' + celldata['Sample']
celldata['new_index'] = celldata['new_index'] + '___' + celldata['Sample']
celldata = celldata.set_index('new_index')

# Add updated metadata to cisTopic object
cistopic_obj.add_cell_data(celldata, split_pattern='-')

# ---------------------- Doublet detection with Scrublet ----------------------
import scrublet as scr

# Initialize scrublet on fragment matrix (transposed, cells as rows)
scrub = scr.Scrublet(cistopic_obj.fragment_matrix.T, expected_doublet_rate=0.1)

# Run doublet prediction
doublet_scores, predicted_doublets = scrub.scrub_doublets()

# Plot distribution of doublet scores (before thresholding)
scrub.plot_histogram()
plt.savefig('scrublet_1.png')

# Call doublets with explicit threshold (0.4) and replot
scrub.call_doublets(threshold=0.4)
scrub.plot_histogram()
plt.savefig('scrublet_2.png')

# Store Scrublet results in dataframe (scores + predictions per cell)
scrublet = pd.DataFrame(
    [scrub.doublet_scores_obs_, scrub.predicted_doublets_],
    columns = cistopic_obj.cell_names,
    index = ['Doublet_scores_fragments', 'Predicted_doublets_fragments']
).T

# Add Scrublet results as metadata to cisTopic object
cistopic_obj.add_cell_data(scrublet, split_pattern='-')

# Count number of doublets detected
sum(cistopic_obj.cell_data.Predicted_doublets_fragments == True)

# ---------------------- Save object before removing doublets ----------------------
pickle.dump(
    cistopic_obj,
    open(os.path.join(work_dir, "cistopic_obj.pkl"), "wb")
)

# ---------------------- Remove doublets ----------------------
# Extract list of singlet cell barcodes (exclude predicted doublets)
singlets = cistopic_obj.cell_data[
    cistopic_obj.cell_data.Predicted_doublets_fragments == False
].index.tolist()

# Subset cisTopic object to singlets only (creates a new object)
cistopic_obj_noDBL = cistopic_obj.subset(singlets, copy=True, split_pattern='-')
print(cistopic_obj_noDBL)

# Save final cisTopic object without doublets
pickle.dump(
    cistopic_obj_noDBL,
    open(os.path.join(work_dir, "cistopic_obj_noDB.pkl"), "wb")
)

#-------------------------------------------------------------------------
# choose the number of topics 

os.environ['MALLET_MEMORY'] = '200G'
from pycisTopic.lda_models import run_cgs_models_mallet
# Configure path Mallet
mallet_path="/home/astavropoulou/Mallet-202108-bin/Mallet-202108/bin/mallet"
# Run models
models=run_cgs_models_mallet(
    cistopic_obj_noDBL,
    n_topics=[2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50],
    n_cpu=12,
    n_iter=500,
    random_state=555,
    alpha=50,
    alpha_by_topic=True,
    eta=0.1,
    eta_by_topic=False,
    tmp_path=os.path.join(work_dir, "temp/mallet"), # make this in temp file 
    save_path=os.path.join(work_dir, "mallet"),
    mallet_path=mallet_path,
)

pickle.dump(
    models,
    open(os.path.join(work_dir, "models.pkl"), "wb")
)


#------------------------------------model selection 
from pycisTopic.lda_models import evaluate_models
model = evaluate_models(
    models,
    select_model = 40,
    return_model = True,
    figsize=(15,7)
)

plt.savefig("model_selection.png",dpi=400)

cistopic_obj_noDBL.add_LDA_model(model)
pickle.dump(
    cistopic_obj_noDBL,
    open(os.path.join(work_dir, "cistopic_obj_noDB.pkl"), "wb")
)

#---------------------- Data exploration ---------------------------------
cistopic_obj_noDBL = pickle.load(open(os.path.join('/home/astavropoulou/SCENICPLUS_runs/scenicplus_D7_D14_D0_fibroblasts', 'cistopic_obj_noDB.pkl'), 'rb'))

#-------------------------------cluster the cells and visualization 
from pycisTopic.clust_vis import (
    find_clusters,
    run_umap,
    run_tsne,
    plot_metadata,
    plot_topic,
    cell_topic_heatmap
)

find_clusters(
    cistopic_obj_noDBL,
    target  = 'cell',
    k = 10,
    res = [0.6, 1.2, 3],
    prefix = 'pycisTopic_',
    scale = True,
    split_pattern = '-'
)

run_umap(
    cistopic_obj_noDBL,
    target  = 'cell', scale=True)

run_tsne(
    cistopic_obj_noDBL,
    target  = 'cell', scale=True)

plot_metadata(
    cistopic_obj_noDBL,
    reduction_name='UMAP',
    variables=['GEX:annotation_scenicPlus','Sample', 'pycisTopic_leiden_10_0.6', 'pycisTopic_leiden_10_1.2', 'pycisTopic_leiden_10_3'],
    target='cell', num_columns=5,
    text_size=10,
    dot_size=5)

plt.savefig('umap.png')

plot_metadata(
    cistopic_obj_noDBL,
    reduction_name='tSNE',
    variables=['log10_unique_fragments_count', 'tss_enrichment', 'Doublet_scores_fragments', 'fraction_of_fragments_in_peaks'],
    target='cell', num_columns=4,
    text_size=10,
    dot_size=5)

plt.savefig('continuous_variables.png')

plot_topic(
    cistopic_obj_noDBL,
    reduction_name = 'UMAP',
    target = 'cell',
    num_columns=5
)

plt.savefig('cistopic_contribution.png')

cell_topic_heatmap(
    cistopic_obj_noDBL,
    variables = ['GEX:annotation_scenicPlus'],
    scale = False,
    legend_loc_x = 1.0,
    legend_loc_y = -1.2,
    legend_dist_y = -1,
    figsize = (10, 10)
)

plt.savefig('topic_heatmap.png')

# ------------------------------------------ topic binirization 
from pycisTopic.topic_binarization import binarize_topics
region_bin_topics_top_3k = binarize_topics(
    cistopic_obj_noDBL, method='ntop', ntop = 3_000,
    plot=True, num_columns=5
)
plt.savefig('binarization.png')

region_bin_topics_otsu = binarize_topics(
    cistopic_obj_noDBL, method='otsu',
    plot=True, num_columns=5
)
plt.savefig('binarization_otsu.png')

binarized_cell_topic = binarize_topics(
    cistopic_obj_noDBL,
    target='cell',
    method='li',
    plot=True,
    num_columns=5, nbins=100)

plt.savefig('binarization_li.png')

from pycisTopic.topic_qc import compute_topic_metrics, plot_topic_qc, topic_annotation
import matplotlib.pyplot as plt
from pycisTopic.utils import fig2img

topic_qc_metrics = compute_topic_metrics(cistopic_obj_noDBL)

fig_dict={}
fig_dict['CoherenceVSAssignments']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Log10_Assignments', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['AssignmentsVSCells_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Log10_Assignments', var_y='Cells_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSCells_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Cells_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSRegions_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Regions_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSMarginal_dist']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Marginal_topic_dist', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSGini_index']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Gini_index', var_color='Gini_index', plot=False, return_fig=True)

# Plot topic stats in one figure
fig=plt.figure(figsize=(40, 43))
i = 1
for fig_ in fig_dict.keys():
    plt.subplot(2, 3, i)
    img=fig2img(fig_dict[fig_]) #To convert figures to png to plot together, see .utils.py. This converts the figure to png.
    plt.imshow(img)
    plt.axis('off')
    i += 1
    plt.savefig('binarization_metrics'+fig_+'.png')

plt.subplots_adjust(wspace=0, hspace=-0.70)
plt.show()
plt.savefig('binarization_metrics.png')

topic_annot = topic_annotation(
    cistopic_obj_noDBL,
    annot_var='GEX:annotation_scenicPlus',
    binarized_cell_topic=binarized_cell_topic,
    general_topic_thr = 0.2
)

# ---------------------------------- compute DARs
from pycisTopic.diff_features import (
    impute_accessibility,
    normalize_scores,
    find_highly_variable_features,
    find_diff_features
)
import numpy as np

imputed_acc_obj = impute_accessibility(
    cistopic_obj_noDBL,
    selected_cells=None,
    selected_regions=None,
    scale_factor=10**6
)

normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)

variable_regions = find_highly_variable_features(
    normalized_imputed_acc_obj,
    min_disp = 0.05,
    min_mean = 0.0125,
    max_mean = 3,
    max_disp = np.inf,
    n_bins=20,
    n_top_features=None,
    plot=True
)
plt.savefig('variable_regions.png')
len(variable_regions)

markers_dict= find_diff_features(
    cistopic_obj_noDBL,
    imputed_acc_obj,
    variable='GEX:annotation_scenicPlus',
    var_features=variable_regions,
    contrasts=None,
    adjpval_thr=0.05,
    log2fc_thr=np.log2(1.5),
    n_cpu=5,
    #_temp_dir=tmp_dir,
    split_pattern = '-'
)

from pycisTopic.clust_vis import plot_imputed_features

plot_imputed_features(
    cistopic_obj_noDBL,
    reduction_name='UMAP',
    imputed_data=imputed_acc_obj,
    features=[markers_dict[x].index.tolist()[0] for x in ['cnt_Trophocytes', 'cnt_PDGFRalo','cnt_Telocytes','D7_Trophocytes', 'D7_PDGFRalo','D7_Telocytes','D14_Trophocytes', 'D14_PDGFRalo','D14_Telocytes']],
    scale=False,
    num_columns=4
)
plt.savefig('markers.png')

print("Number of DARs found:")
print("---------------------")
for x in markers_dict:
    print(f"  {x}: {len(markers_dict[x])}")

os.makedirs(os.path.join(work_dir, "region_sets"), exist_ok = True)
os.makedirs(os.path.join(work_dir, "region_sets", "Topics_otsu"), exist_ok = True)
os.makedirs(os.path.join(work_dir, "region_sets", "Topics_top_3k"), exist_ok = True)
os.makedirs(os.path.join(work_dir, "region_sets", "DARs_cell_type"), exist_ok = True)

from pycisTopic.utils import region_names_to_coordinates

for topic in region_bin_topics_otsu:
    region_names_to_coordinates(
        region_bin_topics_otsu[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(work_dir, "region_sets", "Topics_otsu", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )

for topic in region_bin_topics_top_3k:
    region_names_to_coordinates(
        region_bin_topics_top_3k[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(work_dir, "region_sets", "Topics_top_3k", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )

for cell_type in markers_dict:
    region_names_to_coordinates(
        markers_dict[cell_type].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(work_dir, "region_sets", "DARs_cell_type", f"{cell_type}.bed"),
        sep = "\t",
        header = False, index = False
    )

#----------------------------------------------------------------------------------



#-----------------------------------  Start scenicplus 
#make the same variable for the cistopic object 
import pickle
cistopic_obj_noDBL = pickle.load(open(os.path.join(work_dir, 'cistopic_obj_noDB.pkl'), 'rb'))

cistopic_obj_noDBL.cell_data['GEX:annotation_scenicPlus'][cistopic_obj_noDBL.cell_data['GEX:annotation_scenicPlus'] == 'cnt_Trophocytes'] = "D0_Trophocytes"
cistopic_obj_noDBL.cell_data['GEX:annotation_scenicPlus'][cistopic_obj_noDBL.cell_data['GEX:annotation_scenicPlus'] == 'cnt_PDGFRalo'] = "D0_CD81_stroma"
cistopic_obj_noDBL.cell_data['GEX:annotation_scenicPlus'][cistopic_obj_noDBL.cell_data['GEX:annotation_scenicPlus'] == 'cnt_Telocytes'] = "D0_SEMFs"

cistopic_obj_noDBL.cell_data['GEX:annotation_scenicPlus'][cistopic_obj_noDBL.cell_data['GEX:annotation_scenicPlus'] == 'D7_Trophocytes'] = "D7_Trophocytes"
cistopic_obj_noDBL.cell_data['GEX:annotation_scenicPlus'][cistopic_obj_noDBL.cell_data['GEX:annotation_scenicPlus'] == 'D7_PDGFRalo'] = "D7_CD81_stroma"
cistopic_obj_noDBL.cell_data['GEX:annotation_scenicPlus'][cistopic_obj_noDBL.cell_data['GEX:annotation_scenicPlus'] == 'D7_Telocytes'] = "D7_SEMFs"

cistopic_obj_noDBL.cell_data['GEX:annotation_scenicPlus'][cistopic_obj_noDBL.cell_data['GEX:annotation_scenicPlus'] == 'D14_Trophocytes'] = "D14_Trophocytes"
cistopic_obj_noDBL.cell_data['GEX:annotation_scenicPlus'][cistopic_obj_noDBL.cell_data['GEX:annotation_scenicPlus'] == 'D14_PDGFRalo'] = "D14_CD81_stroma"
cistopic_obj_noDBL.cell_data['GEX:annotation_scenicPlus'][cistopic_obj_noDBL.cell_data['GEX:annotation_scenicPlus'] == 'D14_Telocytes'] = "D14_SEMFs"

# this is the final scATAC file as input to SCENIC+
pickle.dump(
    cistopic_obj_noDBL,
    open(os.path.join(work_dir, "cistopic_obj_noDB_to_scenicplus.pkl"), "wb")
)

#--------------------------------- In order to create my CisTargerDatabase CB
# ----------------------------------in Bash -----------------------------
#!/bin/bash
# REGION_BED="/home/astavropoulou/scenicplus_D7_D14_D0_fibroblasts/consensus_peak_calling/consensus_regions.bed"
# GENOME_FASTA="/home/astavropoulou/scenicplus_D7_D14_D0_fibroblasts/mouse_genome_mm10_UCSC/mm10.fa"
# CHROMSIZES="/home/astavropoulou/scenicplus_D7_D14_D0_fibroblasts/mouse_genome_mm10_UCSC/mm10.chrom.sizes"
# DATABASE_PREFIX="D7_D14_cnt"
# SCRIPT_DIR="/home/astavropoulou/scenicplus_D7_D14_D0_fibroblasts/create_cisTarget_databases"
# 
# ${SCRIPT_DIR}/create_fasta_with_padded_bg_from_bed.sh \
# 	        ${GENOME_FASTA} \
# 		        ${CHROMSIZES} \
# 			        ${REGION_BED} \
# 				        mm10.D7_D14_Cnt.with_1kb_bg_padding.fa \
# 					        1000 \
# 						        yes
#                                 
        
        
#-------------------------------------------------------------
#!/bin/bash
# ./cbust --help
# 
# OUT_DIR=""${PWD}""
# CBDIR="${OUT_DIR}/aertslab_motif_colleciton/v10nr_clust_public/singletons"
# FASTA_FILE="${OUT_DIR}/mm10.D7_D14_Cnt.with_1kb_bg_padding.fa "
# MOTIF_LIST="${OUT_DIR}/motifs.txt"
# 
# REGION_BED="/home/astavropoulou/scenicplus_D7_D14_D0_fibroblasts/consensus_peak_calling/consensus_regions.bed"
# GENOME_FASTA="/home/astavropoulou/scenicplus_D7_D14_D0_fibroblasts/mouse_genome_mm10_UCSC/mm10.fa"
# CHROMSIZES="/home/astavropoulou/scenicplus_D7_D14_D0_fibroblasts/mouse_genome_mm10_UCSC/mm10.chrom.sizes"
# DATABASE_PREFIX="D7_D14_cnt"
# SCRIPT_DIR="/home/astavropoulou/scenicplus_D7_D14_D0_fibroblasts/create_cisTarget_databases"
# 
# "${SCRIPT_DIR}/create_cistarget_motif_databases.py" \
# 	    -f ${FASTA_FILE} \
# 	        -M ${CBDIR} \
# 		    -m ${MOTIF_LIST} \
# 		        -o ${OUT_DIR}/${DATABASE_PREFIX} \
# 			    --bgpadding 1000 \
# 			        -t 20


#-------------------------------- Downstream analysis after running scenicplus -------------------------------------
# -------------------- Meta information: number of pseudo-bulk cells per category 

# Load SCENIC+ mudata object (main result file)
scplus_mdata = mudata.read("/home/labuser/Athanasia/PhD_things/Prepare_paper_D7_D14/scenicplus_output_from_my_PC_the_run_in_server/scplusmdata.h5mu")

# Export regulon metadata tables for inspection (extended and direct)
scplus_mdata.uns["extended_e_regulon_metadata"].to_csv('extended_e_regulon_metadata.csv', sep= '\t')
scplus_mdata.uns["direct_e_regulon_metadata"].to_csv('direct_e_regulon_metadata.csv', sep='\t')

# --------------------- Focus only on DIRECT regulons ------------------------------
# Filter regulons to keep only:
#   - Direct (not extended)
#   - Positive regulation (regulation == 1)
#   - correlation with target genes (rho_R2G > 0)
scplus_mdata.uns['direct_e_regulon_metadata'] = scplus_mdata.uns['direct_e_regulon_metadata'][scplus_mdata.uns['direct_e_regulon_metadata']['is_extended'] == False]
scplus_mdata.uns['direct_e_regulon_metadata'] = scplus_mdata.uns['direct_e_regulon_metadata'][scplus_mdata.uns['direct_e_regulon_metadata']['regulation'] == 1]
scplus_mdata.uns['direct_e_regulon_metadata'] = scplus_mdata.uns['direct_e_regulon_metadata'][scplus_mdata.uns['direct_e_regulon_metadata']['rho_R2G'] > 0]

# Save filtered regulons sorted by importance or triplet rank
scplus_mdata.uns['direct_e_regulon_metadata'].sort_values(['eRegulon_name','importance_TF2G'],ascending=False).to_csv('Pos_pos_regulons_sorted_TF2G.csv', sep= '\t')
scplus_mdata.uns['direct_e_regulon_metadata'].sort_values(['eRegulon_name','triplet_rank'],ascending=True).to_csv('Pos_pos_regulons_sorted_triplet_rank.csv', sep= '\t')

# Add a simplified "celltype_timepoint" annotation from cell barcodes (remove trailing numbers)
scplus_mdata.obs['Celltype_Timepoint'] = scplus_mdata.obs.index.str.replace(r'.(\d*$)','').copy()

# -------------------------------------- Filter regulons based on mean AUC --------------------------------------
# Convert mudata into a scenicplus object (needed for downstream analysis)
scplus_obj = mudata_to_scenicplus(
    mdata = scplus_mdata,
    path_to_cistarget_h5 = "/home/labuser/Athanasia/PhD_things/Prepare_paper_D7_D14/scenicplus_output_from_my_PC_the_run_in_server/ctx_results.hdf5",
    path_to_dem_h5 = "/home/labuser/Athanasia/PhD_things/Prepare_paper_D7_D14/scenicplus_output_from_my_PC_the_run_in_server/dem_results.hdf5"
)

# Get lists of regulons, gene signatures, and region signatures
regulons = list(scplus_mdata.uns['direct_e_regulon_metadata']['eRegulon_name'].unique())
pos_regulons_names = scplus_mdata.uns['direct_e_regulon_metadata']['Gene_signature_name'].unique()
pos_regions_names = scplus_mdata.uns['direct_e_regulon_metadata']['Region_signature_name'].unique()

# Keep only the AUC values for the positive regulons we selected
scplus_obj.uns['eRegulon_AUC']['Gene_based'] = scplus_obj.uns['eRegulon_AUC']['Gene_based'][pos_regulons_names]
scplus_obj.uns['eRegulon_AUC']['Region_based'] = scplus_obj.uns['eRegulon_AUC']['Region_based'][pos_regions_names]

# Add simplified annotation to the AUC tables
scplus_obj.uns['eRegulon_AUC']['Gene_based']['celltype_timepoint'] = scplus_obj.uns['eRegulon_AUC']['Gene_based'].index.str.replace(r'.(\d*$)','').copy()
scplus_obj.uns['eRegulon_AUC']['Region_based']['celltype_timepoint'] = scplus_obj.uns['eRegulon_AUC']['Region_based'].index.str.replace(r'.(\d*$)','').copy()

# Count how many meta-cells per category exist
scplus_obj.uns['eRegulon_AUC']['Gene_based'].groupby("celltype_timepoint").size()

# Compute mean enrichment values for each celltype_timepoint
matrix_gene_based = scplus_obj.uns['eRegulon_AUC']['Gene_based'].groupby(['celltype_timepoint']).mean()[pos_regulons_names]
matrix_region_based = scplus_obj.uns['eRegulon_AUC']['Region_based'].groupby(['celltype_timepoint']).mean()[pos_regions_names]

# --------------------- Filtering regulons based on enrichment ---------------------
# Idea: keep regulons that are highly enriched in at least one category
# 1. For each regulon, get the maximum mean enrichment across categories
# 2. Keep only regulons above the median (50% quantile) of these max values
input_m = matrix_gene_based
input_m.max().sort_values(ascending=False)

# Keep only regulons with max mean enrichment > median
filtered_m = input_m[(input_m.max()[input_m.max() > np.quantile(input_m.max(), 0.5)]).index]
filtered_regulons = filtered_m.columns
len(filtered_regulons)

# --------------------- Further filtering of regulons ---------------------
# Keep only regulons present in the filtered list
scplus_mdata_tmp = scplus_mdata.copy()
scplus_mdata_tmp.uns['direct_e_regulon_metadata'] = scplus_mdata_tmp.uns['direct_e_regulon_metadata'][scplus_mdata_tmp.uns['direct_e_regulon_metadata']['Gene_signature_name'].isin(filtered_regulons)]

# Save filtered regulon metadata
scplus_mdata_tmp.uns['direct_e_regulon_metadata'].sort_values(['eRegulon_name','triplet_rank'],ascending=True).to_csv(
    'direct_e_regulon_metadata_filtered_basedOn_max_mean_enrichment_across_categories_GeneBased.csv',
    sep= '\t'
)

# Final manual curation: remove non-informative regulons
scplus_mdata_tmp.uns['direct_e_regulon_metadata'] = scplus_mdata_tmp.uns['direct_e_regulon_metadata'][~scplus_mdata_tmp.uns['direct_e_regulon_metadata']['TF'].isin(('Mypop','Foxk1','Tbx2','Atf5','Egr2'))]

# Save curated regulons
scplus_mdata_tmp.uns['direct_e_regulon_metadata'].sort_values(['eRegulon_name','triplet_rank'],ascending=True).to_csv(
    '/home/astavropoulou/D0_D7_D14_PAPER_FINAL_FILES_RUNS/analysis_and_plots/5_Scenicplus_D0_D7_D14/direct_e_regulon_metadata_filtered_basedOn_max_mean_enrichment_across_categories_GeneBased.csv',
    sep= '\t'
)
