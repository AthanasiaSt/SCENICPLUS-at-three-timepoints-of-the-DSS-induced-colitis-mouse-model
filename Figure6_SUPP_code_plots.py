#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 22 18:24:53 2025

@author: Athanasia Stavropoulou
"""

# Code for Figure6 SUPP in the manuscript
import scanpy as sc

# dataset including seurat fibrobrast analysis processed with velocyto/ not included in the processed zenodo files
adata = sc.read_h5ad('..../velocyto/scVelo_inputs_outputs/D0_D7_D14_from_Seurat_with_spliced_unspliced.h5ad')

#-------general preprocessing after subsetting 

adata.layers["counts"] = adata.X.copy()

# Normalizing to median total counts
sc.pp.normalize_total(adata)
# Logarithmize the data
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="orig.ident")
sc.pl.highly_variable_genes(adata)
sc.tl.pca(adata)

#----------------------------------integration with harmony 
sc.external.pp.harmony_integrate(adata, 'group')

sc.pl.pca(
  adata,
  color=["orig.ident", "celltypes_origIdent_older"],
  ncols=2,
  size=2,
)

sc.pp.neighbors(
  adata,
  n_neighbors=30,
  n_pcs=10,
  use_rep="X_pca_harmony",
  key_added="harmony_neighbors",
)

sc.tl.umap(
  adata,
  min_dist=0.3,
  neighbors_key="harmony_neighbors"
)

sc.pl.umap(adata, color = "orig.ident")
sc.pl.umap(adata, color = "celltypes_origIdent_older")

# Define your custom colors as Python dicts
colours_conditions = {
  'D0': '#569B9A',
  'D0_kinchen': '#A1CEC5',
  'D7_kinchen': '#8e0f62',
  'D14': '#CD5808'
}

colours_celltypes = {
  'Trophocytes': '#8E7692',
  'PDGFRalo': '#a5af37',
  'Telocytes': '#416522'
}

sc.settings.set_figure_params(dpi=1000, dpi_save=1000, figsize=(6,6))

# UMAP colored by condition
sc.pl.umap(
  adata,
  color="orig.ident",   # column in adata.obs
  palette=colours_conditions,
  size=40,
  save="_condition.pdf"
)

# UMAP colored by celltype
sc.pl.umap(
  adata,
  color="celltypes_origIdent_older",
  palette=colours_celltypes,
  size=40,
  save="_celltype.pdf"
)


sc.tl.diffmap(adata)

sc.pl.scatter(
  adata,
  basis="diffmap",
  color=['orig.ident'],
  components=[1, 2],
  size=40,
  palette=[colours_conditions[x] for x in adata.obs['orig.ident'].cat.categories],
  save="_DPT_pseudotime_all_samples.pdf"
)

sc.pl.scatter(
  adata,
  basis="diffmap",
  color=['celltypes_origIdent_older'],
  components=[1, 2], size=40,
  palette=[colours_celltypes[x] for x in adata.obs['celltypes_origIdent_older'].cat.categories],
  save="_DPT_pseudotime_all_samples_celltypes.pdf"
)


# Setting root cell as described above
root_ixs = adata.obsm["X_diffmap"][:, 2].argmax()
adata.uns["iroot"] = root_ixs

sc.tl.dpt(adata)

sc.pl.scatter(
  adata,
  basis="umap",
  color=["dpt_pseudotime"],
  color_map="gnuplot2", 
  size=40,
  save="_DPT_pseudotime.pdf"
)