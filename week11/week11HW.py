#!/usr/bin/env python

import scanpy as sc
import matplotlib.pyplot as plt

# Read the 10x dataset filtered down to just the highly-variable genes
adata = sc.read_h5ad("variable_data.h5")
adata.uns['log1p']['base'] = None # This is needed due to a bug in scanpy 

sc.pp.neighbors(adata, n_neighbors = 10, n_pcs = 40)

#Exercise 1: Clustering
#Step 1.2

sc.tl.leiden(adata)

#Step 1.3

fig, axes = plt.subplots(ncols = 2, figsize = (12, 5))

#Plot UMAP
sc.tl.umap(adata, maxiter = 900)
sc.pl.umap(adata, color = 'leiden', ax = axes[0], title = 'UMAP', show = False)

#Plot t-SNE
sc.tl.tsne(adata)
sc.pl.tsne(adata, color = 'leiden', ax = axes[1], title = 't-SNE', show = False)

#Show the figure
plt.savefig('t-SNE_UMAP_plots.png', dpi = 300, bbox_inches = 'tight')
plt.show()

#Exercise 2: Identifying cluster marker genes
#Step2.1

wilcoxon_adata = sc.tl.rank_genes_groups(adata, groupby = 'leiden', method = 'wilcoxon', use_raw = True, copy = True)

logreg_adata = sc.tl.rank_genes_groups(adata, groupby = 'leiden', method = 'logreg', use_raw = True, copy = True, max_iter = 900)

#Step2.2

fig, ax = plt.subplots(figsize = (12, 8), ncols = 5, nrows = 5, sharex = True, sharey = False)

sc.pl.rank_genes_groups(wilcoxon_adata, n_genes = 25, ax = ax, title = 'Wilcoxon Gene Expression', use_raw = True, show = False)

fig.suptitle('Top 25 Genes (Wilcoxon Gene Expression)', fontsize=16)
plt.subplots_adjust(top = 0.9)  
plt.savefig('Wilcoxon_based_rank_genes_groups.png', dpi = 300, bbox_inches = 'tight')
plt.close()

fig1, ax1 = plt.subplots(figsize = (12, 8), ncols = 5, nrows = 5, sharex = True, sharey = False)

sc.pl.rank_genes_groups(logreg_adata, n_genes = 25, ax = ax, title = 'Logistic Regression', use_raw = True, show = False)

fig1.suptitle('Top 25 Genes (Logistic Regression)', fontsize=16)
plt.subplots_adjust(top = 0.9)  

plt.savefig('Logreg_based_rank_genes_groups.png', dpi = 300)
plt.close()

#Exercise 3: Identifying cluster cell types

#Step 3.1

leiden = adata.obs['leiden']
umap = adata.obsm['X_umap']
tsne = adata.obsm['X_tsne']
adata = sc.read_h5ad('filtered_data.h5')
adata.obs['leiden'] = leiden
adata.obsm['X_umap'] = umap
adata.obsm['X_tsne'] = tsne

adata.write('filtered_clustered_data.h5')

adata = sc.read_h5ad("filtered_clustered_data.h5")
adata.uns['log1p']['base'] = None # This is needed due to a bug in scanpy 

sc.tl.rank_genes_groups(adata, 'leiden')

#Step 3.2

marker_genes = ['RPS12\n LDHB','S100A9\n LYZ','CD74\n CD79A','CCL5\n NKG7','LST1\n AIF1','NKG7\n GNLY','HLA-DPA1\n HLA-DPB2','PF4\n SDPR']
adata.rename_categories('leiden', marker_genes)

sc.pl.tsne(adata, color = 'leiden', legend_loc = 'on data', show = False)

plt.suptitle('TSNE Plot', fontsize = 16)
plt.subplots_adjust(top = 0.9)

plt.savefig('tsne_top_gene.png', dpi = 300)

#Step 3.3

top_marker_genes = ['Not enough info','Monocytes','B Cell','CD8+T','Macrophage','NK','Antigen presenting Cell','HSC']
adata.rename_categories('leiden', top_marker_genes)

sc.pl.tsne(adata, color = 'leiden', legend_loc = 'on data', show = False)

plt.suptitle('TSNE Plot', fontsize = 16)
plt.subplots_adjust(top = 0.9)

plt.savefig('tsne_with_labels.png', dpi = 300)
plt.show()



