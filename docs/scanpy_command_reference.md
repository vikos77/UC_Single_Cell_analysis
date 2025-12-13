# Scanpy Command Reference

Quick reference guide for common Scanpy operations used in this project.

## Data Loading & Saving

```python
import scanpy as sc

# Read AnnData object
adata = sc.read_h5ad("file.h5ad")

# Save AnnData object
adata.write_h5ad("file.h5ad")
```

## Preprocessing

```python
# Quality control metrics
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

# Filtering
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Normalization
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Highly variable genes
sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor='seurat_v3')

# Scaling
sc.pp.scale(adata, max_value=10)
```

## Dimensionality Reduction

```python
# PCA
sc.tl.pca(adata, n_comps=50)

# Neighborhood graph
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)

# UMAP
sc.tl.umap(adata, min_dist=0.5)
```

## Clustering

```python
# Leiden clustering
sc.tl.leiden(adata, resolution=1.0)
```

## Differential Expression

```python
# Find marker genes
sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')

# Get results
markers = sc.get.rank_genes_groups_df(adata, group='0')
```

## Visualization

```python
# UMAP plots
sc.pl.umap(adata, color='leiden')
sc.pl.umap(adata, color='condition')
sc.pl.umap(adata, color=['gene1', 'gene2'])

# Violin plots
sc.pl.violin(adata, keys='n_genes', groupby='sample')

# Dot plot
sc.pl.dotplot(adata, var_names=marker_genes, groupby='leiden')

# Heatmap
sc.pl.heatmap(adata, var_names=marker_genes, groupby='leiden')

# Ranking plot
sc.pl.rank_genes_groups(adata, n_genes=20)
```

## Gene Scoring

```python
# Score gene sets
sc.tl.score_genes(adata, gene_list=['gene1', 'gene2'], score_name='score_name')
```

## Settings

```python
# Configure Scanpy
sc.settings.verbosity = 2  # 0=errors, 1=warnings, 2=info, 3=hints
sc.settings.set_figure_params(dpi=100, facecolor='white')
sc.settings.figdir = 'figures/'
```

## Common Attributes

```python
# AnnData structure
adata.X           # Expression matrix
adata.obs         # Cell metadata (rows)
adata.var         # Gene metadata (columns)
adata.uns         # Unstructured metadata
adata.obsm        # Multi-dimensional cell annotations (PCA, UMAP)
adata.varm        # Multi-dimensional gene annotations
adata.layers      # Alternative expression matrices

# Dimensions
adata.n_obs       # Number of cells
adata.n_vars      # Number of genes
adata.obs_names   # Cell barcodes
adata.var_names   # Gene names
```

## Useful Functions

```python
# Subsetting
adata_subset = adata[adata.obs['condition'] == 'UC', :]  # Select cells
adata_subset = adata[:, adata.var['highly_variable']]    # Select genes

# Copying
adata_copy = adata.copy()

# Concatenation
adata_merged = adata1.concatenate(adata2, batch_key='batch')
```

## Documentation
- Official docs: https://scanpy.readthedocs.io
- Tutorials: https://scanpy-tutorials.readthedocs.io
