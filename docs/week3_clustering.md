# Week 3: Preprocessing, Dimensionality Reduction & Clustering

## Overview
Transformation of raw counts into analysis-ready data and unsupervised cell clustering.

## Preprocessing Steps
1. **Normalization**: Library size normalization (10,000 counts per cell)
2. **Log transformation**: log1p transformation
3. **HVG Selection**: 3,000 highly variable genes (Seurat v3 method)
4. **Scaling**: Z-score scaling with regression of MT% and total counts

## Dimensionality Reduction
- **PCA**: 50 components computed, 30 used downstream
- **Neighborhood graph**: k=15 neighbors
- **UMAP**: min_dist=0.5 for visualization

## Clustering
- **Algorithm**: Leiden clustering
- **Resolution**: 1.0
- **Clusters identified**: 17

## Related Scripts
- `scripts/04_preprocessing.ipynb`
- `scripts/05_dimensionality_reduction.ipynb`

## Key Findings
- Clear separation of cells by condition in UMAP
- 17 distinct cell populations identified
- Inflammation-specific clusters emerge

## Notes
[Add your notes and observations here]
