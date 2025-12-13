# Week 4: Cell Type Annotation

## Overview
Assignment of biological cell type identities to clusters based on marker gene expression.

## Annotation Strategy
1. **Literature-based markers**: Parikh et al. (2019), Smillie et al. (2019)
2. **Differential expression**: Wilcoxon rank-sum test for cluster markers
3. **Cell type scoring**: Using `sc.tl.score_genes()`
4. **Manual validation**: Cross-reference with condition enrichment

## Cell Type Categories

### Absorptive Lineage
- Mature Enterocytes
- Mature Colonocytes
- Enterocytes (UC-noninfl)
- Colonocytes (UC-noninfl)
- Immature Enterocytes
- Inflammatory Enterocytes
- BEST4+ Enterocytes

### Secretory Lineage
- Goblet cells
- Inflammatory Goblet

### Stem/Progenitor Lineage
- Early Progenitors
- Transit-amplifying
- Inflammatory Stem/TA
- Inflammatory Progenitors

### Specialized Cells
- Stressed Epithelial
- Antigen-presenting Epithelial

### Non-epithelial (Contamination)
- T lymphocytes
- Mast cells

## Key Markers
- **LCN2**: Dominant inflammation marker
- **BEST4**: Specialized enterocyte subset
- **MUC2**: Goblet cells
- **OLFM4**: Stem cells
- **MKI67**: Proliferating cells

## Related Scripts
- `scripts/06_cell_type_annotation.ipynb`
- `scripts/06b_apply_annotations.ipynb`

## Results
- 17 cell types identified
- 4 inflammation-specific populations (>99% UC inflamed)
- 2,006 inflammatory cells (18% of dataset)

## Notes
[Add your notes and observations here]
