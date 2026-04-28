# Single-Cell RNA-seq Analysis Pipeline  
### Using Scanpy, AnnData, and scverse

---

![Python](https://img.shields.io/badge/Python-3.9+-blue)
![Library](https://img.shields.io/badge/Scanpy-scRNA--seq-orange)
![Status](https://img.shields.io/badge/Status-Completed-success)
![License](https://img.shields.io/badge/License-MIT-green)

---

## Table of Contents

- [Overview](#overview)
- [Project Structure](#project-structure)
- [Workflow](#workflow)
- [AnnData Fundamentals](#1-anndata-fundamentals)
- [Quality Control](#2-data-loading--quality-control)
- [Filtering](#3-filtering)
- [Normalization](#4-normalization)
- [Feature Selection](#5-highly-variable-genes)
- [Dimensionality Reduction](#6-dimensionality-reduction)
- [Clustering](#7-clustering)
- [Annotation](#8-cell-type-annotation)
- [Differential Expression](#9-differential-expression)
- [Advanced Operations](#10-advanced-anndata-operations)
- [Installation](#installation)
- [Usage](#how-to-run)
- [Results](#results)
- [Future Work](#future-work)
- [Conclusion](#conclusion)

---

## Overview

This repository provides a complete and structured workflow for **single-cell RNA sequencing (scRNA-seq) analysis** using the **Scanpy** and **AnnData** ecosystem.

It integrates:
- Data structure understanding
- End-to-end analysis pipeline
- Advanced data handling techniques

The workflow progresses from **raw count data to biological insights** through reproducible computational steps.

---

## Project Structure

```
.
├── Getting_started_with_anndata.ipynb
├── scRNA_scanpy_tutorial.ipynb
├── scverse_tutorial.ipynb
├── images/
└── README.md
```

---

## Workflow

```
Raw Data
   ↓
AnnData Object
   ↓
Quality Control
   ↓
Filtering
   ↓
Normalization
   ↓
Feature Selection
   ↓
PCA / UMAP
   ↓
Clustering
   ↓
Annotation
   ↓
Differential Expression
```

---

# 1. AnnData Fundamentals

## Objective
Understand how scRNA-seq data is structured using `AnnData`.

## Code

```python
import numpy as np
import scanpy as sc
from scipy.sparse import csr_matrix

counts = csr_matrix(np.random.poisson(1, (100, 50)))
adata = sc.AnnData(counts)

adata.obs["cell_type"] = ["T-cell"] * 50 + ["B-cell"] * 50
adata.var["gene_symbol"] = [f"Gene{i}" for i in range(50)]
```

---

# 2. Data Loading & Quality Control

## Code

```python
sc.pp.calculate_qc_metrics(adata, inplace=True)
```

## Visualization

```python
sc.pl.violin(adata, ["total_counts", "n_genes_by_counts"], multi_panel=True)
```


---

# 3. Filtering

```python
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
```

### Doublet Detection

```python
sc.pp.scrublet(adata)
```

---

# 4. Normalization

```python
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
```

---

# 5. Highly Variable Genes

```python
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata = adata[:, adata.var.highly_variable]
```


---

# 6. Dimensionality Reduction

## PCA

```python
sc.tl.pca(adata)
sc.pl.pca(adata)
```


## UMAP

```python
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pl.umap(adata, color="cell_type")
```


---

# 7. Clustering

```python
sc.tl.leiden(adata, resolution=0.5)
sc.pl.umap(adata, color="leiden")
```


---

# 8. Cell Type Annotation

```python
marker_genes = {
    "T-cell": ["CD3D"],
    "B-cell": ["MS4A1"]
}

sc.pl.dotplot(adata, marker_genes, groupby="leiden")
```


---

# 9. Differential Expression

```python
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
sc.pl.rank_genes_groups(adata)
```


---

# 10. Advanced AnnData Operations

## Layers

```python
adata.layers["raw"] = adata.X.copy()
```

## Subsetting

```python
adata_subset = adata[adata.obs["cell_type"] == "T-cell"]
```

---

# Installation

```bash
pip install scanpy anndata matplotlib igraph pooch
```

---

# How to Run

```bash
git clone <your-repo-url>
cd your-repo
jupyter notebook
```

Run notebooks in order:

1. Getting_started_with_anndata.ipynb  
2. scRNA_scanpy_tutorial.ipynb  
3. scverse_tutorial.ipynb  

---

# Results

- High-quality cells retained after QC  
- Clear clustering of cell populations  
- UMAP visualization of cell distributions  
- Identification of marker genes  
- Biological interpretation of clusters  

---

# Conclusion

This project demonstrates a complete and reproducible pipeline for scRNA-seq analysis, combining efficient data structures with modern computational techniques to extract meaningful biological insights.

---
