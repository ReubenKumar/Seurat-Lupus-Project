# 🧬 Cutaneous Lupus Erythematosus (CLE) Single-Cell RNA-seq Analysis

🚧 **Project in progress**

This project analyzes **10x Genomics 3' single-cell RNA-seq data** from patients with **Cutaneous Lupus Erythematosus (CLE)** — a chronic inflammatory skin disease. The dataset includes skin tissue samples from three groups:
- **Healthy skin**
- **Non-inflamed CLE-affected skin**
- **Inflamed CLE-affected skin**

---

## 🎯 Objectives

- **Identify distinct cell types** present within the samples using the Seurat R package
- **Discover differentially expressed genes (DEGs)** between the different skin conditions and examine cell counts between conditions.

---

## 📄 Contents

- **`Lupus.md`** — Contains the RMarkdown script and results, including:
  - Code for processing and analysing the scRNA-seq data
  - Plots and visualisations for cell type identification and clustering
  - Supporting commentary and interpretation of the results

---

## 📌 Notes

- The analysis workflow uses the **Seurat** R package for single-cell data processing, normalisation, clustering, and visualisation.

---

## 📦 Requirements

- R (version ≥ 4.x)
- Seurat (version ≥ 5.x)
- dplyr, ggplot2, and other tidyverse packages

---

## 📊 Example Output

![UMAP plot showing clustered cell types](figure/umap_clusters.png)

---


