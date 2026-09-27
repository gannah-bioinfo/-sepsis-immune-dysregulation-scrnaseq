# Deciphering Immune Dysregulation in Sepsis Through Single-Cell Transcriptomics

**Author:** Gannah Sayed Saleh

## Research Question

How does sepsis alter immune-cell composition and transcriptional states, and which immune populations and pathways are associated with disease severity?

## Overview

This project analyzes whole-blood single-cell RNA-seq (scRNA-seq) data from two independent public cohorts to characterize immune dysregulation in sepsis:

- **GSE163668** — 46 individuals (28 sepsis/COVID-19 patients, 18 healthy controls). Includes both individually indexed and genetically demultiplexed pooled 10x libraries.
- **GSE157789** — 30 individuals across four clinical groups (Severe COVID-ARDS, Severe non-COVID/bacterial ICU, Healthy Control, COVID Recovered), used as an independent cross-cohort validation and to compare COVID- vs. bacterial-driven sepsis.

## Pipeline

Implemented in R (Seurat v5):

1. **Preprocessing** — data download (`GEOquery`), demultiplexing of pooled libraries (freemuxlet cluster assignments), merging into per-cohort Seurat objects.
2. **Quality control** — filtering on gene counts and mitochondrial content, doublet removal (`scDblFinder`).
3. **Normalization & feature selection** — log-normalization, top 2,000 variable genes.
4. **Dimensionality reduction** — PCA (15 PCs).
5. **Integration** — batch correction across donors with `Harmony`.
6. **Clustering & annotation** — graph-based clustering (UMAP + Louvain), marker-gene-based cell-type annotation (`FindAllMarkers`).
7. **Statistical analysis** — per-cell and patient-level composition testing (χ², Wilcoxon with FDR correction), differential expression (`FindMarkers`), pathway enrichment (`clusterProfiler`, GO Biological Process).

## Key Findings

- Sepsis drives a myeloid-dominant compositional shift: neutrophils expand from <1% (controls) to ~37% (patients), with concurrent lymphopenia. Confirmed significant both per-cell (χ² = 14,404, p < 2.2×10⁻¹⁶) and at the patient level (n = 46, Wilcoxon FDR < 0.05 for 8/14 populations).
- Monocytes and neutrophils show distinct transcriptional/pathway programs (antiviral/interferon vs. detoxification/metabolic stress).
- Neutrophil proportion scales with clinical severity (24% mild/moderate → 43% severe).
- Cross-cohort validation (GSE157789) shows neutrophil expansion also tracks with disease stage, and reveals an etiology-specific signature: interferon-driven in COVID-ARDS neutrophils vs. bacterial/LPS-driven (IL1B, CXCL8) in non-COVID ARDS monocytes.

## Repository Structure

```
scripts/
  01_GSE163668_pipeline.R   Full analysis pipeline for the primary cohort (GSE163668)
reports/
  Sepsis_scRNAseq_Report.docx   Full scientific report with figures
  Sepsis_OnePager.docx          One-page paper-style summary
```

## Data Availability

Raw data: NCBI GEO accessions [GSE163668](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163668) and [GSE157789](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157789).
