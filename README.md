# Transcriptome of *Anopheles gambiae* oenocytes

Code and data to reproduce the analyses of differential gene expression from **Isolation and transcriptomic analysis of An. gambiae oenocytes enable the delineation of the cuticular hydrocarbon biosynthetic pathway** (Grigoraki et al. 2019).

## Contents

This repository contains the gene expression data and the necessary `R` code to run the main analysis.

Folders:

* `data_expression/`: sample-specific gene expression data, as produced by Salmon. It's formatted for easy loading into `DESeq2`
* `data_genome/`: genome annotations (GO, Pfam, gene names, etc.) required for functional enrichments
* `data_metadata/`: list of samples (12 in total) and sample classification. Sample classification codes are used in `DESeq2` to define specific comparisons between sample groups.
* `helper_scripts`: custom scripts for basic plots (heatmaps, volcano plots) and functional gene enrichment. Required by the main scripts.

### Sample classification

See:

| Code | Category |  Description |
|-----| -- |-------|
| Fem | sex | female samples   |
| Mal | sex | male samples   |
| Oecy | tissue | oenocyte-enriched extracts   |
| Bulk | tissue | carcass |

See:

| sample id | sample group |sex | tissue
| --- | --- | --- | --- |
| s01ACF | OecyF | Fem | Oecy
| s02ACF | OecyF | Fem | Oecy
| s03ACF | OecyF | Fem | Oecy
| s04ACM | OecyM | Mal | Oecy
| s05ACM | OecyM | Mal | Oecy
| s06ACM | OecyM | Mal | Oecy
| s07BCF | BulkF | Fem | Bulk
| s08BCF | BulkF | Fem | Bulk
| s09BCF | BulkF | Fem | Bulk
| s10BCM | BulkM | Mal | Bulk
| s11BCM | BulkM | Mal | Bulk
| s12BCM | BulkM | Mal | Bulk

## Requirements

Analyses for the paper were run on `R` 3.6.1. 

The following libraries are required to run the main script:

```R
library(DESeq2)
library(tximport)
library(ape)
library(gplots)
library(stringr)
library(plyr)
library(tidyr)
library(topGO)
library(readr)
library(pheatmap)
library(cluster)
library(VennDiagram)
```

If you reuse this code in any way, don't forget to cite the original papers for each of these libraries, too. It's easy and makes everyone happy:

```R
> citation("DESeq2")

  Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2 Genome Biology 15(12):550
  (2014)

A BibTeX entry for LaTeX users is

  @Article{,
    title = {Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2},
    author = {Michael I. Love and Wolfgang Huber and Simon Anders},
    year = {2014},
    journal = {Genome Biology},
    doi = {10.1186/s13059-014-0550-8},
    volume = {15},
    issue = {12},
    pages = {550},
  }
```