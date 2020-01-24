# Transcriptome of *Anopheles gambiae* oenocytes

Code and data to reproduce the analyses of differential gene expression from **Isolation and transcriptomic analysis of *An. gambiae* oenocytes enable the delineation of the cuticular hydrocarbon biosynthetic pathway** (Grigoraki et al. 2019).

## Quick how-to

To re-do the main analyses of the paper, just clone this repository and run the following `R` script:

```bash
git clone git@github.com:xgrau/oenocytes-agam.git    # or download it
Rscript s01_differential_expression_deseq2.R         # or open the script and tinker with it
```

The script will:

* load input data,
* perform the differential expression analysis between various groups of samples (grouped by sex and tissue; see below).
* perform functional enrichment analyses
* create figures and tables and various outputs (`results/` folder)

All necessary input data and gene annotations are included in this package (see **Contents** section below).

`R` libraries required to run these analyses are listed below (**Requirements** section).

### Contents

This repository contains the gene expression data and the necessary `R` code to run the main analysis.

Folders:

* `data_expression/`: sample-specific gene expression data, as produced by Salmon. It's formatted for easy loading into `DESeq2`
* `data_genome/`: genome annotations (GO, Pfam, gene names, etc.) required for functional enrichments
* `data_metadata/`: list of samples (12 in total) and sample classification. Sample classification codes are used in `DESeq2` to define specific comparisons between sample groups.
* `helper_scripts`: custom scripts for basic plots (heatmaps, volcano plots) and functional gene enrichment. Required by the main scripts.
* `results/`: results from the analysis.

### Output

The script will produce a number of files with tables and plots, all of them stored in the `results/` folder. Each file has a self-explanatory name, following this convention:

* `de_all`: differential expression analyses comparing all samples (heatmaps, PCAs, PCoAs, Venn diagrams with overlaps between sample groups)
* `de_co.Bulk.Fem-Mal`: files starting with this prefix contain information from a specific comparison (`de_co`), defined by the last three words in the name: in this case, it's a a female-to-male  comparison (`Fem-Mal`) of the carcass (`Bulk` tissue). For each comparison, we include:
  * `de_co.Bulk.Fem-Mal_volcano.pdf`: a volcano plot comparing expression in the bulk tissue, from females and males.
  * `de_co.Bulk.Fem-Mal_difpos.csv`: genes overexpressed in `Fem` (first term in filename).
  * `de_co.Bulk.Fem-Mal_difpos.hygeo.domain`: functional enrichments of Pfam domains in this list of genes, using the hypergeometric test.
  * `de_co.Bulk.Fem-Mal_difpos.topgo.fisherelim` functional enrichments of GOs in this list of genes, using Fisher's test and Elim weighting.
  * `de_co.Bulk.Fem-Mal_difpos.csv`: same, for genes overexpressed in `Mal` (second term in filename).
* `session.deseq_difexp.Bulk.Fem-Mal.csv` a very large table with all the differential expression statistics produced by DESeq, per each gene.
* `subset.FAdecarboxyl_p450`: heatmaps and tables of gene expression for subsets of genes involved in the FA synthesis pathway; including FA decarboxylases, desaturases, elongases, reductases and synthases.

### Sample classification

Definitions:

| Code | Category |  Description |
|-----| -- |-------|
| Fem | sex | female samples   |
| Mal | sex | male samples   |
| Oecy | celltype | oenocyte-enriched extracts   |
| Bulk | celltype | carcass |

List of samples:

| sample id | sample group |sex | celltype
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

In total, we perform four different comparisons, which are defined at the beginning of the script:

```R
# define comparisons:
# 1   -> category within which comparisons are made
# 2&3 -> groups to compare
# 4   -> extra category, to restrict comparisons within it
#        ALL if no restrictions are intended
# 5   -> group within the extra category that will be used;
#        ALL if no restrictions are intended
# 6   -> formula for DEseq comparisons
complis = list(
  convec=c("celltype","Oecy","Bulk","sex","Fem","~celltype"), # comparison of oecy vs. bulk in females
  convec=c("celltype","Oecy","Bulk","sex","Mal","~celltype"), # ...
  convec=c("sex","Fem","Mal","celltype","Bulk","~sex"),
  convec=c("sex","Fem","Mal","celltype","Oecy","~sex")
)
```

## Read mapping

The process of mapping of reads to the predicted has been performed with `Salmon` v0.10.2, as specified in the Methods section of the paper.

If you want to repeat the analysis yourself, these are the relevant commands (adjusting filenames accordingly):

```bash
# salmon index
salmon index -t transcripts.fasta -i transcripts.salmon.index --type quasi -k 31 1> log_index.log 2> &1
# run this command for each sample separately (A, B, etc.)
salmon quant -i Anogam_long.cds_mcherry.salmon.index -l A -p 10 -1 sampleA_1.fastq.gz -2 sampleA_2.fastq.gz -o sampleA_salmon_out 1> log_quant.log 2> &1
# output in sampleA_salmon_out, sampleB_salmon_out, etc. folders can be loaded to DESeq2 using the main script
```

The required files are not provided in this repository, but they can be found in public repositores:

* *Anopheles gambiae* transcripts: from Vectorbase, annotation version `AgamP4.9`.
* raw reads have been deposited to ENA, under the `XXX` accession.

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

If you reuse this code in any way, don't forget to cite the original papers for each of these libraries, too. It's easy and makes everyone happy. For example:

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
