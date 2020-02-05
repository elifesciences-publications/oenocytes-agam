# # load libraries
# library(GenomicFeatures)
# library(ape)
# library(stringr)
# library(plyr)
# library(tidyr)
library(topGO)
# library(Biostrings)
# library(Rsamtools)
# library(Rfast)
# library(cluster)
# library(gplots)

### Define input ####

# input files
si       = "Anogam"                                  # Species name (requires Spi_long.annot.gff and Spi_gDNA.fasta)
ioi_fn   = "data_as_suppa/Anogam_suppa.2020-02-05.ioi_isoforms_isoform.psi"   # IOE suppa file
tpm_fn   = "data_as_suppa/Anogam_suppa.2020-02-05.isoforms.tpm"   # IOE suppa file
outcode  = "results_as/"                             # Output code SUPPA analysis (eg *14mai18*, used for retrieving suppa groupings!)
samgrup  = "data_metadata/samples_classification.csv"                  # sample list (first col) & grouping (2nd, 3rd, etc.)
grupnoms    = c("group","sex","tissue")                 # group names (2nd, 3rd, 4th... cols)
annottab    = "data_genome/Anogam_long.pep_Pfamscan.seqs"               # pfam annots (requires .blat file too)
gomapfile   = "data_genome/Anogam_long.pep_eggnog_diamond.emapper.annotations.GO" # gomap (with transcripts)
tx2genedict = "data_genome/Anogam_tx2ge.csv"                         # transcript-to-gene mappings
genename_fn = "data_genome/Anogam_genes.description.csv"             # gene names

# input variables
pthr = 0.05 # pval threshold for suppa comparisons
k_list = 3:10 # list of k values to try for kmeans clustering

# load functions
source("helper_scripts/dimred_18jul18.R")
source("helper_scripts/pval_from_CDF_8nov18.R")
source("helper_scripts/volcano_plot.R")
source("helper_scripts/geneSetAnalysis.R")
col.fun = colorRampPalette(interpolate="l",c("aliceblue","deepskyblue","dodgerblue4"))


#### Load data ####

sf = read.table(samgrup,header = F)
colnames(sf) = c("sample",grupnoms)

# tx2gene
tx2gene = read.table(tx2genedict)
colnames(tx2gene) = c("transcript","gene")
# pfam annotations
panno = read.table(file = annottab)
colnames(panno) = c("transcript","pstart","pend","pfamid","domain","domseq")
panno = merge(panno, tx2gene, by.x="transcript",by.y="transcript")
pannu = subset(panno, select=c("transcript","pfamid","domain","gene"))
pannu = pannu[!duplicated(pannu), ]
# load GOmappings
gomap = readMappings(gomapfile)
gomap_gene_names = sapply(base::strsplit( names(gomap), split="-"), "[" , 1 )
names(gomap) = gomap_gene_names
# gene names
gname = read.table(genename_fn)
colnames(gname) = c("gene","gene_name")

# matrix of TPMs transcript
tpm = read.table(tpm_fn)

# load isoforom PSIs
ioi = read.table(ioi_fn, header = T)
list_genes = sapply(base::strsplit( rownames(ioi), split=";"), "[" , 1 )
list_transcripts = sapply(base::strsplit( rownames(ioi), split=";"), "[" , 2 )
rownames(ioi) = list_transcripts

# # remove genes that only have one isoform
# list_genes_with_iso = list_genes[table(list_genes) > 1]
# list_transcripts_with_iso = as.vector(tx2gene[tx2gene$gene %in% list_genes_with_iso, "transcript"])
# ioi_with_iso = ioi[ rownames(ioi) %in% list_transcripts_with_iso ,]
# list_transcripts_with_iso = rownames(ioi_with_iso)
# list_genes_with_iso = sapply(base::strsplit( list_transcripts_with_iso, split="-"), "[" , 1 )
# 
# # restrict TPMs to genes with observed isoforms
# tpm_with_iso = tpm[rownames(tpm) %in% rownames(ioi_with_iso),]


# remove invariant isoforms, and isoforms where all entries are nan 
# This is required for proper clustering, but not for differential splicing calculations
# (where we want to include constitutive isoforms too!)
is_invariant_ioi = as.vector(apply(ioi, 1, sd, na.rm = T) == 0 | is.na(apply(ioi, 1, sd, na.rm = T)))
is_nan_ioi = as.vector(apply(is.na(ioi), 1, sum, na.rm = T) > 0)
ioi_clean = ioi[!is_invariant_ioi & !is_nan_ioi,]
list_genes_ioi_clean = sapply(base::strsplit( rownames(ioi_clean), split="-"), "[" , 1 )

# restrict annotations to genes with observed isoforms
# If we don't do that, functional enrichments will be naturally biased towards types
# of genes that have more isoforms (long neural genes?).
pannu_with_iso = pannu[pannu$gene %in% list_genes_ioi_clean,]
gomap_with_iso = gomap[names(gomap) %in% list_genes_ioi_clean]

# clean tpms too?
# tpm_clean = tpm[rownames(tpm) %in% rownames(ioi_clean),]

# now, remove isoforms that do not have alternative isoforms (a priori, totally pointless)
# list_transcripts_clean = rownames(ioi_clean)
# list_genes_clean = list_genes[list_transcripts %in% list_transcripts_clean]
# list_genes_multi_isoform = list_genes_clean[list_genes_clean %in% names(sum(table(list_genes_clean) == 1))]



#### General clustering ####
# PCA, PCoA, clustering
ioi_clus = dimredfun(matriu = as.matrix(ioi_clean) , outputname = sprintf("%s/isoform_psi_clus", outcode), 
                     varname = "PSI", cols_are = "samples", rows_are = "isoforms", isbidi = F,
                     cols_dist_method = "pearson", rows_dist_method = "pearson", clus_method = "ward.D2")

# kmeans clustering
ioi_kclu = kmeansfun(
  matrix=ioi_clean,  
  outputprefix = sprintf("%s/isoform_psi_kmeans", outcode), 
  k_list = c(3:10), ylab = "PSI",ylims = c(0,1))


#### Differential splicing ####

# differentially included isoforms, with pval
ioi_diff_OeCa_F = pval_from_CDF(
  matrix_frequencies = ioi, matrix_tpms = tpm,
  name_group_i = "OecyF",name_group_j = "BulkF",
  samples_group_i = as.vector(sf[sf$group=="OecyF",]$sample) ,samples_group_j = as.vector(sf[sf$group=="BulkF",]$sample),
  comparison_code = "OeCa F",
  pval_correction = "BH",
  ecdf_area = 1000)

ioi_diff_OeCa_M = pval_from_CDF(
  matrix_frequencies = ioi, matrix_tpms = tpm,
  name_group_i = "OecyM",name_group_j = "BulkM",
  samples_group_i = as.vector(sf[sf$group=="OecyM",]$sample) ,samples_group_j = as.vector(sf[sf$group=="BulkM",]$sample),
  comparison_code = "OeCa M",
  pval_correction = "BH",
  ecdf_area = 1000)

# plot tpm v. dPSI values
pdf(file=sprintf("%s/evaluate_dpsi_tpm_pCDF.pdf", outcode),height=5,width=8)
ioi_diff_OeCa_F$plot()
ioi_diff_OeCa_M$plot()
dev.off()

# volcano plots
ioi_diff_OeCa_F_volcano = volcanoexp(
  table = ioi_diff_OeCa_F$result, 
  plotname = "OeCa_F", 
  fileprefix = sprintf("%s/diff", outcode), pthreshold = 0.05, 
  fc_varname = "freq_diff", fcp = "Oe", fcn = "Ca", pval_varname = "pval", xlims = c(-1,1), ylims = c(0,-5))

ioi_diff_OeCa_M_volcano = volcanoexp(
  table = ioi_diff_OeCa_M$result, 
  plotname = "OeCa_M", 
  fileprefix = sprintf("%s/diff", outcode), pthreshold = 0.05, 
  fc_varname = "freq_diff", fcp = "Oe", fcn = "Ca", pval_varname = "pval", xlims = c(-1,1), ylims = c(0,-5))

# save tables
write.table(ioi_diff_OeCa_F$result,file=sprintf("%s/diff.OeCa_F_results.csv", outcode), sep="\t", quote = F, row.names = F)
write.table(ioi_diff_OeCa_M$result,file=sprintf("%s/diff.OeCa_M_results.csv", outcode), sep="\t", quote = F, row.names = F)

# list of genes with isoforms that are included in Oe
ioi_diff_OeCa_F_genes = as.vector(unique(tx2gene[tx2gene$transcript %in% ioi_diff_OeCa_F_volcano$genes_sp,"gene"]))
ioi_diff_OeCa_M_genes = as.vector(unique(tx2gene[tx2gene$transcript %in% ioi_diff_OeCa_M_volcano$genes_sp,"gene"]))




#### Overlaps ####

pdf(file=sprintf("%s/overlap_dpsi_OeCa_FM.pdf", outcode),height=4,width=4)

# overlap at gene level: which genes are differentially spliced in males and females, between Oe and Ca?
overlap_OeCa_FM_genes = venn.two(
  list1 = ioi_diff_OeCa_F_genes, list2 = ioi_diff_OeCa_M_genes, catname1 = "Female", catname2 = "Male", 
  main = "OeCa dPSI genes in females and males", col1 = "cyan4", col2 =  "orange")

# overlap at isoform level: which isoforms are differentially INCLUDED IN OENOCYTES in males and females? (aka dPSI>0)
overlap_OeCa_FM_isoforms = venn.two(
  list1 = ioi_diff_OeCa_F_volcano$genes_sp, list2 = ioi_diff_OeCa_M_volcano$genes_sp, catname1 = "Female", catname2 = "Male", 
  main = "OeCa dPSI>0 isoforms in females and males (Oe)", col1 = "cyan4", col2 =  "orange")

# overlap at isoform level: which isoforms are differentially INCLUDED IN CARCASS in males and females? (aka dPSI<0); MIRROR IMAGE OF PREVIOUS PLOT
overlap_OeCa_FM_isoforms = venn.two(
  list1 = ioi_diff_OeCa_F_volcano$genes_sn, list2 = ioi_diff_OeCa_M_volcano$genes_sn, catname1 = "Female", catname2 = "Male", 
  main = "OeCa dPSI<0 isoforms in females and males (Ca)", col1 = "cyan4", col2 =  "orange")

dev.off()




#### Functional enrichment ####

# pfam enrichments
hygeofun(list_interest=ioi_diff_OeCa_F_genes, 
         annotation=pannu,gene_col="gene",ano_col="domain",
         outputname=sprintf("%s/diff", outcode),
         name_geneset="OeCa_F",topnum = 40, padjbool = F)

hygeofun(list_interest=ioi_diff_OeCa_M_genes, 
         annotation=pannu,gene_col="gene",ano_col="domain",
         outputname=sprintf("%s/diff", outcode),
         name_geneset="OeCa_M",topnum = 40, padjbool = F)

# GO enrichments
suppressMessages(topgofun(
  list_interest=ioi_diff_OeCa_F_genes,gomap=gomap,
  ontologyset=c("BP","MF","CC"),tg_test="fisher",tg_algorithm="elim",
  outputname=sprintf("%s/diff", outcode),
  name_geneset="OeCa_F",topnum=40))

suppressMessages(topgofun(
  list_interest=ioi_diff_OeCa_M_genes,gomap=gomap,
  ontologyset=c("BP","MF","CC"),tg_test="fisher",tg_algorithm="elim",
  outputname=sprintf("%s/diff", outcode),
  name_geneset="OeCa_M",topnum=40))



#### Lists of genes ####


pdf(file=sprintf("%s/fa_genes_diffspl.pdf", outcode),height=6,width=8)


# FA elongases: one
interestlist = as.character(pannu[ pannu$domain == "ELO" ,]$gene)
interestlist = unique(c(ioi_diff_OeCa_F_genes[ioi_diff_OeCa_F_genes %in% interestlist] , ioi_diff_OeCa_M_genes[ioi_diff_OeCa_M_genes %in% interestlist]))

for (gene in interestlist) {
  
  # psi data main heatmap
  txis = as.vector(tx2gene[tx2gene$gene == gene, "transcript"])
  tpsi = ioi[rownames(ioi) %in% txis, ]
  tpsi = tpsi[ order(row.names(tpsi)),]
  
  # retrieve pval for diff splicing and lateral annotation
  tann = data.frame(
    "OeCa_M" = as.numeric(ioi_diff_OeCa_M$result[ioi_diff_OeCa_M$result$event %in% rownames(tpsi), "pval"] < 0.05),
    "OeCa_F" = as.numeric(ioi_diff_OeCa_F$result[ioi_diff_OeCa_F$result$event %in% rownames(tpsi), "pval"] < 0.05),
    row.names = rownames(ioi_diff_OeCa_M$result[ioi_diff_OeCa_M$result$event %in% rownames(tpsi),])
  )
  
  # gene name
  gnom = as.character(gname[gname$gene == gene, "gene_name"])
  
  # clean names for nicer plots
  rownames(tpsi) = gsub("Anogam_","", rownames(tpsi))
  rownames(tann) = gsub("Anogam_","", rownames(tann))
  
  # plot heatmap
  pheatmap(tpsi, color = col.fun(20), breaks = seq(0,1,length.out = 20), 
           cellwidth = 18, cellheight = 12, na_col = "dodgerblue4",
           border_color = "white", cluster_cols=F, cluster_rows=F,display_numbers = T, number_color = "red",
           gaps_col = seq(0,12, by=3), annotation_row = tann, 
           main=sprintf("PSI\n%s\n%s", gene, gnom))
}


# FA desaturases: three
interestlist = as.character(pannu[ pannu$domain == "FA_desaturase" ,]$gene)
interestlist = c(interestlist, "Anogam_AGAP003050") # added manually because of Pfam misannotation
interestlist = unique(c(ioi_diff_OeCa_F_genes[ioi_diff_OeCa_F_genes %in% interestlist] , ioi_diff_OeCa_M_genes[ioi_diff_OeCa_M_genes %in% interestlist]))

for (gene in interestlist) {
  
  # psi data main heatmap
  txis = as.vector(tx2gene[tx2gene$gene == gene, "transcript"])
  tpsi = ioi[rownames(ioi) %in% txis, ]
  tpsi = tpsi[ order(row.names(tpsi)),]
  
  # retrieve pval for diff splicing and lateral annotation
  tann = data.frame(
    "OeCa_M" = as.numeric(ioi_diff_OeCa_M$result[ioi_diff_OeCa_M$result$event %in% rownames(tpsi), "pval"] < 0.05),
    "OeCa_F" = as.numeric(ioi_diff_OeCa_F$result[ioi_diff_OeCa_F$result$event %in% rownames(tpsi), "pval"] < 0.05),
    row.names = rownames(ioi_diff_OeCa_M$result[ioi_diff_OeCa_M$result$event %in% rownames(tpsi),])
  )
  
  # gene name
  gnom = as.character(gname[gname$gene == gene, "gene_name"])
  
  # clean names for nicer plots
  rownames(tpsi) = gsub("Anogam_","", rownames(tpsi))
  rownames(tann) = gsub("Anogam_","", rownames(tann))
  
  # plot heatmap
  pheatmap(tpsi, color = col.fun(20), breaks = seq(0,1,length.out = 20), 
           cellwidth = 18, cellheight = 12, na_col = "dodgerblue4",
           border_color = "white", cluster_cols=F, cluster_rows=F,display_numbers = T, number_color = "red",
           gaps_col = seq(0,12, by=3), annotation_row = tann, 
           main=sprintf("PSI\n%s\n%s", gene, gnom))
}



# FA decarboxylases: three (one not present in males)
interestlist = as.character(pannu[ pannu$domain == "p450" ,]$gene)
interestlist = unique(c(ioi_diff_OeCa_F_genes[ioi_diff_OeCa_F_genes %in% interestlist] , ioi_diff_OeCa_M_genes[ioi_diff_OeCa_M_genes %in% interestlist]))

for (gene in interestlist) {
  
  # psi data main heatmap
  txis = as.vector(tx2gene[tx2gene$gene == gene, "transcript"])
  tpsi = ioi[rownames(ioi) %in% txis, ]
  tpsi = tpsi[ order(row.names(tpsi)),]
  
  # retrieve pval for diff splicing and lateral annotation
  tann = data.frame(
    "OeCa_M" = as.numeric(ioi_diff_OeCa_M$result[ioi_diff_OeCa_M$result$event %in% rownames(tpsi), "pval"] < 0.05),
    "OeCa_F" = as.numeric(ioi_diff_OeCa_F$result[ioi_diff_OeCa_F$result$event %in% rownames(tpsi), "pval"] < 0.05),
    row.names = rownames(ioi_diff_OeCa_M$result[ioi_diff_OeCa_M$result$event %in% rownames(tpsi),])
  )
  
  # gene name
  gnom = as.character(gname[gname$gene == gene, "gene_name"])
  
  # clean names for nicer plots
  rownames(tpsi) = gsub("Anogam_","", rownames(tpsi))
  rownames(tann) = gsub("Anogam_","", rownames(tann))
  
  # plot heatmap
  pheatmap(tpsi, color = col.fun(20), breaks = seq(0,1,length.out = 20), 
           cellwidth = 18, cellheight = 12, na_col = "dodgerblue4",
           border_color = "white", cluster_cols=F, cluster_rows=F,display_numbers = T, number_color = "red",
           gaps_col = seq(0,12, by=3), annotation_row = tann, 
           main=sprintf("PSI\n%s\n%s", gene, gnom))
}



# FA synthases: ignore them, there is not diff splicing
interestlist = as.character(pannu[ pannu$domain == "ketoacyl-synt" ,]$gene)
interestlist = unique(c(ioi_diff_OeCa_F_genes[ioi_diff_OeCa_F_genes %in% interestlist] , ioi_diff_OeCa_M_genes[ioi_diff_OeCa_M_genes %in% interestlist]))

# FA reductases: none
interestlist = as.character(pannu[ pannu$domain == "NAD_binding_4" ,]$gene)
interestlist = unique(c(ioi_diff_OeCa_F_genes[ioi_diff_OeCa_F_genes %in% interestlist] , ioi_diff_OeCa_M_genes[ioi_diff_OeCa_M_genes %in% interestlist]))


dev.off()