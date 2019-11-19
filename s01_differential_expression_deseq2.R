#### Input vars ####

# genome annotations
si      = "Anogam"
gfffile = "data_genome/Anogam_long.annot.gff"
txgdict = "data_genome/Anogam_long.cds.csv"
nmgdict = "data_genome/Anogam_genes.description.csv"
gomapfi = "data_genome/Anogam_long.pep_eggnog_diamond.emapper.annotations.GO"
panfile = "data_genome/Anogam_long.pep_Pfamscan.seqs"
panform = "pfamscan"

# where to store output?
outcode = "results_de/" # (folder + prefix)


# salmon expression files
saloutf = "data_expression/"

# define comparisons: metadata for each sample
samgrup = "data_metadata/samples_classification.csv" # can contain >1 grouping (eg location+treatment)

# sample classification
samclas = c("group","sex","celltype")

# global formula DEseq (comparison-specific formulae indicated below)
expdiss = "~celltype+sex" 

# define comparisons in the following list of vectors:
# 1   -> category within which comparisons are made
# 2&3 -> groups to compare
# 4   -> extra category, to restrict comparisons within it
#        ALL if no restrictions are intended
# 5   -> group within the extra category that will be used;
#        ALL if no restrictions are intended
# 6   -> formula for DEseq comparisons
complis = list(
  convec=c("celltype","Oecy","Bulk","sex","Fem","~celltype"),
  convec=c("celltype","Oecy","Bulk","sex","Mal","~celltype"),
  convec=c("sex","Fem","Mal","celltype","Bulk","~sex"),
  convec=c("sex","Fem","Mal","celltype","Oecy","~sex")
)

# pairs of comparisons (use index from complis) that will be
# compared in the Venn diagrams to see overlap in DE results
# BEWARE: MUST USE SAME CONVENTION FOR OVER/UNDEREXPRESSION!
# Same group of samples must be "first" in both comparisons.
# If no comparison is to be done, just say c(1,1)
complis_pairs = list(
  complis[c(1,2)]
)

# input variables
pval_threshold  = 0.001
alpha_threshold = 0.001
fc_threshold    = 2      # FC threshold used to define over/under expression
# FC = 1  -> log2FC = 0
# FC = 2  -> log2FC = 1
# FC = 10 -> log2FC = 3.27...
logfc_threshold = log2(fc_threshold)


# load libraries
library(ape)
library(gplots)
library(stringr)
library(plyr)
library(tidyr)
library(topGO)
library(tximport)
library(readr)
library(DESeq2)
library(pheatmap)
library(cluster)
library(VennDiagram)

# load functions
source("helper_scripts/geneSetAnalysis.R")
source("helper_scripts/dimred_18jul18.R")
source("helper_scripts/slidingfunctions_18nov18.R")
source("helper_scripts/volcano_plot.R")
col.fun = colorRampPalette(interpolate="l",c("aliceblue","deepskyblue","dodgerblue4"))
graphics.off()




#### Load input ####

# sample classification
sf           = read.table(samgrup,header = F)
colnames(sf) = c("sample",samclas)
rownames(sf) = sf$sample

# transcript to gene dictionary
di           = read.table(txgdict)
colnames(di) = c("long_transcript_id","gene_id","transcript_length")
di           = di[,1:2]

# gene name dictionary
dign           = read.table(nmgdict,sep = "\t")
colnames(dign) = c("gene_id","gene_name")


# functional mappings: GO, pfam
gomap = readMappings(gomapfi)
panno = read.table(file = panfile)
if (panform == "simple") {
  colnames(panno) = c("transcript","pstart","pend","domain")
} else if (panform == "pfamscan") {
  colnames(panno) = c("transcript","pstart","pend","pfamid","domain","domseq")
}
panno = merge(panno,di, by.x="transcript",by.y="long_transcript_id")
pannu = subset(panno, select=c("transcript","pfamid","domain","gene_id"))
pannu = pannu[!duplicated(pannu), ]


# load gff
message("# Loading input: GFF")
exon_name     = "CDS"
gene_name     = "gene"
mRNA_name     = "mRNA"
gf            = read.gff(gfffile)
gf$attributes = gsub("ID=","",gf$attributes)
gf$attributes = gsub("Parent=","",gf$attributes)
gf$attributes = gsub(";.*","",gf$attributes)
gf            = subset(gf, type %in% c(gene_name))
levels(gf$type)[levels(gf$type)==exon_name] = "exon"
levels(gf$type)[levels(gf$type)==mRNA_name] = "gene"
gf            = gf[order(gf[,"seqid"],gf[,"start"]),]
gf$source     = "R"
gi            = GRanges(
  gf$seqid,
  IRanges(start=gf$start, end=gf$end , names=gf$attributes),
  strand = gf$strand)



#### DESEQ all-to-all ####

# deseq input data
ex_lif        = file.path(saloutf,paste(si,"_",sf$sample,".salmon_long.out",sep=""), "quant.sf")
names(ex_lif) = sf$sample
ex_inp        = tximport(files = ex_lif, type = "salmon", tx2gene = di)
ex_dat        = DESeqDataSetFromTximport(txi=ex_inp,colData=sf,design=formula(expdiss))
ex_dat        = DESeq(ex_dat,test="Wald") # if tximport is used, a normMatrix with avgTxLength is used as normalization factor (which superseeds size factors)

# matrix of normalised counts 
# uses normalization factors: sample size, effective len...
ex_dnc = na.omit(DESeq2::counts(ex_dat,normalized=T))
ex_dnc = log(ex_dnc,10)
ex_dnc = abs(ex_dnc)
ex_dnc[!is.finite(ex_dnc)] = 0
ex_dnc = ex_dnc[apply(ex_dnc,1,sd) != 0,]

# dimensionality reduction
# prepare data: scale along rows
ex_dnc_st = ex_dnc
ex_dnc_st = t(apply(ex_dnc_st, 1, scale))
colnames(ex_dnc_st) = colnames(ex_dnc)

ex_dnc_dimred = dimredfun(
  matriu = ex_dnc_st, outputname = paste(outcode,"de_all",sep=""), varname = "standardized log norm counts", 
  cols_are = "Samples",rows_are = "Genes", isbidi = F,
  cols_dist_method = "pearson", rows_dist_method = "pearson",
  clus_method = "ward.D2")


#### DEseq specific comparisons ####

for (com in complis) {
  
  # define columns with groups
  if (com[5] == "ALL") {
    cli = as.character(sf[sf[,com[1]] %in% com[2],]$sample)
    clj = as.character(sf[sf[,com[1]] %in% com[3],]$sample)
  } else {
    cli = as.character(sf[sf[,com[1]] %in% com[2] & sf[,com[4]] %in% com[5],]$sample)
    clj = as.character(sf[sf[,com[1]] %in% com[3] & sf[,com[4]] %in% com[5],]$sample)
  }
  tici   = paste(com[5],".",com[2],"-",com[3],sep="")
  
  # prepare data for DEseq
  message(paste("# DEseq comparison",tici))
  ex_lif_i = ex_lif[names(ex_lif) %in% c(cli,clj)]
  sf_i     = sf[sf$sample %in% c(cli,clj),]
  sf_i     = droplevels(sf_i)
  ex_inp_i = suppressMessages(tximport(files = ex_lif_i, type = "salmon", tx2gene = di))
  ex_dat_i = suppressMessages(DESeqDataSetFromTximport(txi = ex_inp_i, colData = sf_i, design = formula(com[6])))
  ex_dat_i[[com[1]]] = factor(ex_dat_i[[com[1]]], levels=c(com[3],com[2]))
  ex_dat_i = suppressMessages(DESeq(ex_dat_i, test="Wald"))
  
  # obtain effect size (LFCs), sval and pval
  ex_res_i               = results(ex_dat_i,alpha=alpha_threshold,contrast=com[1:3],name=tici) # unshrunken LFC with pval, FDR-adjusted at alpha
  ex_res_s               = lfcShrink(ex_dat_i,coef=2, type="apeglm", quiet=T, svalue = T)      # shrink LFC using apeglm
  ex_res_i               = as.data.frame(ex_res_i)
  ex_res_i$svalue_shrink = ex_res_s$svalue                                                     # s-value: probability of false signs, ie probability that the effect is actually LFC>0 or LFC<0 (as opposed to pval of probability that effect is 0 or not 0)
  ex_res_i$log2FC_shrink = ex_res_s$log2FoldChange                                             # shrunken LFC values
  
  # matrix of normalized counts
  ex_dnc_i = ex_dnc[,c(cli,clj)]
  ex_dnc_i = ex_dnc_i[apply(ex_dnc_i,1,sd) != 0,]
  ex_yes = order(rowVars(ex_dnc_i),decreasing=T)
  
  # volcano plot & list of over/under expressed genes
  message(paste("# Volcano",tici))
  ex_lis_i = volcanoexp(
    table=ex_res_i,plotname=tici,fileprefix=paste(outcode,"de_co",sep=""),
    pthreshold=alpha_threshold,fc_varname="log2FC_shrink",fcp=com[2],fcn=com[3],pval_varname="padj",
    xlims = c(-10,10),ylims = c(0,-100),minfold = logfc_threshold)
  
  # convert lists of genes to lists of transcripts
  ex_lis_i$tx_sp = as.vector(di[di$gene_id %in% ex_lis_i$genes_sp,]$long_transcript_id)
  ex_lis_i$tx_sn = as.vector(di[di$gene_id %in% ex_lis_i$genes_sn,]$long_transcript_id)
  
  # top genes
  message(paste("# Top genes",tici))
  pdf(file=paste(outcode,"de_co.",tici,"_topgenes.pdf",sep=""),height=12,width=8)
  par(mfrow=c(6,4))
  
  ex_res_i_sortp   = ex_res_i[ex_res_i$padj<pval_threshold,]
  ex_res_i_sortp.p = ex_res_i_sortp[order(ex_res_i_sortp$log2FC_shrink,decreasing = T),]
  ex_res_i_sortp.n = ex_res_i_sortp[order(ex_res_i_sortp$log2FC_shrink),]
  for (geni in 1:24) {
    geni_i = rownames(ex_res_i_sortp.p[geni,])
    geni.n = as.character(dign[dign$gene_id == geni_i,]$gene_name)
    DESeq2::plotCounts(
      ex_dat_i, gene=geni_i, intgroup=com[1], col="slategray",xlab="",
      cex_sub=0.7,cex_axis=0.7,cex_main=0.7,cex_lab=0.7,
      sub=paste("FC=",signif(2^(ex_res_i_sortp.p[geni,]$log2FC_shrink),4)
                ,"| p=",signif(ex_res_i_sortp.p[geni,]$padj,4),"\n",geni.n))
  }
  for (geni in 1:24) {
    geni_i = rownames(ex_res_i_sortp.n[geni,])
    geni.n = as.character(dign[dign$gene_id == geni_i,]$gene_name)
    DESeq2::plotCounts(
      ex_dat_i, gene=geni_i, intgroup=com[1], col="slategray",xlab="",
      cex_sub=0.7,cex_axis=0.7,cex_main=0.7,cex_lab=0.7,
      sub=paste("FC=",signif(2^-(ex_res_i_sortp.n[geni,]$log2FC_shrink),4)
                ,"(^-1) | p=",signif(ex_res_i_sortp.n[geni,]$padj,4),"\n",geni.n))
  }
  dev.off()
  
  # functional enrichment
  message(paste("# Functional enrichment domains",tici))
  hygeofun(list_interest=ex_lis_i$tx_sp,annotation=pannu,gene_col="transcript",ano_col="domain",
           outputname=paste(outcode,"de_co",sep=""),
           name_geneset=paste(tici,"difpos",sep="_"),topnum = 30)
  hygeofun(list_interest=ex_lis_i$tx_sn,annotation=pannu,gene_col="transcript",ano_col="domain",
           outputname=paste(outcode,"de_co",sep=""),
           name_geneset=paste(tici,"difneg",sep="_"),topnum = 30)
  
  message(paste("# Functional enrichment topgo pos",tici))
  suppressMessages(topgofun(list_interest=ex_lis_i$tx_sp,gomap=gomap,
                            ontologyset=c("BP","MF","CC"),tg_test="fisher",tg_algorithm="elim",
                            outputname=paste(outcode,"de_co",sep=""),
                            name_geneset=paste(tici,"difpos",sep="_"),topnum=30))
  message(paste("# Functional enrichment topgo neg",tici))
  suppressMessages(topgofun(list_interest=ex_lis_i$tx_sn,gomap=gomap,
                            ontologyset=c("BP","MF","CC"),tg_test="fisher",tg_algorithm="elim",
                            outputname=paste(outcode,"de_co",sep=""),
                            name_geneset=paste(tici,"difneg",sep="_"),topnum=30))
  
  # save comparison (can get reloaded later)
  save(list=c("ex_dat_i","ex_lis_i","ex_dnc_i","ex_res_i","sf_i"), 
       file=paste(outcode,"de_co.",tici,".data.RData",sep=""))
  
}



#### DE gene overlaps ####

ex_lis_joint = list()

pdf(file=paste(outcode,"de_all.gene_overlaps.pdf",sep=""),height=4,width=4)
for (cop in complis_pairs) {
  
  ex_lis_joini = list("pai","paj")
  # define comparison pair
  pai      = cop[[1]]
  paj      = cop[[2]]
  pai.tici = paste(pai[5],".",pai[2],"-",pai[3],sep="")
  paj.tici = paste(paj[5],".",paj[2],"-",paj[3],sep="")
  paji     = paste(pai.tici,paj.tici,sep=" ")
  
  # load expression data for comparison
  load(paste(outcode,"de_co.",pai.tici,".data.RData",sep=""))
  ex_lis_joini[[1]]$genes_sp = ex_lis_i$genes_sp
  ex_lis_joini[[1]]$genes_sn = ex_lis_i$genes_sn
  load(paste(outcode,"de_co.",paj.tici,".data.RData",sep=""))
  ex_lis_joini[[2]]$genes_sp = ex_lis_i$genes_sp
  ex_lis_joini[[2]]$genes_sn = ex_lis_i$genes_sn
  
  # set analysis for OVEREXPRESSED genes
  ex_lis_joint[[paji]]$overexpressed = 
    venn.two(
      list1=ex_lis_joini[[1]]$genes_sp,
      list2=ex_lis_joini[[2]]$genes_sp,
      catname1=pai.tici,
      catname2=paj.tici,
      main = paste("overexp: ",pai.tici,"&",paj.tici),
      col1 = "green3",
      col2 = "blue3"
    )
  
  # set analysis for OVEREXPRESSED genes
  ex_lis_joint[[paji]]$underexpressed = 
    venn.two(
      list1=ex_lis_joini[[1]]$genes_sn,
      list2=ex_lis_joini[[2]]$genes_sn,
      catname1=pai.tici,
      catname2=paj.tici,
      main = paste("underexp: ",pai.tici,"&",paj.tici),
      col1 = "green3",
      col2 = "blue3"
    )
  
}
dev.off()


#### DE gene overlaps functions ####

pdf(file=paste(outcode,"de_all.gene_overlaps_function.pdf",sep=""), height=1.8*length(complis_pairs[[1]]),width=12)
par(mfrow=c(length(complis_pairs[[1]]),4))
for (cop in complis_pairs) {
  
  # define comparison pair
  pai      = cop[[1]]
  paj      = cop[[2]]
  pai.tici = paste(pai[5],".",pai[2],"-",pai[3],sep="")
  paj.tici = paste(paj[5],".",paj[2],"-",paj[3],sep="")
  paji     = paste(pai.tici,paj.tici,sep=" ")
  
  # overexp shared
  print("Functions shared overexpressed")
  tx_list_fun =  di[di$gene_id %in% ex_lis_joint[[paji]]$overexpressed$list_intersect, "long_transcript_id"]
  hygeofun(tx_list_fun,pannu,"transcript","domain",NA,
           paste("shared overexp, n =",length(tx_list_fun),"genes |"),
           10,printfile = F)
  suppressMessages(
    topgofun(tx_list_fun,gomap,NA,
             paste("shared overexp, n =",length(tx_list_fun),"genes |"),
             c("BP","MF","CC"), "fisher", "elim",10, printfile = F))
  
  # underexp shared
  print("Functions shared underexpressed")
  tx_list_fun =  di[di$gene_id %in% ex_lis_joint[[paji]]$underexpressed$list_intersect, "long_transcript_id"]
  hygeofun(tx_list_fun,pannu,"transcript","domain",NA,
           paste("shared underexp, n =",length(tx_list_fun),"genes |"),
           10,printfile = F)
  suppressMessages(
    topgofun(tx_list_fun,gomap,NA,
             paste("shared overexp, n =",length(tx_list_fun),"genes |"),
             c("BP","MF","CC"), "fisher", "elim",10,printfile = F))
  
}
dev.off()


### FAS heatmaps ###

# function for plotting expression levels across samples with significance of DE
# (It inherits variables from this script, it's not entirely standalone)
source("helper_scripts/interest_genes_heatmap.R")


# FA synthesis genes (defined from Pfam domains)
interestlist = as.character(pannu[ pannu$domain == "ketoacyl-synt" ,]$gene_id)
interest_genes_heatmap(interestlist = interestlist, interestname = "FAsynthase_ketoacyl-synt",source_matrix = ex_dnc)

interestlist = as.character(pannu[ pannu$domain == "ELO" ,]$gene_id)
interest_genes_heatmap(interestlist = interestlist, interestname = "FAelongase_ELO",source_matrix = ex_dnc)

interestlist = as.character(pannu[ pannu$domain == "FA_desaturase" ,]$gene_id)
interestlist = c(interestlist, "Anogam_AGAP003050") # added manually because of Pfam misannotation
interest_genes_heatmap(interestlist = interestlist, interestname = "FAdesaturase_FA_desaturase",source_matrix = ex_dnc)

interestlist = as.character(pannu[ pannu$domain == "NAD_binding_4" ,]$gene_id)
interest_genes_heatmap(interestlist = interestlist, interestname = "FAreductase_NAD_binding_4",source_matrix = ex_dnc)

interestlist = as.character(pannu[ pannu$domain == "p450" ,]$gene_id)
interest_genes_heatmap(interestlist = interestlist, interestname = "FAdecarboxyl_p450",source_matrix = ex_dnc)



#### Save ####

message("# Save")

ex_res_tot = data.frame()
for (com in complis) {
  
  # define columns with groups
  if (com[5] == "ALL") {
    cli = as.character(sf[sf[,com[1]] %in% com[2],]$sample)
    clj = as.character(sf[sf[,com[1]] %in% com[3],]$sample)
  } else {
    cli = as.character(sf[sf[,com[1]] %in% com[2] & sf[,com[4]] %in% com[5],]$sample)
    clj = as.character(sf[sf[,com[1]] %in% com[3] & sf[,com[4]] %in% com[5],]$sample)
  }
  
  # load expression data for comparison
  tici = paste(com[5],".",com[2],"-",com[3],sep="")
  load(paste(outcode,"de_co.",tici,".data.RData",sep=""))
  
  ex_dnc_f = na.omit(DESeq2::counts(ex_dat,normalized=T))
  
  # add columns
  ex_res_i$gene          = rownames(ex_res_i)
  ex_res_i$comparison    = tici
  ex_res_i$samples_i     = com[2]
  ex_res_i$normcounts_i  = rowMeans(ex_dnc_f[,cli],na.rm = T)
  ex_res_i$samples_j     = com[3]
  ex_res_i$normcounts_j  = rowMeans(ex_dnc_f[,clj],na.rm = T)
  ex_res_i               = merge(ex_res_i,dign,by.x="gene",by.y="gene_id")
  ex_res_i$is_signif_pos = ex_res_i$gene %in% ex_lis_i$genes_sp
  ex_res_i$is_signif_neg = ex_res_i$gene %in% ex_lis_i$genes_sn
  ex_res_i$is_signif_any = ex_res_i$is_signif_pos | ex_res_i$is_signif_neg
  
  write.table(ex_res_i,file=paste(outcode,"session.deseq_difexp",tici,"csv",sep=""),
              quote = T,row.names = F,sep = "\t")
  
  ex_res_tot             = rbind(ex_res_tot,ex_res_i)
  
}

write.table(ex_res_tot,file=paste(outcode,"session.deseq_difexp.csv",sep=""),
            quote = T,row.names = F,sep = "\t")

save.image(paste(outcode,"session.deseq_all.RData",sep=""))

message("\n\n### FI! ###\n\n")



stop("FINISHED WITHOUT PROBLEMS")





# list of genes expressed in each sample group

expr_threshold=10

for (exi in c(0,1,10,20,50,100)){
for (coi in as.character(unique(sf$group))) {
  
  cli = as.character(sf[sf$group %in% coi,]$sample)
  ex_dnc_cli = na.omit(DESeq2::counts(ex_dat,normalized=T))
  ex_dnc_cli = ex_dnc_cli[,cli]
  ex_dnc_cli = ex_dnc_cli[rowMeans(ex_dnc_cli)>exi,] 
  print(paste(coi,exi,nrow(ex_dnc_cli)))
  
}
}


for (coi in c("OecyF","OecyM")) {
  
  cli = as.character(sf[sf$group %in% coi,]$sample)
  
  ex_dnc_cli = na.omit(DESeq2::counts(ex_dat,normalized=T))
  ex_dnc_cli = ex_dnc_cli[,cli]
  ex_dnc_cli = ex_dnc_cli[rowMeans(ex_dnc_cli)>expr_threshold,]
  print(paste(coi,nrow(ex_dnc_cli),nrow(ex_dnc)))
  interestname = coi
  interestlist = rownames(ex_dnc_cli)
  interestlist_tx = as.character(di[di$gene_id %in% interestlist, "long_transcript_id"])
  
  write.table(interestlist,
              file=paste("genesubset.samples",coi,"min",expr_threshold,"csv",sep="."),
              quote = F, row.names = F, col.names = F)
  
  hygeofun(list_interest=interestlist,annotation=pannu,gene_col="gene_id",ano_col="domain",
           outputname="genesubset.samples",
           name_geneset=paste(interestname,sep="_"),topnum = 30)
  suppressMessages(topgofun(list_interest=interestlist_tx,
                            gomap=gomap,
                            ontologyset=c("BP","MF","CC"),tg_test="fisher",tg_algorithm="elim",
                            outputname="genesubset.samples",
                            name_geneset=paste(interestname,sep="_"),topnum=30))
  
}


pdf(file=paste("genesubset.samples.overlap_min10_OecyMF.pdf",sep="."),height=4,width=4)
list1=read.table("genesubset.samples.OecyF.min.10.csv",header = F)
list2=read.table("genesubset.samples.OecyM.min.10.csv",header = F)
bondia = venn.two(
  list1 = list1$V1,eulerbool = T,
  list2 = list2$V1,
  catname1 = "OecyF",catname2 = "OecyM",main = "genes expressed in Oecy",
  print.mode = T)
dev.off()

