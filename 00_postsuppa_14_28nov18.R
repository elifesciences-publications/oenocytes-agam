# load libraries
library(GenomicFeatures)
library(ape)
library(stringr)
library(plyr)
library(tidyr)
library(topGO)
library(Biostrings)
library(Rsamtools)
library(Rfast)
library(cluster)
library(gplots)

### Define input ####

# input files
setwd("/home/xavi/Documents/Transvar/suppa_Anogam_Oenocytes/")

si          = "Anogam"                                  # Species name (requires Spi_long.annot.gff and Spi_gDNA.fasta)
suppaioef   = "~/dades/Transcriptomes/Anogam_Oenocytes_6nov18/Anogam_suppa.8nov18.ioe" # IOE suppa file
suppaioif   = "~/dades/Transcriptomes/Anogam_Oenocytes_6nov18/Anogam_suppa.8nov18.ioi" # IOE suppa file
exprtpmf    = "~/dades/Transcriptomes/Anogam_Oenocytes_6nov18/salmon/Anogam_salmon.8nov18_long.tpm"
outcode     = "bll"                                    # Output code SUPPA analysis (eg *14mai18*, used for retrieving suppa groupings!)
samgrup     = "samples_expanded.class"                 # sample list (first col) & grouping (2nd, 3rd, etc.)
grupnoms    = c("group","sex","tissue")              # group names (2nd, 3rd, 4th... cols)
genomefile  = "~/dades/Genomes/Anogam_gDNA.fasta"
gfffile     = "~/dades/Genomes/Anogam_long.annot.gff"
annottab    = "~/dades/Genomes/anotacions/Anogam_long.pep_Pfamscan.seqs"               # pfam annots (requires .blat file too)
gomapfile   = "~/dades/Genomes/anotacions/Anogam_eggnog_dipNOG.emapper.annotations.GO" # gomap (with transcripts)
tx2genedict = "~/dades/Genomes/Anogam_long.cds.csv"                                    # transcript-to-gene mappings
list_evtype = c("SE","MX","RI")
list_analisi= c("LOADIN","DOMMAP","GENERAL","VEPCALL","VEPPOST","QUANT","CLUST",
                "ENRICH","SAVE","NEWCOMP")
annotmap    = paste(annottab,".blat",sep="")
outprefix   = basename(suppaioef)

complis = list(
  convec=c("tissue","Oecy","Bulk","sex","Fem","~tissue"),
  convec=c("tissue","Oecy","Bulk","sex","Mal","~tissue"),
  convec=c("sex","Fem","Mal","tissue","Oecy","~sex")
)


# input variables
pthr        = 0.01         # pval threshold for suppa comparisons
k_list      = c(3,4,5,6,7) # list of k values to try for kmeans clustering

# load functions
source("~/Dropbox/Scripts/bin/fastaFromDf.R")
source("~/Dropbox/Scripts/bin/geneSetAnalysis.R")
source("~/Dropbox/Scripts/bin/dimred_18jul18.R")
source("~/Dropbox/Scripts/bin/pval_from_CDF_8nov18.R")
col.fun = colorRampPalette(interpolate="l",c("aliceblue","deepskyblue","dodgerblue4"))

# plotting function
topdiffgenes = function(table,plotname,ntopgenes,dpsi_pos,dpsi_neg) {
  table = head(table,ntopgenes)
  table = droplevels(table)
  # this block simplifies VEP effects for clearer plotting
  csi = c("inframe_deletion","inframe_insertion","frameshift_variant","stop_gained","stop_lost")
  table = table[grepl(pattern = paste(csi,collapse="|"), x = table$event_consequences),]
  csj = unique(unlist(base::strsplit(paste(table$event_consequences,sep=" "),split = ",")))
  csj = subset(csj, !(csj %in% csi))
  csj = c(csj,",")
  table$event_consequences = gsub(pattern = paste(csj,collapse="|"), replacement = "",x= table$event_consequences)
  table$event_consequences = as.factor(paste(table$event_consequences," dom",table$hasdomain,sep=""))
  # colors  & spacing
  virc  = rainbow(nlevels(table$event_consequences),v = 0.8)
  par(mar=c(5,20,4,2)+0.1)
  # plot
  bp = barplot(
    height = table$dPSI,names.arg = table$event_id,las=2,col=virc[table$event_consequences],border=NA,
    cex.names=0.6,cex.axis=0.7,cex.sub=0.7,cex.lab=0.7,main = paste(plotname),sub  = paste("dPSI<0:",dpsi_neg,"| dPSI>0:",dpsi_pos),
    xlim = c(-1,1),horiz = T,xlab = "dPSI")
  legend("topleft", legend = levels(table$event_consequences),fill=virc,cex=0.6,border=NA)
  text(bp,x=1, pos=2,cex=0.6,col="red",labels = paste("p=",formatC(table$pval,format="e",digits=2),sep=""))
}




#### Load input ####
if ("LOADIN" %in% list_analisi) {
  
  # load sample list
  message("# Sample list & groups")
  sf           = read.table(samgrup,header = F)
  colnames(sf) = c("sample",grupnoms)
  
  # load GOmappings, dictionary
  message("# Loading input: GO, dict")
  gomap        = readMappings(gomapfile)
  di           = read.table(tx2genedict)
  colnames(di) = c("long_transcript_id","gene_id","long_transcript_len")
  
  # matrix of TPMs per event (translated from TPMs per transcript!)
  message("# Loading input: expression (TPM)")
  expr_tpm                    = read.table(exprtpmf)
  expr_tpm$long_transcript_id = rownames(expr_tpm)
  
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
  
  
  # Format IOE to iot
  message("# Loading input: IOE")
  ioe          = read.table(suppaioef,header = T)
  iot          = ioe[ioe$seqname!="seqname",]
  iot          = merge(iot,di, by="gene_id")
  iot$event_id = as.vector(iot$event_id)
  iot$iow      = as.vector(base::strsplit(gsub(":|-",";",iot$event_id), split=";"))
  iot$ev_type  = as.factor(sapply(iot$iow, "[", 2))
  # remove unneeded AS types (AS/AL do not agree with reference!)
  iot          = iot[iot$ev_type %in% list_evtype,]
  iot          = droplevels(iot)
  
  # iot event coords
  iot$ev_start  = NA
  iot$ev_end    = NA
  iot$ev_startB = NA
  iot$ev_endB   = NA
  iot$ev_strand = NA
  
  # strand: always in the last field of iow string
  # "-" is encoded as EMPTY; can be corrected later
  iot[iot$ev_type == "SE",]$ev_strand = as.character(sapply(iot[iot$ev_type == "SE",]$iow, "[", 8))
  iot[iot$ev_type == "RI",]$ev_strand = as.character(sapply(iot[iot$ev_type == "RI",]$iow, "[", 8))
  iot[iot$ev_type == "A5",]$ev_strand = as.character(sapply(iot[iot$ev_type == "A5",]$iow, "[", 8))
  iot[iot$ev_type == "A3",]$ev_strand = as.character(sapply(iot[iot$ev_type == "A3",]$iow, "[", 8))
  iot[iot$ev_type == "AF",]$ev_strand = as.character(sapply(iot[iot$ev_type == "AF",]$iow, "[", 10))
  iot[iot$ev_type == "AL",]$ev_strand = as.character(sapply(iot[iot$ev_type == "AL",]$iow, "[", 10))
  iot[iot$ev_type == "MX",]$ev_strand = as.character(sapply(iot[iot$ev_type == "MX",]$iow, "[", 12))
  iot[iot$ev_strand=='',]$ev_strand = "-"
  
  # Coordinates of event
  # SE: e1-S2:E2-s3
  iot[iot$ev_type == "SE",]$ev_start  = as.numeric(sapply(iot[iot$ev_type == "SE",]$iow, "[", 5))
  iot[iot$ev_type == "SE",]$ev_end    = as.numeric(sapply(iot[iot$ev_type == "SE",]$iow, "[", 6))
  # RI: s1-E1:S2-e2
  iot[iot$ev_type == "RI",]$ev_start  = as.numeric(sapply(iot[iot$ev_type == "RI",]$iow, "[", 5))
  iot[iot$ev_type == "RI",]$ev_end    = as.numeric(sapply(iot[iot$ev_type == "RI",]$iow, "[", 6))
  
  # A5+: E2-s3:E1-s3
  iot[iot$ev_type == "A5" & iot$ev_strand == "+",]$ev_start  = as.numeric(sapply(iot[iot$ev_type == "A5" & iot$ev_strand == "+",]$iow, "[", 6))
  iot[iot$ev_type == "A5" & iot$ev_strand == "+",]$ev_end    = as.numeric(sapply(iot[iot$ev_type == "A5" & iot$ev_strand == "+",]$iow, "[", 4))
  # A5-: E2-s3:E1-s3
  iot[iot$ev_type == "A5" & iot$ev_strand == "-",]$ev_start  = as.numeric(sapply(iot[iot$ev_type == "A5" & iot$ev_strand == "-",]$iow, "[", 5))
  iot[iot$ev_type == "A5" & iot$ev_strand == "-",]$ev_end    = as.numeric(sapply(iot[iot$ev_type == "A5" & iot$ev_strand == "-",]$iow, "[", 7))
  
  # A3+: e1-S2:e2-S3
  iot[iot$ev_type == "A3" & iot$ev_strand == "+",]$ev_start  = as.numeric(sapply(iot[iot$ev_type == "A3" & iot$ev_strand == "+",]$iow, "[", 5))
  iot[iot$ev_type == "A3" & iot$ev_strand == "+",]$ev_end    = as.numeric(sapply(iot[iot$ev_type == "A3" & iot$ev_strand == "+",]$iow, "[", 7))
  # A3-: e1-S2:e2-S3
  iot[iot$ev_type == "A3" & iot$ev_strand == "-",]$ev_start  = as.numeric(sapply(iot[iot$ev_type == "A3" & iot$ev_strand == "-",]$iow, "[", 6))
  iot[iot$ev_type == "A3" & iot$ev_strand == "-",]$ev_end    = as.numeric(sapply(iot[iot$ev_type == "A3" & iot$ev_strand == "-",]$iow, "[", 4))
  
  # AF: s1:e1-s3:s2:e2-s3
  iot[iot$ev_type == "AF",]$ev_start  = as.numeric(sapply(iot[iot$ev_type == "AF",]$iow, "[", 4)) 
  iot[iot$ev_type == "AF",]$ev_end    = as.numeric(sapply(iot[iot$ev_type == "AF",]$iow, "[", 5))
  iot[iot$ev_type == "AF",]$ev_startB = as.numeric(sapply(iot[iot$ev_type == "AF",]$iow, "[", 7)) 
  iot[iot$ev_type == "AF",]$ev_endB   = as.numeric(sapply(iot[iot$ev_type == "AF",]$iow, "[", 8))
  # AL: e1-s2:e2:e2-s3:e3
  iot[iot$ev_type == "AL",]$ev_start  = as.numeric(sapply(iot[iot$ev_type == "AL",]$iow, "[", 5)) 
  iot[iot$ev_type == "AL",]$ev_end    = as.numeric(sapply(iot[iot$ev_type == "AL",]$iow, "[", 6))
  iot[iot$ev_type == "AL",]$ev_startB = as.numeric(sapply(iot[iot$ev_type == "AL",]$iow, "[", 8)) 
  iot[iot$ev_type == "AL",]$ev_endB   = as.numeric(sapply(iot[iot$ev_type == "AL",]$iow, "[", 9))
  # MX: e1-s2:e2-s4:e1-s3:e3-s4
  iot[iot$ev_type == "MX",]$ev_start  = as.numeric(sapply(iot[iot$ev_type == "MX",]$iow, "[", 5)) 
  iot[iot$ev_type == "MX",]$ev_end    = as.numeric(sapply(iot[iot$ev_type == "MX",]$iow, "[", 6))
  iot[iot$ev_type == "MX",]$ev_startB = as.numeric(sapply(iot[iot$ev_type == "MX",]$iow, "[", 9))
  iot[iot$ev_type == "MX",]$ev_endB   = as.numeric(sapply(iot[iot$ev_type == "MX",]$iow, "[", 10))
  
  iot$ismultiexonic = !is.na(iot$ev_endB) 
  iot = subset(iot, select = -c(iow))
  
  # format iot 3n event data
  iot$ev_length  = iot$ev_end - iot$ev_start + 1
  iot$ev_lengthB = iot$ev_endB - iot$ev_startB + 1
  
  iot$is3n  = iot$ev_length%%3L==0L
  iot$is3nB = iot$ev_lengthB%%3L==0L
  
  # format iot to GRanges
  iog = GRanges(iot$seqname,IRanges(
    start=iot$ev_start,end=iot$ev_end,names=iot$event_id),
    strand=iot$ev_strand)
  
  ioeB = iot[iot$ismultiexonic==T,]
  iogB = GRanges(ioeB$seqname,IRanges(
    start=ioeB$ev_startB,end=ioeB$ev_endB,names=ioeB$event_id),
    strand=ioeB$ev_strand)
  
  # Load & format pfam annotations
  message("# Loading input: Pfam & Blat")
  pt            = read.table(file = annottab)
  colnames(pt)  = c("transcript","pstart","pend","pfamid","domain","domseq")
  pt$seqname    = paste(pt$transcript,pt$domain,pt$pstart,pt$pend,sep="|")
  pt$tstart     = pt$pstart * 3 - 2
  pt$tend       = pt$pend   * 3
  pannu         = subset(pt, select=c("transcript","pfamid","domain"))
  pannu         = pannu[!duplicated(pannu), ]
  
  am            = read.table(annotmap,header = F)
  colnames(am)  = c("annot_id","chrom","pident","alilen","mismatch","gapop",
                    "an_start","an_end","ch_start","ch_end","eval","bitscore")
  am            = am[am$pident==100,]
  am$string     = as.vector(base::strsplit(as.vector(am$annot_id), split='|', fixed=T))
  am$domain     = as.factor(sapply(am$string, "[", 4))
  am$do_start   = as.factor(sapply(am$string, "[", 2))
  am$do_end     = as.factor(sapply(am$string, "[", 3))
  am$transcript = as.factor(sapply(am$string, "[", 1))
  am$strand     = ifelse((am$ch_end - am$ch_start)>0,"+","-")
  am$ch_start_S = ifelse((am$ch_end - am$ch_start)>0,am$ch_start,am$ch_end)
  am$ch_end_S   = ifelse((am$ch_end - am$ch_start)>0,am$ch_end,am$ch_start)
  
  pi = GRanges(
    am$chrom,
    IRanges(start=am$ch_start_S, end=am$ch_end_S,
            names=paste(am$domain,am$do_start,am$do_end,sep="|")),
    strand=am$strand
  )
  
} else { message("# Skip LOADIN!") }



#### Map pfam to gff ####
if ("DOMMAP" %in% list_analisi) {
  
  # long GFF mapped to pfam annotations (overlap-i)
  message("# Domain overlaps to GFF-I")
  ovi = findOverlapPairs(subject=gi,query=pi,type="within",
                         select="all",ignore.strand=F)
  oti = unique.data.frame(data.frame(seqnames=seqnames(ovi@first),
                                     doch_start=start(ovi@first),
                                     doch_end=end(ovi@first),
                                     domain=names(ovi@first@ranges),
                                     gene=names(ovi@second@ranges),
                                     doch_strand=strand(ovi@first))
  )
  oti$ismultiexonic = !Biobase::isUnique(oti$domain)
  write.table(oti,file=paste(annottab,".GFFoverlap.bed",sep=""),
              row.names = F,col.names = F,quote=F,sep = "\t")
  ogi = GRanges(seqnames = oti$seqnames,
                IRanges(start=oti$doch_start,
                        end  =oti$doch_end,
                        names=oti$domain),
                strand=oti$doch_strand)
  values(ogi) = oti$ismultiexonic
  
  # iot events mapped to pfam annotations (overlap-j)
  message("# Domain overlaps to IOT-J")
  ovj = findOverlapPairs(subject=iog,query=ogi,
                         select="all",ignore.strand=F)
  otj = data.frame(
    seqnames  = seqnames(ovj@second),
    ev_start  = start(ovj@second)-1,
    ev_end    = end(ovj@second),
    ev_strand = strand(ovj@second),
    event_id  = names(ovj@second@ranges),
    do_start  = start(ovj@first)-1,
    do_end    = end(ovj@first),
    do_strand = strand(ovj@first),
    domain    = names(ovj@first@ranges),
    X         = ovj@first@elementMetadata@listData
  )
  colnames(otj)[10] <- "multiex_domain"
  write.table(otj,file=paste(annottab,".ASoverlap.bed",sep=""),
              row.names = F,col.names = F,quote=F,sep = "\t")
  otj = subset(otj, select = -c(ev_start,ev_end,ev_strand,seqnames,do_strand))
  
  # iot events-B mapped to pfam annotations (overlap-k)
  message("# Domain overlaps to IOT-K")
  ovk = findOverlapPairs(subject=iogB,query=ogi,
                         select="all",ignore.strand=F)
  otk = data.frame(
    seqnames  = seqnames(ovk@second),
    ev_start  = start(ovk@second)-1,
    ev_end    = end(ovk@second),
    ev_strand = strand(ovk@second),
    event_id  = names(ovk@second@ranges),
    do_startB = start(ovk@first)-1,
    do_endB   = end(ovk@first),
    do_strand = strand(ovk@first),
    domainB   = names(ovk@first@ranges),
    X         = ovk@first@elementMetadata@listData
  )
  colnames(otk)[10] <- "multiex_domainB"
  write.table(otk,file=paste(annottab,".ASoverlap.bed",sep=""),
              row.names = F,col.names = F,quote=F,sep = "\t")
  otk = subset(otk, select = -c(ev_start,ev_end,ev_strand,seqnames,do_strand))
  
  # merge with main file
  otl = merge(otj,otk, by="event_id",all.x=T,all.y=T)
  otl = otl[!duplicated(otl$event_id), ]
  iot = merge(iot,otl, by="event_id",all.x=T,all.y=F)
  iot = iot[!duplicated(iot$event_id), ]
  
  iot$hasdomain       = !is.na(iot$domain) 
  iot$hasdomainB      = !is.na(iot$domainB) 
  iot$hasdomainAny    = iot$hasdomain | iot$hasdomainB
  iot$hasdomcomplete  = !iot$multiex_domain & iot$hasdomain
  iot$hasdomcompleteB = !iot$multiex_domainB & iot$hasdomainB
  iot$hasdomcompleteB = (iot$do_startB - iot$ev_startB >= 0) & (iot$do_endB - iot$ev_endB >= 0)
  
} else { message("# Skip DOMMAP!") }

# remove uneeded objects
remove(list=c("ioeB","iogB","ogi","oti","otj","otk","otl","ovi","ovj","ovk"))



#### General plots ####

if ("GENERAL" %in% list_analisi) {
  
  pdf(file=paste(outcode,".general.counts.pdf",sep=""),height=4,width=4.5)
  
  # plot number of events affected by each AS type
  counts = table(iot$ev_type)
  bp=barplot(counts,las=1,col="slategray3",border=NA,
             main="# AS events, per type",
             sub=paste("total AS events =",nrow(iot)),
             ylab="# events")
  text(x=bp,y = 0,pos=3,labels=paste("n =",c(counts)),cex=0.7)
  
  # plot number of genes affected by each AS type
  uniquegeneev = unique(iot[,c("ev_type","gene_id")])
  counts = table(uniquegeneev$ev_type)
  barplot(counts,las=1,col="slategray3",border=NA,
          main="# genes affected by AS, per type",
          sub=paste("total genes affected by AS =",length(unique(iot$gene_id))),
          ylab="# genes")
  text(x=bp,y = 0,pos=3,labels=paste("n =",c(counts)),cex=0.7)
  
  # plot: events/gene
  uniquegeneev = unique(iot[,c("ev_type","gene_id")])
  counts = table(iot$ev_type)/table(uniquegeneev$ev_type)
  barplot(counts,las=1,col="slategray3",border=NA,
          main="# events/gene, per type",
          sub=paste("total genes affected by AS =",length(unique(iot$gene_id))),
          ylab="# events/gene")
  text(x=bp,y = 0,pos=3,labels=paste("n =",c(signif(counts,3))),cex=0.7)
  
  dev.off()
  
} else { message("# Skip GENERAl!") }



#### VEP: events as pseudoindels ####

if ("VEPCALL" %in% list_analisi) {
  
  # Prepare pseudoindels
  message("# Prepare VEP pseudoindels")
  fa = FaFile(genomefile)
  fi = getSeq(fa)
  
  psi = data.frame(
    seqname=iot$seqname,
    in_start=NA,
    in_end=NA,
    string=NA,
    in_strand=iot$ev_strand,
    event_id=iot$event_id
  )
  
  # SE: codified as deletion (XXX/-) that starts and ends where the exon should be
  acceptable = iot$ev_type == "SE"
  psi[acceptable,]$in_start  = iot[acceptable,]$ev_start
  psi[acceptable,]$in_end    = iot[acceptable,]$ev_end
  psi[acceptable,]$string    = paste(
    as.vector(subseq(x=    fi[psi[acceptable,]$seqname],
                     start=psi[acceptable,]$in_start,
                     end=  psi[acceptable,]$in_end)),
    "/-",
    sep=""
  )
  psi[acceptable,]$string    = paste(psi[acceptable,]$string,"/-",sep="")
  
  # RI, A3, A5: codified as insertions (-/XXX) that start 
  # at ev_start-1, and end at ev_start (reversed format!)
  acceptable = 
    iot$ev_type == "RI" | 
    iot$ev_type == "A3" | 
    iot$ev_type == "A5"
  
  psi[acceptable,]$in_start  = iot[acceptable,]$ev_start+2
  psi[acceptable,]$in_end    = iot[acceptable,]$ev_start+1
  psi[acceptable,]$string    = paste(
    "-/",
    as.vector(subseq(x =   fi[psi[acceptable,]$seqname],
                     start=iot[acceptable,]$ev_start,
                     end=  iot[acceptable,]$ev_end-1)),
    sep=""
  )
  
  # MX, AF, AL: codified as complex events (YY/XXX) that start 
  # at ev_start-1, and end at ev_start (like insertions - reversed!)
  acceptable = 
    iot$ev_type == "MX" | 
    iot$ev_type == "AF" | 
    iot$ev_type == "AL"
  
  psi[acceptable,]$in_start  = iot[acceptable,]$ev_start+1
  psi[acceptable,]$in_end    = iot[acceptable,]$ev_start
  psi[acceptable,]$string    = paste(
    "-/",
    as.vector(subseq(x =   fi[psi[acceptable,]$seqname],
                     start=iot[acceptable,]$ev_start,
                     end=  iot[acceptable,]$ev_end)),
    "/",
    as.vector(subseq(x =   fi[psi[acceptable,]$seqname],
                     start=iot[acceptable,]$ev_startB,
                     end=  iot[acceptable,]$ev_endB)),
    sep=""
  )
  
  # format VEP pseudoindel file
  psi = na.omit(psi)
  write.table(psi,file=paste(outcode,".vep.csv",sep=""),
              row.names = F,col.names = F, quote = F, sep = "\t")
  
  # format annotation & genome
  message("# Prepare VEP input")
  system(paste(
    "bash ~/Dropbox/Scripts/bin/formatGTF_fromlongff_VEP.sh",
    si,
    gfffile,
    genomefile)
  )
  
  # Call VEP
  message("# VEP call")
  system(paste(
    "/home/xavi/Programes/ensembl-vep/vep",
    " -i ",paste(outcode,".vep.csv",sep=""),
    " -o ",paste(outcode,".vep.out",sep=""),
    " --gtf ",si,"_vep_annot.gtf.gz",
    " --fasta ",si,"_vep_gDNA.fasta.gz ",
    " --plugin Downstream",
    " --plugin SpliceRegion",
    " --plugin MaxEntScan,/home/xavi/Programes/ensembl-vep/plugins_meus/maxentscan/",
    " --force_overwrite",
    " --fork 10",
    sep="")
  )
  
  # Remove VEP predictions that don't make sense when SUPPA output is compared
  # to the reference sequence annotation:
  # - transcript_ablation : only appears when a nested gene is incorrectly mixed with the larger gene, and then skipped
  # - downstream_gene_variant, 
  #   upstream_gene_variant, 
  #   intergenic_variant : entirely out of the CDS
  # - intron_variant : skipping of a middle exon that doesn't exist in the reference
  vepout = paste(outcode,".vep.out",sep="")
  system(
    paste(
      "grep -v \"transcript_ablation\\|stop_lost\\|downstream_gene_variant\\|intergenic_variant\\|upstream_gene_variant\" ",
      vepout,"| sed \"s/intron_variant,//\" | grep -v \"intron_variant\"",
      ">",paste(vepout,"-TMP",sep="")))
  system(
    paste(
      "mv",paste(vepout,"-TMP",sep=""),vepout
    ))
  
  message("# Load VEP output")
  psj = read.table(file=paste(outcode,".vep.out",sep=""),
                   comment.char="#")
  
  colnames(psj) = c("event_id","coords_indel","indel_seq",
                    "indel_gene","indel_transcript","indel_trtype",
                    "event_consequences","indel_cdna_pos","indel_cds_pos",
                    "indel_pep_pos","indel_pep_seq","indel_codon_seq",
                    "indel_variation","indel_metadata")
  
  psj = psj[c("event_id","indel_seq","indel_gene","event_consequences",
              "indel_pep_pos","indel_cds_pos","indel_metadata")]
  
  iot = merge(iot, psj, by.x = "event_id", by.y = "event_id")
  iot = iot[!duplicated(iot$event_id), ]
  
  # is inframe?
  csi = c("inframe_deletion","inframe_insertion","frameshift_variant")
  csj = unique(unlist(base::strsplit(paste(iot$event_consequences,sep=" "),split = ",")))
  csj = subset(csj, !(csj %in% csi))
  csj = c(csj,",")
  iot$frame_consequences = gsub(pattern = paste(csj,collapse="|"), 
                                replacement = "",
                                x=iot$event_consequences)
  iot$is_inframe         = grepl(pattern="inframe",iot$frame_consequences)
  
  # remove uneeded heavy objects
  remove(list = c("fi","fa"))
  
} else { message("# Skip VEPCALL!") }



#### VEP: postanalysis ####

if ("VEPPOST" %in% list_analisi) {
  
  # plot function
  propbarplot = function(title, table) {
    table = t(table(table$ev_type, table$event_consequences))
    iotr = rowSums(table)
    iots = colSums(table)
    # table = prop.table(table, 2)
    iotS = sum(iots)
    cols = rainbow(length(rownames(table)),v = 0.8)
    
    if (iotS>0){
      bp = barplot(table, main=paste(title,", n=",iotS,sep=""),
                   col=cols,
                   border = "white",
                   las=2,horiz = T) 
      text(bp, x=0,pos = 4, col="white",
           labels = paste("n =",iots))
      legend("topright", legend = paste(rownames(table),", n=",iotr,sep=""),
             cex = 0.4,fill = cols)
    }
  }
  
  # VEP effects
  pdf(file=paste(outcode,".vep.evconseq.pdf",sep=""),height=4.5,width=12)
  
  # VEP effects: main effects
  message(paste("# VEP consequences: high-level"))
  csi = c("inframe_deletion","inframe_insertion","frameshift_variant","stop_gained","stop_lost")
  ioti = iot[grepl(pattern = paste(csi,collapse="|"), x = iot$event_consequences),]
  csj = unique(unlist(base::strsplit(paste(ioti$event_consequences,sep=" "),split = ",")))
  csj = subset(csj, !(csj %in% csi))
  csj = c(csj,",")
  ioti$event_consequences = gsub(pattern = paste(csj,collapse="|"), 
                                 replacement = "",
                                 x= ioti$event_consequences)
  propbarplot("Main effects", ioti)
  
  # VEP effects: effects on frame
  message(paste("# VEP consequences: reading frame"))
  csi = c("inframe_deletion","inframe_insertion","frameshift_variant")
  ioti = iot[grepl(pattern = paste(csi,collapse="|"), x = iot$event_consequences),]
  csj = unique(unlist(base::strsplit(paste(ioti$event_consequences,sep=" "),split = ",")))
  csj = subset(csj, !(csj %in% csi))
  csj = c(csj,",")
  ioti$event_consequences = gsub(pattern = paste(csj,collapse="|"), 
                                 replacement = "",
                                 x= ioti$event_consequences)
  propbarplot("Main effects", ioti)
  
  # VEP effects: effects with stop
  message(paste("# VEP consequences: stops"))
  csi = c("stop_gained","stop_lost")
  ioti = iot[grepl(pattern = paste(csi,collapse="|"), x = iot$event_consequences),]
  csj = unique(unlist(base::strsplit(paste(ioti$event_consequences,sep=" "),split = ",")))
  csj = subset(csj, !(csj %in% csi))
  csj = c(csj,",")
  ioti$event_consequences = gsub(pattern = paste(csj,collapse="|"), 
                                 replacement = "",
                                 x= ioti$event_consequences)
  propbarplot("Stop?", ioti)
  
  # VEP effects: all, per event type
  message("# VEP consequences on events")
  for (csi in list_evtype) {
    ioti = iot[iot$ev_type == csi,]
    ioti$event_consequences = droplevels(ioti$event_consequences)
    propbarplot(csi, ioti)
  }
  
  # VEP effects: subeffects
  message("# VEP consequences within main effects, per evtype")
  for (csi in c("inframe_deletion",
                "inframe_insertion",
                "frameshift_variant",
                "stop_gained",
                "stop_lost")) {
    ioti = iot[grepl(pattern = csi, x = iot$event_consequences),]
    ioti$event_consequences = droplevels(ioti$event_consequences)
    propbarplot(csi, ioti)
  }
  
  dev.off()
  
} else { message("# Skip VEPPOST!") }



# Save session
if ("SAVE" %in% list_analisi) {
  message("# Save")
  save.image(paste(outcode,"session.iot","RData",sep="."))
} else { message("# Skip SAVE!") }





#### Functional enrichment ####

if ("ENRICH" %in% list_analisi) {
  
  # inframe SE
  message("# functional enrichments: SE inframe")
  geneset = unique(as.vector(iot[iot$ev_type=="SE" & iot$is_inframe == T,]$long_transcript_id))
  if (length(geneset)>0){
    hygeofun(geneset,pannu,"transcript","domain",paste(outcode,"fun",sep="."),"SEinf",20)
    suppressMessages(topgofun(geneset,gomap,paste(outcode,"fun",sep="."),"SEinf", c("BP","MF","CC"), "fisher", "elim",20))
  }
  
  # RI
  message("# functional enrichments: RI all")
  geneset = unique(as.vector(iot[iot$ev_type=="RI" ,]$long_transcript_id))
  if (length(geneset)>0){
    hygeofun(geneset,pannu,"transcript","domain",paste(outcode,"fun",sep="."),"RIall",20)
    suppressMessages(topgofun(geneset,gomap,paste(outcode,"fun",sep="."),"RIall", c("BP","MF","CC"), "fisher", "elim",20))
  }
  
} else { message("# Skip ENRICH!") }


#### Quantitative profile ####

if ("QUANT" %in% list_analisi) {
  
  # Load & format PSI table
  message("# Loading IOE-PSI")
  ioepsii               = read.table(file=paste(suppaioef,"_events.psi",sep=""),header = T)
  ioepsii               = subset(ioepsii,select=c(as.vector(sf$sample)))
  ioepsii_col           = colnames(ioepsii)
  ioepsii$event_id      = rownames(ioepsii)
  ioepsii               = ioepsii[ioepsii$event_id %in% iot$event_id,]
  ioepsii_iot           = merge(ioepsii,iot, by="event_id")
  rownames(ioepsii_iot) = ioepsii_iot$event_id
  
  pdf(file=paste(outcode,".general.quant.pdf",sep=""),height=4,width=16)
  par(mfrow=c(1,4))
  for (tii in list_evtype) {
    
    # prepare matrix
    psii = ioepsii_iot[ioepsii_iot$ev_type == tii, ioepsii_col]
    psii = na.omit(psii)
    psii = as.matrix(psii)
    
    # events with intermediate PSI (10-90%)
    bp=barplot(colSums((psii>0.1 & psii<0.9)*1),las=2, main=paste(tii,"events with PSI 10-90%"),ylab="# events", col="slategray3")
    text(x=bp,y=0, pos=3,cex=0.7,col="red",labels = colSums((psii>0.1 & psii<0.9)*1))
    # median PSI
    bp=barplot(colMedians(psii),las=2, main=paste(tii,"median PSI"),ylab="PSI", col="slategray3", names=colnames(psii), ylim=c(0,1))
    text(x=bp,y=1, pos=1,cex=0.7,col="red",labels = signif(colMedians(psii),3))
    # PSI distribution boxplot
    boxplot(psii,las=2, main=paste(tii,"PSI distribution"),ylab="PSI", col="slategray3")
    # PSI distribution ecdf
    plot(0,pch=NA,xlim=c(0,1),ylim=c(0,1),las=1, main=paste(tii,"PSI CDF"),xlab="PSI",ylab="Fraction", col="slategray3")
    colors = rainbow(ncol(psii)+ncol(psii),v = 0.8)[1:ncol(psii)]
    names(colors) = colnames(psii)
    for (sai in colnames(psii)) {
      psii_ecdf = ecdf(psii[,sai])
      lines(x=seq(0,1,0.01), y=psii_ecdf(v=seq(0,1,0.01)),col=colors[sai])
    }
    legend("topleft",legend = colnames(psii), col = colors, lty=1, bty = "n", cex = 0.7)
    
  }
  
  dev.off()
  
}



#### Clustering & dimred ####

if ("CLUST" %in% list_analisi) {
  
  ioepsii_clusterings   = list() # list of outputs (will be filled later)
  
  # cycle through relevant event types
  for (tii in c("SE","RI")) {
    
    if (tii %in% c("SE","MX")) {
      psii = ioepsii_iot[ioepsii_iot$ev_type == tii & ioepsii_iot$is_inframe, ioepsii_col]
    } else {
      psii = ioepsii_iot[ioepsii_iot$ev_type == tii, ioepsii_col]
    }
    
    # prepare matrix
    psii    = na.omit(psii)
    psii    = psii[apply(psii,1,sd) != 0,]
    psii    = as.matrix(psii)
    psii_st = t(apply(psii, 1, scale))
    colnames(psii_st) = colnames(psii)
    
    # dimensionality reduction
    message("# dimred ",tii)
    psiipsii_listt = dimredfun(
      matriu = psii_st, outputname = paste(outcode,".clu.psi_",tii,sep=""),varname="st PSI",
      cols_are = "sample",          rows_are = paste(tii,"events"),
      cols_dist_method = "pearson", rows_dist_method = "pearson",
      isbidi = F,clus_method = "ward.D2")
    
    # k-means clustering: profiles
    message("# kmeans clustering ",tii)
    ioepsii_clusterings[[tii]] = kmeansfun(
      matrix=psii,
      outputprefix = paste(outcode,".clu.psi_",tii,sep=""), 
      k_list = k_list, ylab = "PSI",ylims = c(0,1))
    
    # k-means clustering: functional enrichments in each cluster
    pdf(file=paste(outcode,".clu.psi_",tii,".kmeans.fun.pdf",sep=""), height=2.2*max(k_list),width=12)
    for (ki in k_list) {
      par(mfrow=c(max(k_list),4))
      kc     = paste("k",ki,sep="")
      message("# k-means funcs ", kc)
      for (kj in 1:ki) {
        psii_tabi = ioepsii_clusterings[[tii]]$table
        psii_tabi = psii_tabi[psii_tabi[,kc] == kj,]
        psii_gens = as.vector(iot[iot$event_id %in% rownames(psii_tabi),"long_transcript_id"])
        hygeofun(psii_gens,pannu,"transcript","domain",NA,paste(kc,"k =",kj,tii,sep=" "),20,
                 printfile = F)
        suppressMessages(
          topgofun(psii_gens,gomap,paste(outcode,"fun",sep="."),
                   paste(kc,"k =",kj,tii,sep=" "), c("BP","MF","CC"), "fisher", "elim",20,
                   printfile = F))
      }
    } 
    dev.off()
    
  } # end of event type loop
  
} else { message("# Skip CLUST!") }




#### IOT group comparisons (new pval) ####

pval_threshold=0.1

if ("NEWCOMP" %in% list_analisi) {
  
  # comparisons
  message("# compare groups: loop dPSI & pval")
  ioepsii_diff = data.frame()
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
    message("  pval CDF ",tici," ")
    
    
    for (tii in c("SE","RI")) {
    
      message("# calculate pval CDF ",tici,": ",tii)
      # matrix of frequencies per event (PSI)
      tac           = merge(expr_tpm,subset(iot, select=c("event_id","long_transcript_id","ev_type")),all.y=T, by="long_transcript_id")
      rownames(tac) = tac$event_id
      tac           = tac[tac$ev_type==tii,]
      tac           = tac[-which(names(tac) %in% c("event_id","ev_type","long_transcript_id"))]
      tac           = tac[,c(cli,clj)] 
      # matrix of frequencies per event (PSI): taf
      taf           = merge(ioepsii,subset(iot, select=c("event_id","long_transcript_id","ev_type")),all.y=T, by="event_id")
      rownames(taf) = taf$event_id
      taf           = taf[taf$ev_type==tii,]
      taf           = taf[,-which(names(taf) %in% c("event_id","ev_type","long_transcript_id"))]
      taf           = taf[,c(cli,clj)] 
      taf_st        = na.omit(t(apply(taf, 1, scale)))
      colnames(taf_st) = colnames(taf)
      
      # calculate pvalues, with correction
      ioepsii_diff_o = pval_from_CDF(
        matrix_frequencies = taf, matrix_tpms = tac,
        name_group_i = com[2],name_group_j = com[3],
        samples_group_i = cli,samples_group_j = clj,
        comparison_code = tici,
        pval_correction = "BH",
        ecdf_area = 1000)
      
      # format function output
      ioepsii_diff_i = ioepsii_diff_o$result
      colnames(ioepsii_diff_i) = c("event_id","dPSI","mean_logtpm","pval","pval.corr",
                                   "comparison","group.i","group.j","mean_PSI.i","mean_PSI.j")
      ioepsii_diff_i$ev_type   = tii
      ioepsii_diff_i           = na.omit(ioepsii_diff_i)
      ioepsii_diff             = rbind(ioepsii_diff,ioepsii_diff_i)
    
      # CDF pval: does it make sense?
      pdf(file=paste(outcode,".dpsi_co.",tici,"_",tii,"_dPSI_volcano_pvalCDF.pdf",sep=""),height=4.5,width=6)
      ioepsii_diff_o$plot()
      dev.off()
      
      # volcano plot
      message("  volcano ",tii)
      rownames(ioepsii_diff_i) = ioepsii_diff_i$event_id
      ioepsii_diff_j = subset(ioepsii_diff_i, select = c("dPSI","pval"))
      ioepsii_dili_i = volcanoexp(
        table = ioepsii_diff_j,
        plotname = paste(tici,tii,"dPSI",sep="_"),
        fileprefix = paste(outcode,"dpsi_co",sep="."),
        pthreshold = pval_threshold,
        fc_varname = "dPSI",minfold = 0.1,
        fcp = com[2], fcn = com[3],
        pval_varname = "pval",
        xlims=c(-1,1),ylims = c(0,-4))
      
      # functional enrichments
      message("  funcs dPSI>0 ", tii)
      psii_gens = as.vector(iot[iot$event_id %in% ioepsii_dili_i$genes_sp,"long_transcript_id"])
      hygeofun(psii_gens,pannu,"transcript","domain",paste(outcode,"dpsi_co",sep="."),paste(tici,tii,"dPSI_difpos",sep="_"),30,list_population=unique(as.character(iot$long_transcript_id)),padjbool = F)
      suppressMessages(topgofun(psii_gens,gomap,paste(outcode,"dpsi_co",sep="."),paste(tici,tii,"dPSI_difpos",sep="_"),c("BP","MF","CC"),"fisher","elim",30))
      message("  funcs dPSI<0 ", tii)
      psii_gens = as.vector(iot[iot$event_id %in% ioepsii_dili_i$genes_sn,"long_transcript_id"])
      hygeofun(psii_gens,pannu,"transcript","domain",paste(outcode,"dpsi_co",sep="."),paste(tici,tii,"dPSI_difneg",sep="_"),30,list_population=unique(as.character(iot$long_transcript_id)), padjbool = F)
      suppressMessages(topgofun(psii_gens,gomap,paste(outcode,"dpsi_co",sep="."),paste(tici,tii,"dPSI_difneg",sep="_"),c("BP","MF","CC"),"fisher","elim",30))
      
    }
    
  }
  
} else { message("# Skip NEWCOMP!") }

message("# FET!")


# Save session again
if ("SAVE" %in% list_analisi) {
  message("# Save")
  save.image(paste(outcode,"session.iot","RData",sep="."))
} else { message("# Skip SAVE!") }


stop("ARA!")












# adhoc plots
# for genes of interest related to FA metabolism (desaturases, synthases, etc)

library(pheatmap)

# function:
interest_genes_heatmap = function(interestlist, interestname, source_matrix, diff_matrix, evtype = "SE", inframe_only = T )  {
  
  original_len = length(interestlist)
  
  if (inframe_only) {
    interestmatr = source_matrix[source_matrix$long_transcript_id  %in% interestlist & source_matrix$ev_type == evtype & source_matrix$is_inframe == T,]  
  } else {
    interestmatr = source_matrix[source_matrix$long_transcript_id  %in% interestlist & source_matrix$ev_type == evtype,]  
  }
  interestmatr = subset(interestmatr, select=as.vector(sf$sample))
  interestmatr[is.na(interestmatr)] <- 0
  
  interestlist = rownames(interestmatr)
  
  interestmatr.ann = data.frame(row.names = interestlist)
  interestmatr.exp = matrix(nrow = 0,ncol = 7)
  
  ann_colors = list()
  
  for (com in complis) {
    tici     = paste(com[5],".",com[2],"-",com[3],sep="")
    ticidiff = diff_matrix[diff_matrix$comparison == tici,]
    ticidiff$isp = 
    ticidiff$isn = ticidiff$pval < pval_threshold & ticidiff$dPSI < 0
    isposi = (interestlist %in% ticidiff[ticidiff$pval < pval_threshold & ticidiff$dPSI > 0,]$event_id ) * 1
    isnegi = (interestlist %in% ticidiff[ticidiff$pval < pval_threshold & ticidiff$dPSI < 0,]$event_id ) * -1
    interestmatr.ann[tici] = as.factor(isposi + isnegi)
    ann_colors[[tici]]     = c("-1" = "magenta3","0"="slategray1","1"="green3")
  }
  
  if (nrow(interestmatr) > 1) {
    pdf(file=paste("genesubset.de_adhoc",interestname,evtype,"pdf",sep="."),height=10+length(interestlist)/4,width=10)
    pheatmap(interestmatr, color = col.fun(21), breaks = seq(0,1,length.out = 20),
             border_color = "white", annotation_row = interestmatr.ann,cellheight = 8,cellwidth = 8,cluster_cols=F, 
             scale="none",annotation_colors = ann_colors,
             gaps_col = seq(min(table(as.numeric(sf$group))),to=length(sf$group), by=min(table(as.numeric(sf$group)))),
             main=paste("DE: ",interestname," (expressed: ",length(interestlist)," out of ",original_len,")",sep=""))
    dev.off()
  } else {
    print("No n'hi ha prous!")
  }
  
}

# SE
interest_genes_heatmap(
  interestlist  = pannu[ pannu$domain == "ketoacyl-synt" ,]$transcript, 
  interestname  = "FAsynthase_ketoacyl-synt", 
  source_matrix = ioepsii_iot, diff_matrix = ioepsii_diff, inframe_only = F, evtype = "SE" )

interest_genes_heatmap(
  interestlist  = pannu[ pannu$domain == "ELO" ,]$transcript, 
  interestname  = "FAelongase_ELO", 
  source_matrix = ioepsii_iot, diff_matrix = ioepsii_diff, inframe_only = F, evtype = "SE" )

interest_genes_heatmap(
  interestlist  = pannu[ pannu$domain == "FA_desaturase" ,]$transcript, 
  interestname  = "FAdesaturase_FA_desaturase", 
  source_matrix = ioepsii_iot, diff_matrix = ioepsii_diff, inframe_only = F, evtype = "SE")

interest_genes_heatmap(
  interestlist  = pannu[ pannu$domain == "NAD_binding_4" ,]$transcript, 
  interestname  = "FAreductase_NAD_binding_4", 
  source_matrix = ioepsii_iot, diff_matrix = ioepsii_diff, inframe_only = F, evtype = "SE" )

interest_genes_heatmap(
  interestlist  = pannu[ pannu$domain == "p450" ,]$transcript, 
  interestname  = "FAdecarboxyl_p450", 
  source_matrix = ioepsii_iot, diff_matrix = ioepsii_diff, inframe_only = F, evtype = "SE" )



# RI
interest_genes_heatmap(
  interestlist  = pannu[ pannu$domain == "ketoacyl-synt" ,]$transcript, 
  interestname  = "FAsynthase_ketoacyl-synt", 
  source_matrix = ioepsii_iot, diff_matrix = ioepsii_diff, inframe_only = F, evtype = "RI" )

interest_genes_heatmap(
  interestlist  = pannu[ pannu$domain == "ELO" ,]$transcript, 
  interestname  = "FAelongase_ELO", 
  source_matrix = ioepsii_iot, diff_matrix = ioepsii_diff, inframe_only = F, evtype = "RI" )

interest_genes_heatmap(
  interestlist  = pannu[ pannu$domain == "FA_desaturase" ,]$transcript, 
  interestname  = "FAdesaturase_FA_desaturase", 
  source_matrix = ioepsii_iot, diff_matrix = ioepsii_diff, inframe_only = F, evtype = "RI")

interest_genes_heatmap(
  interestlist  = pannu[ pannu$domain == "NAD_binding_4" ,]$transcript, 
  interestname  = "FAreductase_NAD_binding_4", 
  source_matrix = ioepsii_iot, diff_matrix = ioepsii_diff, inframe_only = F, evtype = "RI" )

interest_genes_heatmap(
  interestlist  = pannu[ pannu$domain == "p450" ,]$transcript, 
  interestname  = "FAdecarboxyl_p450", 
  source_matrix = ioepsii_iot, diff_matrix = ioepsii_diff, inframe_only = F, evtype = "RI" )





























stress_group="RFE"
n_quantiles=4

pdf(file=paste(outcode,".general.quant_fraction.pdf",sep=""),height=4.5,width=3*(n_quantiles+2))
for (tii in c("SE","RI")) {
 
  # separate genes by fraction of transcript length
  plotquantPSI(taula = ioepsii_iot, variable = "long_transcript_len",tii=tii, stress_group=stress_group)
   
  # separate genes by fraction of transcript length
  expr_tpm_i          = expr_tpm[,sf$sample]
  expr_tpm_m          = data.frame(row.names = rownames(expr_tpm_i))
  expr_tpm_m$gene     = rownames(expr_tpm_i)
  expr_tpm_m$mean_tpm = rowMeans(expr_tpm_i)
  expr_tpm_a          = merge(ioepsii_iot, expr_tpm_m, by.x="long_transcript_id",by.y="gene", all=T)
  plotquantPSI(taula = expr_tpm_a, variable = "mean_tpm",tii=tii, stress_group=stress_group)
  
}

dev.off()






plotquantPSI = function(taula,variable,tii,stress_group,n_quantiles=4) {
  
  taula=taula[taula$ev_type == tii,]
  variable=variable
  colors = rainbow(n_quantiles*2,v = 0.8)[1:n_quantiles]
  par(mfrow=c(1,n_quantiles+2))
  otherlevels = levels(sf$group)[!levels(sf$group) %in% stress_group]
  
  # plot correlation between variable and dpsi, in all samples
  coli = 1
  clj  = as.character(sf[sf$group %in% stress_group,]$sample)
  
  taula$quantiles = cut(taula[,variable], quantile(taula[,variable], probs=0:10/10, na.rm=T),na.rm=T, include.lowest=TRUE, labels=FALSE)
  boxplot(rowMeans(taula[,as.character(sf$sample)]) ~ taula$quantiles ,xlab=paste("quantiles",variable),ylab="PSI",ylim=c(0,1),
       main=paste(tii," correlation PSI  ~ ",variable," all",sep=""))
  cortest = cor.test(x=taula[,variable],y=rowMeans(taula[,clj]),method = "sp")
  text(x=0.5,pos=4,y=coli/20,labels = paste(stress_group, "p =",signif(cortest$p.value,3),"/ rho =", signif(cortest$estimate,3)),col=colors[coli],cex=0.7)
  for (sagi in otherlevels) {
    coli = coli+1
    cli = as.character(sf[sf$group %in% sagi,]$sample)
    cortest = cor.test(taula[,variable],rowMeans(taula[,cli]))
    text(x=0.5,pos=4,y=coli/20,labels = paste(sagi, "p =",signif(cortest$p.value,3),"/ rho =", signif(cortest$estimate,3)),col=colors[coli],cex=0.7)
  }
  

  # plot ecdfs of PSI in all genes
  coli = 1
  psij = taula[, clj]
  psij = as.vector(as.matrix(psij[!is.na(psij)]))
  psij_ecdf = ecdf(psij)
  plot(0,pch=NA,xlim=c(0,1),ylim=c(0,1),las=1,
       main=paste(tii," PSI CDF ",variable," all",sep=""),
       xlab="PSI",ylab="Fraction",sub=paste("KS one-sided (l/g), n =",length(psij)))
  lines(x=seq(0,1,0.01), y=psij_ecdf(v=seq(0,1,0.01)),col=colors[coli])
  for (sagi in otherlevels) {
    coli = coli+1
    cli = as.character(sf[sf$group %in% sagi,]$sample)
    psii = taula[, cli]
    psii = as.vector(as.matrix(psii[!is.na(psii)]))
    psii_ecdf = ecdf(psii)
    lines(x=seq(0,1,0.01), y=psii_ecdf(v=seq(0,1,0.01)),col=colors[coli])
    kstest_l = ks.test(psii,psij,alternative = "less")
    kstest_g = ks.test(psii,psij,alternative = "greater")
    text(x=0.5,y=coli/20,labels = paste(stress_group,sagi, "p =",signif(kstest_l$p.value,3),"/",signif(kstest_g$p.value,3)),col=colors[coli],cex=0.7)
  }
  legend("topleft",legend = levels(sf$group), col = colors, pch=1, cex = 0.7) 
  
  # plot ecdfs of PSI in genes, fractioned according to *variable* quantiles
  fra_top = as.vector(quantile(taula[,variable], probs=c(seq(0,1,length.out = n_quantiles+1)), na.rm=T))
  for (fra_toi in 1:(length(fra_top)-1)) {
    
    coli = 1
    clj  = as.character(sf[sf$group %in% stress_group,]$sample)
    psij = taula[taula[,variable] > fra_top[fra_toi] & taula[,variable] < fra_top[fra_toi+1], clj]
    psij = as.vector(as.matrix(psij[!is.na(psij)]))
    psij_ecdf = ecdf(psij)
    plot(0,pch=NA,xlim=c(0,1),ylim=c(0,1),las=1,
         main=paste(tii," PSI CDF ",variable," Q",fra_toi,": ",signif(fra_top[fra_toi],3),"-",signif(fra_top[fra_toi+1],3),sep=""),
         xlab="PSI",ylab="Fraction",sub=paste("KS one-sided (l/g), n =",length(psij)))
    lines(x=seq(0,1,0.01), y=psij_ecdf(v=seq(0,1,0.01)),col=colors[coli])
    for (sagi in otherlevels) {
      coli = coli+1
      cli = as.character(sf[sf$group %in% sagi,]$sample)
      psii = taula[taula[,variable] > fra_top[fra_toi] & taula[,variable] < fra_top[fra_toi+1], cli]
      psii = as.vector(as.matrix(psii[!is.na(psii)]))
      psii_ecdf = ecdf(psii)
      lines(x=seq(0,1,0.01), y=psii_ecdf(v=seq(0,1,0.01)),col=colors[coli])
      kstest_l = ks.test(psii,psij,alternative = "less")
      kstest_g = ks.test(psii,psij,alternative = "greater")
      text(x=0.5,y=coli/20,labels = paste(stress_group,sagi, "p =",signif(kstest_l$p.value,3),"/",signif(kstest_g$p.value,3)),col=colors[coli],cex=0.7)
    }
    legend("topleft",legend = levels(sf$group), col = colors, pch=1, cex = 0.7)
    
  }
}




