# library(ape)
# gff_fn = "/home/xavi/Documents/pirimeth-resistance-agam/data_genome/Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.9.gff3"
# gff_fn = "/home/xavi/Documents/oenocytes-agam/data_genome/Anogam_long.annot.gff"
# gff_df = read.gff(gff_fn)


# plot one gene
plot_gene = function(gff_df, gene_id, pos_start=NA, pos_end=NA, flanking_bp=100,
                     gene_label = "gene", transcript_label = "mRNA", exon_label = "CDS", 
                     utr5_label = "five_prime_UTR", utr3_label = "three_prime_UTR") {

  # subset gff
  gff_gene = gff_df[ gff_df$type == gene_label,]
  gff_txis = gff_df[ gff_df$type == transcript_label,]
  gff_exon = gff_df[ gff_df$type == exon_label,]
  gff_utr5 = gff_df[ gff_df$type == utr5_label,]
  gff_utr3 = gff_df[ gff_df$type == utr3_label,]
  
  # obtain transcript list
  gene = gff_gene[ grep(sprintf("ID=%s", gene_id), gff_gene$attributes) , ]
  txis = gff_txis[ grep(sprintf("Parent=%s", gene_id), gff_txis$attributes) , ]
  txis_ids = gsub("ID=","",txis$attributes)
  txis_ids = gsub(";.*","",txis_ids)
  
  # genewise info
  chrom  = as.character(gene$seqid)
  strand = as.character(gene$strand)
  if (is.na(pos_start)) { pos_start = gene$start - flanking_bp }
  if (is.na(pos_end))   { pos_end = gene$end + flanking_bp }
  
  plot_y=1
  plot(NA,NA, xlim=c(pos_start, pos_end), ylim = c(0,10), xlab = "bp", ylab="", sub=chrom, yaxt='n', bty="n",
       main=gene_id)
  for (txi in txis_ids) {
    exon = gff_exon[ grep(sprintf("Parent=%s", txi), gff_exon$attributes) , ]  
    utr5 = gff_utr5[ grep(sprintf("Parent=%s", txi), gff_utr5$attributes) , ]  
    utr3 = gff_utr3[ grep(sprintf("Parent=%s", txi), gff_utr3$attributes) , ]
    txii = gff_txis[ grep(sprintf("ID=%s", txi), gff_txis$attributes) , ]
    segments(x0=txii$start, x1=txii$end, y0=plot_y, y1=plot_y, col = "darkblue", lwd = 1.2 ,lend = 1, frame.plot = F)
    text(x = txii$start, y=plot_y+0.5, labels = sprintf("%s (%s)", txi, strand) ,col="darkblue", cex=0.8, adj=0)
    if (nrow(exon) > 0 ) { rect(xleft=exon$start, xright = exon$end, ybottom=plot_y-0.25, ytop = plot_y+0.25, col = "blue", border = NA) }
    if (nrow(utr5) > 0 ) { rect(xleft=utr5$start, xright = utr5$end, ybottom=plot_y-0.15, ytop = plot_y+0.15, col = "lightskyblue", border = NA) }
    if (nrow(utr3) > 0 ) { rect(xleft=utr3$start, xright = utr3$end, ybottom=plot_y-0.15, ytop = plot_y+0.15, col = "lightskyblue", border = NA) }
    plot_y = plot_y + 1.5
  }
  
}



