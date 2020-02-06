source("helper_scripts/genes_from_gff.R")
library("ape")


# load data
gff_fn = "data_genome/Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.9.gff3"
gff_df = read.gff(gff_fn)
outcode = "results_as/"

gene_list = c("AGAP004373",
              "AGAP001713",
              "AGAP003051", 
              "AGAP004572", 
              "AGAP001076", 
              "AGAP002195",
              "AGAP013511"
              )

pdf(file=sprintf("%s/fa_genes_coordinates.pdf", outcode),height=5,width=8)
for (gene in gene_list) {
  
  plot_gene(gff_df = gff_df, gene_id = gene)
  
}
dev.off()