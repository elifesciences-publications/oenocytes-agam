### Define input ####

# input files
library(ape)
library(phytools)

phy_list  = c("FAelongase/FAelongase.00.lt.iqt.treefile",
              "FAdesaturase/FAdesaturase.00.lt.iqt.treefile",
              "FAreductase/FAreductase.00.lt.iqt.treefile",
              "FAsynthase/FAsynthase.00.lt.iqt.treefile")

pdf(file="gene_phylogenies.pdf",height=12,width=9)
for (phi in phy_list) {
  
  phy = ape::read.tree(phi)
  phy = phytools::midpoint.root(phy)
  phy = ladderize(phy, right = FALSE)
  
  # paint branches
  phy = phytools::paintBranches(
    phy,edge=sapply(phy$tip.label[grepl("Anogam_",phy$tip.label)],match,phy$tip.label),
    state="Anogam", anc.state = "none")
  phy = phytools::paintBranches(
    phy,edge=sapply(phy$tip.label[grepl("Aedaeg_",phy$tip.label)],match,phy$tip.label),
    state="Aedaeg", anc.state = "none")
  phy = phytools::paintBranches(
    phy,edge=sapply(phy$tip.label[grepl("Dromel_",phy$tip.label)],match,phy$tip.label),
    state="Dromel", anc.state = "none")
  phy_col_d = phy$mapped.edge
  phy_col_d = phy_col_d[,c("Anogam","Aedaeg","Dromel","none")]
  phy_col_v = apply(phy_col_d, 1, function(x) which(x>0))
  phy_col_c = c("firebrick2","darkorange1","springgreen4","slategray")[phy_col_v]
  
  # paint tips
  phy_tips_label = sapply(base::strsplit(phy$tip.label, split = "_"), "[", 1)
  phy_tips_color = phy_tips_label
  phy_tips_color[phy_tips_color=="Anogam"] = "firebrick2"
  phy_tips_color[phy_tips_color=="Aedaeg"] = "darkorange1"
  phy_tips_color[phy_tips_color=="Dromel"] = "springgreen4"
  
  # plot
  ape::plot.phylo(
    main=phi,
    phy,
    use.edge.length=T, 
    show.tip.label=T, 
    edge.color = phy_col_c, 
    tip.color = phy_tips_color,
    font = 1, 
    align.tip.label = T,
    show.node.label = T,
    cex=0.8, edge.width = 1.2,
    underscore = T
    )
  ape::add.scale.bar(x=0, y=0, lcol="slategray", cex=0.8)
  
}
dev.off()



# separate larger page for p450

phy_list = c("FAdecarboxy/FAdecarboxy.00.lt.iqt.treefile")

pdf(file="gene_phylogenies_p450.pdf",height=48,width=9)
for (phi in phy_list) {
  
  phy = ape::read.tree(phi)
  phy = phytools::midpoint.root(phy)
  phy = ladderize(phy, right = FALSE)
  
  # paint branches
  phy = phytools::paintBranches(
    phy,edge=sapply(phy$tip.label[grepl("Anogam_",phy$tip.label)],match,phy$tip.label),
    state="Anogam", anc.state = "none")
  phy = phytools::paintBranches(
    phy,edge=sapply(phy$tip.label[grepl("Aedaeg_",phy$tip.label)],match,phy$tip.label),
    state="Aedaeg", anc.state = "none")
  phy = phytools::paintBranches(
    phy,edge=sapply(phy$tip.label[grepl("Dromel_",phy$tip.label)],match,phy$tip.label),
    state="Dromel", anc.state = "none")
  phy_col_d = phy$mapped.edge
  phy_col_d = phy_col_d[,c("Anogam","Aedaeg","Dromel","none")]
  phy_col_v = apply(phy_col_d, 1, function(x) which(x>0))
  phy_col_c = c("firebrick2","darkorange1","springgreen4","slategray")[phy_col_v]
  
  # paint tips
  phy_tips_label = sapply(base::strsplit(phy$tip.label, split = "_"), "[", 1)
  phy_tips_color = phy_tips_label
  phy_tips_color[phy_tips_color=="Anogam"] = "firebrick2"
  phy_tips_color[phy_tips_color=="Aedaeg"] = "darkorange1"
  phy_tips_color[phy_tips_color=="Dromel"] = "springgreen4"
  
  # plot
  ape::plot.phylo(
    main=phi,
    phy,
    use.edge.length=T, 
    show.tip.label=T, 
    edge.color = phy_col_c, 
    tip.color = phy_tips_color,
    font = 1, 
    align.tip.label = T,
    show.node.label = T,
    cex=0.8, edge.width = 1.2,
    underscore = T
  )
  ape::add.scale.bar(x=0, y=0, lcol="slategray", cex=0.8)
  
}
dev.off()




