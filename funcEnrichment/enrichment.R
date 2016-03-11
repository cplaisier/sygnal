library(topGO)
# Read in GO mappings to affymetrix probe ids
library('org.Hs.eg.db')

# Automated for all clusters
d1 = read.csv('../output/cluster.members.genes.txt',header=F)
biclustMembership = list()
for(j in 1:length(d1[,1])) {
    biclustMembership[[j]] = strsplit(as.character(d1[j,]),split=' ')[[1]][-1]
}
xx <- annFUN.org("BP", mapping = "org.Hs.eg.db", ID = "symbol")
geneNames <- intersect(unique(unlist(biclustMembership)), unique(unlist(xx)))
tmp1 <- biclustMembership[[1]]
geneList <- factor(as.integer(geneNames %in% tmp1))
names(geneList) <- geneNames

# Make Biological Process GOData object
GOdata.BP <- new("topGOdata", ontology='BP', allGenes = geneList, annot = annFUN.org, mapping = 'org.Hs.eg.db', ID = 'symbol')
m1.BP <- matrix(nrow = length(GOdata.BP@graph@nodes), ncol = length(biclustMembership), dimnames = list(GOdata.BP@graph@nodes, 1:length(biclustMembership)))
m2.BP <- matrix(nrow = length(biclustMembership), ncol=2, dimnames = list(1:length(biclustMembership),c('Top10.Terms.BP','BH.sig.GO.Ids.BP')))
for( biclust in (1:length(biclustMembership)) ) {
    # Expand gene list and change factor in GOdata
    biclustGenes <- biclustMembership[[biclust]]
    GOdata.BP@allScores <- factor(as.integer(geneNames %in% biclustGenes))
    # Biological process
    r1.BP = runTest(GOdata.BP, algorithm = 'classic', statistic = 'fisher')
    m1.BP[,biclust] = r1.BP@score
    m2.BP[biclust,1] = gsub(',','_',paste(GenTable(GOdata.BP, r1.BP)[,2],collapse=';'))
    m2.BP[biclust,2] = paste(names(which(p.adjust(r1.BP@score,method='BH')<=0.05)),collapse=';')
}

write.csv(m2.BP,'../output/biclusterEnrichment_GOBP.csv') 

