library(data.table)
library(apcluster)
library(reshape)
library(stringr)
library(RedeR)
library(ggplot2)

setwd('/Users/urmovosa/Documents/move_to_mac/CRP/')

to.read = file("reconstituted_genesets_150901.binary.dat", "rb")

rows <- read.table('reconstituted_genesets_150901.binary.rows.txt', head = F)
cols <- read.table('reconstituted_genesets_150901.binary.columns.txt', head = F)

# how many indices 
length(cols$V1) * length(rows$V1)

and <- matrix(nrow = length(rows$V1), ncol = length(cols$V1), data = NA)
and <- as.data.table(and)

rownames(and) <- rows$V1
colnames(and) <- as.character(cols$V1)

a1 <- readBin(to.read, numeric(), n = 219217416 + 3,  endian = "big")[-1]

# tulem <- data.table(data = NA, nrow = length(rows), ncol = length(cols))
# for (i in 1:length(rows$V1)){
#   ab1 <- i
#   ab2 <- i + length(cols$V1)
#   ab3 <- as.data.table(a1[ab1:ab2])
#   
#   tulem <- cbind(tulem, ab3)
#   print(i)
# }

a2 <- matrix(a1, ncol = 19987)
colnames(a2) <- rows$V1
rownames(a2) <- cols$V1
#a2 <- t(a2)

# find significant pathways
genesetenrichment1 <- fread('/Users/urmovosa/Documents/move_to_mac/CRP/CRP_data/combined_HapMap_1kG_5e8_MAF001_halfN_perm5000_rep50_unique_genesetenrichment.txt')

genesetenrichment2 <- fread('/Users/urmovosa/Documents/move_to_mac/CRP/CRP_data/combined_HapMap_1kG_5e8_MAF001_halfN_perm10000_rep500_unique_1_genesetenrichment.txt')
genesetenrichment1[str_detect(genesetenrichment1$`Original gene set ID`, 'KEGG'), ]
genesetenrichment1[str_detect(genesetenrichment1$`Original gene set ID`, 'GO:'), ]
genesetenrichment1[str_detect(genesetenrichment1$`Original gene set ID`, 'REACTOME'), ]
genesetenrichment1[str_detect(genesetenrichment1$`Original gene set ID`, 'MP:'), ]

genesetenrichment1 <- genesetenrichment1[genesetenrichment1$`False discovery rate` %in% c('<0.01', '<0.05'), ]
genesetenrichment1 <- genesetenrichment2[genesetenrichment2$`False discovery rate` %in% c('<0.01', '<0.05'), ]
table(genesetenrichment1$`Original gene set ID` %in% genesetenrichment2$`Original gene set ID`)
table(genesetenrichment2$`Original gene set ID` %in% genesetenrichment1$`Original gene set ID`)

# find significantly prioritized genes:
geneprioritization <- fread('/Users/urmovosa/Documents/move_to_mac/CRP/CRP_data/combined_HapMap_1kG_5e8_MAF001_halfN_perm10000_rep500_unique_1_geneprioritization.txt')
geneloci <- fread('/Users/urmovosa/Documents/move_to_mac/CRP/CRP_data/combined_HapMap_1kG_5e8_MAF001_halfN_perm10000_rep500_unique_1_loci.txt')
geneprioritization$Gene <- paste(geneprioritization$`Gene symbol`, geneprioritization$`Ensembl gene ID`, geneprioritization$`Chromosome and position`, sep = ' ')
geneprioritization$abi <- str_replace(geneprioritization$`Chromosome and position`, ':.*', '')
geneprioritization$abi2 <- str_replace(geneprioritization$`Chromosome and position`, '.*:', '')
geneprioritization$abi2 <- str_replace(geneprioritization$abi2, '-.*', '')
geneprioritization$abi <- str_replace(geneprioritization$abi, 'chr', '')
geneprioritization$abi <- paste(geneprioritization$abi, geneprioritization$abi2, sep = ' ')

geneprioritization <- geneprioritization[geneprioritization$`False discovery rate` %in% c('<=0.01', '<0.05'), ]
geneprioritization <- geneprioritization[order(geneprioritization$abi, 
                                               geneprioritization$`Nominal P value`), ]

geneprioritization$Gene <- factor(geneprioritization$Gene, levels = rev(as.character(geneprioritization$Gene)))

geneprioritization$varv <- 'proov'
for(i in seq(from = 1, to = length(unique(geneprioritization$abi)), by = 3)){
  geneprioritization[geneprioritization$abi == unique(geneprioritization$abi)[i], ]$varv <- 'red'
}
for(i in seq(from = 2, to = length(unique(geneprioritization$abi)), by = 3)){
  geneprioritization[geneprioritization$abi == unique(geneprioritization$abi)[i], ]$varv <- 'lightblue'
}
for(i in seq(from = 3, to = length(unique(geneprioritization$abi)), by = 3)){
  geneprioritization[geneprioritization$abi == unique(geneprioritization$abi)[i], ]$varv <- 'black'
}


p <- ggplot(geneprioritization, aes(x = Gene, y = -log10(`Nominal P value`), shape = `False discovery rate`)) + 
  geom_point(colour = as.character(geneprioritization$varv)) + theme_bw() + coord_flip() +
  theme(legend.position="none") + ylab(expression(paste(-log[10], (P)))) + scale_y_continuous(limits = c(0, 8)) + 
  scale_color_manual(values = as.character(geneprioritization$varv)) + ylab('Locus')

write.table(geneprioritization[order(geneprioritization$`Nominal P value`), 1:10, with = F], 'DEPICT_genenprioritization_perm10000_rep500.txt', 
            sep = '\t', row.names = F, quote = F)



pdf('DEPICT_prioritized_loci.pdf', width = 5, height = 7)
p
dev.off()



sets <- as.character(genesetenrichment1$`Original gene set ID`)

# construct character vector

paths <- as.character(cols$V1)
sig_paths <- paths %in% sets
sig_paths_names <- paths[paths %in% sets]

a3 <- a2[rownames(a2) %in% sig_paths_names, ]
genesetenrichment1$setname <- 'set'
genesetenrichment1$`Original gene set ID` <- as.character(genesetenrichment1$`Original gene set ID`)
genesetenrichment1$`Original gene set description` <- as.character(genesetenrichment1$`Original gene set description`)
genesetenrichment1$`Original gene set description` <- str_replace(genesetenrichment1$`Original gene set description`, 'KEGG_', '')
genesetenrichment1$`Original gene set description` <- str_replace(genesetenrichment1$`Original gene set description`, 'REACTOME_', '')
genesetenrichment1$`Original gene set description` <- str_replace_all(genesetenrichment1$`Original gene set description`, '_', ' ')
genesetenrichment1$`Original gene set description` <- str_replace_all(genesetenrichment1$`Original gene set description`, 'PPI ', '')
genesetenrichment1[str_detect(genesetenrichment1$`Original gene set ID`, 'KEGG'), ]$`Original gene set description` <- tolower(genesetenrichment1[str_detect(genesetenrichment1$`Original gene set ID`, 'KEGG'), ]$`Original gene set description`)
genesetenrichment1[str_detect(genesetenrichment1$`Original gene set ID`, 'REACTOME'), ]$`Original gene set description` <- tolower(genesetenrichment1[str_detect(genesetenrichment1$`Original gene set ID`, 'REACTOME'), ]$`Original gene set description`)
genesetenrichment1$`Original gene set description` <- str_replace_all(genesetenrichment1$`Original gene set description`, ':', '-')

genesetenrichment1[str_detect(genesetenrichment1$`Original gene set ID`, 'MP:'), ]$setname <- paste('MP: ', genesetenrichment1[str_detect(genesetenrichment1$`Original gene set ID`, 'MP:'), ]$`Original gene set description`, sep = '')
genesetenrichment1[str_detect(genesetenrichment1$`Original gene set ID`, 'KEGG'), ]$setname <- paste('KEGG: ', genesetenrichment1[str_detect(genesetenrichment1$`Original gene set ID`, 'KEGG'), ]$`Original gene set description`, sep = '')
genesetenrichment1[str_detect(genesetenrichment1$`Original gene set ID`, 'GO:'), ]$setname <- paste('GO: ', genesetenrichment1[str_detect(genesetenrichment1$`Original gene set ID`, 'GO:'), ]$`Original gene set description`, sep = '')
genesetenrichment1[str_detect(genesetenrichment1$`Original gene set ID`, 'REACTOME'), ]$setname <- paste('REACTOME: ', genesetenrichment1[str_detect(genesetenrichment1$`Original gene set ID`, 'REACTOME'), ]$`Original gene set description`, sep = '')
genesetenrichment1[str_detect(genesetenrichment1$`Original gene set ID`, 'ENSG'), ]$setname <- paste('PPI: ', genesetenrichment1[str_detect(genesetenrichment1$`Original gene set ID`, 'ENSG'), ]$`Original gene set description`, sep = '')
write.table(genesetenrichment1, '/Users/urmovosa/Documents/move_to_mac/CRP/CRP_data/DEPICT_genesetenrichment_per10000_rep500.txt', sep = '\t', quote = F, row.names = F)

geneset_abi <- genesetenrichment1[, c(1, ncol(genesetenrichment1), 3), with = F]
geneset_abi <- geneset_abi[order(geneset_abi$`Original gene set ID`), ]
a3 <- a3[order(rownames(a3)), ]
rownames(a3) <- geneset_abi$setname

# apcluster

cor_matrix <- corSimMat(a3)
pr_range <- preferenceRange(cor_matrix)
cor_cluster <- apcluster(cor_matrix, details = T, maxits = 10000, seed = 4321)
#cor_cluster <- apcluster(cor_matrix, details = T, q = 0.5)
#cor_cluster <- apclusterK(cor_matrix, details = T, K = 50)
plot(cor_cluster)
heatmap(cor_cluster, cor_matrix)
plot(cor_cluster, a3)



promAgg <- aggExCluster(cor_matrix, cor_cluster)
plot(promAgg)
promAgg <-cutree(promAgg, 40)
heatmap(promAgg, cor_matrix)

# for RedeR
cor_mat2 <- as.data.frame(cor_matrix)
cor_mat2$ID1 <- row.names(cor_mat2)
cor_mat2 <- melt(cor_mat2)
colnames(cor_mat2)[2] <- 'ID2'

cor_mat2 <- as.data.table(cor_mat2)
abi1 <- cor_mat2
abi2 <- cor_mat2
abi1$comb <- paste(abi1$ID1, abi1$ID2)
abi2$comb <- paste(abi2$ID2, abi2$ID1)
abi <- rbind(abi1, abi2)
setkey(abi, comb)

abi <- as.data.frame(abi)

# cor_mat2 <- melt(cor_mat2, id.vars = c('ID1', 'ID2'), measure.vars = 'value')
# cor_mat2 <- cor_mat2[!cor_mat2$comb2 == cor_mat2$comb1, ]

cor_mat2 <- as.data.frame(cor_mat2)

abi <- abi[!abi$comb == cor_mat2$comb]

#cor_mat2 <- cor_mat2[cor_mat2$value < 1, ]
cor_mat2 <- cor_mat2[!cor_mat2$ID1 == cor_mat2$ID2, ]

pathways <- data.frame(pathway = unique(c(cor_mat2$ID1, cor_mat2$ID2)))
path_abi <- geneset_abi[geneset_abi$setname %in% pathways$pathway, 1:3, with = F]
path_abi$logP <- -log10(path_abi$`Nominal P value`)
path_abi <- path_abi[, -1, with = F]


pal <- colorRampPalette(c("lightblue", "red"),
                        space = "rgb")

cor_mat3 <- cor_mat2[cor_mat2$value < 1 & cor_mat2$value > 0.3, ]
g <- graph.data.frame(cor_mat3, directed = F, vertices = path_abi)
#g <- att.setv(g, from = "logP", to = "nestAlias", breaks = seq(0, 6, 0.2), pal = 2)
g <- att.setv(g, from = "logP", to = "nodeColor", breaks = seq(0, 6, 0.2), pal = 2)
g <- att.setv(g, from = "logP", to = "nodeLineColor", breaks = seq(0, 6, 0.2), pal = 2)
#g <- att.setv(g, from = "setname", to = "nodeAlias", breaks = seq(0, 6, 0.2), pal = 1)
g <- att.sete(g, from = "value", to = "edgeWidth", xlim = c(1, 5, 0.1))



# plot
rdp <- RedPort() 
calld(rdp)
V(g)$nodeFontSize <- 7
V(g)$nodeLineColor <- 'black'

V(g)$nodeFontSize <- 10
V(g)$nodeLineColor <- 'black'
E(g)$edgeColor <- 'lightgrey'

addGraph(rdp, g)
relax(rdp)


cluster_relevance <- data.frame(exemplar = names(cor_cluster@exemplars), min_p = NA, max_p = NA, nr_of_nodes = NA, mean_p = NA, median_p = NA)

for (i in 1:length(names(cor_cluster@exemplars))){
  
  abi <- path_abi[path_abi$setname %in% names(cor_cluster@clusters[[i]]), ]
  cluster_relevance$min_p[i] <- min(abi$`Nominal P value`)
  cluster_relevance$max_p[i] <- max(abi$`Nominal P value`)
  cluster_relevance$nr_of_nodes[i] <- length(abi$`Nominal P value`)
  cluster_relevance$mean_p[i] <- mean(abi$`Nominal P value`)
  cluster_relevance$median_p[i] <- median(abi$`Nominal P value`)
  
}

cluster_relevance <- cluster_relevance[order(cluster_relevance$nr_of_nodes, decreasing = T), ]

for(i in 1:length(names(cor_cluster@exemplars))){
  nestNodes(rdp, nodes = names(cor_cluster@clusters[[i]]), theme='tm4', isAnchor = F, gatt = list(isNest = T, 
                                                                                                  nestAlias = names(cor_cluster@exemplars)[i], 
                                                                                                  nestLineColor = 'black', 
                                                                                                  nestColor = 'white', 
                                                                                                  nodeFontSize = 30, 
                                                                                                  nestFontColor = 'red', 
                                                                                                  nestFontX = -8,
                                                                                                  nestImage = 'transparent',
                                                                                                  nestLineWidth = 3))
  
  
  
  print(i)
}

relax(rdp)
mergeOutEdges(rdp, rescale = F)
#mergeOutEdges(rdp, rescale = T)
#mergeOutEdges(rdp, rescale = T)

#V(g)$nodeAlias <- names(cor_cluster@exemplars)
updateGraph(rdp)

relax(rdp)


updateContainerSize(rdp)  
updateGraph(rdp)

scl <- g$legNodeColor$scale
leg <- g$legNodeColor$legend
addLegend.color(rdp, colvec = scl, labvec=leg, title = expression(paste(-log[10], "(P)")))

resetd(rdp)


# construct graph for exemplars from AP clustering results:

cor_mat4 <- cor_mat2[cor_mat2$value < 1, ]
cor_mat4 <- cor_mat4[cor_mat4$ID1 %in% names(cor_cluster@exemplars) & cor_mat4$ID2 %in% names(cor_cluster@exemplars), ]
cor_mat4 <- cor_mat4[cor_mat4$value > 0.2, ]
cor_mat4 <- head(cor_mat4, 1062)

pathways <- data.frame(pathway = unique(c(cor_mat4$ID1, cor_mat4$ID2)))
path_abi <- geneset_abi[geneset_abi$setname %in% pathways$pathway, 1:3, with = F]
path_abi$logP <- -log10(path_abi$`Nominal P value`)
path_abi <- path_abi[, -1, with = F]

g2<- graph.data.frame(cor_mat4, directed = F, vertices = path_abi)
#g <- att.setv(g, from = "logP", to = "nestAlias", breaks = seq(0, 6, 0.2), pal = 2)
g2 <- att.setv(g2, from = "logP", to = "nodeColor", breaks = seq(0, 6, 0.2), pal = 2)
g2<- att.setv(g2, from = "logP", to = "nodeLineColor", breaks = seq(0, 6, 0.2), pal = 2)
#g <- att.setv(g, from = "setname", to = "nodeAlias", breaks = seq(0, 6, 0.2), pal = 1)
g2 <- att.sete(g2, from = "value", to = "edgeWidth", xlim = c(1, 5, 0.1))
g2 <- att.sete(g2, from = "value", to = "edgeWeight", xlim = c(1, 5, 0.1))
g2 <- att.sete(g2, from = "value", to = "edgeColor", xlim = c(1, 5, 0.1))

rdp <- RedPort() 
calld(rdp)
V(g2)$nodeFontSize <- 4
V(g2)$nodeLineColor <- 'black'
E(g2)$edgeColor <- 'lightgrey'
coords <- layout_nicely(g2)
addGraph(rdp, g2, coords)

resetd(rdp)


# draw each cluster separately:

to_cytoscape <- data.frame()

for(i in 1:length(cor_cluster@exemplars)){
  
  if (nrow(cor_mat2[cor_mat2$ID1 %in% names(cor_cluster@clusters[[i]]) & cor_mat2$ID2 %in% names(cor_cluster@clusters[[i]]), ])) {
    cor_mat4 <- cor_mat2[cor_mat2$ID1 %in% names(cor_cluster@clusters[[i]]) & cor_mat2$ID2 %in% names(cor_cluster@clusters[[i]]), ] } else {
      cor_mat4 <- data.frame(ID1 = names(cor_cluster@clusters[[i]])[1], ID2 = names(cor_cluster@clusters[[i]])[2], value = NA)}
  cor_mat4$cluster <- NA
  cor_mat4$cluster <- names(cor_cluster@exemplars)[i]
  to_cytoscape <- rbind(to_cytoscape, cor_mat4)
  
}

# filter cluster on R > 0.3

#to_cytoscape <- to_cytoscape[to_cytoscape$value > 0.3, ]

# P-values for subnetworks
geneset_abi$setname <- str_replace(geneset_abi$setname, 'KEGG', 'KE')
geneset_abi$setname <- str_replace(geneset_abi$setname, 'REACTOME', 'RE')
#geneset_abi$setname <- str_replace(geneset_abi$setname, 'PPI', 'PI')
geneset_abi$logP <- -log10(geneset_abi$`Nominal P value`)

to_cytoscape$ID1 <- str_replace(to_cytoscape$ID1, 'KEGG', 'KE')
to_cytoscape$ID1  <- str_replace(to_cytoscape$ID1, 'REACTOME', 'RE')
to_cytoscape$ID2 <- str_replace(to_cytoscape$ID2, 'KEGG', 'KE')
to_cytoscape$ID2  <- str_replace(to_cytoscape$ID2, 'REACTOME', 'RE')
to_cytoscape <- to_cytoscape[!is.na(to_cytoscape$ID2), ]

write.table(geneset_abi, 'cytoscape_subnetworks_nodes.txt', sep = '\t', quote = F, row.names = F)
write.table(to_cytoscape, 'cytoscape_subnetworks_edges.txt', sep = '\t', quote = F, row.names = F)


# extract exemplars:


to_cytoscape2 <- data.frame()





cor_mat4 <- cor_mat2[cor_mat2$ID1 %in% names(cor_cluster@exemplars) & cor_mat2$ID2 %in% names(cor_cluster@exemplars), ]
to_cytoscape2 <- cor_mat4
to_cytoscape2$ID1 <- str_replace(to_cytoscape2$ID1, 'KEGG', 'KE')
to_cytoscape2$ID1  <- str_replace(to_cytoscape2$ID1, 'REACTOME', 'RE')
to_cytoscape2$ID2 <- str_replace(to_cytoscape2$ID2, 'KEGG', 'KE')
to_cytoscape2$ID2  <- str_replace(to_cytoscape2$ID2, 'REACTOME', 'RE')

#write.table(geneset_abi, 'cytoscape_subnetworks_nodes.txt', sep = '\t', quote = F, row.names = F)
write.table(to_cytoscape2, 'cytoscape_exemplars_edges.txt', sep = '\t', quote = F, row.names = F)


