library(readr)
library(shape)
library(biomaRt)
library(data.table)
library(stringr)
library(gridBase)
library(Gviz)
library(plotrix)
library(zoo)
library(dplyr)
library(Cairo)

setwd('/Users/urmovosa/')

GWAS_results <- fread('Documents/move_to_mac/ALS/ALS_MLMA_summary_statistics/mlma_p_for_depict_check.txt')

GWAS_results_p10e8 <- GWAS_results[GWAS_results$p < 5e-8, ]


# what is missing 
chr3_proxy <- read.table('Documents/move_to_mac/ALS/ALS_MLMA_summary_statistics/chr3_80_proxySearch.results.csv', head = T, sep = '\t')
chr8_proxy <- read.table('Documents/move_to_mac/ALS/ALS_MLMA_summary_statistics/chr8_80_proxySearch.results.csv', head = T, sep = '\t')
chr9_proxy <- read.table('Documents/move_to_mac/ALS/ALS_MLMA_summary_statistics/chr9_80_proxySearch.results.csv', head = T, sep = '\t')
chr14_proxy <- read.table('Documents/move_to_mac/ALS/ALS_MLMA_summary_statistics/chr14_80_proxySearch.results.csv', head = T, sep = '\t')
chr17_proxy <- read.table('Documents/move_to_mac/ALS/ALS_MLMA_summary_statistics/chr17_80_proxySearch.results.csv', head = T, sep = '\t')
chr19_proxy <- read.table('Documents/move_to_mac/ALS/ALS_MLMA_summary_statistics/chr19_80_proxySearch.results.csv', head = T, sep = '\t')
chr19_proxy <- rbind(chr3_proxy,
                     chr8_proxy,
                     chr9_proxy,
                     chr14_proxy,
                     chr17_proxy,
                     chr19_proxy)

#LD info
chr3_block <- fread('Documents/move_to_mac/ALS/ALS_MLMA_summary_statistics/chr3_blockannotation.results.csv', head = T, sep = '\t')
chr8_block <- fread('Documents/move_to_mac/ALS/ALS_MLMA_summary_statistics/chr8_blockannotation.results.csv', head = T, sep = '\t')
chr9_block <- fread('Documents/move_to_mac/ALS/ALS_MLMA_summary_statistics/chr9_blockannotation.results.csv', head = T, sep = '\t')
chr14_block <- fread('Documents/move_to_mac/ALS/ALS_MLMA_summary_statistics/chr14_blockannotation.results.csv', head = T, sep = '\t')
chr17_block <- fread('Documents/move_to_mac/ALS/ALS_MLMA_summary_statistics/chr17_blockannotation.results.csv', head = T, sep = '\t')
chr19_block <- fread('Documents/move_to_mac/ALS/ALS_MLMA_summary_statistics/chr19_blockannotation.results.csv', head = T, sep = '\t')
all_block <- rbind(chr3_block,
                   chr8_block,
                   chr9_block,
                   chr14_block,
                   chr17_block,
                   chr19_block)


chr19_proxy$CHR <- as.character(chr19_proxy$CHR)
rm(chr3_proxy, chr8_proxy, chr9_proxy, chr14_proxy, chr17_proxy)

# read in brain eQTLs
low_confidence_eQTL <- read_tsv('Documents/move_to_mac/ALS/eQLT_data_mining_analyses/low_confidence_eQTLs.txt')
# remove error
err <- read.delim('Documents/move_to_mac/ALS/eQLT_data_mining_analyses/err_gibbs.txt', head = F)
low_confidence_eQTL <- low_confidence_eQTL[!low_confidence_eQTL$Study == 'Gibbs' & !low_confidence_eQTL$SNP %in% err$V1, ]
rm(err)




high_confidence_eQTL <- read.table('Documents/move_to_mac/ALS/eQLT_data_mining_analyses/Hu_meta_analysis_eQTLs.txt')
braineac_eQTLs <- fread('Documents/move_to_mac/ALS/BRAINEAC_results/braineac_FDR_0.05.txt', sep = '\t')

colnames(low_confidence_eQTL)[6] <- 'Group'

low_confidence_eQTL$type <- 'transcript_level'
high_confidence_eQTL$type <- 'transcript_level'

braineac_eQTLs$Group <- 'BRAINEAC'

braineac_eQTLs <- data.frame(SNP = braineac_eQTLs$rsid,
                             HUGO = braineac_eQTLs$SYMBOL,
                             P = braineac_eQTLs$adj_p,
                             Group = paste('BRAINEAC', braineac_eQTLs$variable, sep = '_'),
                             type = braineac_eQTLs$type)

# Deelen et al. eQTLs

Deelen <- fread('Documents/move_to_mac/ALS/eQLT_data_mining_analyses/Deelen/Deelen_eQTLs.txt')
Deelen <- data.frame(SNP = Deelen$SNPName, HUGO = Deelen$HGNCName, P = Deelen$FDR, Group = 'Deelen', type = 'transcript_level')

# merge together
brain_eQTLs <- unique(rbind(high_confidence_eQTL[, c(2, 6, 7, 8, 13)],
                            low_confidence_eQTL[, c(4, 3, 5, 6, 12)],
                            Deelen,
                            braineac_eQTLs))

brain_eQTLs <- brain_eQTLs[!is.na(brain_eQTLs$HUGO),]

### add Myers data from seeQTL
Myers <- read.delim('Documents/move_to_mac/ALS/eQLT_data_mining_analyses/Myers/Myers_seeQTL_eQTLs.txt', head = T)
brain_eQTLs <- rbind(brain_eQTLs, Myers)

# define window
win <- 500000


# locus:
i <- 19

for(i in c(3, 8, 9, 14, 17, 19, 21)){
  
  # pdf(paste('Documents/move_to_mac/ALS/ALS_MLMA_summary_statistics/regional_plots_for_presentation', i, '.pdf', sep = ''), width = 17/2, height = 10/1.5, useDingbats = F)
  svg(paste('Documents/move_to_mac/ALS/ALS_MLMA_summary_statistics/regional_plots_Chromm_2', i, '.svg', sep = ''), width = 17/2, height = 6.5, family = 'Arial')
  #CairoFontMatch()
  
  par(mfrow = c(4, 1), mar = c(0.25, 6, 0.25, 3), oma = c(3, 1, 4, 6))
  layout(matrix(c(1, 2, 3, 4)), heights = c(2.25, 0.6, 0.6, 0.6))
  
  # chromosome
  chr <- paste('chr', i, sep = '')
  
  abi  <-  GWAS_results_p10e8[as.character(GWAS_results_p10e8$CHR) == i, ]
  abi <- abi[abi$p == min(abi$p), ]
  
  abi2 <- GWAS_results[GWAS_results$CHR == i & (GWAS_results$position < abi$position + win &
                                                  GWAS_results$position > abi$position - win), ]
  
  all_block_a <- all_block[all_block$CHR == i, ]
  #abi2 <- abi2[abi2$p < 0.05]
  
  write.table(abi2[, c(1, 4), with = F], 
              paste('Documents/move_to_mac/ALS/ALS_MLMA_summary_statistics/GWAS_results_for_snipa_chr', i, '.txt', sep = ''), 
              sep = '\t', quote = F, row.names = F, col.names = F)
  
  colnames(chr19_proxy)[2] <- 'SNP'
  abi3 <- merge(abi2, chr19_proxy, by = 'SNP')
  abi3[is.na(abi3$R2), ]$R2 <- 0
  
  abi3$col <- rep(NA, times = nrow(abi3))
  abi3 <- as.data.frame(abi3)
  abi3[abi3$R2 >= 0 & abi3$R2 <= 0.2, ]$col <- 'khaki1'
  abi3[abi3$R2 > 0.2 & abi3$R2 <= 0.5, ]$col <- 'yellow'
  abi3[abi3$R2 > 0.5 & abi3$R2 <= 0.8, ]$col <- 'orange'
  abi3[abi3$R2 > 0.8, ]$col <- 'red'
  abi3$shape <- 1
  
  
  if (length(abi3[!abi3$EQTLGENES %in% '-', ]$shape) > 0) {
    abi3[!abi3$EQTLGENES %in% '-', ]$shape <- 2 }
  
  # graafik
  gene_coord <- c()
  for (j in c(1, -1)){
    ab <- paste(i, abi$position - win, abi$position + win, j, sep = ':')
    
    ensembl <- useMart(host='feb2014.archive.ensembl.org', biomart = 'ENSEMBL_MART_ENSEMBL', dataset = "hsapiens_gene_ensembl")
    #ensembl54 <- useMart(host='may2009.archive.ensembl.org', biomart = 'ENSEMBL_MART_ENSEMBL', dataset = "hsapiens_gene_ensembl")
    
    gene <- getBM(attributes = c('external_gene_id', 'start_position', 'end_position', 'strand'), 
                  filters = 'chromosomal_region', 
                  values = as.character(ab), 
                  mart = ensembl)
    
    gene_coord <- rbind(gene_coord, gene)
  }
  
  gene_coord$start <- gene_coord$start_position
  gene_coord$end <- gene_coord$end_position
  
  gene_coord[gene_coord$strand == '-1', ]$start <- gene_coord[gene_coord$strand == '-1', ]$end_position
  gene_coord[gene_coord$strand == '-1', ]$end <- gene_coord[gene_coord$strand == '-1', ]$start_position
  gene_coord <- gene_coord[, -c(2, 3)]
  gene_coord <- gene_coord[order(gene_coord$start), ]
  gene_coord$y <- rep(seq(from = -max(-log10(abi2$p))/12, to = -max(-log10(abi2$p))/3, by = (-max(-log10(abi2$p))/3)/7), 
                      length = nrow(gene_coord))
  
  
  ###  overlap with brain-eQTL
  abi3$brain_eQTL <- 'no'
  abi3$brain_eQTLGene <- '-'
  abi3$brain_exon_eQTL <- 'no'
  abi3$brain_exon_eQTLGene <- '-'
  
  for (l in 1:nrow(abi3)){
    if (abi3$SNP[l] %in% brain_eQTLs$SNP){
      if (length(as.character(brain_eQTLs[brain_eQTLs$SNP %in% as.character(abi3$SNP[l]), ]$HUGO)) == 1) {
        
        abi3[l, ]$shape <- 24
        abi3[l, ]$brain_eQTL <- 'yes'
        abi3[l, ]$brain_eQTLGene <- as.character(brain_eQTLs[brain_eQTLs$SNP %in% as.character(abi3$SNP[l]), ]$HUGO)
        #print(abi3$SNP[l])
        
      }
      if (length(as.character(brain_eQTLs[brain_eQTLs$SNP %in% as.character(abi3$SNP[l]), ]$HUGO)) > 1) {
        abi3[l, ]$shape <- 24
        abi3[l, ]$brain_eQTL <- 'yes'
        abi3[l, ]$brain_eQTLGene <- unique(as.character(brain_eQTLs[brain_eQTLs$SNP %in% as.character(abi3$SNP[l]), ]$HUGO))
        #print(as.character(brain_eQTLs[brain_eQTLs$SNP %in% as.character(abi3$SNP[l]), ]$HUGO))
      }
      
      if (length(as.character(brain_eQTLs[brain_eQTLs$type %in% 'exon_level' & brain_eQTLs$SNP %in% as.character(abi3$SNP[l]), ]$HUGO)) == 1) {
        abi3[l, ]$shape <- 24
        abi3[l, ]$brain_exon_eQTL <- 'yes'
        abi3[l, ]$brain_exon_eQTLGene <- unique(as.character(brain_eQTLs[brain_eQTLs$SNP %in% as.character(abi3$SNP[l]), ]$HUGO))
        print(as.character(brain_eQTLs[brain_eQTLs$SNP %in% as.character(abi3$SNP[l]), ]$HUGO))
      }
      
      if (length(as.character(brain_eQTLs[brain_eQTLs$type %in% 'exon_level' & brain_eQTLs$SNP %in% as.character(abi3$SNP[l]), ]$HUGO)) > 1) {
        abi3[l, ]$shape <- 24
        abi3[l, ]$brain_exon_eQTL <- 'yes'
        abi3[l, ]$brain_exon_eQTLGene <- unique(as.character(brain_eQTLs[brain_eQTLs$SNP %in% as.character(abi3$SNP[l]), ]$HUGO))
        print(as.character(brain_eQTLs[brain_eQTLs$SNP %in% as.character(abi3$SNP[l]), ]$HUGO))
      }
      
    }
  }  
  
  
  ### NB! Maybe would be better to remove SNPs not present in the 1000G v1
  
  plot(abi2[!abi2$SNP %in% abi3$SNP, ]$position/1000000, -log10(abi2[!abi2$SNP %in% abi3$SNP, ]$p), 
       ylim = c(-max(-log10(abi2$p))/3, round(max(-log10(abi2$p)) + (max(-log10(abi2$p)* 0.1)), 0)),
       xlim = c((abi$position - win)/1000000, (abi$position + win)/1000000),
       ylab = '',
       xlab = '',
       col = 'lightgrey',
       pch = 21, 
       cex = 0.7,
       axes = F)
  
  axis(2, at = seq(0, round(max(-log10(abi2$p)) + (max(-log10(abi2$p)* 0.1)), 0), by = 2), las = 2)
  axis(4, at = seq(0, round(max(-log10(abi2$p)) + (max(-log10(abi2$p)* 0.1)), 0), by = round(max(-log10(abi2$p)) + (max(-log10(abi2$p)* 0.1)), 0)/5), labels = seq(0, 100, by = 20), las = 2)
  mtext(expression(paste(-log[10], (P))), side = 2, line = 2, cex = 0.75)
  mtext('Recombination rate (cM/Mb)', side = 4, line = 2, cex = 0.75)
  mtext(unique(abi3$QRSID), side = 3, outer = TRUE, adj = 0.5, line = 0, cex = 1.5, font = 2)
  box()
  
  
  #eqtl effect
  points(abi3$position/1000000, -log10(abi3$p), col = abi3$col, pch = abi3$shape, bg = abi3$col)
  
  
  lines(all_block_a$POS1/1000000, rescale(all_block_a$CMMB, c(0, max(-log10(abi2$p) * max(all_block_a$CMMB))/100)), col = 'cyan4', lwd = 0.75)
  lines(c((abi$position - win)/1000000 - 0.04, (abi$position + win)/1000000 + 0.04) , c(-log10(5e-8), -log10(5e-8)), col = 'firebrick', lty = 2, lwd = 2)
  
  points(abi3[abi3$brain_eQTL == 'yes', ]$position/1000000, -log10(abi3[abi3$brain_eQTL == 'yes', ]$p), 
         col = abi3[abi3$brain_eQTL == 'yes', ]$col, pch = abi3[abi3$brain_eQTL == 'yes', ]$shape, 
         bg = 'black')
  
  points(abi3[abi3$brain_exon_eQTL == 'yes', ]$position/1000000, -log10(abi3[abi3$brain_exon_eQTL == 'yes', ]$p), 
         col = abi3[abi3$brain_exon_eQTL == 'yes', ]$col, pch = abi3[abi3$brain_exon_eQTL == 'yes', ]$shape, 
         bg = 'black')
  
  
  # mark eqtl genes
  gene_coord$col = rep('lightblue', times = nrow(gene_coord))
  
  # check this part, may contain bug!!!
  
  abi4 <- unique(abi3[abi3$CISEQTL == 'yes' & !abi3$EQTLGENES == '-', c(1, 25)])
  
  # ordinary eqtl
  ord_eqtl <- c()
  for (x in 1:nrow(abi4)){
    
    eq <- unlist(str_split(abi4[x, ]$EQTLGENES, ', '))
    kl <- data.frame(comb = paste(as.character(abi4[x, ]$SNP), eq, sep = ','))
    ord_eqtl <- rbind(ord_eqtl, kl)
  }
  
  library(tidyr)
  
  ord_eqtl <- separate(ord_eqtl, comb, c('SNP', 'gene'), sep = ',')
  
  # brain eqtl
  
  brain_eqtl <- brain_eQTLs[as.character(brain_eQTLs$SNP) %in% as.character(abi3$SNP), c(1, 2)]
  brain_exon_eqtl <- brain_eQTLs[brain_eQTLs$type %in% 'exon_level' & as.character(brain_eQTLs$SNP) %in% as.character(abi3$SNP), c(1, 2)]
  eqtl_gene <- merge(ord_eqtl, brain_eqtl, by = 'SNP', all = T)
  eqtl_gene <- merge(eqtl_gene, brain_exon_eqtl, by = 'SNP', all = T)
  
  if (nrow(eqtl_gene[!eqtl_gene$SNP == 'NA', ]) > 0){
    
    colnames(eqtl_gene) <- c('SNP', 'ord_eQTL', 'brain_eqtl', 'brain_exon_eqtl')
    eqtl_gene$ord_eQTL <- as.character(eqtl_gene$ord_eQTL)
    eqtl_gene$brain_eqtl <- as.character(eqtl_gene$brain_eqtl)
    eqtl_gene$brain_exon_eqtl <- as.character(eqtl_gene$brain_exon_eqtl)
    
    if (nrow(eqtl_gene[is.na(eqtl_gene$ord_eQTL), ]) > 0){
      eqtl_gene[is.na(eqtl_gene$ord_eQTL), ]$ord_eQTL <- '-'}
    if (nrow(eqtl_gene[is.na(eqtl_gene$brain_eqtl), ]) > 0){
      eqtl_gene[is.na(eqtl_gene$brain_eqtl), ]$brain_eqtl <- '-'}
    eqtl_gene <- unique(eqtl_gene)
  }
  
  for (k in 1:nrow(gene_coord)){
    
    if (nrow(eqtl_gene) > 0){
      
      
      if (gene_coord[k, ]$external_gene_id %in% eqtl_gene$ord_eQTL){gene_coord[k, ]$col = 'seagreen'}
      if (gene_coord[k, ]$external_gene_id %in% eqtl_gene$brain_eqtl){gene_coord[k, ]$col = 'black'}
      if (gene_coord[k, ]$external_gene_id %in% eqtl_gene$brain_exon_eqtl){gene_coord[k, ]$col = 'black'}
    }
  }
  
  
  
  
  
  korgus <- max(-log10(abi2$p)) + abs(max(gene_coord[gene_coord$strand == '-1', ]$y))
  
  gene_coord$adj  <- 0
  gene_coord[gene_coord$strand == '-1', ]$adj  <- 1
  
  ebaol <- gene_coord[gene_coord$col == 'lightblue', ]
  ol <- gene_coord[gene_coord$col %in% c('orange', 'black', 'seagreen', 'navy'), ]
  
  
  
  text(ebaol$start/1000000, ebaol$y + korgus/50, 
       labels = as.character(ebaol$external_gene_id), 
       cex = 0.5, font = 3, col = ebaol$col, adj = 0.5)
  
  Arrows(x0 = ebaol$start/1000000, y0 = ebaol$y, x1 = ebaol$end/1000000, 
         y1 = ebaol$y, cex = 0.5, arr.length = 7/70, col = ebaol$col, lwd = 1, arr.type = 'triangle')  
  
  if (nrow(ol) > 0){
    text(ol$start/1000000, ol$y + korgus/50, 
         labels = as.character(ol$external_gene_id), 
         cex = 0.6, font = 4, col = ol$col, adj = 0.5)
    
    Arrows(x0 = ol$start/1000000, y0 = ol$y, x1 = ol$end/1000000, 
           y1 = ol$y, cex = 0.5, arr.length = 7/70, col = ol$col, lwd = 1, arr.type = 'triangle')  
    
  }
  
  
  legend('topright', legend = c(expression(paste(R^2>0.8)),
                                expression(paste(R^2>0.5)),
                                expression(paste(R^2>0.2)),
                                expression(paste(R^2>0.1)),
                                expression(paste(italic('cis'), '-eQTL')),
                                expression(paste('Brain ', italic('cis'), '-eQTL'))),
         col = c('red', 'orange', 'yellow', 'khaki1', 'black', 'black'),
         cex = 1, pch = c(1, 1, 1, 1, 2, 24, 24), bty = 'o', pt.bg = c(rep('black', times = 6)), bg = 'white')
  
  #lines(c((abi$position - win)/1000000, (abi$position - win)/1000000), c(-100, 100))
  #lines(c((abi$position + win)/1000000, (abi$position + win)/1000000), c(-100, 100))
  
  
  # RNA-seq data
  load(file = paste('Documents/move_to_mac/ALS/brain_epigenome_data/chr', i, '_female_RNA_seq_pos.RData', sep = ''))
  load(file = paste('Documents/move_to_mac/ALS/brain_epigenome_data/chr', i, '_female_RNA_seq_neg.RData', sep = ''))
  chr <- paste('chr', i, sep = '')
  
  
  and <- as.data.frame(female_H3K4me1@range)
  and <- as.data.table(and)
  and$score <- as.numeric(female_H3K4me1@data)
  rm(female_H3K4me1)
  
  abi_and <- seq(from = abi$position - 510000, to = abi$position + 510000)
  abi_and <- data.table(start = abi_and)
  and <- merge(and, abi_and, by = 'start', all = T)
  and[is.na(and$seqnames), ]$score <- 0
  
  # down
  and2 <- as.data.frame(female_H3K4me1_neg@range)
  and2 <- as.data.table(and2)
  and2$score <- as.numeric(female_H3K4me1_neg@data)
  rm(female_H3K4me1_neg)
  #and2$roll <- rollmean(and2$score, 3)
  
  abi_and2 <- seq(from = abi$position - 510000, to = abi$position + 510000)
  abi_and2 <- data.table(start = abi_and2)
  and2 <- merge(and2, abi_and2, by = 'start', all = T)
  and2[is.na(and2$seqnames), ]$score <- 0
  
  and <- and %>%
    mutate(rolling_average_score = rollmean(x = score, 3*1000, align = "right", fill = 0))
  and2 <- and2 %>%
    mutate(rolling_average_score = rollmean(x = score, 3*1000, align = "right", fill = 0))
  and[nrow(and),]$rolling_average_score <- 0
  and2[nrow(and2),]$rolling_average_score <- 0
  
  plot(NULL, xlim = c((abi$position - win)/1000000, (abi$position + win)/1000000), ylim = c(-40, 40), axes = F, ylab = '', xlab = '')
  polygon(and$start/1000000, and$rolling_average_score, col = 'navy', border = NA)
  polygon(and2$start/1000000, and2$rolling_average_score, col = 'firebrick', border = NA)
  axis(2)
  
  mtext('RNA-seq', side = 2, line = 2, cex = 0.75)
  #lines(c((abi$position - win)/1000000, (abi$position - win)/1000000), c(-50, 50))
  #lines(c((abi$position + win)/1000000, (abi$position + win)/1000000), c(-50, 50))
  box()
  
  # add CHOMM data
  # Download states
  
  female_brain <- fread("Documents/move_to_mac/ALS/brain_epigenome_data/core_marks/E082_15_coreMarks_stateno.bed", sep = '\t')
  female_brain$tissue <- 'Fetal Brain Female'
  male_brain <- fread("Documents/move_to_mac/ALS/brain_epigenome_data/core_marks/E081_15_coreMarks_stateno.bed", sep = '\t')
  male_brain$tissue <- 'Fetal Brain Male'
  angular_gyrus <- fread("Documents/move_to_mac/ALS/brain_epigenome_data/core_marks/E067_15_coreMarks_stateno.bed", sep = '\t')
  angular_gyrus$tissue <- 'Brain Angular Gyrus'
  anterior_caudate <- fread("Documents/move_to_mac/ALS/brain_epigenome_data/core_marks/E068_15_coreMarks_stateno.bed", sep = '\t')
  anterior_caudate$tissue <- 'Brain Anterior Caudate'
  cingular_gyrus <- fread("Documents/move_to_mac/ALS/brain_epigenome_data/core_marks/E069_15_coreMarks_stateno.bed", sep = '\t')
  cingular_gyrus$tissue <- 'Brain Cingulate Gyrus'
  germinal_matrix <- fread("Documents/move_to_mac/ALS/brain_epigenome_data/core_marks/E070_15_coreMarks_stateno.bed", sep = '\t')
  germinal_matrix$tissue <- 'Brain Germinal Matrix'
  hippocampus_middle <- fread("Documents/move_to_mac/ALS/brain_epigenome_data/core_marks/E071_15_coreMarks_stateno.bed", sep = '\t')
  hippocampus_middle$tissue <- 'Brain Hippocampus Middle'
  temporal_lobe <- fread("Documents/move_to_mac/ALS/brain_epigenome_data/core_marks/E072_15_coreMarks_stateno.bed", sep = '\t')
  temporal_lobe$tissue <- 'Brain Inferior Temporal Lobe'
  prefrontal_cortex <- fread("Documents/move_to_mac/ALS/brain_epigenome_data/core_marks/E073_15_coreMarks_stateno.bed", sep = '\t')
  prefrontal_cortex$tissue <- 'Brain Dorsolateral Temporal Cortex'
  substantia_nigra <- fread("Documents/move_to_mac/ALS/brain_epigenome_data/core_marks/E074_15_coreMarks_stateno.bed", sep = '\t')
  substantia_nigra$tissue <- 'Brain Substantia Nigra'
  
  epigenome <- rbind(female_brain, male_brain, angular_gyrus, anterior_caudate, cingular_gyrus,
                     germinal_matrix, hippocampus_middle, temporal_lobe, 
                     prefrontal_cortex, substantia_nigra)
  
  epigenome$V2 <- as.numeric(epigenome$V2)
  epigenome$V3 <- as.numeric(epigenome$V3)
  epigenome$col <- NA
  epigenome <- as.data.frame(epigenome)
  
  epigenome[epigenome$V4 == '1', ]$col <- 'red'
  epigenome[epigenome$V4 == '2', ]$col <- 'orangered'
  epigenome[epigenome$V4 == '3', ]$col <- 'limegreen'
  epigenome[epigenome$V4 == '4', ]$col <- 'green'
  epigenome[epigenome$V4 == '5', ]$col <- 'darkgreen'
  epigenome[epigenome$V4 == '6', ]$col <- 'greenyellow'
  epigenome[epigenome$V4 == '7', ]$col <- 'yellow'
  epigenome[epigenome$V4 == '8', ]$col <- 'mediumaquamarine'
  epigenome[epigenome$V4 == '9', ]$col <- 'paleturquoise'
  epigenome[epigenome$V4 == '10', ]$col <- 'indianred'
  epigenome[epigenome$V4 == '11', ]$col <- 'darksalmon'
  epigenome[epigenome$V4 == '12', ]$col <- 'darkkhaki'
  epigenome[epigenome$V4 == '13', ]$col <- 'dimgray'
  epigenome[epigenome$V4 == '14', ]$col <- 'gainsboro'
  epigenome[epigenome$V4 == '15', ]$col <- 'white'
  
  epigenome <- as.data.frame(epigenome[epigenome$V1 == chr & 
                                         as.numeric(epigenome$V2) > abi$position - win - 100000 & 
                                         as.numeric(epigenome$V3) < abi$position + win + 100000, ])
  
  plot(NULL, xlim = c((abi$position - win)/1000000, (abi$position + win)/1000000), ylim = c(0, 10), axes = F, xlab = paste(chr, ' position(Mb)(hg19)', sep = ''), ylab = '')
  mtext(paste(chr, ' position (Mb)(hg19)', sep = ''), side = 1,  line = 3)
  
  j <- 0
  
  for (i in rev(unique(epigenome$tissue))){
    
    rect(epigenome[epigenome$tissue == i, ]$V2/1000000, rep(j, nrow(epigenome[epigenome$tissue == i, ])), 
         epigenome[epigenome$tissue == i, ]$V3/1000000, rep(j + 1, nrow(epigenome[epigenome$tissue == i, ])), 
         col = epigenome[epigenome$tissue == i, ]$col, border = NA)
    j <- j+1
  }
  
  axis(1)
  axis(2, at = seq(from = 0.5, to = 9.5, by = 1), labels = rev(unique(epigenome$tissue)), las = 2, cex.axis = 0.75)
  #lines(c((abi$position - 500000)/1000000, (abi$position - 500000)/1000000), c(-100, 100))
  #lines(c((abi$position + 500000)/1000000, (abi$position + 500000)/1000000), c(-100, 100))
  box()
  
  dev.off()
  
}



# DNAse data 
load(file = 'Documents/move_to_mac/ALS/brain_epigenome_data/chr8_male_DNase.RData')
load(file = 'Documents/move_to_mac/ALS/brain_epigenome_data/chr8_female_DNase.RData')

svg('Documents/move_to_mac/ALS/brain_epigenome_data/chr8_DNAse.svg', width = 17/2, height = 2)
plotTracks(list(male_H3K4me1, female_H3K4me1),
           stackHeight = 0.3,
           from = abi$position - 500000,
           to = abi$position + 500000,
           window = 300)
dev.off()

# legend for ChroMM
epi <- data.frame(chromm = c('Active TSS',
                             'Flanking Active TSS',
                             "Transcr. at gene 5\' and 3\'",
                             'Strong transcription',
                             'Weak transcription',
                             'Genic enhancers',
                             'Enhancers',
                             'ZNF genes & repeats',
                             'Heterochromatin',
                             'Bivalent/Poised TSS',
                             'Flanking Bivalent TSS/Enh',
                             'Bivalent Enhancer',
                             'Repressed PolyComb',
                             'Weak Repressed PolyComb',
                             'Quiescent/Low'),
                  varv = c('red',
                           'orangered',
                           'limegreen',
                           'green',
                           'darkgreen',
                           'greenyellow',
                           'yellow',
                           'mediumaquamarine',
                           'paleturquoise',
                           'indianred',
                           'darksalmon',
                           'darkkhaki',
                           'dimgray',
                           'gainsboro',
                           'white'
                  ))

svg('Documents/move_to_mac/ALS/brain_epigenome_data/legend.svg', width = 7, height = 7, family = 'Arial')
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "" )
legend("bottomright", legend = as.character(epi$chromm), 
       fill = as.character(epi$varv), col = as.character(epi$varv))
dev.off()