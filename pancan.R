#start
rm(list=ls())
setwd('D:/work/TCGA_Pancancer')
library(fdrtool)
library(ggplot2)
library(plyr)
library(VennDiagram)
library(xtable)

#def
read.expr <- function() {
  expr <- read.table('EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena', T, sep='\t', 
                     stringsAsFactors = F)
  expr <- expr[grep('^[A-Za-z].+$', expr$sample), ]
  expr <- expr[!duplicated(expr$sample), ]
  row.names(expr) <- expr$sample
  expr <- expr[, -1]
  return(expr)
}

read.clin <- function() {
  clinical<- read.table('Survival_SupplementalTable_S1_20171025_xena_sp', 
                        T, sep='\t', stringsAsFactors = F)
  clinical$sample <- gsub('-', '.', clinical$sample)
  clinical$X_PATIENT <- gsub('-', '.', clinical$X_PATIENT)
  clinical$histological_grade[clinical$histological_grade=='G1'|
                              clinical$histological_grade=='G2'] <- 'High Grade'
  clinical$histological_grade[clinical$histological_grade=='G3'|
                              clinical$histological_grade=='G4'] <- 'Low Grade'
  clinical <- clinical[clinical$histological_grade=='High Grade'|
                       clinical$histological_grade=='Low Grade',]
  clinical <- clinical[substring(clinical$sample,14,15)=='01',]
  write.table(clinical, 'clinical_grade.txt', sep='\t', quote=F, row.names=F)
  return(clinical)
}

cancer.design <- function(cancer) {
  sample <- Clin$sample[Clin$cancer.type.abbreviation == cancer]
  sample <- intersect(sample, colnames(Expr))
  grade <- Clin$histological_grade[match(sample, Clin$sample)]
  des <- data.frame(sample, grade)
  write.table(des, paste('./design/design.', cancer, '.txt', sep=''),
              quote = F, row.names = F, sep = '\t')
  return(des)
}

diff <- function(cancer, lfc_cutoff=1, alpha=0.05) {
  high <- as.character(Des$sample[Des$grade == 'High Grade'])
  low <- as.character(Des$sample[Des$grade == 'Low Grade'])
  p <- rep(1, nrow(Expr))
  lfc <- 1:nrow(Expr)
  gene <- 1:nrow(Expr)
  for(i in 1:nrow(Expr)){
    a <- t(Expr[i, na.omit(match(high, colnames(Expr)))])
    b <- t(Expr[i, na.omit(match(low, colnames(Expr)))])
    if(sd(a) != 0 & sd(b) != 0) {
      t <- t.test(a, b, na.rm=T)
      p[i] <- t$p.value
      gene[i] <- rownames(Expr)[i]
      lfc[i] <- mean(b) - mean(a)
    }
  }
  p.adj <- fdrtool(p, statistic="pvalue")
  results <- data.frame(gene, lfc, p, p.adj$qval)
  with(results, {
    results$trend <<- as.factor(ifelse(p.adj.qval <= alpha & abs(lfc) >= lfc_cutoff, 
                                       ifelse(lfc >= lfc_cutoff ,'UP','DOWN'),'NOT'))
  })
  deg <- results[(abs(results$lfc)>=lfc_cutoff & results$p.adj.qval <= alpha),]
  write.table(results, paste('./results/results_', cancer, '.txt', sep=''), 
              sep='\t', quote = F, row.names = F)
  write.table(deg, paste('./results/deg_', cancer, '.txt', sep=''), 
              sep='\t', quote = F, row.names = F)
  return(results)
}

find.deg <- function(cancer) {
  deg <- read.table(paste('./results/deg_', cancer, '.txt', sep=''), T, sep='\t', 
                    stringsAsFactors = F)
  return(deg)
}

volcano <- function(cancer, lfc_cutoff=1.1, alpha=0.05) {
  this_title <<- paste0('The volcano plot of ',cancer,
                        '\nCutoff for logFC is ',round(lfc_cutoff, 1),
                        '\nUp-regulated gene number: ',nrow(Res[Res$trend =='UP',]) ,
                        '\nDown-regulated gene number: ',nrow(Res[Res$trend =='DOWN',]))
  Deg <- arrange(Deg, abs(Deg$lfc), decreasing = T)
  label <- Deg$gene[1:10]
  Res <- Res[Res$p != 1,]
  Res$label <- ''
  Res$label[Res$gene %in% label] <- as.character(label)
  pdf(paste('./volcano/volcano_', cancer, '.pdf', sep=''))
  ggplot(data=Res, aes(x=lfc, y=-log10(p.adj.qval), color=trend)) +
  geom_point(alpha=0.4, size=1.75) +geom_text(aes(label = label), size=3,vjust=-0.5, alpha=1) +
  theme_set(theme_set(theme_bw(base_size=13)))+
  xlab("log2 fold change") + ylab("-log10 p-value") +
  ggtitle( this_title ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
  scale_colour_manual(values = c('blue','grey','red')) +
  geom_vline(xintercept = -lfc_cutoff, linetype = "dashed", color = "grey", size = 1) +
  geom_vline(xintercept = lfc_cutoff, linetype = "dashed", color = "grey", size = 1)+
  geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "grey", size = 1)
}

#main <- function(cancer, lfc_cutoff=1.5, alpha=0.05) {
#  Expr <- read.expr()
#  Clin <- read.clin()
#  Expr.this <- extract.expr('HNSC')
#  Clin.this <- extract.clin('HNSC')
#  design <- pair.design()
#  res <- pair.t(cancer, design, lfc_cutoff, alpha)
#  Deg <- find.deg(cancer)
#  volcano(cancer, lfc_cutoff, alpha)
#  com <- venn(cancer)
#}

#main
Expr <- read.expr()
Clin <- read.clin()
main <- function(cancer, lfc_cutoff=1, alpha=0.05) {
  Des <<- cancer.design(cancer)
  Res <<- diff(cancer, lfc_cutoff, alpha)
  Deg <<- find.deg(cancer)
  volcano(cancer, lfc_cutoff, alpha)
  dev.off()
}
cancer.list <- c("BLCA","CESC","CHOL","ESCA","HNSC","KIRC","LGG","LIHC","OV","PAAD","STAD","UCEC")
sapply(cancer.list, main)
