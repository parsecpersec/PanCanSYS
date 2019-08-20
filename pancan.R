#start
rm(list=ls())
setwd('D:/work/TCGA_Pancancer/')
library(fdrtool)
library(gplots)
library(ggplot2)
library(plyr)
library(Rtsne)

read_cancer <- 

HNSC <- read.table('./HNSC/HiSeqV2', T, sep='\t', stringsAsFactors = F, row.names = 'sample')
ESCA <- read.table('./ESCA/HiSeqV2', T, sep='\t', stringsAsFactors = F, row.names = 'sample')
STAD <- read.table('./STAD/HiSeqV2', T, sep='\t', stringsAsFactors = F, row.names = 'sample')
clinical_HNSC <- read.table('HNSC_clinicalMatrix.txt', T, sep='\t', stringsAsFactors = F)
clinical_ESCA <- read.table('ESCA_clinicalMatrix.txt', T, sep='\t', stringsAsFactors = F)
clinical_STAD <- read.table('STAD_clinicalMatrix.txt', T, sep='\t', stringsAsFactors = F)

#replace char
clinical_HNSC$sampleID <- gsub('-', '.', clinical_HNSC$sampleID)
clinical_ESCA$sampleID <- gsub('-', '.', clinical_ESCA$sampleID)
clinical_STAD$sampleID <- gsub('-', '.', clinical_STAD$sampleID)

#find tumor and grade
for (i in 1:length(clinical_HNSC[,1])) {
  num <- as.numeric(substring(clinical_HNSC[i,1],14,15))
  hisgrade <- clinical_HNSC$neoplasm_histologic_grade[i]
  if (num %in% seq(1,9)) {clinical_HNSC$Is_tumor[i] <- 'Tumor'}
  if (num %in% seq(10,29)) {clinical_HNSC$Is_tumor[i] <- 'Normal'}
  if (hisgrade == 'G1' | hisgrade == 'G2') {clinical_HNSC$grade[i] <- 'high'}
  if (hisgrade == 'G3' | hisgrade == 'G4') {clinical_HNSC$grade[i] <- 'low'}
}
for (i in 1:length(clinical_ESCA[,1])) {
  num <- as.numeric(substring(clinical_ESCA[i,1],14,15))
  hisgrade <- clinical_ESCA$neoplasm_histologic_grade[i]
  if (num %in% seq(1,9)) {clinical_ESCA$Is_tumor[i] <- 'Tumor'}
  if (num %in% seq(10,29)) {clinical_ESCA$Is_tumor[i] <- 'Normal'}
  if (hisgrade == 'G1' | hisgrade == 'G2') {clinical_ESCA$grade[i] <- 'high'}
  if (hisgrade == 'G3' | hisgrade == 'G4') {clinical_ESCA$grade[i] <- 'low'}
}
for (i in 1:length(clinical_STAD[,1])) {
  num <- as.numeric(substring(clinical_STAD[i,1],14,15))
  hisgrade <- clinical_STAD$neoplasm_histologic_grade[i]
  if (num %in% seq(1,9)) {clinical_STAD$Is_tumor[i] <- 'Tumor'}
  if (num %in% seq(10,29)) {clinical_STAD$Is_tumor[i] <- 'Normal'}
  if (hisgrade == 'G1' | hisgrade == 'G2') {clinical_STAD$grade[i] <- 'high'}
  if (hisgrade == 'G3' | hisgrade == 'G4') {clinical_STAD$grade[i] <- 'low'}
}
clinical_HNSC$Is_tumor <- as.factor(clinical_HNSC$Is_tumor)
clinical_ESCA$Is_tumor <- as.factor(clinical_ESCA$Is_tumor)
clinical_STAD$Is_tumor <- as.factor(clinical_STAD$Is_tumor)
clinical_HNSC$grade <- as.factor(clinical_HNSC$grade)
clinical_ESCA$grade <- as.factor(clinical_ESCA$grade)
clinical_STAD$grade <- as.factor(clinical_STAD$grade)

#delete normal
tumor_HNSC <- clinical_HNSC[clinical_HNSC$Is_tumor == 'Tumor',]
tumor_ESCA <- clinical_ESCA[clinical_ESCA$Is_tumor == 'Tumor',]
tumor_STAD <- clinical_STAD[clinical_STAD$Is_tumor == 'Tumor',]
write.table(tumor_HNSC, 'tumor_HNSC.txt', sep='\t', row.names = F)
write.table(tumor_ESCA, 'tumor_ESCA.txt', sep='\t', row.names = F)
write.table(tumor_STAD, 'tumor_STAD.txt', sep='\t', row.names = F)

#HNSC design
order <- match(colnames(HNSC), tumor_HNSC$sampleID)
order2 <- (1:length(order))[!is.na(order)]
grade <- tumor_HNSC$grade[order]
logged_HNSC <- HNSC[, order2]
design_HNSC <- data.frame(colnames(logged_HNSC), grade[order2])
colnames(design_HNSC) <- c('sample', 'grade')

#HNSC t test
high <- as.character(design_HNSC$sample[design_HNSC$grade == 'high'])
low <- as.character(design_HNSC$sample[design_HNSC$grade == 'low'])
p <- rep(1, nrow(logged_HNSC))
lfc <- 1:nrow(logged_HNSC)
gene <- 1:nrow(logged_HNSC)
for(i in 1:nrow(logged_HNSC)){
  a <- t(logged_HNSC[i, match(high, colnames(logged_HNSC))])
  b <- t(logged_HNSC[i, match(low, colnames(logged_HNSC))])
  if(sd(a) != 0 & sd(b) != 0) {
    t <- t.test(a, b)
    p[i] <- t$p.value
    gene[i] <- rownames(logged_HNSC)[i]
    lfc[i] <- mean(a) - mean(b)
  }
}
p.adj <- fdrtool(p, statistic="pvalue")
compared_HNSC <- data.frame(gene, lfc, p, p.adj$qval)
lfc_cutoff=1.5
with(compared_HNSC, {
  compared_HNSC$trend <<- as.factor(ifelse(p.adj.qval <= 0.05 & abs(lfc) >= lfc_cutoff, 
                                           ifelse(lfc >= lfc_cutoff ,'UP','DOWN'),'NOT'))
})
write.table(compared_HNSC, 'compared_HNSC.txt', sep='\t', row.names = F)

#HNSC DEG
deg_HNSC <- compared_HNSC[(abs(compared_HNSC$lfc)>=lfc_cutoff & compared_HNSC$p.adj.qval <= 0.05),]
write.table(deg_HNSC, 'deg_HNSC.txt', sep='\t', row.names = F)
exp_of_deg_HNSC <- logged_HNSC[match(deg_HNSC$gene, rownames(logged_HNSC)),]
exp_of_deg_HNSC <- exp_of_deg_HNSC[1:20, ]


#HNSC Heatmap
x_lab <- paste(rep('sample', ncol(exp_of_deg_HNSC)), as.character(1:ncol(exp_of_deg_HNSC)), sep='')
heatmap.2(as.matrix(exp_of_deg_HNSC),col = redblue(75), 
          scale = "row", dendrogram = 'col',
          key = TRUE, keysize = 0.1, symkey = FALSE, density.info = "none", 
          trace = "none", cexRow = 0.5, labCol = x_lab,
          xlab="sample", ylab="gene", main = "Heatmap of HNSC",
          lmat=rbind(c(0,3), c(2,1), c(0,4)), lhei=c(1.5, 5, 2)
)

#HNSC Volcano
with(compared_HNSC, {
  this_title <<- paste0('Cutoff for logFC is ',round(lfc_cutoff, 1),
                        '\nUp-regulated gene number: ',nrow(compared_HNSC[trend =='UP',]) ,
                        '\nDown-regulated gene number: ',nrow(compared_HNSC[trend =='DOWN',]))
})

deg_HNSC <- arrange(deg_HNSC, abs(deg_HNSC$lfc), decreasing = T)
label <- deg_HNSC$gene[1:10]
compared_HNSC <- compared_HNSC[compared_HNSC$p != 1,]
compared_HNSC$label <- ''
compared_HNSC$label[compared_HNSC$gene %in% label] <- as.character(label)

with(compared_HNSC, {
  g <<- ggplot(data=compared_HNSC, aes(x=lfc, y=-log10(p.adj.qval), color=trend)) +
    geom_point(alpha=0.4, size=1.75) +geom_text(aes(label = label), size = 3,vjust=-0.5, alpha=1) +
    theme_set(theme_set(theme_bw(base_size=13)))+
    xlab("log2 fold change") + ylab("-log10 p-value") +
    ggtitle( this_title ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
    scale_colour_manual(values = c('blue','grey','red'))
})

g=g+geom_vline(xintercept = -lfc_cutoff, linetype = "dashed", color = "grey", size = 1) +
  geom_vline(xintercept = lfc_cutoff, linetype = "dashed", color = "grey", size = 1)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey", size = 1)
g

#HNSC t-SNE

iris_unique <- unique(exp_of_deg_HNSC) # Remove duplicates
iris_matrix <- as.matrix(iris_unique)
iris <- t(normalize_input(iris_matrix))
set.seed(22)
tsne_out <- Rtsne(t(iris_matrix), pca=F, perplexity=30) # Run TSNE
tsne_plot <- data.frame(x = tsne_out$Y[,1], y = tsne_out$Y[,2], grade = design_HNSC$grade)
g <- ggplot(data=tsne_plot, aes(x=x, y=y, color=grade)) +
  geom_point(alpha=0.4, size=1.75) +
  theme_set(theme_set(theme_bw(base_size=13)))+
  xlab("t-SNE 1") + ylab("t-SNE 2") +
  ggtitle('t-SNE of HNSC') + theme(plot.title = element_text(size=15,hjust = 0.5))+
  scale_colour_manual(values = c('blue','red'))
g


#GO annotation
library(Hmisc)
DrawGOBubblePlot <- function(dat, dat_name, top.number = 10, col="blue"){
  dat1 = dat[c(1:top.number),c(2,3,4,5)]
  colnames(dat1)[3] = 'GeneRatio'
  dat1$Term = substr(dat1$Term,12,200)
  dat1$Term = capitalize(dat1$Term)
  dat1$Term = factor(dat1$Term,levels=dat1$Term[length(dat1$Term):1])
  dat1$PValue = -log10(dat1$PValue)
  p = ggplot(dat1,aes(GeneRatio,Term)) +
    geom_point(aes(size=Count,colour=PValue)) +
    scale_colour_gradient(low=col,high="red") + 
    labs(colour=expression(-log[10]("P Value")),size="Gene counts",  
         x="Gene Ratio",y="",title=paste('Biological process (', dat_name, ')', sep='')) +
    theme_bw() +
    scale_x_continuous(limits = c(0,max(dat1$GeneRatio) * 1.2)) 
  return(p)
}
go_HNSC <- read.table('GO_HNSC.txt', T, sep='\t', stringsAsFactors = F)
DrawGOBubblePlot(go_HNSC, 'HNSC', 10, "blue")
go_ESCA <- read.table('GO_ESCA.txt', T, sep='\t', stringsAsFactors = F)
DrawGOBubblePlot(go_ESCA, 'ESCA', 10, "blue")
go_STAD <- read.table('GO_STAD.txt', T, sep='\t', stringsAsFactors = F)
DrawGOBubblePlot(go_STAD, 'STAD', 9, "blue")

#Venn
venn <- venn.diagram(list('HNSC'=deg_HNSC$gene, 
                          'ESCA'=deg_ESCA$gene, 
                          'STAD'=deg_STAD$gene), 'venn diagram.png')
common <- intersect(intersect(deg_HNSC$gene, deg_ESCA$gene), deg_STAD$gene)
com1 <- intersect(deg_HNSC$gene, deg_ESCA$gene)
com2 <- intersect(deg_HNSC$gene, deg_STAD$gene)
com3 <- intersect(deg_STAD$gene, deg_ESCA$gene)
common <- c(common, as.character(deg_HNSC$trend[deg_HNSC$gene==common]),
            as.character(deg_ESCA$trend[deg_ESCA$gene==common]),
            as.character(deg_STAD$trend[deg_STAD$gene==common]))
common <- data.frame(common[1], common[2], common[3], common[4])
colnames(common) <- c('gene', 'HNSC', 'ESCA', 'STAD')
common1 <- data.frame(com1, deg_HNSC[match(com1, deg_HNSC$gene),5],
                      deg_ESCA[match(com1, deg_ESCA$gene), 5])
colnames(common1) <- c('gene', 'HNSC', 'ESCA')
common2 <- data.frame(com2, deg_HNSC[match(com2, deg_HNSC$gene),5],
                      deg_STAD[match(com2, deg_STAD$gene), 5])
colnames(common2) <- c('gene', 'HNSC', 'STAD')
common3 <- data.frame(com3, deg_STAD[match(com3, deg_STAD$gene),5],
                      deg_ESCA[match(com3, deg_ESCA$gene), 5])
colnames(common3) <- c('gene', 'STAD', 'ESCA')
comlist <- rbind(com1, com2)
write.table(common, 'common.txt', sep='\t', row.names = F, quote = F)
write.table(common1, 'common1.txt', sep='\t', row.names = F, quote = F)
write.table(common2, 'common2.txt', sep='\t', row.names = F, quote = F)
write.table(common3, 'common3.txt', sep='\t', row.names = F, quote = F)
