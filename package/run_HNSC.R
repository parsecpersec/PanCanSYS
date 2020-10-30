#### start ####
# module 0
setwd('D:/work/view_gene_package/')
library(ggplot2)
library(xtable)
library(survminer)
library(survival)
library(beeswarm)
library(clusterProfiler)
library(org.Hs.eg.db)
library(Hmisc)
library(plyr)
library(tsne)

#### set global varibles ####
'BMI1' %in% colnames(HNSC)

targene <- 'BMI1' # must exist in expdata
disp.name <- 'BMI1' # visualization in tables and figures
Time <- Sys.Date()
dir.name <- paste0(Time, '_', targene)
dir.create(paste0('./results_HNSC/', dir.name))
dest <- paste0('./results_HNSC/', dir.name, '/')

#### load ####
load('./input/pair.RData')
mRNA <- exp
rm(corm,pairclin,pairclin2,pairrna,pairrna2,pairrna2normal,exp,clin,clinpair)
Expr <- read.csv('./input/pair_tumor_and_normal.csv', T, stringsAsFactors=F, row.names='X')
myclin <- read.csv('./input/clinical_HNSC.csv', T, stringsAsFactors = F)
sams <- intersect(rownames(mRNA), myclin$sample)
HNSC <- mRNA[sams,]
myclin <- myclin[match(sams, myclin$sample),]
# node
myclin$node <- NA
myclin$node[grep('N0', myclin$tnm)] <- 'No'
myclin$node[grep('N[1-3]', myclin$tnm)] <- 'Yes'

#### 1.1 boxplot paired 43+43 ####
tumor <- myclin$sample
normal <- gsub('[.]01$', '.11', tumor)
normal <- intersect(normal, rownames(Expr))
tumor <- gsub('[.]11$', '.01', normal)
Exprc <- t(Expr)
tm <- Exprc[targene, tumor]
nm <- Exprc[targene, normal]
t1 <- t.test(tm, nm, paired = T)
pairclin <- data.frame(Group = c(rep('Tumor', length(tm)), rep('Normal', length(nm))), 
                       exp = c(tm, nm), stringsAsFactors = F)
pv <- t1$p.value
if(pv < 0.001) {ptext <- 'P < 0.001'}else {ptext <- paste0('P = ', round(pv, 3))}
bswm <- beeswarm(exp ~ Group, data = pairclin, spacing = 1.5)
colnames(bswm)[6] <- 'Group'
tiff(paste0(dest, 'boxplot of ', targene, ' (paired).tiff'), 
     height = 2500, width = 2500, res = 300)
print(ggplot(data = bswm, aes(x = Group, y = y, color = Group)) + 
        geom_boxplot(alpha=0, color='black', width = 0.5) + 
        geom_point(aes(x = x, shape = Group), size = 5, alpha = 0.6) +
        ggtitle(paste0('Differential Expression of ', disp.name, ' in HNSCC (TCGA)')) +
        theme_bw() + ylab(disp.name) + xlab(ptext) + scale_shape_manual(values = c(19, 15)) + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.title.x=element_text(vjust=0,size = 20,face='italic'),
              axis.title.y=element_text(size=25),
              title = element_text(size = 20), plot.title = element_text(hjust = 0.5), 
              axis.text.x = element_text(size = 25), axis.text.y=element_text(size=15),
              legend.text=element_text(size=15), legend.position='none') + 
        scale_color_manual(values = c('blue', 'red')))
dev.off()
rm(Exprc, pairclin, bswm, t1, pv, ptext)

#### 1.2 boxplot non-paired 520+43 ####
tumor <- myclin$sample
normal <- gsub('[.]01$', '.11', tumor)
normal <- intersect(normal, rownames(Expr))
Exprc <- t(mRNA)
tm <- Exprc[targene, tumor]
nm <- Exprc[targene, normal]
t1 <- t.test(tm, nm, paired = F)
pairclin <- data.frame(Group = c(rep('Tumor', length(tm)), rep('Normal', length(nm))), 
                       exp = c(tm, nm), stringsAsFactors = F)

pv <- t1$p.value
if(pv < 0.001) {ptext <- 'P < 0.001'}else {ptext <- paste0('P = ', round(pv, 3))}
bswm <- beeswarm(exp ~ Group, data = pairclin, spacing = 1.2)
colnames(bswm)[6] <- 'Group'
tiff(paste0(dest, 'boxplot of ', targene, ' (non-paired).tiff'), 
     height = 2500, width = 2500, res = 300)
print(ggplot(data = bswm, aes(x = Group, y = y, color = Group)) + 
        geom_boxplot(alpha=0, color='black', width = 0.5) + 
        geom_point(aes(x = x, shape = Group), size = 5, alpha = 0.6) +
        ggtitle(paste0('Differential Expression of ', disp.name, ' in HNSCC (TCGA)')) +
        theme_bw() + ylab(disp.name) + xlab(ptext) + scale_shape_manual(values = c(19, 15)) + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.title.x=element_text(vjust=0,size = 20,face='italic'),
              axis.title.y=element_text(size=25),
              title = element_text(size = 20), plot.title = element_text(hjust = 0.5), 
              axis.text.x = element_text(size = 25), axis.text.y=element_text(size=15),
              legend.text=element_text(size=15), legend.position='none') + 
        scale_color_manual(values = c('blue', 'red')))
dev.off()
rm(Exprc, pairclin, bswm, t1, pv, ptext)

#### 1.3 boxplot of clin-path traits and expression ####


# module 2
#### 2.1 KM curve ####
myclin$main.exp <- HNSC[, targene]
myclin$maingroup <- ifelse(myclin$main.exp >= quantile(myclin$main.exp, 0.5), 'High', 'Low')
# KM
myclin$OS.5.years[myclin$OS.5.years=='Alive'] <- 0
myclin$OS.5.years[myclin$OS.5.years=='Dead'] <- 1
myclin$OS.5.years <- as.numeric(myclin$OS.5.years)
# main group
sur.result <- survfit(Surv(time=myclin$OS.time.5.years, event=myclin$OS.5.years, 
                           type = 'right')~maingroup, data=myclin)
tiff(paste0(dest, 'KM curve of ', targene, '.tiff'), 
     height = 2500, width = 2500, res = 300)
print(ggsurvplot(sur.result, data=myclin, conf.int=F, pval=T, pval.method=T,
                 risk.table='absolute', 
                 legend.title=disp.name,
                 legend.labs=c('High', 'Low'),
                 palette=c('red', 'blue'),  
                 title=paste0('Kaplan-Meier Curve of HNSCC Overall Survival'), 
                 xlab='Time (months)',
                 risk.table.height=.30, 
                 ggtheme=theme_bw()+theme(plot.title=element_text(size=20,hjust=0.5),
                                          axis.title.x=element_text(size=18),
                                          axis.title.y=element_text(size=20), axis.text=element_text(size=12), 
                                          axis.text.y=element_text(size=12), legend.text=element_text(size=18),
                                          panel.background=element_blank(), axis.line=element_line(colour='black'),
                                          legend.title=element_text(size=18), text=element_text(face='plain'),
                                          panel.grid.major = element_blank(), panel.grid.minor = element_blank())))
dev.off()

#### 2.2 cox main ####
myclin$maingroup <- factor(myclin$maingroup, levels = c('Low', 'High'))
cox <- coxph(Surv(time=myclin$OS.time.5.years, event=myclin$OS.5.years, 
                  type = 'right')~radiotherapy+chemotherapy+
               gender+age+stage+grade+node+maingroup, data=myclin)  # select cov
cox.sum1 <- xtable(summary(cox)$coefficient)
cox.sum2 <- xtable(summary(cox)$conf.int)
cox.sum <- cbind(cox.sum1, cox.sum2)
cox.sum <- cox.sum[c(1:8),]
rownames(cox.sum) <- c('Radiotherapy', 'Chemotherapy', 'Male', 'Age',
                       'Stage III-IV', 'G3-G4', 'Nodal Metastasis', 'High Expression')
colnames(cox.sum)[c(2,8:9)] <- c('HR','lower','upper')
cox.sum$var <- rownames(cox.sum)
write.csv(cox.sum, paste0(dest, '/cox_', targene, '.csv'), quote = F, row.names = F)
pd <- position_dodge(0.1)
cox.sum$Factors <- 'Covariate'
cox.sum$Factors[cox.sum$var == 'High Expression'] <- disp.name
cox.sum$Factors <- factor(cox.sum$Factors, levels=c('Covariate', disp.name))
cox.sum$var <- factor(cox.sum$var, levels=rev(c('Male', 'Age', 'Stage III-IV', 'G3-G4', 'Nodal Metastasis', 
                                                'Radiotherapy', 'Chemotherapy', 'High Expression')))
tiff(paste0(dest, 'Cox HR of ', targene, '.tiff'), 
     height = 2500, width = 2500, res = 300)
print(ggplot(cox.sum, aes(x=var, y=HR, color = Factors)) + 
        ggtitle(paste0('Cox Propotional Hazard Model of HNSCC Overall Survival')) +
        geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2, position=pd, size=1.5) +
        geom_line(position=pd, size = 1.5) + xlab('Variables') + ylab('Hazard Ratio (HR)') +
        geom_hline(yintercept = 1, linetype = 'dashed', color = "blue", size = 1.5) +
        theme_bw() + coord_flip() +
        theme(axis.text.x = element_text(angle = 0, size = 18),
              axis.text.y = element_text(size = 15), axis.line=element_line(colour='black'),
              plot.title = element_text(size = 21, hjust=0.5), axis.title=element_text(size=20),
              legend.text=element_text(size=15), legend.title=element_text(size=18), 
              panel.background=element_blank(), legend.key=element_rect(),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
        geom_point(position=pd, size=3) + 
        scale_color_manual(values = c('black', 'red')))
dev.off()
rm(cox.sum1,cox.sum2,pd)

#
##
###
#### 3.1 coexpression ####
# design
library(limma)
Group <- factor(myclin$maingroup, levels=c('Low', 'High'))
design <- model.matrix(~0+Group)
colnames(design) <- c('Low', 'High')
fit <- lmFit(as.data.frame(t(HNSC), stringsAsFactors=F), design)
contrast.matrix <- makeContrasts(High-Low,
                                 levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
Deg <- topTable(fit2,adjust.method="BH",coef=1,p.value=1,
                lfc=0,number=Inf,sort.by = 'logFC')
write.csv(Deg, paste0(dest, 'coexp1_', targene, '.csv'), quote=F, row.names=T)

# cor
corm <- cor(HNSC[,targene], HNSC, method='pearson', use='complete.obs')
corm <- as.data.frame(t(corm), stringsAsFactors=F)
colnames(corm) <- 'Value'
corm$Abs <- abs(corm$Value)
corm$Gene <- rownames(corm)
corm <- corm[corm$Gene!=targene,]
corm$Trend <- ifelse(corm$Value > 0, 'Positive', 'Negative')
write.csv(corm, paste0(dest, 'coexp2_', targene, '.csv'), quote=F, row.names=T)

#### coexp ####
selectDeg <- function(lfc=1, pv=0.05, mode='ALL') {
  if(mode == 'ALL') {
    deg <- Deg[abs(Deg$logFC) >= lfc & Deg$adj.P.Val <= pv,]
  }
  if(mode == 'UP') {
    deg <- Deg[Deg$logFC >= lfc & Deg$adj.P.Val <= pv,]
  }
  if(mode == 'DOWN') {
    deg <- Deg[Deg$logFC <= -lfc & Deg$adj.P.Val <= pv,]
  }
  return(rownames(deg))
}
set1 <- selectDeg(lfc=1, pv=0.05, mode='ALL')
set1up <- selectDeg(lfc=1, pv=0.05, mode='UP')
set1down <- selectDeg(lfc=1, pv=0.05, mode='DOWN')

selectCor <- function(cutoff=0.7, mode='ALL') {
  if(mode == 'ALL') {
    coexp <- corm[corm$Abs >= cutoff,]
  }
  if(mode == 'UP') {
    coexp <- corm[corm$Value >= cutoff,]
  }
  if(mode == 'DOWN') {
    coexp <- corm[corm$Value <= -cutoff,]
  }
  coexp <- na.omit(coexp)
  return(coexp$Gene)
}
set2 <- selectCor(cutoff=0.6, mode='ALL')
set2up <- selectCor(cutoff=0.6, mode='UP')
set2down <- selectCor(cutoff=0.6, mode='DOWN')

set <- intersect(set1, set2)
setup <- intersect(set1up, set2up)

#### 3.2 cluster ####
gens <- c('ALDH1A1', 'BMI1', 'CD44', 'KLF4', 'MET', 'NANOG', 'POU5F1', 'PROM1', 'SOX2')
gens <- c('MET', 'NANOG', 'SOX2')
gens <- setup
prof <- HNSC[,gens]
clst1 <- kmeans(prof, centers=2, iter.max=10)
myclin$clst1 <- clst1$cluster
table(myclin$maingroup, myclin$clst1)

sur.result <- survfit(Surv(time=OS.time.5.years, event=OS.5.years, 
                           type = 'right')~clst1, data=myclin)
print(ggsurvplot(sur.result, data=myclin, conf.int=F, pval=T, pval.method=T,
                 risk.table='absolute', 
                 legend.title='Group',
                 legend.labs=c('1', '2'),
                 palette=c('blue', 'red'),  
                 title=paste0('Kaplan-Meier Curve of HNSCC Overall Survival'), 
                 xlab='Time (months)',
                 risk.table.height=.30, 
                 ggtheme=theme_bw()+theme(plot.title=element_text(size=20,hjust=0.5),
                                          axis.title.x=element_text(size=18),
                                          axis.title.y=element_text(size=20), axis.text=element_text(size=12), 
                                          axis.text.y=element_text(size=12), legend.text=element_text(size=18),
                                          panel.background=element_blank(), axis.line=element_line(colour='black'),
                                          legend.title=element_text(size=18), text=element_text(face='plain'),
                                          panel.grid.major = element_blank(), panel.grid.minor = element_blank())))

myclin1 <- myclin[myclin$clst1 == 1,]
myclin2 <- myclin[myclin$clst1 == 2,]
myclin1 <- myclin[myclin$maingroup == 'High',]
myclin2 <- myclin[myclin$maingroup == 'Low',]
sur.result <- survfit(Surv(time=OS.time.5.years, event=OS.5.years, 
                           type = 'right')~maingroup, data=myclin1)
print(ggsurvplot(sur.result, data=myclin1, conf.int=F, pval=T, pval.method=T,
                 risk.table='absolute', 
                 legend.title='SOX2',
                 legend.labs=c('High', 'Low'),
                 palette=c('blue', 'red'),  
                 title=paste0('Kaplan-Meier Curve of HNSCC Overall Survival'), 
                 xlab='Time (months)',
                 risk.table.height=.30, 
                 ggtheme=theme_bw()+theme(plot.title=element_text(size=20,hjust=0.5),
                                          axis.title.x=element_text(size=18),
                                          axis.title.y=element_text(size=20), axis.text=element_text(size=12), 
                                          axis.text.y=element_text(size=12), legend.text=element_text(size=18),
                                          panel.background=element_blank(), axis.line=element_line(colour='black'),
                                          legend.title=element_text(size=18), text=element_text(face='plain'),
                                          panel.grid.major = element_blank(), panel.grid.minor = element_blank())))
myclin2$grade <- factor(myclin2$grade, levels = c('G1-G2', 'G3-G4'))
cox <- coxph(Surv(time=OS.time.5.years, event=OS.5.years, 
                  type = 'right')~radiotherapy+chemotherapy+
               gender+age+stage+grade+node, data=myclin2)  # select cov
cox.sum1 <- xtable(summary(cox)$coefficient)
cox.sum2 <- xtable(summary(cox)$conf.int)
cox.sum <- cbind(cox.sum1, cox.sum2)
cox.sum <- cox.sum[c(1:7),]
rownames(cox.sum) <- c('Radiotherapy', 'Chemotherapy', 'Male', 'Age',
                       'Stage III-IV', 'G3-G4', 'Nodal Metastasis')
colnames(cox.sum)[c(2,8:9)] <- c('HR','lower','upper')
cox.sum$var <- rownames(cox.sum)
# write.csv(cox.sum, paste0(dest, '/cox_', targene, '.csv'), quote = F, row.names = F)
pd <- position_dodge(0.1)
cox.sum$Factors <- 'Covariate'
cox.sum$Factors[cox.sum$var == 'G3-G4'] <- 'Grade'
cox.sum$var <- factor(cox.sum$var, levels=rev(c('Male', 'Age', 'Stage III-IV', 'G3-G4', 'Nodal Metastasis', 
                                                'Radiotherapy', 'Chemotherapy')))
print(ggplot(cox.sum, aes(x=var, y=HR, color = Factors)) + 
        ggtitle(paste0('Cox Propotional Hazard Model of HNSCC Overall Survival')) +
        geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2, position=pd, size=1.5) +
        geom_line(position=pd, size = 1.5) + xlab('Variables') + ylab('Hazard Ratio (HR)') +
        geom_hline(yintercept = 1, linetype = 'dashed', color = "blue", size = 1.5) +
        theme_bw() + coord_flip() +
        theme(axis.text.x = element_text(angle = 0, size = 18),
              axis.text.y = element_text(size = 15), axis.line=element_line(colour='black'),
              plot.title = element_text(size = 21, hjust=0.5), axis.title=element_text(size=20),
              legend.text=element_text(size=15), legend.title=element_text(size=18), 
              panel.background=element_blank(), legend.key=element_rect(),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
        geom_point(position=pd, size=3) + 
        scale_color_manual(values = c('black', 'red')))

TSNE <- tsne(X=prof, k=2)
TSNE <- as.data.frame(TSNE)
TSNE$clst1 <- factor(clst1$cluster)
ggplot(data=TSNE, aes(x=V1, y=V2, color=clst1)) + geom_point()

#### 4 pathway ####
GOKEGG <- function(geneset=NULL, output, filename, trend) {
  if(output == 'GO') {
    ego <- enrichGO(gene          = bitr(geneset,
                                         fromType = 'SYMBOL', toType = 'ENTREZID',
                                         OrgDb = 'org.Hs.eg.db')$ENTREZID,
                    # universe      = names(geneList),
                    OrgDb         = org.Hs.eg.db,
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)
    GO_ <- as.data.frame(ego)
#    write.table(GO_, paste0('./output_OSCC/', filename,'_GO_',trend,'.txt'), 
#                sep = '\t', row.names = F, quote = F)
    return(GO_)}
  if(output == 'KEGG') {
    kk <- enrichKEGG(gene          = bitr(geneset,
                                          fromType = 'SYMBOL', toType = 'ENTREZID',
                                          OrgDb = 'org.Hs.eg.db')$ENTREZID,
                     organism      = 'hsa',
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,)
    KEGG_ <- as.data.frame(setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID"))
#    write.table(KEGG_, paste0('./output_OSCC/', filename,'_KEGG_',trend,'.txt'), 
#                sep = '\t', row.names = F, quote = F)
    return(KEGG_)}
}

DrawEnrich <- function(dat, type, top.number = 5, col="blue", trend){
  if(type == 'BP') {tit <- 'Biological Process of '}
  if(type == 'CC') {tit <- 'Cellular Component of '}
  if(type == 'MF') {tit <- 'Molecular Function of '}
  if(type == 'KEGG') {tit <- 'KEGG Pathway of '}
  if(type == 'GO') {tit <- 'Gene Ontology Enrichment of '}
  dat1 = dat[c(1:top.number),]
  dat1$Description = capitalize(dat1$Description)
  dat1$Description = factor(dat1$Description,levels=dat1$Description[length(dat1$Description):1])
  dat1$PValue = -log10(dat1$p.adjust)
  dat1$GeneRatio <- dat1$Count / as.numeric(gsub('^.*/', '', dat1$GeneRatio))
  dat1$Description <- gsub(', ', '\n', dat1$Description)
  p = ggplot(dat1,aes(GeneRatio, Description)) +
    geom_point(aes(size=Count,colour=PValue)) +
    scale_colour_gradient(low=col,high="red") + 
    labs(colour=expression(-log[10]("P Value")),size="Gene counts",  
         x="Gene Ratio",y="",title=paste0(tit, 'Top ', trend, ' Correlated Genes')) +
    theme_bw() + theme(axis.text=element_text(size=14),plot.title=element_text(size=18,hjust=0.5),
                       legend.text=element_text(size=12)) +
    scale_x_continuous(limits = c(0,max(dat1$GeneRatio) * 1.2)) 
  return(p)
}

GO <- GOKEGG(geneset=setup, output='GO', filename=disp.name, trend='pos50')
KEGG <- GOKEGG(geneset=setup, output='KEGG', filename=disp.name, trend='pos50')
DrawEnrich(GO, 'GO', min(nrow(GO), 5), 'blue', '50 Positively')
DrawEnrich(KEGG, 'KEGG', min(nrow(KEGG), 5), 'blue', '50 Positively')

GO <- GOKEGG(geneset=setdown, output='GO', filename='TCGA', trend='DOWN')
KEGG <- GOKEGG(geneset=setdown, output='KEGG', filename='TCGA', trend='DOWN')
DrawEnrich(GO, 'GO', 5, 'blue', 'Down')
DrawEnrich(KEGG, 'KEGG', 5, 'blue', 'Down')
