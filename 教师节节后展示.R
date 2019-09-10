#start
rm(list=ls())
setwd('D:/work/TCGAbiolinks') # 工作路径

library(TCGAbiolinks)
library(fdrtool)
library(ggplot2)
library(plyr)
library(SIMLR)
library(survival)
library(survminer)
library(xtable)
library(Hmisc)

# 下载
cancer.list <- c('HNSC') # 需要的癌症种类（头颈癌）
for(cancer in cancer.list) {
  project.title <- paste('TCGA', cancer, sep = '-')
  
  # 下载完整临床数据
  
  query <- GDCquery(project = project.title,
                    data.category = 'Clinical',
                    file.type = 'xml')
  GDCdownload(query)
  clinical.list <- c('drug', 'admin', 'follow_up', 'radiation', 'patient', 
                     'stage_event', 'new_tumor_event')
  for(clinical.file in clinical.list) {
    clin <- GDCprepare_clinic(query, clinical.file)
    if(clinical.file == 'drug') {clinical <- clin}
    else {clinical <- merge(clin, clinical)}
    write.table(clin, paste(cancer, '_clinical_', clinical.file, '.txt', sep = ''), 
                quote = F, sep = '\t', row.names = F)
  }
  write.table(clinical, paste(cancer, '_clinical.txt', sep = ''), 
              quote = F, sep = '\t', row.names = F)
  
  # 下载rna-seq的counts数据
  query <- GDCquery(project = project.title, 
                    data.category = 'Transcriptome Profiling', 
                    data.type = 'Gene Expression Quantification', 
                    workflow.type = 'HTSeq - FPKM')
  
  GDCdownload(query, method = 'api', files.per.chunk = 2)
  expdat <<- GDCprepare(query = query)
  count_matrix=assay(expdat)
  write.table(count_matrix, paste(cancer, '_counts.txt', sep = ''), sep = '\t', quote = F, 
              row.names = T)
}
# 下载

# def
read.expr <- function(cancer) {
  expr <- read.table(paste(cancer, '_HiSeqV2', sep=''), T, sep='\t', 
                     stringsAsFactors = F, row.names = 'sample')
  return(expr)
}

# 数据清洗
read.clin <- function(cancer) {
  clinical <- read.table(paste(cancer,'_clinical.txt', sep=''), 
                         T, sep='\t', stringsAsFactors = F)
  clinical$sampleID <- gsub('-', '.', clinical$sampleID)
  clinical[clinical == ''] <- NA
  clinical[clinical == '[Discrepancy]'] <- NA
  # survival
  clinical$OS <- as.numeric(clinical$OS)
  clinical$OS.time <- as.numeric(clinical$OS.time)
  clinical$OS.time.month <- as.numeric(clinical$OS.time/30)
  clinical$OS.5.years <- clinical$OS
  clinical$OS.5.years[clinical$OS.time.month>60] <- 0
  clinical$OS.time.5.years <- clinical$OS.time.month
  clinical$OS.time.5.years[clinical$OS.time.month>60] <- 60
  # stage and grade
  clinical$stage <- NA
  clinical$stage[substring(clinical$clinical_stage,1,8)=='Stage II'|
                 substring(clinical$clinical_stage,1,7)=='Stage I'] <- '1-2'
  clinical$stage[substring(clinical$clinical_stage,1,9)=='Stage III'|
                 substring(clinical$clinical_stage,1,8)=='Stage IV'] <- '3-4'
  clinical$grade[clinical$neoplasm_histologic_grade=='GX'] <- NA
  clinical$grade[clinical$neoplasm_histologic_grade=='G1'|
                   clinical$neoplasm_histologic_grade=='G2'] <- 'high'
  clinical$grade[clinical$neoplasm_histologic_grade=='G3'|
                   clinical$neoplasm_histologic_grade=='G4'] <- 'low'
  return(clinical)
}

# 匹配与临床资料ID相同的基因表达值
design <- function(cancer) {
  sample <- Clin$sampleID[!is.na(Clin$grade)]
  sample <- intersect(sample, colnames(Expr))
  grade <- Clin$grade[match(sample, Clin$sampleID)]
  expr <- Expr[, sample]
  des <- data.frame(sample, grade)
  colnames(des) <- c('sample', 'grade')
  write.table(des, paste('./saves/design_of_', cancer, '.txt', sep=''), 
              row.names=F, quote=F, sep='\t')
  return(expr)
}

# 差异性表达基因(DEG)的筛选
diff.gene <- function(cancer, des, expr, lfc_cutoff=1.5, alpha=0.05) {
  high <- as.character(des$sample[des$grade == 'high'])
  low <- as.character(des$sample[des$grade == 'low'])
  p <- rep(1, nrow(expr))
  lfc <- 1:nrow(expr)
  gene <- 1:nrow(expr)
  for(i in 1:nrow(expr)){
    a <- t(expr[i, na.omit(match(high, colnames(expr)))]) # expression of high grade
    b <- t(expr[i, na.omit(match(low, colnames(expr)))]) # expression of low grade
    if(all(!is.na(a)) & all(!is.na(b))) {
      if(sd(a, na.rm = T) != 0 & sd(b, na.rm = T) != 0) {
        t <- t.test(a, b, na.rm = T)
        p[i] <- t$p.value
        gene[i] <- rownames(Expr)[i]
        lfc[i] <- mean(b, na.rm = T) - mean(a, na.rm = T)
        # log2 fold change = log2(low + 1) - log2(high + 1), where '+ 1' is a 
        # widely used modification
      }
    }
  }
  p.adj <- fdrtool(p, statistic="pvalue", plot = F) # Benjamini-Hochberg method
  results <- data.frame(gene, lfc, p, p.adj$qval)
  with(results, {
    results$trend <<- as.factor(ifelse(p.adj.qval <= alpha & abs(lfc) >= lfc_cutoff, 
                                       ifelse(lfc >= lfc_cutoff ,'UP','DOWN'),'NOT'))
  })
  deg <- results[(abs(results$lfc)>=lfc_cutoff & results$p.adj.qval <= alpha),]
  write.table(results, paste('./results/results_of_', cancer, '.txt', sep=''), 
              sep='\t', quote = F, row.names = F)
  write.table(deg, paste('./results/deg_of_', cancer, '.txt', sep=''), 
              sep='\t', quote = F, row.names = F)
  return(results)
}

find.deg <- function(cancer) {
  deg <- read.table(paste('./results/deg_of_', cancer, '.txt', sep=''), 
                    T, sep='\t', stringsAsFactors = F)
  return(deg)
}

# 差异性表达基因(DEG)结果可视化
volcano <- function(cancer, res, lfc_cutoff=1.5, alpha=0.05) {
  this_title <<- paste0('The volcano plot of ',cancer,
                        '\nCutoff for logFC is ',round(lfc_cutoff, 1),
                        '\nUp-regulated gene number: ',nrow(res[res$trend =='UP',]) ,
                        '\nDown-regulated gene number: ',nrow(res[res$trend =='DOWN',]))
  # find top 10 differently expressed genes and mark them
  Deg <- arrange(Deg, abs(Deg$lfc), decreasing = T)
  label <- Deg$gene[1:10]
  res <- res[res$p != 1,]
  res$label <- ''
  res$label[res$gene %in% label] <- as.character(label)
  g = ggplot(data=res, aes(x=lfc, y=-log10(p.adj.qval), color=trend)) +
    geom_point(alpha=0.4, size=1.75) +geom_text(aes(label = label), size=3,vjust=-0.5, alpha=1) +
    theme_set(theme_set(theme_bw(base_size=13)))+
    xlab("log2 fold change") + ylab("-log10 p-value") +
    ggtitle( this_title ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
    scale_colour_manual(values = c('blue','grey','red'))
  g=g+geom_vline(xintercept = -lfc_cutoff, linetype = "dashed", color = "grey", size = 1) +
    geom_vline(xintercept = lfc_cutoff, linetype = "dashed", color = "grey", size = 1)+
    geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "grey", size = 1)
  return(g)
}

# SIMLR
simlr <- function(genes='all', ncluster=2) {
  if(genes == 'all') {
    Exp.of.deg <- Expr
  }
  if(genes == 'deg') {
    Exp.of.deg <- Expr[match(Deg$gene, rownames(Expr)),]
  }
  Exp.of.deg <- na.omit(Exp.of.deg)
  sim <- SIMLR(X = Exp.of.deg, c = ncluster) # run SIMLR 
  sim_plot <- data.frame(x = sim$ydata[,1], y = sim$ydata[,2], grade = Des$grade, 
                         cluster = sim$y[1])
  Des$cluster <<- as.factor(unlist(sim$y[1])) # 记录分型
  Clin$cluster <<- NA
  Clin$cluster[match(Des$sample, Clin$sampleID)] <<- Des$cluster
  g <- ggplot(data = sim_plot, aes(x = x, y = y, color = grade)) +
    geom_point(alpha = 0, size = 1.75) + geom_text(aes(label = cluster), alpha = 1) +
    theme_set(theme_set(theme_bw(base_size = 13)))+
    xlab("SIMLR component 1") + ylab("SIMLR component 2") +
    ggtitle('SIMLR Visualization') + theme(plot.title = element_text(size = 15,hjust = 0.5))+
    scale_colour_manual(values = c('blue','red'))
  return(g)
}

# 以分型为组行生存分析
prog <- function(cancer, group='radiotherapy') {
  sur <<- Surv(time = Clin$OS.time.5.years, event = Clin$OS.5.years, type = 'right')
  if(group == 'radiotherapy') {
    sur.result <<- survfit(sur ~ cluster + radiation_therapy, data = Clin) # 比较放疗及分型
    g <- ggsurvplot(sur.result, data = Clin, conf.int = F, pval = T, pval.method = T,
                    risk.table = 'absolute', 
                    legend.labs = c(paste('cluster 1 \nwithout ', group, sep = ''),
                                    paste('cluster 1 \nwith ', group, sep = ''),
                                    paste('cluster 2 \nwithout ', group, sep = ''),
                                    paste('cluster 2 \nwith ', group, sep = '')
                    ), 
                    legend.title = "group",
                    palette = c("blue", "red", 'dark green', 'purple'),  
                    title = paste("Kaplan-Meier Curve for ", cancer, " Survival", sep = ''), 
                    xlab = 'Time(months)',
                    risk.table.height = .30)
  }
  return(g)
}

cox.multi <- function(clst) {
  cox <- coxph(sur ~ radiation_therapy + 
                 age_at_initial_pathologic_diagnosis + gender + stage, data = Clin, 
               subset = Clin$cluster == clst)
  cox.sum1 <- xtable(summary(cox)$coefficient)
  cox.sum2 <- xtable(summary(cox)$conf.int)
  cox.sum <- cbind(cox.sum1, cox.sum2)
  rownames(cox.sum) <- c('Radiotherapy', 'Age', 'Male', 'Stage 3-4')
  colnames(cox.sum)[c(2,8:9)] <- c('HR','lower','upper')
  cox.sum$var <- rownames(cox.sum)
  return(cox.sum)
}

cox.plot <- function(subtype) {
  pd <- position_dodge(0.1)
  cox.summary <- cox.multi(subtype)
  g <- ggplot(cox.summary, aes(x = var, y = HR, color = 'orange')) + 
    geom_point(size = 1.5) + ggtitle('Cox Propotional Hazard Regression Model in cluster 1') +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.6, position = pd) +
    geom_line(position = pd, size = 0.75) + xlab('Variables') + ylab('Hazard Ratio (HR)') +
    geom_hline(yintercept = 1, linetype = 'dashed', color = "blue", size = 0.1) +
    theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.6)) +
    geom_point(position = pd)
  return(g)
}

# main
Expr <- read.expr('HNSC')
Clin <- read.clin('HNSC')
Expr <- design('HNSC')
Des <- read.table(paste('./saves/design_of_', cancer, '.txt', sep=''), 
                  T, sep='\t', stringsAsFactors = F)
Res <- diff.gene('HNSC', Des, Expr, 1, 0.05)
Deg <- find.deg('HNSC')
volcano('HNSC', Res, 1, 0.05)
simlr(genes = 'all', ncluster = 2) # genes can be 'all' or 'deg'
prog('HNSC', 'radiotherapy')
cox.plot(subtype=1) # subtype is an integer
cox.plot(subtype=2)

# KEGG pathway analysis
DrawGOBubblePlot <- function(dat, dat_name, top.number = 10, col="blue"){
  dat1 = dat[c(1:top.number),c(2,3,4,5)]
  colnames(dat1)[3] = 'GeneRatio'
  dat1$Term = substr(dat1$Term,10,200)
  dat1$Term = capitalize(dat1$Term)
  dat1$Term = factor(dat1$Term,levels=dat1$Term[length(dat1$Term):1])
  dat1$PValue = -log10(dat1$PValue)
  p = ggplot(dat1,aes(GeneRatio,Term)) +
    geom_point(aes(size=Count,colour=PValue)) +
    scale_colour_gradient(low=col,high="red") + 
    labs(colour=expression(-log[10]("P Value")),size="Gene counts",  
         x="Gene Ratio",y="",title=paste('KEGG Pathway (', dat_name, ')', sep='')) +
    theme_bw() +
    scale_x_continuous(limits = c(0,max(dat1$GeneRatio) * 1.2)) 
  return(p)
}
# This file should be generated in advance
KEGG <- read.table('./results/KEGG.txt', T, sep='\t', stringsAsFactors = F)
DrawGOBubblePlot(GO, 'OSCC', 10, "blue")
