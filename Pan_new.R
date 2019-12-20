#start
setwd('D:/work/Pancancer_new/')
library(ggplot2)
library(xtable)
library(survminer)
library(survival)
library(limma)
library(VennDiagram)
# library(TCGAbiolinks)

# exp
read.expr <- function() {
  expr <- read.table('../HNSC_mining/reads/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena',
                     T, sep='\t', stringsAsFactors = F)
  expr <- expr[grep('^[A-Za-z].+$', expr$sample), ]
  expr <- expr[!duplicated(expr$sample), ]
  row.names(expr) <- expr$sample
  expr <- expr[, -1]
  return(expr)
}
Expr <- read.expr()

"
# download clin
cancer.list <- c('BLCA','CESC','HNSC','KIRC','LGG','LIHC','OV','STAD','UCEC')
cancer.list <- c('UCEC')
for(cancer in cancer.list) {
  project.title <- paste('TCGA', cancer, sep = '-')
  query <- GDCquery(project = project.title,
                    data.category = 'Clinical',
                    file.type = 'xml')
  GDCdownload(query)
  clinical.list <- c('admin', 'drug', 'follow_up', 'new_tumor_event', 'patient',
                     'radiation', 'stage_event')
  for(clinical.file in clinical.list) {
    clin <- GDCprepare_clinic(query, clinical.file)
    write.table(clin, paste0('./raw_clins/', cancer, '_clinical_', clinical.file, '.txt'), 
                quote = F, row.names = F, sep = '\t')
  }
}
"

# 
read.clin <- function(cancer) {
  clin.radiation <- read.table(paste0('./raw_clins/', cancer, '_clinical_radiation.txt'), 
                               T, sep = '\t', stringsAsFactors = F, quote = "")
  clin.drug <- read.table(paste0('./raw_clins/', cancer, '_clinical_drug.txt'), 
                          T, sep = '\t', stringsAsFactors = F, quote = "", fill = T)
  clin.followup <- read.table(paste0('./raw_clins/', cancer, '_clinical_follow_up.txt'), 
                              T, sep = '\t', stringsAsFactors = F, quote = "")
  clin.patient <- read.table(paste0('./raw_clins/', cancer, '_clinical_patient.txt'), 
                             T, sep = '\t', stringsAsFactors = F, quote = "")
  clinical <- unique(clin.patient)
  clinical[clinical == ''] <- NA
  if(cancer == 'OSCC') {
    clinical <- clinical[clinical$anatomic_neoplasm_subdivision %in%
                           c('Base of tongue', 'Oral Tongue', 'Alveolar Ridge',
                             'Buccal Mucosa', 'Floor of mouth', 'Oral Cavity', 'Lip',
                             'Hard Palate'),] # 选部位
  }
  clinical$sample <- paste0(clinical$bcr_patient_barcode, '-01') # 使之与tumor组织匹配
  clinical$sample <- gsub('-', '.', clinical$sample) # replace char
  # grade
  clinical$grade <- clinical$neoplasm_histologic_grade
  clinical$grade[clinical$neoplasm_histologic_grade=='G1'|
                 clinical$neoplasm_histologic_grade=='G2'|
                 clinical$neoplasm_histologic_grade=='GB'] <- 'Low Grade'
  # GB 是 OV 特有的
  clinical$grade[clinical$neoplasm_histologic_grade=='G3'|
                 clinical$neoplasm_histologic_grade=='G4'] <- 'High Grade'
  clinical$grade[clinical$neoplasm_histologic_grade=='GX'] <- NA
  clinical <- clinical[!is.na(clinical$grade),]
  # stage
  num_clin <- length(na.omit(clinical$stage_event_clinical_stage))
  num_path <- length(na.omit(clinical$stage_event_pathologic_stage))
  if(cancer != 'LGG') {
    if(num_clin >= num_path) {
      clinical$stage <- clinical$stage_event_clinical_stage
      clinical$stage[grepl('^Stage I{1,2}[A-C]?[0-9]?$', 
                           clinical$stage_event_clinical_stage)] <- 'Early Stage'
      clinical$stage[grepl('^Stage I(II|V)[A-C]?[0-9]?$', 
                           clinical$stage_event_clinical_stage)] <- 'Late Stage'
    }else {
      clinical$stage <- clinical$stage_event_pathologic_stage
      clinical$stage[grepl('^Stage I{1,2}[A-C]?[0-9]?$', 
                           clinical$stage_event_pathologic_stage)] <- 'Early Stage'
      clinical$stage[grepl('^Stage I(II|V)[A-C]?[0-9]?$', 
                           clinical$stage_event_pathologic_stage)] <- 'Late Stage'
    }
  }
  
  "
  # exposure
  clinical$smoke <- NA
  clinical$smoke[clinical$tobacco_smoking_history==1] <- 'NO'
  clinical$smoke[!is.na(clinical$tobacco_smoking_history) & 
                          clinical$tobacco_smoking_history!=1] <- 'YES'
  clinical$drink <- NA
  clinical$drink <- clinical$alcohol_history_documented
  "
  # therapy
  clinical$chemotherapy <- 'NO'
  drugs <- unique(clin.drug$bcr_patient_barcode)
  clinical$chemotherapy[clinical$bcr_patient_barcode %in% drugs] <- 'YES'
  clinical$radiotherapy <- 'NO'
  radiation <- unique(clin.radiation$bcr_patient_barcode)
  clinical$radiotherapy[clinical$bcr_patient_barcode %in% radiation] <- 'YES'
  # survival
  clin.followup[clin.followup == ''] <- NA
  clin.followup <- clin.followup[!is.na(clin.followup$vital_status),]
  clin.followup$OS.time <- clin.followup$days_to_death
  clin.followup$OS.time[is.na(clin.followup$OS.time)] <- 
    clin.followup$days_to_last_followup[is.na(clin.followup$OS.time)] # 合并alive与dead的时间
  bcrcode <- unique(clin.followup$bcr_patient_barcode)
  blank <- rep('', length(bcrcode))
  followup <<- data.frame(blank, blank, blank, stringsAsFactors = F)
  colnames(followup) <<- c('code', 'status', 'time')
  for(i in 1:length(bcrcode)) {  # 选择较后的随访（时间和状态）
    time <- max(clin.followup$OS.time[clin.followup$bcr_patient_barcode == bcrcode[i]])
    if(all(clin.followup$vital_status[clin.followup$bcr_patient_barcode == bcrcode[i]] 
           == 'Alive')) {status <- 'Alive'}
    else {status <- 'Dead'}
    followup$code[i] <<- bcrcode[i]
    followup$status[i] <<- status
    followup$time[i] <<- time
  }
  clinical$OS <- NA
  clinical$OS <- followup$status[match(clinical$bcr_patient_barcode, followup$code)]
  clinical$OS.time <- NA
  clinical$OS.time <- followup$time[match(clinical$bcr_patient_barcode, followup$code)]
  clinical$OS.time.month <- as.numeric(clinical$OS.time)/30
  clinical$OS.5.years <- clinical$OS
  clinical$OS.5.years[clinical$OS.time.month>60] <- 'Alive'
  clinical$OS.time.5.years <- clinical$OS.time.month
  clinical$OS.time.5.years[clinical$OS.time.month>60] <- 60
  if(cancer == 'LGG') {
    ext <- c('sample', 'bcr_patient_barcode', 'gender', 
             'age_at_initial_pathologic_diagnosis',
             'grade', 'neoplasm_histologic_grade',
             'chemotherapy', 'radiotherapy', 'OS', 'OS.time',
             'OS.time.month', 'OS.time.5.years', 'OS.5.years')
  }else {
    ext <- c('sample', 'bcr_patient_barcode', 'gender', 
             'age_at_initial_pathologic_diagnosis',
             'grade', 'neoplasm_histologic_grade',
             'stage', 'stage_event_clinical_stage', 'stage_event_pathologic_stage',
             'chemotherapy', 'radiotherapy', 'OS', 'OS.time',
             'OS.time.month', 'OS.time.5.years', 'OS.5.years')
  }
  clinical <- clinical[, ext]
  colnames(clinical)[4] <- 'age'
  clinical <- clinical[!duplicated(clinical$sample),]
  
  write.table(clinical, paste0('./clean_clins/clinical_grade_', cancer, '.txt'), 
              sep='\t', quote=F, row.names=F)
  return(clinical)
}
cancer.list <- c('BLCA','CESC','HNSC','KIRC','LGG','LIHC','OV','STAD','UCEC')
read.clin(cancer.list[7])
for(cancer in cancer.list) {
  Clin <- read.clin(cancer)
}
rm(Clin)

# post clinical
blca <- read.table('./clean_clins/clinical_grade_BLCA.txt', T, stringsAsFactors = F, sep = '\t')
blca$age <- as.numeric(blca$age)
blca$age[blca$age > 200] <- NA
blca <- blca[blca$grade != 'Bladder - NOS',]
write.table(blca, './clean_clins/clinical_grade_BLCA.txt', row.names = F, quote = F, sep = '\t')
rm(blca)

load.Clin <- function(cancer) {
  clin <- read.table(paste0('./clean_clins/clinical_grade_', cancer, '.txt'), 
                     T, stringsAsFactors = F, sep = '\t')
  return(clin)
}





# deg grade
get.deg <- function(cancer, lfc_cutoff=1, alpha=0.05) {
  results <- topTable(fit2,adjust.method="BH",coef=1,p.value=1,
                      lfc=log(1,2),number=Inf,sort.by = 'logFC')
  diff <- topTable(fit2,adjust.method="BH",coef=1,p.value=alpha,
                   lfc=lfc_cutoff,number=Inf,sort.by = 'logFC')
  results$trend <- as.factor(ifelse(results$adj.P.Val <= alpha & abs(results$logFC)>=lfc_cutoff, 
                                    ifelse(results$logFC >= lfc_cutoff ,'UP','DOWN'),'NOT'))
  diff$trend <- as.factor(ifelse(diff$adj.P.Val <= alpha & abs(diff$logFC) >= lfc_cutoff, 
                                 ifelse(diff$logFC >= lfc_cutoff ,'UP','DOWN'),'NOT'))
  write.csv(results, paste0('./deg1/results_',cancer,'.csv'), row.names = T, quote = F)
  write.csv(diff, paste0('./deg1/degs_',cancer,'.csv'), row.names = T, quote = F)
}

deg1 <- function(cancer) {
  sams <- intersect(colnames(Expr), Clin$sample)
  Clin <- Clin[Clin$sample %in% sams, ]
  Exprc <- Expr[,Clin$sample]
  # 使用limma要注意“顺序必须是一样的”
  # design
  Group <- factor(Clin$grade, levels=c('Low Grade', 'High Grade'))
  design <- model.matrix(~0+Group)
  colnames(design) <- c('LowGrade','HighGrade')
  # diff
  fit <- lmFit(Exprc, design)
  contrast.matrix <- makeContrasts(HighGrade-LowGrade, levels=design)
  fit2 <<- contrasts.fit(fit, contrast.matrix)
  fit2 <<- eBayes(fit2)
  get.deg(cancer, 1, 0.05)
  Res <<- read.csv(paste0('./deg1/results_',cancer,'.csv'), T, stringsAsFactors = F)
  Deg <<- read.csv(paste0('./deg1/degs_',cancer,'.csv'), T, stringsAsFactors = F)
}
cancer.list <- c('BLCA','CESC','HNSC','KIRC','LGG','LIHC','OV','STAD','UCEC')
Clin <- load.Clin(cancer.list[7])
# OV wu DEG
deg1(cancer.list[7])



# deg tumor normal
get.deg2 <- function(cancer, lfc_cutoff=1, alpha=0.05) {
  results <- topTable(fit2,adjust.method="BH",coef=1,p.value=1,
                      lfc=log(1,2),number=Inf,sort.by = 'logFC')
  diff <- topTable(fit2,adjust.method="BH",coef=1,p.value=alpha,
                   lfc=lfc_cutoff,number=Inf,sort.by = 'logFC')
  results$trend <- as.factor(ifelse(results$adj.P.Val <= alpha & abs(results$logFC)>=lfc_cutoff, 
                                    ifelse(results$logFC >= lfc_cutoff ,'UP','DOWN'),'NOT'))
  diff$trend <- as.factor(ifelse(diff$adj.P.Val <= alpha & abs(diff$logFC) >= lfc_cutoff, 
                                 ifelse(diff$logFC >= lfc_cutoff ,'UP','DOWN'),'NOT'))
  write.csv(results, paste0('./deg2/results_',cancer,'.csv'), row.names = T, quote = F)
  write.csv(diff, paste0('./deg2/degs_',cancer,'.csv'), row.names = T, quote = F)
}
deg2 <- function(cancer) {
  sams <- intersect(colnames(Expr), Clin$sample)
  Clin <- Clin[Clin$sample %in% sams, ]
  Exprc <- Expr[,Clin$sample]
  sams_nom <- gsub('[.]01$', '.11', sams)
  sams_ava_nom <- intersect(sams_nom, colnames(Expr))
  ava <- c(sams_ava_nom, gsub('[.]11$', '.01', sams_ava_nom))
  Exprc <- Expr[, ava]
  Group <- c(rep('Normal', length(sams_ava_nom)), rep('Tumor', length(sams_ava_nom)))
  Group <- factor(Group, levels = c('Tumor', 'Normal'))
  design <- model.matrix(~0+Group)
  colnames(design) <- c('Tumor','Normal')
  # diff
  fit <- lmFit(Exprc, design)
  contrast.matrix <- makeContrasts(Tumor-Normal, levels=design)
  fit2 <<- contrasts.fit(fit, contrast.matrix)
  fit2 <<- eBayes(fit2)
  get.deg2(cancer, 1, 0.05)
  Res <<- read.csv(paste0('./deg2/results_',cancer,'.csv'), T, stringsAsFactors = F)
  Deg <<- read.csv(paste0('./deg2/degs_',cancer,'.csv'), T, stringsAsFactors = F)
}
cancer.list <- c('BLCA','CESC','HNSC','KIRC','LGG','LIHC','OV','STAD','UCEC')
Clin <- load.Clin(cancer.list[5])
# LGG and OV no normal comparison
deg2(cancer.list[5])


# deg prognosis
cox.gene <- function(gene, cut=0.5, cancer) {
  sams <- intersect(colnames(Expr), Clin$sample)
  Clin <- Clin[Clin$sample %in% sams, ]
  Exprc <- Expr[,Clin$sample]
  Clin$exp <- as.numeric(c(Exprc[gene,]))
  Clin$group <- ifelse(Clin$exp >= quantile(Clin$exp, cut, na.rm = T), 'High', 'Low')
  Clin$group <- factor(Clin$group, levels = c('Low', 'High'))
  Clin$OS.5.years[Clin$OS.5.years=='Alive'] <- 0
  Clin$OS.5.years[Clin$OS.5.years=='Dead'] <- 1
  Clin$OS.5.years <- as.numeric(Clin$OS.5.years)
  if(cancer %in% c('CESC', 'OV', 'UCEC')) {
    cox <- coxph(Surv(Clin$OS.time.5.years, Clin$OS.5.years, type = 'right') ~ 
                   radiotherapy+chemotherapy+stage+grade+age+group, data = Clin)
    return(c(xtable(cox)$p[6], xtable(cox)$`exp(coef)`[6]))
  }else if(cancer %in% c('LGG')) {
    cox <- coxph(Surv(Clin$OS.time.5.years, Clin$OS.5.years, type = 'right') ~ 
                   radiotherapy+chemotherapy+grade+age+gender+group, data = Clin)
    return(c(xtable(cox)$p[6], xtable(cox)$`exp(coef)`[6]))
  }else {
    cox <- coxph(Surv(Clin$OS.time.5.years, Clin$OS.5.years, type = 'right') ~ 
                   radiotherapy+chemotherapy+stage+grade+age+gender+group, data = Clin)
    return(c(xtable(cox)$p[7], xtable(cox)$`exp(coef)`[7]))
  }
}
gg <- rownames(Expr)
deg3 <- function(cancer) {
  phr <- sapply(X = gg, FUN = cox.gene, cancer=cancer)
  Prog <- as.data.frame(gg)
  Prog$cox.P <- phr[1,]
  Prog$HR <- phr[2,]
  colnames(Prog)[1] <- 'gene'
  write.csv(Prog, paste0('./deg3/gene_cox_', cancer, '.csv'), quote = F, row.names = F)
  Sig <- Prog[!is.na(Prog$cox.P) & Prog$cox.P < 0.05,]
  write.csv(Sig, paste0('./deg3/sig_cox_', cancer, '.csv'), quote = F, row.names = F)
}
cancer.list <- c('BLCA','CESC','HNSC','KIRC','LGG','LIHC','OV','STAD','UCEC')
Clin <- load.Clin(cancer.list[3])
# 'OV' and "STAD' removed genes that had NA
deg3(cancer.list[3])

gen <- rownames(Expr)
gg <- rownames(na.omit(Exprc[!any(NA %in% c(Exprc[gen, ])), ]))

rm(gen, gg)

# create list
genList <- function(cancer) {
  res1 <- read.csv(paste0('./deg1/results_',cancer,'.csv'), T, stringsAsFactors = F)
  res2 <- read.csv(paste0('./deg2/results_',cancer,'.csv'), T, stringsAsFactors = F)
  res3 <- read.csv(paste0('./deg3/gene_cox_',cancer,'.csv'), T, stringsAsFactors = F)
  res <- res1[, c(1,2,5,6,8)]
  colnames(res)[2:5] <- paste0(colnames(res)[2:5], '_1')
  res <- merge(res, res2[, c(1,2,5,6,8)], by.x = 'X')
  colnames(res)[6:9] <- paste0(colnames(res)[6:9], '_2')
  colnames(res3)[1] <- 'X'
  res <- merge(res, res3, by.x = 'X', all.x = T)
  write.csv(res, paste0('./deg_sum/degs_', cancer, '.csv'), row.names = F, quote = F)
}
cancer.list <- c('BLCA','CESC','HNSC','KIRC','LGG','LIHC','OV','STAD','UCEC')
cancer.list <- c('BLCA','CESC','HNSC','KIRC','LIHC','STAD','UCEC')
# OV and LGG no normal comparsion
genList(cancer.list[3])

# Venn
drawVenn <- function(cancer) {
  deg1 <- read.csv(paste0('./deg1/degs_',cancer,'.csv'), T, stringsAsFactors = F)$X
  deg2 <- read.csv(paste0('./deg2/degs_',cancer,'.csv'), T, stringsAsFactors = F)$X
  deg3 <- read.csv(paste0('./deg3/sig_cox_',cancer,'.csv'), T, stringsAsFactors = F)$gene
  venn.diagram(x = list(By_grade = deg1, 
                        By_tumor = deg2,
                        By_prognosis_value = deg3), 
               fill = c('yellow', 'purple', 'red'),
               col = 'transparent',
               main = paste0('Common DEGs Identified Using Different Comparison Methods in ',
                             cancer),
               filename = paste0('./deg_sum/Venn_', cancer, '.tiff'), 
               imagetype = 'tiff', cat.pos = c(0, 0, 180))
}
cancer.list <- c('BLCA','CESC','HNSC','KIRC','LIHC','STAD','UCEC')
drawVenn(cancer.list[7])

saveList <- function(cancer) {
  deg1 <- read.csv(paste0('./deg1/degs_',cancer,'.csv'), T, stringsAsFactors = F)$X
  deg2 <- read.csv(paste0('./deg2/degs_',cancer,'.csv'), T, stringsAsFactors = F)$X
  deg3 <- read.csv(paste0('./deg3/sig_cox_',cancer,'.csv'), T, stringsAsFactors = F)$gene
  coms <- intersect(deg1, intersect(deg2, deg3))
  res <- read.csv(paste0('./deg_sum/degs_',cancer,'.csv'),T,stringsAsFactors=F,row.names='X')
  rescom <- res[rownames(res) %in% coms, ]
  write.csv(rescom, paste0('./deg_sum/coms_', cancer, '.csv'), row.names = T, quote = F)
}
cancer.list <- c('BLCA','CESC','HNSC','KIRC','LIHC','STAD','UCEC')
saveList(cancer.list[7])

# cluster and cox plot
set.seed(42)
km.clst <- function(cancer) {
  deg <- read.csv(paste0('./deg_sum/coms_',cancer,'.csv'),T,stringsAsFactors=F)$X
  sams <- intersect(colnames(Expr), Clin$sample)
  Clin <<- Clin[Clin$sample %in% sams, ]
  Exprc <- Expr[deg, Clin$sample]
#  Exprc <- na.omit(Exprc)
  km <- kmeans(x=t(Exprc), center=2, iter.max=10)
  Clin$cluster <<- as.numeric(km$cluster)
  print(table(Clin$grade, Clin$cluster))
}

cox.km <- function(kmclst, cancer, def_high_grade) {
  Clin$OS.5.years[Clin$OS.5.years=='Alive'] <- 0
  Clin$OS.5.years[Clin$OS.5.years=='Dead'] <- 1
  Clin$OS.5.years <- as.numeric(Clin$OS.5.years)
  if(kmclst == def_high_grade) {lab <- 'High-like Grade'}else {lab <- 'Low-like Grade'}
  if(cancer %in% c('CESC', 'OV', 'UCEC')) {
    cox <- coxph(Surv(time = Clin$OS.time.5.years, event = Clin$OS.5.years, 
                      type = 'right')~radiotherapy+chemotherapy+
                   age+stage, data=Clin, 
                 subset = Clin$cluster == kmclst)
    cox.sum1 <- xtable(summary(cox)$coefficient)
    cox.sum2 <- xtable(summary(cox)$conf.int)
    cox.sum <- cbind(cox.sum1, cox.sum2)
    cox.sum <- cox.sum[1:4,]
    rownames(cox.sum) <- c('Radiation', 'Chemotherapy', 'Age', 'Stage III-IV')
    colnames(cox.sum)[c(2,8:9)] <- c('HR','lower','upper')
    cox.sum$var <- rownames(cox.sum)
  }else {
    cox <- coxph(Surv(time = Clin$OS.time.5.years, event = Clin$OS.5.years, 
                      type = 'right')~radiotherapy+chemotherapy+
                   gender+age+stage, data=Clin, 
                 subset = Clin$cluster == kmclst)
    cox.sum1 <- xtable(summary(cox)$coefficient)
    cox.sum2 <- xtable(summary(cox)$conf.int)
    cox.sum <- cbind(cox.sum1, cox.sum2)
    cox.sum <- cox.sum[1:5,]
    rownames(cox.sum) <- c('Radiation', 'Chemotherapy', 'Male', 'Age', 'Stage III-IV')
    colnames(cox.sum)[c(2,8:9)] <- c('HR','lower','upper')
    cox.sum$var <- rownames(cox.sum)
  }
  write.csv(cox.sum, paste0('./cluster/cox_', cancer, kmclst, '.csv'), quote = F, row.names = F)
  pd <- position_dodge(0.1)
  g <- ggplot(cox.sum, aes(x=var, y=HR)) + 
    ggtitle(paste0('Cox Propotional Hazard Model of ', cancer, ' in ', lab)) +
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.6, position=pd, color='#F68780',size=1.2) +
    geom_line(position=pd, size = 0.75) + xlab('Variables') + ylab('Hazard Ratio (HR)') +
    geom_hline(yintercept = 1, linetype = 'dashed', color = "blue", size = 1) +
    theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.6, size = 15),
          axis.text.y = element_text(size = 15),
          title = element_text(size = 20)) +
    geom_point(position=pd, color='#F68780', size=3)
  tiff(filename = paste0('./cluster/cox_',cancer,kmclst,'.tiff'),width = 800, height = 600)
  print(g)
  dev.off()
}
cancer.list <- c('BLCA','CESC','HNSC','KIRC','LGG','LIHC','OV','STAD','UCEC')
cancer.list <- c('BLCA','HNSC','KIRC','LIHC','STAD','UCEC')
# CESC no DEGs
Clin <- load.Clin(cancer.list[1])
km.clst(cancer.list[1])
cox.km(1, cancer.list[1], 1)
cox.km(2, cancer.list[1], 1)

Clin <- load.Clin(cancer.list[2])
km.clst(cancer.list[2])
cox.km(1, cancer.list[2], 1)
cox.km(2, cancer.list[2], 1)

Clin <- load.Clin(cancer.list[3])
km.clst(cancer.list[3])
cox.km(1, cancer.list[3], 2)
cox.km(2, cancer.list[3], 2)

Clin <- load.Clin(cancer.list[4])
km.clst(cancer.list[4])
cox.km(1, cancer.list[4], 2)
cox.km(2, cancer.list[4], 2)

# STAD and UCEC remove NA

Clin <- load.Clin(cancer.list[5])
km.clst(cancer.list[5])
cox.km(1, cancer.list[5], 2)
cox.km(2, cancer.list[5], 2)

Clin <- load.Clin(cancer.list[6])
km.clst(cancer.list[6])
cox.km(1, cancer.list[6], 1)
cox.km(2, cancer.list[6], 1)

# com in com
cancer.list <- c('BLCA','HNSC','KIRC','LIHC','STAD','UCEC')
for(c1 in 1:5) {
  for(c2 in (1+c1):6) {
    d1 <- read.csv(paste0('./deg_sum/coms_',cancer.list[c1],'.csv'),T,stringsAsFactors = F)$X
    d2 <- read.csv(paste0('./deg_sum/coms_',cancer.list[c2],'.csv'),T,stringsAsFactors = F)$X
    venn.diagram(x=list(x1=d1, x2=d2), width = 1500, height = 1500, cat.pos = c(180, 180),
                 filename = paste0('./deg_sum/comcom_',cancer.list[c1],'_',
                                   cancer.list[c2],'.tiff'), imagetype = 'tiff',
                 resolution = 300, main = paste0('Common DEGs of ', cancer.list[c1], 
                                                 ' and ', cancer.list[c2]),
                 fill = c('yellow', 'purple'), col = 'transparent',
                 category.names = c(cancer.list[c1], cancer.list[c2]))
    write.csv(intersect(d1,d2), paste0('./deg_sum/comcom_',cancer.list[c1],'_',
                                       cancer.list[c2],'.csv'), quote = F)
  }
}

