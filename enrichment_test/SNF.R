# load library
library(CancerSubtypes)
library(survival)

# clear enviroment
rm(list = ls())



# load the data
# --------------------------------------------
# change the cacer name
# 1 "BREAST"
# 2 "COLON"
# 3 "GLIO"
# 4 "KIDNEY"
# 5 "LUNG"
# 6 "OVARY"
cancername <- "GLIO"
classi <- 3
parentpath <- "H://PyCharm//PyCharmProjects//MvSCN//DSFC201909//data//"
genefile <- paste(parentpath, cancername, "_Gene_Expression.txt", sep="")
methyfile <- paste(parentpath, cancername, "_Methy_Expression.txt", sep="")
mirnafile <- paste(parentpath, cancername, "_Mirna_Expression.txt", sep="")
survfile <- paste(parentpath, cancername, "_Survival.txt", sep="")
# read mRNA
gene <- read.table(genefile, header = T)
dim(gene)
# read Methylatio
methy <- read.table(methyfile, header = T)
dim(methy)
# read miRNA
mirna <- read.table(mirnafile, header = T)
dim(mirna)
# read survival
surv <- read.table(survfile, header = T)
dim(surv)

# -----------------------------------
# ExecuteSNF.CC
# names is not same, check by names(gene) == names(methy)
# change to same first
names(methy) <- names(gene)
names(mirna) <- names(gene)

file_dir = "D://enrichment_test//GBM//SNF//SNF_param_comb//Labels//"

save_dir_2 = "D://enrichment_test//GBM//SNF//surv_pval.txt"

p_value = c()
alpha = 0.3
i=1
for( j in 1:6){
  for (k in 10:30){
      print("------------------------------------------------------------------------------------------------------")
      print(k)
      print(alpha)
      save_dir_1 = paste0(file_dir, i, ".txt")
      datas=list(GeneExp=gene,MethyExp=methy,miRNAExp=mirna)
      resultSNF=ExecuteSNF(datasets=datas, clusterNum=classi, K=k, alpha = alpha, t = 20, plot = FALSE)
      
      surv$Labels<-resultSNF$group
      print(surv$Labels)
      write.table(surv$Labels, save_dir_1, row.names = F, col.names = F)
      

      survresult<-survdiff(Surv(Survival,Death)~Labels, data=surv)
      p.val <- 1 - pchisq(survresult$chisq, length(survresult$n) - 1)
      p_value[i] <- p.val
      print(p_value)
      write.table(p_value, save_dir_2, row.names = F, col.names = F)
      i=i+1
  } 
  alpha = alpha + 0.1
}

