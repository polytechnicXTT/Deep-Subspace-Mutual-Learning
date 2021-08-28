source("H:\\PyCharm\\PycharmProjects\\NEMO_master\\NEMO\\R\\NEMO.R")
library(CancerSubtypes)
library(survival)
library("SNFtool")

# load the data
# --------------------------------------------
# change the cacer name
# 1 "BREAST"
# 2 "COLON"
# 3 "GLIO"
# 4 "KIDNEY"
# 5 "LUNG"
# 6 "OVARY"
cancername <- "LUNG"
classi <- 4


parentpath <- "H:/PyCharm/PycharmProjects/MvSCN/DSFC201909/data/"
genefile <- paste(parentpath, cancername, "_Gene_Expression.txt", sep="")
methyfile <- paste(parentpath, cancername, "_Methy_Expression.txt", sep="")
mirnafile <- paste(parentpath, cancername, "_Mirna_Expression.txt", sep="")
survfile <- paste(parentpath, cancername, "_Survival.txt", sep="")

num = 1
file_dir <- "D://enrichment_test//LSCC//NEMO//"
save_dir_1 = paste0(file_dir, "NEMO_param_comb//Labels//label",num, ".txt")
save_dir_2 = paste0(file_dir, "surv_pval",num,".txt")

gene <- read.table(genefile, header = T)
methy <- read.table(methyfile, header = T)
mirna <- read.table(mirnafile, header = T)
surv <- read.table(survfile, header = T)

names(methy) <- names(gene)
names(mirna) <- names(gene)

gene = as.matrix(gene)
num_sam = ncol(gene)
print(num_sam)
methy = as.matrix(methy)
mirna = as.matrix(mirna)

data_list <- list(gene, methy, mirna)

num_neighbors = num_sam%/%classi
print(num_neighbors)

Labels = nemo.clustering(data_list, classi, num_neighbors)
print(Labels)
write.table(Labels, save_dir_1, row.names = F, col.names = F)


survresult<-survdiff(Surv(Survival, Death)~Labels, data=surv)
p.val <- 1 - pchisq(survresult$chisq, length(survresult$n) - 1)
print(p.val)
write.table(p.val, save_dir_2, row.names = F, col.names = F)