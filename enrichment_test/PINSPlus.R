## ----message=FALSE------------------------------------------------------------
library(PINSPlus)
library(CancerSubtypes)
library(survival)

# load the data
# --------------------------------------------
# change the cacer name
# 1 "BREAST"
# 2 "COLON"
# 3 "GLIO"
# 4 "KIDNEY"
# 5 "LUNG"
# 6 "OVARY"

cancername <- "COLON"
classi <- 3
parentpath <- "H:/PyCharm/PycharmProjects/MvSCN/DSFC201909/data/"
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
# ExecuteCC
# names is not same, check by names(gene) == names(methy)
# change to same first
names(methy) <- names(gene)
names(mirna) <- names(gene)

gene = t(as.matrix(gene))
print(nrow(gene))
methy = t(as.matrix(methy))
mirna = t(as.matrix(mirna))


file_dir = "D://enrichment_test//COAD//PINSPlus//PINSPlus_param_comb//Labels//"
save_dir_1 = paste0(file_dir, "label.txt")
save_dir_2 = "D://enrichment_test//COAD//PINSPlus//surv_pval.txt"


data=list(gene, methy, mirna)
# print(nrow(data[[1]]))


## -----------------------------------------------------------------------------
pins.ret <- PINSPlus::SubtypingOmicsData(dataList=data, k=classi)
Labels <- pins.ret$cluster2
print(Labels)
write.table(Labels, save_dir_1, row.names = F, col.names = F)

# Labels_dir = paste0(file_dir, "label.txt")
# Labels <- read.table(Labels_dir)
# dim(Labels)
survresult<-survdiff(Surv(Survival,Death)~Labels, data=surv)
p.val <- 1 - pchisq(survresult$chisq, length(survresult$n) - 1)
print(p.val)
write.table(p.val, save_dir_2, row.names = F, col.names = F)


