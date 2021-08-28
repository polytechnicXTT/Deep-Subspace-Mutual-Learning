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
survfile <- paste(parentpath, cancername, "_Survival.txt", sep="")

file_dir <- "D://enrichment_test//LSCC//NEMO//partial_clustering//"
graph_dir <- paste0(file_dir, "final_graph//final_graph_0.1.txt")

save_dir_1 = paste0(file_dir, "Labels//1//label_0.1.txt")
save_dir_2 = paste0(file_dir, "surv_pval//1//surv_pval_0.1.txt")

sim_graph <- read.table(graph_dir, header = T)
sim_graph <- as.matrix(sim_graph)

surv <- read.table(survfile, header = T)

Labels = spectralClustering(sim_graph, classi)
print(Labels)
write.table(Labels, save_dir_1, row.names = F, col.names = F)


survresult<-survdiff(Surv(Survival, Death)~Labels, data=surv)
p.val <- 1 - pchisq(survresult$chisq, length(survresult$n) - 1)
print(p.val)
write.table(p.val, save_dir_2, row.names = F, col.names = F)