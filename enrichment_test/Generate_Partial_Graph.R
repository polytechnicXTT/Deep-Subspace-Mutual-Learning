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
cancername <- "GLIO"
classi <- 3


parentpath <- "H:/PyCharm/PycharmProjects/MvSCN/DSFC201909/data/"
genefile <- paste(parentpath, cancername, "_Gene_Expression.txt", sep="")
methyfile <- paste(parentpath, cancername, "_Methy_Expression.txt", sep="")
mirnafile <- paste(parentpath, cancername, "_Mirna_Expression.txt", sep="")
survfile <- paste(parentpath, cancername, "_Survival.txt", sep="")

file_dir <- "D://enrichment_test//GBM//NEMO//partial_clustering\\"
# gene_dir <- paste0(file_dir, "gene_graph.txt")
# methy_dir <- paste0(file_dir, "methy_graph.txt")
# mirna_dir <- paste0(file_dir, "mirna_graph.txt")
mean_2_dir <- paste0(file_dir, "mean_2_graph.txt")
mean_3_dir <- paste0(file_dir, "mean_3_graph.txt")


gene <- read.table(genefile, header = T)
methy <- read.table(methyfile, header = T)
mirna <- read.table(mirnafile, header = T)

names(methy) <- names(gene)
names(mirna) <- names(gene)

num_sam <- ncol(gene)
print(num_sam)
num_neighbors = num_sam%/%classi
print(num_neighbors)


gene <- list(gene)
gene_graph <- nemo.affinity.graph(gene, k=num_neighbors)
# write.table(gene_graph, gene_dir, row.names = T, col.names = T)

methy <- list(methy)
methy_graph <- nemo.affinity.graph(methy, k=num_neighbors)
# write.table(methy_graph, methy_dir, row.names = T, col.names = T)

mirna <- list(mirna)
mirna_graph <- nemo.affinity.graph(mirna, k=num_neighbors)
# write.table(mirna_graph, mirna_dir, row.names = T, col.names = T)

mean_2_graph <- (methy_graph + mirna_graph) / 2
write.table(mean_2_graph, mean_2_dir, row.names = T, col.names = T)

mean_3_graph <- (gene_graph + methy_graph + mirna_graph) / 3
write.table(mean_3_graph, mean_3_dir, row.names = T, col.names = T)




