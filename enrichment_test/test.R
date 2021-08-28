# install.packages("xlsx", repos="https://mran.microsoft.com/snapshot/2019-02-01/")
library("xlsx")


cal.age.enrichment <- function(age, cluster){
  # 计算数值型参数的KW检验
  # cluster <- read.xlsx("D:\\enrichment_test\\lung\\lung_clustering.xlsx", sheetIndex=1)
  # age <- read.xlsx("D:\\enrichment_test\\lung\\lung_age.xlsx", sheetIndex=1)
  age.test.result = kruskal.test(as.numeric(unlist(age)), as.numeric(unlist(cluster)))
  
  return(age.test.result)
}

cal.discrete.enrichment <- function(file_dir){
  # 计算离散型参数的卡方检验
  para_tb = table(as.data.frame(file_dir[,2:3]))
  enrichment.result <- chisq.test(para_tb)
  return(enrichment.result)
}


# age = "D:\\enrichment_test\\lung\\age.xls"
# age_file <- read.xlsx(age, sheetIndex = 1)
# # print(age_file)
# cluster_file_dir = "D:\\enrichment_test\\lung\\cluster\\cluster_test"
# gender_file_dir = "D:\\enrichment_test\\lung\\gender\\gender_test"
# path_M_file_dir = "D:\\enrichment_test\\lung\\pathologic_M\\pathologic_M_test"
# path_N_file_dir = "D:\\enrichment_test\\lung\\pathologic_N\\pathologic_N_test"
# path_T_file_dir = "D:\\enrichment_test\\lung\\pathologic_T\\pathologic_T_test"
# path_stage_file_dir = "D:\\enrichment_test\\lung\\pathologic_stage\\pathologic_stage_test"
# 
# epoch = 20
# 
# for (i in 1:epoch-1){
#   cluster_name = paste0(cluster_file_dir, i, ".xls")
#   print(cluster_name)
#   cluster_file <- read.xlsx(cluster_name, sheetIndex = 1)
#   age_test_res = cal.age.enrichment(age_file, cluster_file)
#   print(age_test_res)
#   p_val = as.matrix(age_test_res)
#   # p_value = as.matrix(p_val[3,])
#   p_value = as.data.frame(p_val[3,])
#   
#   write.table(p_value,"D:\\enrichment_test\\lung\\result.xls",append = T, row.names = F, col.names = F)
# }
# 
# print("_______________________________________________________________________________")

# cluster <- read.xlsx("D:\\enrichment_test\\lung\\cluster\\cluster_test0.xls", sheetIndex=1)
# age <- read.xlsx("D:\\enrichment_test\\lung\\age.xls", sheetIndex=1)
# age.test.result = cal.age.enrichment(age,cluster)
# print(age.test.result)
# p_val = as.matrix(age.test.result)
# p_value = as.matrix(p_val[3,])
# 
# write.table(p_val[3,],"D:\\enrichment_test\\lung\\result.xls",append = T, sep = "\n", col.names = F)

path_M <- read.xlsx("D:\\enrichment_test\\LSCC\\1\\pathologic_M\\pathologic_M_test0.xls", sheetIndex = 1)
path_M.pval = cal.discrete.enrichment(path_M)
print(path_M.pval)

