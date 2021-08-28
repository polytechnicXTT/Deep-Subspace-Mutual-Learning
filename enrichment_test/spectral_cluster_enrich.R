# install.packages("xlsx", repos="https://mran.microsoft.com/snapshot/2019-02-01/")
library("xlsx")


cal.age.enrichment <- function(age, cluster){
  
  # cluster <- read.xlsx("D:\\enrichment_test\\lung\\lung_clustering.xlsx", sheetIndex=1)
  # age <- read.xlsx("D:\\enrichment_test\\lung\\lung_age.xlsx", sheetIndex=1)
  age.test.result = kruskal.test(as.numeric(unlist(age)), as.numeric(unlist(cluster)))
  
  return(age.test.result)
}

cal.discrete.enrichment <- function(file_dir){
  
  para_tb = table(as.data.frame(file_dir[,2:3]))
  # print(para_tb)
  enrichment.result <- chisq.test(para_tb)
  return(enrichment.result)
}

model = "model"
#missing_rate = "0.5"
missing_rate = "view_3"
# file_dir = paste0("D:\\PyCharm\\PyCharm_Projects\\DSML\\cluster_result\\GBM\\", model, "\\")
file_dir = paste0("D:\\PyCharm\\PyCharm_Projects\\DSML\\cluster_result\\GBM\\spectral_cluster_res\\", missing_rate, "\\", model, "\\")
param_dir = paste0(file_dir, model, "_param_comb\\")
enrich_dir = paste0(file_dir, model, "_enrich_surv\\")


age = "D:\\PyCharm\\PyCharm_Projects\\DSML\\data\\GBM\\GLIO_age.xlsx"
age_file <- read.xlsx(age, sheetIndex = 1)
# print(age_file)
# cluster_file_dir = "D:\\enrichment_test\\LSCC\\4\\cluster\\cluster_test"
# gender_file_dir = "D:\\enrichment_test\\LSCC\\4\\gender\\gender_test"
# path_M_file_dir = "D:\\enrichment_test\\LSCC\\4\\pathologic_M\\pathologic_M_test"
# path_N_file_dir = "D:\\enrichment_test\\LSCC\\4\\pathologic_N\\pathologic_N_test"
# path_T_file_dir = "D:\\enrichment_test\\LSCC\\4\\pathologic_T\\pathologic_T_test"
# path_stage_file_dir = "D:\\enrichment_test\\LSCC\\4\\pathologic_stage\\pathologic_stage_test"

for (i in 1:126){
  cluster_file_dir = paste0(param_dir, i, "\\cluster\\cluster_test")
  gender_file_dir = paste0(param_dir, i, "\\gender\\gender_test")
  path_M_file_dir = paste0(param_dir, i, "\\pathologic_M\\pathologic_M_test")
  path_N_file_dir = paste0(param_dir, i, "\\pathologic_N\\pathologic_N_test")
  path_T_file_dir = paste0(param_dir, i,"\\pathologic_T\\pathologic_T_test")
  path_stage_file_dir = paste0(param_dir, i, "\\pathologic_stage\\pathologic_stage_test")

  cluster_name = paste0(cluster_file_dir, ".xls")
  print(cluster_name)
  cluster_file <- read.xlsx(cluster_name, sheetIndex = 1)
  age_test_res = cal.age.enrichment(age_file, cluster_file)
  # print(age_test_res)
  p_val = as.matrix(age_test_res)
  p_value = as.data.frame(p_val[3,])
    
  age_dir = paste0(enrich_dir,"age_res_", i, ".txt")
    
  # write.table(p_value,"D:\\enrichment_test\\LSCC\\4\\age_res.txt",append = T, row.names = F, col.names = F)
    
  write.table(p_value,age_dir,append = T, row.names = F, col.names = F)
  print("_______________________________________________________________________________")
  
  

  gender_name = paste0(gender_file_dir, ".xls")
  print(gender_name)
  gender_file = read.xlsx(gender_name, sheetIndex = 1)
  gender.pval = cal.discrete.enrichment(gender_file)
  # print(gender.pval)
  p_val = as.matrix(gender.pval)
  p_value = as.data.frame(p_val[3,])
    
  gender_dir = paste0(enrich_dir,"gender_res_", i, ".txt")
    
  # write.table(p_value,"D:\\enrichment_test\\LSCC\\4\\gender\\gender_res.xls",append = T, row.names = F, col.names = F)
    
  write.table(p_value,gender_dir,append = T, row.names = F, col.names = F)
  
  print("_______________________________________________________________________________")
  
  
  # path_M_name = paste0(path_M_file_dir, ".xls")
  # print(path_M_name)
  # pathM_file = read.xlsx(path_M_name, sheetIndex = 1)
  # pathM.pval = cal.discrete.enrichment(pathM_file)
  # # print(pathM.pval)
  # p_val = as.matrix(pathM.pval)
  # p_value = as.data.frame(p_val[3,])
  #   
  # pathM_dir = paste0(enrich_dir,"pathM_res_", i, ".txt")
  #   
  # # write.table(p_value,"D:\\enrichment_test\\LSCC\\4\\pathologic_M\\pathM_res.xls",append = T, row.names = F, col.names = F)
  #   
  # write.table(p_value,pathM_dir,append = T, row.names = F, col.names = F)
  # print("_______________________________________________________________________________")
  # 
  # 
  # path_N_name = paste0(path_N_file_dir, ".xls")
  # print(path_N_name)
  # pathN_file = read.xlsx(path_N_name, sheetIndex = 1)
  # pathN.pval = cal.discrete.enrichment(pathN_file)
  # # print(pathN.pval)
  # p_val = as.matrix(pathN.pval)
  # p_value = as.data.frame(p_val[3,])
  #   
  # pathN_dir = paste0(enrich_dir,"pathN_res_", i, ".txt")
  #   
  # # write.table(p_value,"D:\\enrichment_test\\LSCC\\4\\pathologic_N\\pathN_res.xls",append = T, row.names = F, col.names = F)
  #   
  # write.table(p_value, pathN_dir, append = T, row.names = F, col.names = F)
  # 
  # print("_______________________________________________________________________________")
  # 
  # 
  # 
  # 
  # path_T_name = paste0(path_T_file_dir, ".xls")
  # print(path_T_name)
  # pathT_file = read.xlsx(path_T_name, sheetIndex = 1)
  # pathT.pval = cal.discrete.enrichment(pathT_file)
  # # print(pathT.pval)
  # p_val = as.matrix(pathT.pval)
  # p_value = as.data.frame(p_val[3,])
  #   
  # pathT_dir = paste0(enrich_dir,"pathT_res_", i, ".txt")
  #   
  # # write.table(p_value,"D:\\enrichment_test\\LSCC\\4\\pathologic_T\\pathT_res.xls",append = T, row.names = F, col.names = F)
  #   
  # write.table(p_value,pathT_dir,append = T, row.names = F, col.names = F)
  # 
  # print("_______________________________________________________________________________")
  # 
  # 
  # 
  # path_stage_name = paste0(path_stage_file_dir, ".xls")
  # print(path_stage_name)
  # path_stage_file = read.xlsx(path_stage_name, sheetIndex = 1)
  # path_stage.pval = cal.discrete.enrichment(path_stage_file)
  # # print(path_stage.pval)
  # p_val = as.matrix(path_stage.pval)
  # p_value = as.data.frame(p_val[3,])
  #   
  # path_stage_dir = paste0(enrich_dir,"path_stage_res_", i, ".txt")
  #   
  # # write.table(p_value,"D:\\enrichment_test\\LSCC\\4\\pathologic_stage\\path_stage_res.xls",append = T, row.names = F, col.names = F)
  #   
  # write.table(p_value,path_stage_dir,append = T, row.names = F, col.names = F)
  # 
  
}


# cluster <- read.xlsx("D:\\enrichment_test\\lung\\cluster\\cluster_test0.xls", sheetIndex=1)
# age <- read.xlsx("D:\\enrichment_test\\lung\\age.xls", sheetIndex=1)
# age.test.result = cal.age.enrichment(age,cluster)
# print(age.test.result)
# gender <- read.xlsx("D:\\enrichment_test\\lung\\lung_gender.xlsx", sheetIndex = 1)
# gender.pval = cal.discrete.enrichment(gender)
# print(gender.pval)

