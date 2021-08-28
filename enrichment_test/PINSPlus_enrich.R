# install.packages("xlsx", repos="https://mran.microsoft.com/snapshot/2019-02-01/")
library("xlsx")


cal.age.enrichment <- function(age, cluster){
  # 璁＄畻骞撮緞鐨凨W妫€楠屽€?
  # cluster <- read.xlsx("D:\\enrichment_test\\lung\\lung_clustering.xlsx", sheetIndex=1)
  # age <- read.xlsx("D:\\enrichment_test\\lung\\lung_age.xlsx", sheetIndex=1)
  age.test.result = kruskal.test(as.numeric(unlist(age)), as.numeric(unlist(cluster)))
  
  return(age.test.result)
}

cal.discrete.enrichment <- function(file_dir){
  # 璁＄畻绂绘暎鍨嬪弬鏁扮殑瀵岄泦鍒嗘瀽
  para_tb = table(as.data.frame(file_dir[,2:3]))
  print(para_tb)
  enrichment.result <- chisq.test(para_tb)
  return(enrichment.result)
}

file_dir = "D:\\enrichment_test\\COAD\\PINSPlus\\"
param_dir = paste0(file_dir, "PINSPlus_param_comb\\")
enrich_dir = paste0(file_dir, "PINSPlus_enrich_surv\\")

# age enrichment
age = "D:\\enrichment_test\\COAD\\age.xls"
age_file <- read.xlsx(age, sheetIndex = 1)
# print(age_file)
cluster_name = paste0(param_dir, "label_1.xls")
print(cluster_name)
cluster_file <- read.xlsx(cluster_name, sheetIndex = 1)
age_test_res = cal.age.enrichment(age_file, cluster_file)
print(age_test_res)
p_val = as.matrix(age_test_res)
p_value = as.data.frame(p_val[3,])
age_dir = paste0(enrich_dir,"age_res.txt")
write.table(p_value,age_dir,append = T, row.names = F, col.names = F)

print("_______________________________________________________________________________")


# gender enrichment analysis
gender_name = paste0(param_dir, "gender.xls")
print(gender_name)

gender_file = read.xlsx(gender_name, sheetIndex = 1)
gender.pval = cal.discrete.enrichment(gender_file)
print(gender.pval)
p_val = as.matrix(gender.pval)
p_value = as.data.frame(p_val[3,])

gender_dir = paste0(enrich_dir,"gender_res.txt")
write.table(p_value,gender_dir,append = T, row.names = F, col.names = F)

print("_______________________________________________________________________________")

# path_M enrichment analysis
path_M_name = paste0(param_dir, "pathologic_M.xls")
print(path_M_name)

pathM_file = read.xlsx(path_M_name, sheetIndex = 1)
pathM.pval = cal.discrete.enrichment(pathM_file)
print(pathM.pval)
p_val = as.matrix(pathM.pval)
p_value = as.data.frame(p_val[3,])

pathM_dir = paste0(enrich_dir,"pathM_res.txt")
write.table(p_value,pathM_dir,append = T, row.names = F, col.names = F)
print("_______________________________________________________________________________")


# path_N enrichment analysis
path_N_name = paste0(param_dir, "pathologic_N.xls")
print(path_N_name)

pathN_file = read.xlsx(path_N_name, sheetIndex = 1)
pathN.pval = cal.discrete.enrichment(pathN_file)
print(pathN.pval)
p_val = as.matrix(pathN.pval)
p_value = as.data.frame(p_val[3,])

pathN_dir = paste0(enrich_dir,"pathN_res.txt")
write.table(p_value, pathN_dir, append = T, row.names = F, col.names = F)

print("_______________________________________________________________________________")


# path_T enrichment analysis
path_T_name = paste0(param_dir, "pathologic_T.xls")
print(path_T_name)

pathT_file = read.xlsx(path_T_name, sheetIndex = 1)
pathT.pval = cal.discrete.enrichment(pathT_file)
print(pathT.pval)
p_val = as.matrix(pathT.pval)
p_value = as.data.frame(p_val[3,])

pathT_dir = paste0(enrich_dir,"pathT_res.txt")
write.table(p_value,pathT_dir,append = T, row.names = F, col.names = F)

print("_______________________________________________________________________________")

# path_stage enrichment analysis
path_stage_name = paste0(param_dir, "pathologic_stage.xls")
print(path_stage_name)

path_stage_file = read.xlsx(path_stage_name, sheetIndex = 1)
path_stage.pval = cal.discrete.enrichment(path_stage_file)
print(path_stage.pval)
p_val = as.matrix(path_stage.pval)
p_value = as.data.frame(p_val[3,])

path_stage_dir = paste0(enrich_dir,"path_stage_res.txt")

write.table(p_value,path_stage_dir,append = T, row.names = F, col.names = F)
