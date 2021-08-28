library(survival)

epoch = 20

model = "model19"
missing_rate = "0.5"
# missing_rate = "view_3"
# file_dir = paste0("D:\\PyCharm\\PyCharm_Projects\\DSML\\cluster_result\\GBM\\", model, "\\")
file_dir = paste0("D:\\PyCharm\\PyCharm_Projects\\DSML\\cluster_result\\COAD\\partial_res\\", missing_rate, "\\", model, "\\")
param_dir = paste0(file_dir, model, "_param_comb\\")
enrich_dir = paste0(file_dir, model, "_enrich_surv\\")

for (j in 1:126){
  file_dir_1 = paste0(param_dir, j, "\\txt\\colon_test")
  p_value <- c()
  for (i in 0:19)
  {
    file_name = paste0(file_dir_1, i, ".txt")
    print(file_name)
    survdata<-read.table(file = file_name, header = T)
    survresult <- survdiff(Surv(survival,death)~Labels, data=survdata)
    p.val <- 1 - pchisq(survresult$chisq, length(survresult$n) - 1)
    print(survresult)
    print(p.val)
    p_value[i+1] <- p.val
    print(p_value)
    
    surv_pval = format(p_value, scientific = TRUE)
    
    write.table(surv_pval, paste0(enrich_dir,"\\surv_pval_", j, ".txt"), row.names = F, col.names = F)
    
  }
} 

