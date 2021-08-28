Algorithm_name <- c("SNF", "CC", "SNF.CC", "PINS", "LRAcluster", "NEMO", "DSML")

# 1. DSML-our method
mean_dsml_surv <- 2.4
mean_dsml_enrich <- 2.6

# 2. SNF
mean_snf_surv <- 1.5
mean_snf_enrich <- 1.4

# 3. CC
mean_cc_surv <- 1.4
mean_cc_enrich <- 1.8

# 4. SNF.CC
mean_snfcc_surv <- 2.0
mean_snfcc_enrich <- 1.6

# 5. PINS
mean_pins_surv <- 1.6
mean_pins_enrich <- 1.8

# 6. LRAcluster
mean_lracluster_surv <- 1.8
mean_lracluster_enrich <- 1.6

# 7. NEMO
mean_nemo_surv <- 1.6
mean_nemo_enrich <- 1.4

alg_cols <- c("cyan", "magenta", "blue", "green", "indianred4", "orange", "red")
pch = c(16, 17, 18, 2, 3, 8, 15)

num_enrich <- c(mean_snf_enrich, mean_cc_enrich, mean_snfcc_enrich, mean_pins_enrich, mean_lracluster_enrich, mean_nemo_enrich, mean_dsml_enrich)
surv_pval <- c(mean_snf_surv, mean_cc_surv, mean_snfcc_surv, mean_pins_surv, mean_lracluster_surv, mean_nemo_surv, mean_dsml_surv)

available.indices = !is.na(surv_pval)
surv_pval = surv_pval[available.indices]
num_enrich = num_enrich[available.indices]
current.cols = alg_cols[available.indices]
current.pch = pch[available.indices]

xlab = 'Mean number of enriched clinical parameters'
ylab = 'Mean -log10 of logrank p-value' 

plot(num_enrich, surv_pval, xlab=xlab, ylab=ylab, xlim=c(0.0, max(num_enrich+1.0)), 
    ylim=c(0.0, max(surv_pval, surv_significance)+0.2), col = alg_cols, pch=pch, cex.lab = 1.5, cex=2, cex.axis=1.2, cex.main=1.6, lwd=3)

abline(h = seq(0, 5.8, by = 0.2),v = seq(0, 5, by=0.2), col = "lightgray", lty = 3)
abline(h = c(mean_dsml_surv), v = c(mean_dsml_enrich), col="red", lty = 3)

# grid(20, 20, col="lightgray")
# text(mean_dsml_enrich-0.4, mean_dsml_surv, "DSML")
# text(mean_snf_enrich+0.4, mean_snf_surv, "SNF")
# text(num_enrich, surv_pval, Algorithm_name)
legend("bottomright", Algorithm_name, pch=pch, col=alg_cols, cex=0.9)


