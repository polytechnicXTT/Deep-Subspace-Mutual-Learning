Algorithm_name <- c("SNF", "CC", "SNF.CC", "PINS", "LRAcluster", "NEMO", "DSML")

# 1. DSML-our method
# LSCC BIC KRCCC COAD GBM
# dsml_surv <- c(-log10(1.45e-4), -log10(1.79e-4), -log10(4.21e-4), -log10(1.03e-4), -log10(1.36e-2))
# mean_dsml_surv <- mean(dsml_surv)
# dsml_enrich <- c(3, 3, 4, 3, 2)
# mean_dsml_enrich <- mean(dsml_enrich)
mean_dsml_surv <- 3.5
mean_dsml_enrich <- 3.2
print(mean_dsml_surv)
print(mean_dsml_enrich)

# 2. SNF
snf_surv <- c(-log10(7.44e-3), -log10(1.87e-5), -log10(2.15e-1), -log10(3.74e-2), -log10(5.75e-4))
mean_snf_surv <- mean(snf_surv)
snf_enrich <- c(2, 1, 3, 2, 1)
mean_snf_enrich <- mean(snf_enrich)
print(mean_snf_surv)
print(mean_snf_enrich)

# 3.CC
cc_surv <- c(-log10(2.39e-3), -log10(3.57e-2), -log10(9.31e-4), -log10(9.30e-1), -log10(7.44e-1))
mean_cc_surv <- mean(cc_surv)
cc_enrich <- c(1, 3, 2, 1, 1)
mean_cc_enrich <- mean(cc_enrich)
print(mean_cc_surv)
print(mean_cc_enrich)

# 4.SNF.CC
snfcc_surv <- c(-log10(3.40e-3), -log10(3.93e-5), -log10(1.87e-1), -log10(3.68e-2), -log10(4.29e-4))
mean_snfcc_surv <- mean(snfcc_surv)
snfcc_enrich <- c(2, 1, 3, 2, 1)
mean_snfcc_enrich <- mean(snfcc_enrich)
print(mean_snfcc_surv)
print(mean_snfcc_enrich)

# 5.PINS
pins_surv <- c(-log10(2.00e-3), -log10(2.78e-1), -log10(9.60e-1), -log10(1.34e-1), -log10(3.32e-2))
mean_pins_surv <- mean(pins_surv)
pins_enrich <- c(1, 1, 1, 1, 1)
mean_pins_enrich <- mean(pins_enrich)
print(mean_pins_surv)
print(mean_pins_enrich)

# 6.LRAcluster
lracluster_surv <- c(-log10(1.98e-2), -log10(1.66e-1), -log10(5.92e-1), -log10(7.12e-4), -log10(2.03e-2))
mean_lracluster_surv <- mean(lracluster_surv)
lracluster_enrich <- c(2, 1, 2, 1, 1)
mean_lracluster_enrich <- mean(lracluster_enrich)
print(mean_lracluster_surv)
print(mean_lracluster_enrich)

# 7.NEMO
nemo_surv <- c(-log10(1.25e-2), -log10(2.12e-3), -log10(2.43e-2), -log10(1.99e-2), -log10(2.21e-4))
mean_nemo_surv <- mean(nemo_surv)
nemo_enrich <- c(1, 1, 3, 1, 1)
mean_nemo_enrich <- mean(nemo_enrich)
print(mean_nemo_surv)
print(mean_nemo_enrich)

alg_cols <- c("cyan", "magenta", "blue", "green", "indianred4", "orange", "red")
pch = c(16, 17, 18, 2, 3, 8, 15)
# alg_cols <- c(rep("black", 7))
# pch = c(rep(19, 7))


available.indices = !is.na(surv_pval)
surv_pval = surv_pval[available.indices]
num_enrich = num_enrich[available.indices]
current.cols = alg_cols[available.indices]
current.pch = pch[available.indices]

xlab = 'Mean number of enriched clinical parameters'
ylab = 'Mean -log10 of logrank p-value' 

num_enrich <- c(mean_snf_enrich, mean_cc_enrich, mean_snfcc_enrich, mean_pins_enrich, mean_lracluster_enrich, mean_nemo_enrich, mean_dsml_enrich)
surv_pval <- c(mean_snf_surv, mean_cc_surv, mean_snfcc_surv, mean_pins_surv, mean_lracluster_surv, mean_nemo_surv, mean_dsml_surv)


plot(num_enrich, surv_pval, xlab=xlab, ylab=ylab, xlim=c(0.0, max(num_enrich+1.0)), 
    ylim=c(0.0, max(surv_pval, surv_significance)+0.2), col = alg_cols, pch=pch, cex.lab = 1.5, cex=2, cex.axis=1.2, cex.main=1.6, lwd=3)

abline(h = seq(0, 5.8, by = 0.2),v = seq(0, 5, by=0.2), col = "lightgray", lty = 3)
abline(h = c(mean_dsml_surv), v = c(mean_dsml_enrich), col="red", lty = 3)

# grid(20, 20, col="lightgray")
# text(mean_dsml_enrich-0.4, mean_dsml_surv, "DSML")
# text(mean_snf_enrich+0.4, mean_snf_surv, "SNF")
# text(num_enrich, surv_pval, Algorithm_name)
legend("bottomright", Algorithm_name, pch=pch, col=alg_cols, cex=0.9)


