# plot.name = "LSCC.png"

cancer_type = "LSCC"
Algorithm <- c("DSML", "SNF", "CC", "SNF.CC", "PINS", "LRAcluster")

# parameter for LSCC
surv_pval <- c(-log10(1.45e-4), -log10(7.44e-3), -log10(2.39e-3), -log10(3.40e-3), -log10(2.00e-3), -log10(1.98e-2))
num_enrich <- c(3, 2, 1, 2, 1, 2)

# parameter for BIC
# surv_pval <- c(-log10(1.95e-8), -log10(1.87e-5), -log10(3.57e-2), -log10(3.93e-5), -log10(2.78e-1), -log10(1.66e-1))
# num_enrich <- c(5, 1, 3, 1, 1, 1)

# parameter for KRCCC
# surv_pval <- c(-log10(1.43e-7), -log10(2.15e-1), -log10(9.31e-4), -log10(1.87e-1), -log10(9.60e-1), -log10(5.92e-1))
# num_enrich <- c(5, 3, 2, 3, 1, 2)

# parameter for COAD
#surv_pval <- c(-log10(2.55e-8), -log10(3.74e-2), -log10(9.30e-1), -log10(3.68e-2), -log10(1.34e-1), -log10(7.12e-4))
#num_enrich <- c(5, 2, 1, 2, 1, 1)

# parameter for GBM
# surv_pval <- c(-log10(1.36e-2), -log10(5.75e-4), -log10(7.44e-1), -log10(4.29e-4), -log10(3.32e-2), -log10(2.03e-2))
# num_enrich <- c(2, 1, 1, 1, 1, 1)



# png(file.path("D:\\enrichment_test\\LSCC\\", plot.name), width=1200)
# alg_cols = c("royalblue","tomato", "black", "green"£©
alg_cols <- c("red", "cyan", "magenta", "blue", "green", "pink")
# print(alg_cols)
pch = c(15, 16, 17, 18, 2, 3)

available.indices = !is.na(surv_pval)
surv_pval = surv_pval[available.indices]
num_enrich = num_enrich[available.indices]
current.cols = alg_cols[available.indices]
current.pch = pch[available.indices]

surv_significance <- -log10(0.05)


xlab = '-log10(logrank pvalue)'
ylab = 'enriched clinical parameters'
 

plot(surv_pval, num_enrich, main=cancer_type, xlab=xlab, ylab=ylab, xlim=c(0.0, max(surv_pval, surv_significance)+0.2), 
    ylim=c(0.0, max(num_enrich+1.0)), col = alg_cols, pch=pch, cex.lab = 1.5, cex=3, cex.axis=1.2, cex.main=1.6, lwd=3)

abline(v=surv_significance, col="red")

legend("bottomright", Algorithm, pch=pch, col=alg_cols)

