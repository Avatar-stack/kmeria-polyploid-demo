#!/usr/bin/env Rscript
# GWASpoly (additive model, kinship from markers) on GATK ploidy=4 dosages.
# Usage: gwaspoly.R geno.csv pheno.tsv outprefix n_core
suppressMessages(library(GWASpoly))
a <- commandArgs(TRUE); geno <- a[1]; ph <- a[2]; out <- a[3]; nc <- as.integer(a[4])
p <- read.delim(ph, check.names = FALSE)
pf <- paste0(out, ".pheno.csv"); write.csv(p, pf, row.names = FALSE, quote = FALSE)
d <- read.GWASpoly(ploidy = 4, pheno.file = pf, geno.file = geno, format = "numeric",
                   n.traits = ncol(p) - 1, delim = ",")
d <- set.K(d, LOCO = FALSE, n.core = nc)
d <- GWASpoly(d, models = "additive", traits = colnames(p)[-1], n.core = nc)
for (t in colnames(p)[-1]) {
  s <- d@scores[[t]]
  write.table(data.frame(marker = rownames(s), p = 10^(-s[, "additive"])),
              paste0(out, ".", t, ".p"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}
