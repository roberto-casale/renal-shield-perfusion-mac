# holdout_permutation_unpaired.R -- external validation, step 3 (Section 4.6, notebook 3).
#
# For one comparison (HOLDOUT_TAG = "GSE267650" or "GSE267650direct"), limma is fitted on
# log2 CPM of GSE267650. Two statistics are computed for the signature genes:
#   conc    fraction of genes whose log2 fold change has the sign of their direction;
#   mean_t  mean moderated t, with the sign reversed for genes assigned "down".
# The group labels are reassigned at random B times (default 1,000, seed 20260828) and
# p = (1 + r) / (B + 1), r being the number of reassignments with a statistic at least
# as large as the observed one. Writes results/validation_holdout_perm_<tag>.csv
# (row perm = 0 holds the observed values).

.libPaths(c(file.path(getwd(), "env", "rlib"), .libPaths()))
Sys.setenv(RENAL_IRI_ROOT = getwd())
suppressWarnings(suppressMessages({
  source("R/helpers.R"); source("R/pipeline_meta.R"); source("R/pipeline_de.R"); library(limma)
}))

a   <- commandArgs(trailingOnly = TRUE)
tag <- if (exists("HOLDOUT_TAG")) HOLDOUT_TAG else a[1]
B   <- as.integer(if (exists("HOLDOUT_B")) HOLDOUT_B else if (is.na(a[2])) 1000 else a[2])

raw <- read.delim(gzfile("data/holdout/GSE267650_Counts_clean.txt.gz"), check.names = FALSE,
                  stringsAsFactors = FALSE)
cnt <- as.matrix(raw[, -(1:3)]); rownames(cnt) <- raw$symbol
g <- sub("_[0-9]+$", "", colnames(cnt))
cnt <- cnt[, !grepl("^Azin1|^AZIN1", g)]; g <- g[!grepl("^Azin1|^AZIN1", g)]
ord <- order(rownames(cnt), -rowSums(cnt))
cnt <- cnt[ord, ][!duplicated(rownames(cnt)[ord]), ]

CTL <- "WT_Ctl-No-IR"
case <- grep("^WT_IRI_", unique(g), value = TRUE)
REF <- if (tag == "GSE267650direct") "WT_Ischemia_20min_without_reperfusion" else CTL
s <- c(colnames(cnt)[g %in% REF], colnames(cnt)[g %in% case])
y <- cnt[, s]
y <- y[rowSums(y >= 10) >= max(2, floor(0.2 * ncol(y))), ]
lcpm <- log2(t(t(y) / colSums(y)) * 1e6 + 1)
n0 <- sum(g %in% REF); n1 <- length(s) - n0
lab0 <- rep(c(0L, 1L), c(n0, n1))
log_msg(sprintf("%s: %d vs %d, %d genes, B = %d sampled labellings", tag, n0, n1, nrow(lcpm), B))

sig <- read.csv("results/Table2_signature.csv", stringsAsFactors = FALSE)
probe <- data.frame(gene = rownames(lcpm), logFC = 0, AveExpr = rowMeans(lcpm),
                    P.Value = 1, adj.P.Val = 1, organism = "Mus musculus", stringsAsFactors = FALSE)
hp <- map_to_human(probe, "Mus musculus")
hp <- hp[hp$human %in% sig$human, ]
hp <- hp[order(hp$human, -hp$AveExpr), ]; hp <- hp[!duplicated(hp$human), ]
row_of <- match(hp$gene, rownames(lcpm)); dir_of <- sig$direction[match(hp$human, sig$human)]
orient <- ifelse(dir_of == "up", 1, -1)
log_msg(sprintf("signature genes traceable to a row: %d of %d", nrow(hp), nrow(sig)))

stats_of <- function(lab) {
  d <- model.matrix(~ factor(lab))

  tt <- topTable(eBayes(lmFit(lcpm, d), trend = TRUE), coef = 2, number = Inf, sort.by = "none")
  l <- tt$logFC[row_of]; inn <- orient * tt$t[row_of]
  c(conc = mean((dir_of == "up" & l > 0) | (dir_of == "down" & l < 0), na.rm = TRUE),
    mean_t = mean(inn, na.rm = TRUE))
}
obs <- stats_of(lab0)

total <- choose(n0 + n1, n0)
exact <- total <= 5000
if (exact) {
  log_msg(sprintf("permutation group = %d labellings -> ENUMERATING ALL", total))
  idx <- combn(n0 + n1, n0)
  null <- matrix(NA_real_, ncol(idx), 2, dimnames = list(NULL, names(obs)))
  for (i in seq_len(ncol(idx))) {
    l <- rep(1L, n0 + n1); l[idx[, i]] <- 0L
    null[i, ] <- stats_of(l)
  }
} else {
  set.seed(20260828)
  log_msg(sprintf("permutation group = %.3g labellings -> SAMPLING B = %d", total, B))
  null <- matrix(NA_real_, B, 2, dimnames = list(NULL, names(obs)))
  t0 <- Sys.time()
  for (i in seq_len(B)) {
    null[i, ] <- stats_of(sample(lab0))
    if (i %% 200 == 0) cat(sprintf("   %d / %d (%.0f s)
", i, B,
                                   as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
}

cat(sprintf("\n==== label permutation (%s, B = %d) ====\n", tag, B))
for (k in colnames(null)) {
  o <- obs[[k]]; ge <- sum(null[, k] >= o - 1e-12)

  p       <- if (exact) ge / nrow(null) else (1 + ge) / (nrow(null) + 1)
  floor_p <- if (exact) 1 / nrow(null) else 1 / (nrow(null) + 1)
  cat(sprintf("  %-7s observed %8.3f | null median %7.3f, max %7.3f | >= observed: %4d/%d | p = %.4f%s
",
              k, o, median(null[, k]), max(null[, k]), ge, nrow(null), p,
              if (abs(p - floor_p) < 1e-12) "  (the smallest attainable p)" else ""))
}

write.csv(rbind(data.frame(perm = 0L, t(obs)), data.frame(perm = seq_len(nrow(null)), null)),
          sprintf("results/validation_holdout_perm_%s.csv", tag), row.names = FALSE)
