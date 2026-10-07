# holdout_replication.R -- external validation, step 2 (Section 4.6, notebook 3).
#
# For one comparison (HOLDOUT_ACC = "GSE267650" or "GSE267650direct"), the direction of
# each signature gene, assigned by the corpus, is compared with the sign of its log2 fold
# change in GSE267650 (DESeq2, step 1). Prints the concordance and the median oriented
# log2 fold change (log2 fold change multiplied by -1 for genes assigned "down"), and
# writes results/validation_holdout_<acc>.csv.

.libPaths(c(file.path(getwd(), "env", "rlib"), .libPaths()))
Sys.setenv(RENAL_IRI_ROOT = getwd())
suppressWarnings(suppressMessages({
  source("R/helpers.R"); source("R/pipeline_meta.R"); source("R/pipeline_de.R")
}))

acc <- if (exists("HOLDOUT_ACC")) HOLDOUT_ACC else commandArgs(trailingOnly = TRUE)[1]
if (is.na(acc)) stop("set HOLDOUT_ACC, or: Rscript R/validation/holdout_replication.R <GSE_ACCESSION>")
stopifnot(!(acc %in% names(REGISTRY)))   # the validation series is not one of the nine

sig <- read.csv("results/Table2_signature.csv", stringsAsFactors = FALSE)
sig <- sig[, c("human", "direction", "mean_logFC", "FDR_posthoc")]
de_file <- sprintf("results/HOLDOUT_DE_%s.csv", acc)
if (!file.exists(de_file)) stop(de_file, " is missing: run R/validation/holdout_load_GSE267650.R first")
hd <- read.csv(de_file, stringsAsFactors = FALSE)
hd <- collapse_human(map_to_human(hd, unique(hd$organism)[1]))

m <- merge(sig, hd[, c("human", "logFC", "P.Value")], by = "human")
m$agrees <- (m$direction == "up" & m$logFC > 0) | (m$direction == "down" & m$logFC < 0)
oriented <- ifelse(m$direction == "up", 1, -1) * m$logFC
cat(sprintf("%s: %d of %d signature genes measured | concordant %d (%.1f%%) | median oriented log2 fold change %+.2f\n",
            acc, nrow(m), nrow(sig), sum(m$agrees), 100 * mean(m$agrees), median(oriented)))

out <- sprintf("results/validation_holdout_%s.csv", acc)
write.csv(m[order(m$FDR_posthoc), ], out, row.names = FALSE)
log_msg(sprintf("WROTE %s", out))
