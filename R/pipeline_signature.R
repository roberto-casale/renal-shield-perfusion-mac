# pipeline_signature.R -- the signature rule (Section 4.4, notebook 2).
#
# Only genes measured in at least 5 of the 9 datasets are tested. Benjamini-Hochberg
# correction is applied to the aggregation score of each gene in its assigned
# direction, separately for the genes assigned "up" and for those assigned "down". A
# gene enters the signature when its adjusted p-value is below 0.05 and it changes in
# that direction in at least half of the datasets that measure it (at least 2).

suppressWarnings(suppressMessages({
  source(file.path(local({.r <- Sys.getenv("RENAL_IRI_ROOT", ""); if (nzchar(.r) && dir.exists(.r)) .r else { .d <- normalizePath(getwd(), winslash = "/"); while (!file.exists(file.path(.d, "run_all.ps1")) && dirname(.d) != .d) .d <- dirname(.d); .d }}), "R", "helpers.R"))
  library(dplyr)
}))

# Thresholds of the signature.
SIG_RULE <- list(
  MIN_DATASETS = 5,
  FDR_MAX      = 0.05
)

posthoc_signature <- function(cons, min_datasets = SIG_RULE$MIN_DATASETS,
                              fdr_max = SIG_RULE$FDR_MAX) {
  keep <- cons[cons$n_datasets >= min_datasets, ]
  up   <- keep$direction == "up"
  keep$FDR_posthoc <- NA_real_
  keep$FDR_posthoc[up]  <- p.adjust(keep$RRA_up[up],   method = "BH")
  keep$FDR_posthoc[!up] <- p.adjust(keep$RRA_down[!up], method = "BH")
  keep$signature <- keep$FDR_posthoc < fdr_max &
                    keep$concordant >= pmax(2L, ceiling(keep$n_datasets / 2))
  log_msg(sprintf("Post-hoc signature: %d genes (%d tested at coverage >= %d)",
                  sum(keep$signature), nrow(keep), min_datasets))
  keep
}
