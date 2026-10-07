# holdout_load_GSE267650.R -- external validation, step 1 (Section 4.6, notebook 3).
#
# Downloads the GSE267650 counts, keeps the wild-type mice (non-operated controls,
# ischemia-only, and reperfused at six timepoints) and runs DESeq2, with group as the
# only term, for the two comparisons: reperfused vs controls, and reperfused vs
# ischemia-only. Writes results/HOLDOUT_DE_GSE267650.csv,
# results/HOLDOUT_DE_GSE267650direct.csv and results/HOLDOUT_design_GSE267650.csv.

.libPaths(c(file.path(getwd(), "env", "rlib"), .libPaths()))
Sys.setenv(RENAL_IRI_ROOT = getwd())
suppressWarnings(suppressMessages({ source("R/helpers.R"); source("R/pipeline_de.R") }))

hold_dir <- DIRS$holdout
f <- file.path(hold_dir, "GSE267650_Counts_clean.txt.gz")
if (!file.exists(f)) {
  log_msg("GSE267650 supplementary file absent; downloading from GEO")
  robust_fetch("geo_supp_GSE267650_counts", fn = function()
    GEOquery::getGEOSuppFiles("GSE267650", baseDir = hold_dir, filter_regex = "Counts"))

  cand <- list.files(hold_dir, pattern = "GSE267650_Counts.*[.]gz$",
                     recursive = TRUE, full.names = TRUE)
  stopifnot("GEO returned no GSE267650 counts file" = length(cand) >= 1)
  if (!file.exists(f)) file.copy(cand[1], f)
}
stopifnot(file.exists(f))
raw <- read.delim(gzfile(f), check.names = FALSE, stringsAsFactors = FALSE)
sym <- raw$symbol
cnt <- as.matrix(raw[, -(1:3)])
rownames(cnt) <- sym
log_msg(sprintf("GSE267650: %d genes x %d samples", nrow(cnt), ncol(cnt)))

grp_of <- function(x) sub("_[0-9]+$", "", x)
g <- grp_of(colnames(cnt))
wt <- !grepl("^Azin1|^AZIN1", g)
log_msg(sprintf("dropped %d Azin1 samples (genotype confound)", sum(!wt)))
cnt <- cnt[, wt]; g <- g[wt]

CTL  <- "WT_Ctl-No-IR"
ISCH <- "WT_Ischemia_20min_without_reperfusion"
TIMES <- grep("^WT_IRI_", unique(g), value = TRUE)
TIMES <- TIMES[order(as.numeric(sub("^WT_IRI_([0-9]+)h$", "\1", TIMES)))]
log_msg(sprintf("control n=%d | ischemia-only n=%d | reperfusion arms: %s",
                sum(g == CTL), sum(g == ISCH), paste(TIMES, collapse = ", ")))
stopifnot(sum(g == CTL) >= 3, sum(g == ISCH) >= 3, length(TIMES) == 6)

ord <- order(rownames(cnt), -rowSums(cnt))
cnt <- cnt[ord, ][!duplicated(rownames(cnt)[ord]), ]
cnt <- cnt[rownames(cnt) != "" & !is.na(rownames(cnt)), ]
log_msg(sprintf("after collapsing duplicate symbols: %d genes", nrow(cnt)))

one <- function(case, label, tag) {
  s <- c(colnames(cnt)[g == CTL], colnames(cnt)[g %in% case])
  meta <- data.frame(sample = s,
                     group = rep(c("control", "IRI"), c(sum(g == CTL), sum(g %in% case))),
                     stringsAsFactors = FALSE)
  d <- run_de(list(mat = cnt[, s], meta = meta, paired = FALSE, is_counts = TRUE))
  d$accession <- tag; d$organism <- "Mus musculus"; d$modality <- "rnaseq"
  out <- sprintf("results/HOLDOUT_DE_%s.csv", tag)
  write.csv(d, out, row.names = FALSE)
  log_msg(sprintf("%-46s %2d vs %2d -> %s (%d genes, %d flagged)",
                  label, sum(g == CTL), sum(g %in% case), out, nrow(d), sum(d$sig)))
}

one_ref <- function(case, ref, label, tag) {
  s <- c(colnames(cnt)[g %in% ref], colnames(cnt)[g %in% case])
  meta <- data.frame(sample = s,
                     group = rep(c("control", "IRI"), c(sum(g %in% ref), sum(g %in% case))),
                     stringsAsFactors = FALSE)
  d <- run_de(list(mat = cnt[, s], meta = meta, paired = FALSE, is_counts = TRUE))
  d$accession <- tag; d$organism <- "Mus musculus"; d$modality <- "rnaseq"
  write.csv(d, sprintf("results/HOLDOUT_DE_%s.csv", tag), row.names = FALSE)
  log_msg(sprintf("%-46s %2d vs %2d -> %d genes, %d flagged",
                  label, sum(g %in% ref), sum(g %in% case), nrow(d), sum(d$sig)))
}

design <- data.frame(
  arm = c("control", "ischemia_only", "reperfusion", "excluded_Azin1"),
  n = c(sum(g == CTL), sum(g == ISCH), sum(g %in% TIMES), sum(!wt)),
  detail = c("WT_Ctl-No-IR, no surgery",
             "clamp 20 min, harvested with the clamp still closed",
             sprintf("clamp released, %d timepoints: %s", length(TIMES),
                     paste(sub("^WT_IRI_", "", TIMES), collapse = " ")),
             "Azin1-locked and Azin1-uneditable: a genotype contrast is not what is tested"),
  stringsAsFactors = FALSE)
write.csv(design, "results/HOLDOUT_design_GSE267650.csv", row.names = FALSE)
log_msg(sprintf("design written: %s", paste(sprintf("%s=%d", design$arm, design$n),
                                            collapse = " | ")))

one(TIMES, "reperfused vs non-operated controls", "GSE267650")
one_ref(TIMES, ISCH, "reperfused vs ischemia-only", "GSE267650direct")
