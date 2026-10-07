# validation_GSE267650.R -- external validation of the signature (Section 4.6, notebook 3).
#
# GSE267650 (mouse, 20 min bilateral renal pedicle clamping): non-operated controls,
# ischemia-only animals and animals reperfused for 4-72 h; genetically modified animals
# (Azin1) are excluded. Two comparisons: reperfused vs non-operated controls, and
# reperfused vs ischemia-only.
#
# For each comparison, genes with at least 10 counts in at least 20% of the samples
# (rounded down, minimum 2) are kept, counts are converted to log2 counts per million,
# and limma is fitted with group as the only term (eBayes with trend = TRUE). When several
# mouse genes map to the same human gene, the one with the highest mean expression is
# used; this choice does not depend on the group labels. For the signature genes:
#   concordance    fraction whose log2 fold change has the sign of their direction;
#   median         median oriented log2 fold change (multiplied by -1 for genes "down");
#   mean_t         mean moderated t, with the sign reversed for genes "down".
# Significance: the group labels are reassigned at random 1,000 times (seed 20260828) and
# p = (1 + r) / (1,000 + 1), r being the number of reassignments with a mean_t at least as
# large as the observed one.
#
# Writes results/validation_summary.csv (the numbers of Section 2.4),
# results/validation_genes_<comparison>.csv (one row per signature gene),
# results/validation_permutations_<comparison>.csv (row perm = 0 holds the observed values)
# and results/validation_design.csv.

suppressWarnings(suppressMessages({
  source(file.path(local({.r <- Sys.getenv("RENAL_IRI_ROOT", ""); if (nzchar(.r) && dir.exists(.r)) .r else { .d <- normalizePath(getwd(), winslash = "/"); while (!file.exists(file.path(.d, "run_all.ps1")) && dirname(.d) != .d) .d <- dirname(.d); .d }}), "R", "helpers.R"))
  source(file.path(PROJECT, "R", "pipeline_meta.R"))
  library(limma)
}))

# Counts of GSE267650, downloaded once from the GEO supplementary files.
f <- file.path(DIRS$holdout, "GSE267650_Counts_clean.txt.gz")
if (!file.exists(f)) {
  log_msg("GSE267650 supplementary file absent; downloading from GEO")
  robust_fetch("geo_supp_GSE267650_counts", fn = function()
    GEOquery::getGEOSuppFiles("GSE267650", baseDir = DIRS$holdout, filter_regex = "Counts"))
  cand <- list.files(DIRS$holdout, pattern = "GSE267650_Counts.*[.]gz$",
                     recursive = TRUE, full.names = TRUE)
  stopifnot("GEO returned no GSE267650 counts file" = length(cand) >= 1)
  if (!file.exists(f)) file.copy(cand[1], f)
}
raw <- read.delim(gzfile(f), check.names = FALSE, stringsAsFactors = FALSE)
cnt <- as.matrix(raw[, -(1:3)]); rownames(cnt) <- raw$symbol
g <- sub("_[0-9]+$", "", colnames(cnt))
wt <- !grepl("^Azin1|^AZIN1", g)
cnt <- cnt[, wt]; g <- g[wt]
ord <- order(rownames(cnt), -rowSums(cnt))
cnt <- cnt[ord, ][!duplicated(rownames(cnt)[ord]), ]

CTL   <- "WT_Ctl-No-IR"
ISCH  <- "WT_Ischemia_20min_without_reperfusion"
TIMES <- grep("^WT_IRI_", unique(g), value = TRUE)
TIMES <- TIMES[order(as.numeric(sub("^WT_IRI_([0-9]+)h$", "\\1", TIMES)))]
stopifnot(sum(g == CTL) >= 3, sum(g == ISCH) >= 3, length(TIMES) == 6)

# Groups and numbers of animals.
design <- data.frame(
  arm = c("control", "ischemia_only", "reperfusion", "excluded_Azin1"),
  n = c(sum(g == CTL), sum(g == ISCH), sum(g %in% TIMES), sum(!wt)),
  detail = c("WT_Ctl-No-IR, no surgery",
             "clamp 20 min, harvested with the clamp still closed",
             sprintf("clamp released, %d timepoints: %s", length(TIMES),
                     paste(sub("^WT_IRI_", "", TIMES), collapse = " ")),
             "Azin1-locked and Azin1-uneditable: a genotype contrast is not what is tested"),
  stringsAsFactors = FALSE)
save_csv(design, "validation_design.csv")

sig <- read.csv("results/Table2_signature.csv", stringsAsFactors = FALSE)
B <- 1000

# One comparison: reperfused animals against the reference group `ref`.
validate <- function(ref, label) {
  s <- c(colnames(cnt)[g %in% ref], colnames(cnt)[g %in% TIMES])
  y <- cnt[, s]
  y <- y[rowSums(y >= 10) >= max(2, floor(0.2 * ncol(y))), ]
  lcpm <- log2(t(t(y) / colSums(y)) * 1e6 + 1)
  n0 <- sum(g %in% ref); n1 <- length(s) - n0
  lab0 <- rep(c(0L, 1L), c(n0, n1))

  # signature genes: human ortholog of each mouse gene, the most expressed one per human gene
  probe <- data.frame(gene = rownames(lcpm), logFC = 0, AveExpr = rowMeans(lcpm),
                      P.Value = 1, adj.P.Val = 1, organism = "Mus musculus", stringsAsFactors = FALSE)
  hp <- map_to_human(probe, "Mus musculus")
  hp <- hp[hp$human %in% sig$human, ]
  hp <- hp[order(hp$human, -hp$AveExpr), ]; hp <- hp[!duplicated(hp$human), ]
  row_of <- match(hp$gene, rownames(lcpm)); dir_of <- sig$direction[match(hp$human, sig$human)]
  orient <- ifelse(dir_of == "up", 1, -1)

  fit_of <- function(lab) topTable(eBayes(lmFit(lcpm, model.matrix(~ factor(lab))), trend = TRUE),
                                   coef = 2, number = Inf, sort.by = "none")
  stats_of <- function(lab) {
    tt <- fit_of(lab)
    l <- tt$logFC[row_of]; inn <- orient * tt$t[row_of]
    c(conc = mean((dir_of == "up" & l > 0) | (dir_of == "down" & l < 0), na.rm = TRUE),
      mean_t = mean(inn, na.rm = TRUE))
  }

  # observed values, gene by gene
  tt <- fit_of(lab0)
  genes <- data.frame(human = hp$human, mouse_gene = hp$gene, direction = dir_of,
                      logFC = tt$logFC[row_of], t = tt$t[row_of],
                      oriented_logFC = orient * tt$logFC[row_of], stringsAsFactors = FALSE)
  genes$concordant <- genes$oriented_logFC > 0
  save_csv(genes, sprintf("validation_genes_%s.csv", label))
  obs <- stats_of(lab0)

  # permutation test
  set.seed(20260828)
  null <- matrix(NA_real_, B, 2, dimnames = list(NULL, names(obs)))
  for (i in seq_len(B)) null[i, ] <- stats_of(sample(lab0))
  write.csv(rbind(data.frame(perm = 0L, t(obs)), data.frame(perm = seq_len(B), null)),
            file.path(DIRS$results, sprintf("validation_permutations_%s.csv", label)), row.names = FALSE)
  r <- sum(null[, "mean_t"] >= obs[["mean_t"]] - 1e-12)

  data.frame(comparison = label, n_reference = n0, n_reperfused = n1,
             signature_genes_measured = nrow(genes),
             concordant = sum(genes$concordant), concordance = obs[["conc"]],
             median_oriented_logFC = median(genes$oriented_logFC),
             mean_t = obs[["mean_t"]], permutations = B, r = r, p = (1 + r) / (B + 1),
             stringsAsFactors = FALSE)
}

summary_tab <- rbind(validate(CTL, "vs_controls"), validate(ISCH, "vs_ischemia"))
save_csv(summary_tab, "validation_summary.csv")
print(summary_tab, row.names = FALSE, digits = 4)
