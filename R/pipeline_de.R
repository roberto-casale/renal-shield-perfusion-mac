# pipeline_de.R -- the nine corpus series (notebook 1).
#
# Each load_GSE*() downloads one GEO series and assigns every sample to "control" or
# "IRI", as described in Table 1. run_de() then compares IRI with control using the
# model that matches the design (Section 4.2):
#   * microarrays: limma (empirical Bayes moderated t);
#   * GSE98622, processed RNA-seq values without counts (log2(x + 1) by maybe_log2()):
#     limma with trend = TRUE;
#   * GSE126805, RNA-seq counts from pre- and post-reperfusion biopsies of the same
#     kidney: DESeq2 with the kidney as a term in the model;
#   * GSE27274, cortex and medulla of the same rat: tissue as a covariate and the
#     animal as a random block (limma duplicateCorrelation).
# Every gene is kept with its log2 fold change and p-value. The column `sig`
# (BH-adjusted p < 0.05 and |log2 fold change| >= 0.585) is descriptive only.
#
# The file also holds the GEO query of Section 4.1 and the counts drawn in Figure 1.

suppressWarnings(suppressMessages({
  source(file.path(local({.r <- Sys.getenv("RENAL_IRI_ROOT", ""); if (nzchar(.r) && dir.exists(.r)) .r else { .d <- normalizePath(getwd(), winslash = "/"); while (!file.exists(file.path(.d, "run_all.ps1")) && dirname(.d) != .d) .d <- dirname(.d); .d }}), "R", "helpers.R"))
  library(GEOquery); library(Biobase); library(limma)
  library(DESeq2); library(readxl); library(matrixStats)
}))
options(GEOquery.inmemory.gpl = FALSE)

# Puts expression values on the log2 scale when they are still linear.
maybe_log2 <- function(mat) {

  qx <- as.numeric(stats::quantile(mat, c(0, 0.25, 0.5, 0.75, 0.99, 1),
                                   na.rm = TRUE))
  linear <- (qx[6] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
  if (linear) {
    mat[mat < 0] <- NA
    mat <- log2(mat + 1)
    attr(mat, "logged") <- TRUE
  } else attr(mat, "logged") <- FALSE
  mat
}

# One row per gene symbol: when several probes share a symbol, keep the most expressed.
collapse_to_symbol <- function(mat, symbols) {

  keep <- !is.na(symbols) & symbols != "" & symbols != "---"
  mat <- mat[keep, , drop = FALSE]; symbols <- symbols[keep]
  o <- order(rowMeans(mat, na.rm = TRUE), decreasing = TRUE)
  mat <- mat[o, , drop = FALSE]; symbols <- symbols[o]
  dup <- duplicated(symbols)
  mat <- mat[!dup, , drop = FALSE]
  rownames(mat) <- symbols[!dup]
  mat
}

# First gene symbol of a multi-symbol probe annotation.
first_symbol <- function(x) {

  x <- as.character(x)
  x <- sub("\\s*///.*$", "", x)
  has <- grepl("//", x)
  if (any(has)) {
    x[has] <- vapply(strsplit(x[has], "\\s*//\\s*"), function(p)
      if (length(p) >= 2) p[2] else NA_character_, character(1))
  }
  trimws(x)
}

# The gene-symbol column of a GEO platform table, whatever its name.
find_symbol_col <- function(fd) {
  cands <- c("GENE_SYMBOL","Gene Symbol","Gene symbol","GeneSymbol","SYMBOL",
             "gene_symbol","ILMN_Gene","Symbol","GENE")
  for (cc in cands) if (!is.null(fd[[cc]])) return(fd[[cc]])
  sc <- grep("symbol", names(fd), ignore.case = TRUE, value = TRUE)
  if (length(sc)) return(fd[[sc[1]]])
  NULL
}

# GEO series matrix (cached).
get_geo_es <- function(acc) {
  r <- robust_fetch(paste0("geo_probe_", acc), fn = function()
    getGEO(acc, GSEMatrix = TRUE, getGPL = TRUE, destdir = DIRS$geo))
  if (!r$ok) stop("getGEO failed for ", acc)
  r$value[[1]]
}

# GSE30718, human allografts: acute kidney injury vs stable protocol biopsies.
# Repeat biopsies of the same patient are dropped, keeping the earliest.
load_GSE30718 <- function() {
  es <- get_geo_es("GSE30718")
  mat <- exprs(es); fd <- fData(es); pd <- pData(es)
  mat <- maybe_log2(mat)
  mat <- collapse_to_symbol(mat, first_symbol(fd[["Gene Symbol"]]))

  dg <- as.character(pd[["diagnosis:ch1"]])
  grp <- ifelse(grepl("Acute kidney injury", dg), "IRI",
         ifelse(grepl("pristine", dg), "control", NA))
  meta <- data.frame(sample = colnames(mat), group = grp,
                     title = as.character(pd[["title"]]), stringsAsFactors = FALSE)

  pat <- sub("^.*?[ ]([0-9]+)A[0-9].*$", "\\1", meta$title)
  bx  <- suppressWarnings(as.integer(sub("^.*?[ ][0-9]+A([0-9]).*$", "\\1", meta$title)))
  dup <- meta$group %in% "IRI" & !is.na(bx)
  drop <- rep(FALSE, nrow(meta))
  for (pt in unique(pat[dup])) {
    i <- which(dup & pat == pt)
    if (length(i) > 1) drop[i[bx[i] > min(bx[i])]] <- TRUE
  }
  if (any(drop))
    log_msg(sprintf("GSE30718: dropped %d repeat biopsies from %d patients (earliest kept per GEO)",
                    sum(drop), length(unique(pat[drop]))))
  keep <- !is.na(meta$group) & !drop
  list(mat = mat[, keep], meta = meta[keep, ], is_counts = FALSE, paired = FALSE,
       organism = "Homo sapiens", platform = "GPL570", modality = "transcriptomics",
       contrast = "allograft AKI vs stable (pristine) allograft protocol biopsy")
}

# GSE43974, human deceased-donor kidneys: biopsy after reperfusion (T3) vs after
# cold ischemia (T2), brain-dead and cardiac-dead donors, compared as unpaired groups.
load_GSE43974 <- function() {
  es <- get_geo_es("GSE43974")
  mat <- exprs(es); fd <- fData(es); pd <- pData(es)
  mat <- maybe_log2(mat)
  mat <- collapse_to_symbol(mat, first_symbol(fd[["ILMN_Gene"]]))
  ti  <- as.character(pd$title)
  m   <- regmatches(ti, regexec("Kidney_donor_([A-Za-z]+)_T([0-9])_([0-9]+)", ti))
  dtype <- sapply(m, function(z) if (length(z) == 4) z[2] else NA)
  tp    <- as.integer(sapply(m, function(z) if (length(z) == 4) z[3] else NA))

  dec <- dtype %in% c("BD", "DCD")
  grp <- ifelse(dec & tp == 2, "control", ifelse(dec & tp == 3, "IRI", NA))
  meta0 <- data.frame(sample = colnames(mat), group = grp, stringsAsFactors = FALSE)
  meta0 <- meta0[!is.na(meta0$group), ]

  log_msg(sprintf("GSE43974: %d control (T2, post-cold-ischemia) vs %d IRI (T3, post-reperfusion), unpaired, full cohort",
                  sum(meta0$group == "control"), sum(meta0$group == "IRI")))
  meta <- data.frame(sample = meta0$sample, group = meta0$group, stringsAsFactors = FALSE)
  list(mat = mat[, meta$sample], meta = meta, is_counts = FALSE, paired = FALSE,
       organism = "Homo sapiens", platform = "GPL10558", modality = "transcriptomics",
       contrast = "deceased-donor transplant: post-reperfusion (T3) vs post-cold-ischemia (T2), unpaired")
}

# GSE126805, human transplants, RNA-seq counts: post- vs pre-reperfusion biopsy of the
# same kidney.
load_GSE126805 <- function() {
  dir <- file.path(DIRS$geo, "GSE126805")
  f <- file.path(dir, "GSE126805_all.gene_counts.txt.gz")
  if (!file.exists(f))
    robust_fetch("geo_supp_GSE126805_counts", fn = function()
      getGEOSuppFiles("GSE126805", baseDir = DIRS$geo, filter_regex = "gene_counts"))
  d <- read.delim(gzfile(f), check.names = FALSE, stringsAsFactors = FALSE)
  sym <- d$external_gene_name
  scols <- grep("^sample\\.", colnames(d), value = TRUE)

  sym0 <- as.character(sym)
  cnt <- as.matrix(d[, scols]); rownames(cnt) <- make.unique(sym0)

  pre  <- grep("baseline_pre$",  scols, value = TRUE)
  post <- grep("baseline_post$", scols, value = TRUE)
  kid_pre  <- sub("^sample\\.(kidney_[0-9]+)_baseline_pre$",  "\\1", pre)
  kid_post <- sub("^sample\\.(kidney_[0-9]+)_baseline_post$", "\\1", post)
  common <- intersect(kid_pre, kid_post)
  sel_pre  <- pre [match(common, kid_pre)]
  sel_post <- post[match(common, kid_post)]
  mat <- cbind(cnt[, sel_pre, drop = FALSE], cnt[, sel_post, drop = FALSE])

  stopifnot(length(sym0) == nrow(mat))
  mat <- rowsum(mat, group = sym0)
  meta <- data.frame(
    sample  = colnames(mat),
    group   = rep(c("control", "IRI"), c(length(sel_pre), length(sel_post))),
    subject = c(common, common), stringsAsFactors = FALSE)
  list(mat = mat, meta = meta, is_counts = TRUE, paired = TRUE,
       organism = "Homo sapiens", platform = "GPL21290 (RNA-seq)",
       modality = "transcriptomics",
       contrast = "post-reperfusion vs pre-reperfusion biopsy, paired within kidney")
}

# GSE39548, mouse: renal IRI vs naive control, both without a protective manoeuvre.
load_GSE39548 <- function() {
  es <- get_geo_es("GSE39548")
  mat <- exprs(es); fd <- fData(es); pd <- pData(es)
  mat <- maybe_log2(mat)
  sym <- fd[["GENE_SYMBOL"]]; if (is.null(sym)) sym <- fd[["GENE_NAME"]]
  mat <- collapse_to_symbol(mat, first_symbol(sym))
  proc <- as.character(pd[["surgery procedure:ch1"]])

  prot <- as.character(pd[["protection maneuver against iri:ch1"]])
  grp <- ifelse(proc == "none" & prot == "none", "control",
         ifelse(proc == "IRI"  & prot == "none", "IRI", NA))
  meta <- data.frame(sample = colnames(mat), group = grp, stringsAsFactors = FALSE)
  keep <- !is.na(meta$group)
  list(mat = mat[, keep], meta = meta[keep, ], is_counts = FALSE, paired = FALSE,
       organism = "Mus musculus", platform = "GPL7202", modality = "transcriptomics",
       contrast = "renal IRI vs naive control (unprotected)")
}

# GSE98622, mouse, RNA-seq values from the supplementary file: IRI (2-72 h) vs sham.
load_GSE98622 <- function() {
  dir <- file.path(DIRS$geo, "GSE98622")
  f <- list.files(dir, pattern = "xlsx$", full.names = TRUE)
  if (!length(f)) {
    robust_fetch("geo_supp_GSE98622_xlsx", fn = function()
      getGEOSuppFiles("GSE98622", baseDir = DIRS$geo, filter_regex = "xlsx"))
    f <- list.files(dir, pattern = "xlsx$", full.names = TRUE)
  }
  d <- suppressMessages(readxl::read_excel(f[1], sheet = 1))
  d <- as.data.frame(d)
  sym <- d$symbol
  scols <- setdiff(colnames(d), c(colnames(d)[1], "symbol", "type"))
  sym0 <- as.character(sym)
  mat <- as.matrix(d[, scols]); rownames(mat) <- make.unique(sym0)

  iri  <- grep("^IRI(2h|4h|24h|48h|72h)-", scols, value = TRUE)
  sham <- grep("^SHAM(4h|24h)-", scols, value = TRUE)
  sel <- c(sham, iri)
  mat <- mat[, sel, drop = FALSE]
  mat <- maybe_log2(mat)

  stopifnot(length(sym0) == nrow(mat))
  mat <- collapse_to_symbol(mat, sym0)
  meta <- data.frame(sample = sel,
                     group = rep(c("control", "IRI"), c(length(sham), length(iri))),
                     stringsAsFactors = FALSE)

  list(mat = mat, meta = meta, is_counts = FALSE, paired = FALSE,
       rnaseq_continuous = TRUE,
       organism = "Mus musculus", platform = "GPL13112 (RNA-seq)",
       modality = "transcriptomics",
       contrast = "acute renal IRI (2-72h) vs sham surgery")
}

# GSE27274, rat: IRI vs sham; cortex and medulla of the same 12 rats (24 arrays).
load_GSE27274 <- function() {
  es <- get_geo_es("GSE27274")
  mat <- exprs(es); fd <- fData(es); pd <- pData(es)
  mat <- maybe_log2(mat)
  sym <- fd[["ILMN_Gene"]]; if (is.null(sym)) sym <- fd[["Symbol"]]
  mat <- collapse_to_symbol(mat, first_symbol(sym))
  tr <- as.character(pd[["treatment:ch1"]])
  grp <- ifelse(grepl("sham", tr), "control",
         ifelse(grepl("ischemic-reperfusion", tr), "IRI", NA))
  tissue <- as.character(pd[["tissue:ch1"]])

  ttl <- as.character(pd[["title"]])
  rat <- sub("^(.*?)\\s+(cortex|medulla)\\s+([0-9]+)$", "\\1_\\3", ttl)
  meta <- data.frame(sample = colnames(mat), group = grp, tissue = tissue,
                     animal = rat, stringsAsFactors = FALSE)
  keep <- !is.na(meta$group)
  stopifnot(length(unique(meta$animal[keep])) == 12L,
            all(table(meta$animal[keep]) == 2L))
  list(mat = mat[, keep], meta = meta[keep, ], is_counts = FALSE, paired = FALSE,
       covariate = "tissue", block = rat[keep], accession = "GSE27274",
       organism = "Rattus norvegicus", platform = "GPL6101", modality = "transcriptomics",
       contrast = "renal IRI (6/24/120h, cortex+medulla) vs sham; cortex and medulla of the same rat modelled as a random block (12 animals, 24 arrays)")
}

# GSE58438, rat: IRI with saline vehicle vs naive control; drug-treated arms excluded.
load_GSE58438 <- function() {
  es <- get_geo_es("GSE58438")
  mat <- exprs(es); fd <- fData(es); pd <- pData(es)
  mat <- maybe_log2(mat)
  mat <- collapse_to_symbol(mat, first_symbol(fd[["gene_assignment"]]))
  inj <- as.character(pd[["injury:ch1"]])
  trt <- as.character(pd[["treatment:ch1"]])
  grp <- ifelse(inj == "Yes", "IRI", ifelse(inj == "No", "control", NA))
  meta <- data.frame(sample = colnames(mat), group = grp, treatment = trt,
                     stringsAsFactors = FALSE)
  keep <- !is.na(meta$group) & meta$treatment %in% c("Control", "Saline")
  list(mat = mat[, keep], meta = meta[keep, ], is_counts = FALSE, paired = FALSE,
       organism = "Rattus norvegicus", platform = "GPL11534", modality = "transcriptomics",
       contrast = "renal IRI, saline vehicle vs naive no-injury control (VPA and dexamethasone arms excluded)")
}

# GSE34351, mouse: wild-type ischemic kidney vs sham.
load_GSE34351 <- function() {
  es <- get_geo_es("GSE34351")
  mat <- exprs(es); fd <- fData(es); pd <- pData(es)
  mat <- maybe_log2(mat)
  mat <- collapse_to_symbol(mat, first_symbol(find_symbol_col(fd)))
  geno <- as.character(pd[["genotype:ch1"]])
  surg <- as.character(pd[["surgery:ch1"]])
  grp <- ifelse(surg == "Ischemic", "IRI", ifelse(surg == "Sham", "control", NA))
  meta <- data.frame(sample = colnames(mat), group = grp, genotype = geno,
                     stringsAsFactors = FALSE)
  keep <- !is.na(meta$group) & meta$genotype == "WT"
  list(mat = mat[, keep], meta = meta[keep, ], is_counts = FALSE, paired = FALSE,
       organism = "Mus musculus", platform = "GPL1261", modality = "transcriptomics",
       contrast = "wild-type ischemic kidney vs sham")
}

# GSE148420, rat: unilateral IRI (day 3 and day 7) vs sham operation.
load_GSE148420 <- function() {
  es <- get_geo_es("GSE148420")
  mat <- exprs(es); fd <- fData(es); pd <- pData(es)
  mat <- maybe_log2(mat)
  mat <- collapse_to_symbol(mat, first_symbol(find_symbol_col(fd)))
  op <- as.character(pd[["operation:ch1"]])
  grp <- ifelse(grepl("ischemia reperfusion", op), "IRI",
         ifelse(grepl("Sham", op), "control", NA))
  sac <- as.character(pd[["sacrifice:ch1"]])
  meta <- data.frame(sample = colnames(mat), group = grp, sacrifice = sac,
                     stringsAsFactors = FALSE)

  keep <- !is.na(meta$group)
  list(mat = mat[, keep], meta = meta[keep, ], is_counts = FALSE, paired = FALSE,
       organism = "Rattus norvegicus", platform = "GPL14746", modality = "transcriptomics",
       contrast = "unilateral renal IRI (day 3/7) vs sham operation")
}

# Differential expression, IRI vs control, with the model that matches the design.
run_de <- function(obj) {

  obj$meta$group <- factor(obj$meta$group, levels = c("control", "IRI"))
  if (obj$is_counts) {
    cnt <- round(obj$mat)
    keep <- rowSums(cnt >= 10) >= max(2, floor(0.2 * ncol(cnt)))
    cnt <- cnt[keep, , drop = FALSE]
    cd <- obj$meta; rownames(cd) <- cd$sample
    if (isTRUE(obj$paired)) {
      cd$subject <- factor(cd$subject)
      dds <- DESeqDataSetFromMatrix(cnt, cd, design = ~ subject + group)
    } else {
      dds <- DESeqDataSetFromMatrix(cnt, cd, design = ~ group)
    }
    dds <- DESeq(dds, quiet = TRUE)
    res <- as.data.frame(results(dds, name = "group_IRI_vs_control"))
    out <- data.frame(gene = rownames(res), logFC = res$log2FoldChange,
                      AveExpr = log2(res$baseMean + 1),
                      P.Value = res$pvalue, adj.P.Val = res$padj,
                      stringsAsFactors = FALSE)
  } else {
    mat <- obj$mat
    mat <- mat[rowSums(!is.na(mat)) >= ncol(mat) * 0.7, , drop = FALSE]

    mat <- mat[rowSums(mat != 0, na.rm = TRUE) > 0, , drop = FALSE]
    covariate <- obj$covariate %||% NULL
    if (isTRUE(obj$paired)) {
      subject <- factor(obj$meta$subject)
      design <- model.matrix(~ subject + obj$meta$group)
      colnames(design) <- make.names(colnames(design))
      coef <- tail(colnames(design), 1)
      fit <- lmFit(mat, design); fit <- eBayes(fit, trend = isTRUE(obj$rnaseq_continuous))
      res <- topTable(fit, coef = coef, number = Inf, sort.by = "none")
    } else if (!is.null(covariate)) {
      cov <- factor(obj$meta[[covariate]])
      design <- model.matrix(~ cov + obj$meta$group)
      colnames(design) <- make.names(colnames(design))
      coef <- tail(colnames(design), 1)

      blk <- obj$block %||% NULL
      if (!is.null(blk)) {
        blk <- factor(blk)
        dc  <- limma::duplicateCorrelation(mat, design, block = blk)
        fit <- lmFit(mat, design, block = blk, correlation = dc$consensus)
        log_msg(sprintf("  %s: random block on %d units, consensus correlation %.4f",
                        obj$accession %||% "dataset", nlevels(blk), dc$consensus))
      } else {
        fit <- lmFit(mat, design)
      }
      fit <- eBayes(fit, trend = isTRUE(obj$rnaseq_continuous))
      res <- topTable(fit, coef = coef, number = Inf, sort.by = "none")
    } else {
      design <- model.matrix(~ obj$meta$group)
      fit <- lmFit(mat, design); fit <- eBayes(fit, trend = isTRUE(obj$rnaseq_continuous))
      res <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
    }
    out <- data.frame(gene = rownames(res), logFC = res$logFC,
                      AveExpr = res$AveExpr, P.Value = res$P.Value,
                      adj.P.Val = res$adj.P.Val, stringsAsFactors = FALSE)
  }
  out <- out[!is.na(out$P.Value) & !is.na(out$logFC), ]
  out <- out[order(out$P.Value), ]
  out$sig <- (out$adj.P.Val < CONST$P_ADJ_MAX) &
             (abs(out$logFC) >= CONST$LFC_MIN)
  out$sig[is.na(out$sig)] <- FALSE
  out
}

# The GEO query of Section 4.1.
GEO_QUERY <- paste0(
  '(kidney[Title] OR renal[Title]) AND (ischemi*[Title] OR reperfusion[Title]) ',
  'AND ("expression profiling by array"[Filter] OR ',
  '"expression profiling by high throughput sequencing"[Filter]) AND gse[Filter]')

# Number of GEO series returned by the query. The answer is cached, because the
# count can change over time.
geo_identification_search <- function(query = GEO_QUERY,
                                      out = file.path(DIRS$results, "geo_identification.csv")) {

  r <- robust_fetch(paste0("geo_identification_", digest_str(query)), fn = function()
        as.integer(rentrez::entrez_search(db = "gds", term = query, retmax = 0)$count))
  if (!isTRUE(r$ok)) stop("GEO identification search failed: ", r$error)
  n <- as.integer(r$value)
  log_msg(sprintf("GEO identification: %d series match the declared query", n))

  utils::write.csv(data.frame(query = query, n_identified = n,
                              searched_on = substr(r$date, 1, 10)), out, row.names = FALSE)
  n
}

# Which corpus series, and the validation series, the query returns (Figure 1).
geo_corpus_coverage <- function(query = GEO_QUERY,
                                accessions = names(REGISTRY),
                                holdout = "GSE267650",
                                out = file.path(DIRS$results, "geo_corpus_coverage.csv"),
                                out_all = file.path(DIRS$results,
                                                    "geo_identified_series.csv")) {

  r <- robust_fetch(paste0("geo_serieslist_", digest_str(query), "_retmax500_summary"),
                    fn = function() {
                      h <- rentrez::entrez_search(db = "gds", term = query, retmax = 500)
                      s <- rentrez::entrez_summary(db = "gds", id = h$ids)
                      g <- function(f) unname(vapply(s, function(x) {
                        v <- x[[f]]
                        if (is.null(v) || !length(v)) NA_character_ else as.character(v)[1]
                      }, character(1)))
                      data.frame(accession = g("accession"), title = g("title"),
                                 organism = g("taxon"), n_samples = g("n_samples"),
                                 gds_type = g("gdstype"), stringsAsFactors = FALSE)
                    })
  if (!isTRUE(r$ok)) stop("GEO series listing failed: ", r$error)
  all <- r$value
  hits <- all$accession

  all$in_corpus <- all$accession %in% accessions
  all$held_out  <- all$accession %in% holdout
  utils::write.csv(all[order(all$accession), ], out_all, row.names = FALSE)
  log_msg(sprintf("WROTE       results/%s (%d rows)", basename(out_all), nrow(all)))

  tab <- data.frame(
    accession = c(accessions, holdout),
    role = c(rep("corpus", length(accessions)), rep("held-out", length(holdout))),
    returned_by_query = c(accessions, holdout) %in% hits,
    stringsAsFactors = FALSE)
  utils::write.csv(tab, out, row.names = FALSE)
  log_msg(sprintf("corpus coverage: %d of %d analysed series returned by the declared query; %d identified separately",
                  sum(tab$returned_by_query[tab$role == "corpus"]), length(accessions),
                  sum(!tab$returned_by_query[tab$role == "corpus"])))
  tab
}

# The nine corpus series, in the order of Table 1.
REGISTRY <- list(
  GSE30718  = load_GSE30718,
  GSE43974  = load_GSE43974,
  GSE126805 = load_GSE126805,
  GSE39548  = load_GSE39548,
  GSE98622  = load_GSE98622,
  GSE34351  = load_GSE34351,
  GSE27274  = load_GSE27274,
  GSE58438  = load_GSE58438,
  GSE148420 = load_GSE148420
)
