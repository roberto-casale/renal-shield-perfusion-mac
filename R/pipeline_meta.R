# pipeline_meta.R -- from nine per-dataset tables to one rank aggregation (notebook 2).
#
# 1. All genes are expressed as human gene symbols (Section 4.3). Human symbols are
#    checked against org.Hs.eg.db. Rodent symbols are first updated to their current
#    name (org.Mm.eg.db, org.Rn.eg.db) and then mapped to human orthologs with
#    babelgene (default settings: pairs supported by at least 3 orthology databases,
#    best-supported match). Genes without a human ortholog are dropped. When several
#    rodent genes map to the same human gene, the one with the smallest p-value is
#    kept, so that each human gene appears once per dataset.
# 2. Each dataset gives two lists, genes with log2 fold change > 0 and genes with
#    log2 fold change < 0, each ordered by p-value; no threshold is applied.
# 3. Robust Rank Aggregation combines each direction across the nine datasets, with
#    every list normalised by its own length (rankMatrix(full = TRUE)). Each gene is
#    assigned the direction with the smaller aggregation score.

suppressWarnings(suppressMessages({
  source(file.path(local({.r <- Sys.getenv("RENAL_IRI_ROOT", ""); if (nzchar(.r) && dir.exists(.r)) .r else { .d <- normalizePath(getwd(), winslash = "/"); while (!file.exists(file.path(.d, "run_all.ps1")) && dirname(.d) != .d) .d <- dirname(.d); .d }}), "R", "helpers.R"))
  library(babelgene); library(RobustRankAggreg); library(dplyr)
}))

# Official HGNC symbols and unambiguous aliases (org.Hs.eg.db), loaded once.
.hgnc_env <- new.env(parent = emptyenv())

.load_hgnc <- function() {
  if (!is.null(.hgnc_env$official)) return(invisible(TRUE))

  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE) ||
      !requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop("org.Hs.eg.db and AnnotationDbi are REQUIRED for symbol harmonisation. ",
         "Install them (they are in R/00_install.R's bioc_pkgs) and re-run. ",
         "Continuing without them would silently change the signature.")
  }
  suppressMessages({
    off <- AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "SYMBOL")
    a2s <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
             keys = AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "ALIAS"),
             keytype = "ALIAS", columns = "SYMBOL")
  })
  a2s <- a2s[!is.na(a2s$ALIAS) & !is.na(a2s$SYMBOL), ]

  a2s$KEY <- toupper(a2s$ALIAS)
  amb <- names(which(tapply(a2s$SYMBOL, a2s$KEY, function(z) length(unique(z))) > 1))
  a2s <- a2s[!(a2s$KEY %in% amb), ]
  a2s <- a2s[!(a2s$KEY %in% toupper(off)), ]
  .hgnc_env$official <- off
  .hgnc_env$alias <- setNames(a2s$SYMBOL, a2s$KEY)
  invisible(TRUE)
}

# Gene symbols that spreadsheet software converts into dates (e.g. SEPT2 -> 2-Sep).
EXCEL_MONTHS <- c(Mar = "MARCH", Sep = "SEPT", Dec = "DEC", Oct = "OCT",
                  Feb = "FEB", Apr = "APR", Nov = "NOV")

# Official symbols that contain a hyphen and must keep it.
BLOCK_DEHYPHEN <- c("MT-TP")

# Human symbol -> official HGNC symbol: repairs date-converted names, ORF case and
# hyphens, then replaces unambiguous aliases.
harmonise_human_symbol <- function(x) {
  .load_hgnc()
  off <- .hgnc_env$official; ali <- .hgnc_env$alias
  y <- trimws(as.character(x))

  m <- regmatches(y, regexec("^([0-9]{1,2})-(Mar|Sep|Dec|Oct|Feb|Apr|Nov)$", y))
  hit <- lengths(m) == 3
  if (any(hit)) y[hit] <- vapply(m[hit], function(z)
    paste0(EXCEL_MONTHS[[z[3]]], z[2]), character(1))
  m2 <- regmatches(y, regexec("^(Mar|Sep|Dec|Oct|Feb|Apr|Nov)-([0-9]{1,2})$", y))
  hit2 <- lengths(m2) == 3
  if (any(hit2)) y[hit2] <- vapply(m2[hit2], function(z)
    paste0(EXCEL_MONTHS[[z[2]]], z[3]), character(1))

  up <- toupper(y)

  orf <- grepl("^C[0-9XY]+ORF[0-9]+$", up)
  if (any(orf)) {
    cand <- sub("ORF", "orf", up[orf])
    take <- cand %in% off
    up[orf][take] <- cand[take]
  }

  hy <- grepl("-", up) & !(up %in% off) & !(up %in% BLOCK_DEHYPHEN)
  if (any(hy)) {
    cand <- gsub("-", "", up[hy])
    take <- cand %in% off
    up[hy][take] <- cand[take]
  }

  na <- !(up %in% off) & (up %in% names(ali))
  if (any(na)) up[na] <- unname(ali[up[na]])
  up
}

# Current mouse and rat symbols and unambiguous aliases, loaded once per species.
.rodent_env <- new.env(parent = emptyenv())
RODENT_DB <- c("Mus musculus" = "org.Mm.eg.db", "Rattus norvegicus" = "org.Rn.eg.db")

# Outdated rodent symbol -> current symbol, so that babelgene recognises it.
harmonise_rodent_symbol <- function(x, organism) {
  pkg <- RODENT_DB[[organism]]
  if (is.null(pkg)) return(x)
  if (is.null(.rodent_env[[pkg]])) {
    if (!requireNamespace(pkg, quietly = TRUE) ||
        !requireNamespace("AnnotationDbi", quietly = TRUE)) {
      stop(pkg, " is REQUIRED to harmonise ", organism, " symbols before ortholog ",
           "mapping. Without it about 10,000 rodent genes are dropped silently. ",
           "It is listed in R/00_install.R.")
    }
    db <- getExportedValue(pkg, pkg)
    suppressMessages({
      off <- AnnotationDbi::keys(db, keytype = "SYMBOL")
      a2s <- AnnotationDbi::select(db, keys = AnnotationDbi::keys(db, keytype = "ALIAS"),
                                   keytype = "ALIAS", columns = "SYMBOL")
    })
    a2s <- a2s[!is.na(a2s$ALIAS) & !is.na(a2s$SYMBOL), ]
    a2s$KEY <- toupper(a2s$ALIAS)
    amb <- names(which(tapply(a2s$SYMBOL, a2s$KEY, function(z) length(unique(z))) > 1))
    a2s <- a2s[!(a2s$KEY %in% amb) & !(a2s$KEY %in% toupper(off)), ]
    .rodent_env[[pkg]] <- list(off = toupper(off),
                               alias = setNames(a2s$SYMBOL, a2s$KEY))
  }
  m <- .rodent_env[[pkg]]
  y <- trimws(as.character(x)); k <- toupper(y)
  repair <- !(k %in% m$off) & (k %in% names(m$alias))
  y[repair] <- unname(m$alias[k[repair]])
  attr(y, "n_repaired") <- sum(repair)
  y
}

# Adds the human gene symbol of every gene (column `human`). For rodent data the
# suffixes some platforms add (_PREDICTED, _MAPPED) are removed first.
map_to_human <- function(de, organism) {
  de$gene <- as.character(de$gene)
  if (organism == "Homo sapiens") {
    de$human <- harmonise_human_symbol(de$gene)
    return(de)
  }
  sp <- c("Mus musculus" = "mouse", "Rattus norvegicus" = "rat")[organism]

  clean <- trimws(sub("_(PREDICTED|MAPPED)$", "", de$gene, ignore.case = TRUE))
  titlecase <- paste0(toupper(substr(clean, 1, 1)), tolower(substring(clean, 2)))

  updated <- harmonise_rodent_symbol(titlecase, organism)
  if (!is.null(attr(updated, "n_repaired")) && attr(updated, "n_repaired") > 0)
    log_msg(sprintf("  %s: %d outdated rodent symbols updated before ortholog mapping",
                    organism, attr(updated, "n_repaired")))
  lookup <- unique(rbind(
    data.frame(gene = de$gene, variant = clean,     stringsAsFactors = FALSE),
    data.frame(gene = de$gene, variant = titlecase, stringsAsFactors = FALSE),
    data.frame(gene = de$gene, variant = as.character(updated), stringsAsFactors = FALSE)))
  orth <- babelgene::orthologs(genes = unique(lookup$variant),
                               species = sp, human = FALSE)
  map <- unique(orth[, c("symbol", "human_symbol")])
  mm <- merge(lookup, map, by.x = "variant", by.y = "symbol")
  mm <- unique(mm[, c("gene", "human_symbol")])
  merged <- merge(de, mm, by = "gene")
  merged$human <- harmonise_human_symbol(merged$human_symbol)
  merged$human_symbol <- NULL
  merged
}

# One row per human gene, keeping the row with the smallest p-value.
collapse_human <- function(de) {
  de <- de[!is.na(de$human) & de$human != "", ]
  de <- de[order(de$P.Value), ]
  de <- de[!duplicated(de$human), ]
  de
}

# One list per dataset: genes going up (or down), ordered by p-value.
build_ranklists <- function(human_de, direction = c("up", "down")) {
  direction <- match.arg(direction)
  lapply(human_de, function(d) {
    d <- d[!is.na(d$logFC) & !is.na(d$P.Value), ]
    d <- if (direction == "up") d[d$logFC > 0, ] else d[d$logFC < 0, ]
    d <- d[order(d$P.Value), ]
    d$human
  })
}

# Robust Rank Aggregation, each list normalised by its own length.
run_rra <- function(ranklists, universe_n = NULL) {

  ranklists <- ranklists[sapply(ranklists, length) > 0]
  rmat <- RobustRankAggreg::rankMatrix(ranklists, full = TRUE)
  ag <- RobustRankAggreg::aggregateRanks(rmat = rmat, method = "RRA")
  ag$FDR <- p.adjust(ag$Score, method = "BH")
  ag
}

# Both directions aggregated, and each gene summarised across datasets: datasets
# measuring it, how many go up or down, mean log2 fold change, aggregation scores,
# assigned direction.
build_consensus <- function(human_de) {

  up   <- run_rra(build_ranklists(human_de, "up"))
  down <- run_rra(build_ranklists(human_de, "down"))
  names(up)[names(up) == "Score"]   <- "RRA_up"
  names(down)[names(down) == "Score"] <- "RRA_down"

  all_genes <- unique(unlist(lapply(human_de, function(d) d$human)))
  lfc_mat <- sapply(human_de, function(d) d$logFC[match(all_genes, d$human)])
  rownames(lfc_mat) <- all_genes
  n_up   <- rowSums(lfc_mat > 0, na.rm = TRUE)
  n_down <- rowSums(lfc_mat < 0, na.rm = TRUE)
  n_meas <- rowSums(!is.na(lfc_mat))
  mean_lfc <- rowMeans(lfc_mat, na.rm = TRUE)

  cons <- data.frame(human = all_genes,
                     n_datasets = n_meas, n_up = n_up, n_down = n_down,
                     mean_logFC = mean_lfc, stringsAsFactors = FALSE)
  cons$RRA_up   <- up$RRA_up[match(cons$human, up$Name)]
  cons$RRA_down <- down$RRA_down[match(cons$human, down$Name)]
  cons$FDR_up   <- up$FDR[match(cons$human, up$Name)]
  cons$FDR_down <- down$FDR[match(cons$human, down$Name)]
  cons$RRA_up[is.na(cons$RRA_up)]     <- 1
  cons$RRA_down[is.na(cons$RRA_down)] <- 1
  cons$FDR_up[is.na(cons$FDR_up)]     <- 1
  cons$FDR_down[is.na(cons$FDR_down)] <- 1

  cons$direction <- ifelse(cons$RRA_up <= cons$RRA_down, "up", "down")
  cons$meta_FDR  <- ifelse(cons$direction == "up", cons$FDR_up, cons$FDR_down)
  cons$concordant <- ifelse(cons$direction == "up", cons$n_up, cons$n_down)

  cons$signed_score <- ifelse(cons$direction == "up",
                              -log10(cons$RRA_up), log10(cons$RRA_down))

  cons$min_conc <- pmax(2L, ceiling(cons$n_datasets / 2))
  cons$consensus <- cons$meta_FDR < CONST$RRA_FDR_MAX &
                    cons$concordant >= cons$min_conc
  cons <- cons[order(cons$meta_FDR, -abs(cons$signed_score)), ]

  cons
}
