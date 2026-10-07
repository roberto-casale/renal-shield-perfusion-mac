# pipeline_ora.R -- over-representation analysis with WebGestaltR (Section 4.5, notebook 2).
#
# Four databases (KEGG and Reactome pathways, DrugBank and GLAD4U drug-associated gene
# sets); categories of 5 to 2,000 genes; Benjamini-Hochberg FDR < 0.05 within each
# database. The background (reference) gene list is passed by the caller: in notebook 2
# it is the set of tested genes. Each answer is cached under data/cache/.

suppressWarnings(suppressMessages({
  source(file.path(local({.r <- Sys.getenv("RENAL_IRI_ROOT", ""); if (nzchar(.r) && dir.exists(.r)) .r else { .d <- normalizePath(getwd(), winslash = "/"); while (!file.exists(file.path(.d, "run_all.ps1")) && dirname(.d) != .d) .d <- dirname(.d); .d }}), "R", "helpers.R"))
  library(WebGestaltR); library(dplyr)
}))

# Databases queried.
ORA_DBS <- c("pathway_KEGG", "pathway_Reactome", "drug_DrugBank", "drug_GLAD4U")

# Parameters of the analysis.
ORA <- list(MIN_NUM = 5, MAX_NUM = 2000, FDR_THR = 0.05, REF = "genome_protein-coding")

# One database: the enriched terms, with the signature genes in each (`overlap_genes`).
ora_one <- function(genes, database, tag, background = NULL, refresh = FALSE) {
  genes <- unique(toupper(trimws(genes[!is.na(genes) & genes != ""])))
  background <- if (is.null(background)) NULL else
    unique(toupper(trimws(background[!is.na(background) & background != ""])))

  key <- sprintf("webgestalt_ora_%s_%s_%s_%s_%s_n%d_%d_f%g", tag, database,
                 if (is.null(background)) ORA$REF else "customref",
                 if (is.null(background)) "na" else digest_str(paste(sort(background), collapse = ",")),
                 digest_str(paste(sort(genes), collapse = ",")),
                 ORA$MIN_NUM, ORA$MAX_NUM, ORA$FDR_THR)
  r <- robust_fetch(key, refresh = refresh, fn = function() {
    res <- WebGestaltR::WebGestaltR(
      enrichMethod     = "ORA",
      organism         = "hsapiens",
      enrichDatabase   = database,
      interestGene     = genes,
      interestGeneType = "genesymbol",
      referenceGene     = background,
      referenceGeneType = if (is.null(background)) NULL else "genesymbol",
      referenceSet      = if (is.null(background)) ORA$REF else NULL,
      minNum           = ORA$MIN_NUM,
      maxNum           = ORA$MAX_NUM,
      sigMethod        = "fdr",
      fdrMethod        = "BH",
      fdrThr           = ORA$FDR_THR,
      isOutput         = FALSE,
      nThreads         = 2)

    if (is.null(res)) data.frame() else res
  })
  log_resource("WebGestalt", sprintf("ORA %s (%s): %d interest genes",
                                     database, tag, length(genes)))
  if (!r$ok || is.null(r$value) || !nrow(r$value)) {
    return(data.frame(database = character(0), geneSet = character(0),
                      description = character(0), size = integer(0),
                      overlap = integer(0), enrichmentRatio = numeric(0),
                      pValue = numeric(0), FDR = numeric(0),
                      overlap_genes = character(0), gene_set = character(0),
                      stringsAsFactors = FALSE))
  }
  d <- r$value

  ov <- if (!is.null(d$userId)) d$userId else if (!is.null(d$overlapId)) d$overlapId else NA_character_
  data.frame(
    database        = database,
    geneSet         = as.character(d$geneSet),
    description     = as.character(d$description),
    size            = as.integer(d$size),
    overlap         = as.integer(d$overlap),
    enrichmentRatio = as.numeric(d$enrichmentRatio),
    pValue          = as.numeric(d$pValue),
    FDR             = as.numeric(d$FDR),
    overlap_genes   = as.character(ov),
    gene_set        = tag,
    stringsAsFactors = FALSE)
}

# All four databases, combined and ordered by FDR.
run_ora_set <- function(genes, tag, background = NULL) {
  parts <- lapply(ORA_DBS, function(db) ora_one(genes, db, tag, background))
  res <- do.call(rbind, parts)
  if (is.null(res)) res <- data.frame()
  rownames(res) <- NULL
  res <- res[order(res$FDR, -res$enrichmentRatio), ]
  res
}
