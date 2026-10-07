# 00_install.R -- installs the R packages of the analysis into env/rlib, a library local
# to the project. The exact versions used for the article are listed in renv.lock.

project <- local({.r <- Sys.getenv("RENAL_IRI_ROOT", ""); if (nzchar(.r) && dir.exists(.r)) .r else { .d <- normalizePath(getwd(), winslash = "/"); while (!file.exists(file.path(.d, "run_all.ps1")) && dirname(.d) != .d) .d <- dirname(.d); .d }})
rlib    <- file.path(project, "env", "rlib")
dir.create(rlib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(rlib, .libPaths()))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org", lib = rlib)

options(
  repos          = c(CRAN = "https://cloud.r-project.org"),
  pkgType        = "binary",
  install.packages.check.source = "no",
  timeout        = 600,
  Ncpus          = max(1, parallel::detectCores() - 1)
)

message("R: ", R.version.string)
message("Bioconductor target: ", as.character(BiocManager::version()))
message(".libPaths():"); print(.libPaths())

cran_pkgs <- c("RobustRankAggreg", "babelgene", "matrixStats",
               "jsonlite", "rentrez",
               "WebGestaltR", "dplyr", "readxl",
               "ggplot2", "renv", "IRkernel")

bioc_pkgs <- c("GEOquery", "limma", "DESeq2", "Biobase",
               "org.Hs.eg.db", "AnnotationDbi",
               "org.Mm.eg.db", "org.Rn.eg.db")

have <- rownames(installed.packages())
need <- setdiff(c(cran_pkgs, bioc_pkgs), have)
message("Already available: ", paste(intersect(c(cran_pkgs, bioc_pkgs), have), collapse = ", "))
message("Will install: ", paste(need, collapse = ", "))

if (length(need)) {

  BiocManager::install(need, update = FALSE, ask = FALSE, lib = rlib,
                       type = "binary", Ncpus = getOption("Ncpus"))

  still <- setdiff(need, rownames(installed.packages(lib.loc = c(rlib, .libPaths()))))
  ann   <- intersect(still, c("org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db"))
  if (length(ann)) {
    message("No binary for: ", paste(ann, collapse = ", "), " -- retrying from source ",
            "(data packages: nothing is compiled)")
    install.packages(ann, lib = rlib, type = "source",
                     repos = BiocManager::repositories()[["BioCann"]])
  }
}

required <- c(

  "GEOquery","limma","DESeq2","Biobase",

  "RobustRankAggreg","babelgene","matrixStats",

  "Matrix",

  "WebGestaltR",

  "org.Hs.eg.db","org.Mm.eg.db","org.Rn.eg.db","AnnotationDbi",

  "jsonlite","rentrez","dplyr","readxl",

  "ggplot2","renv","IRkernel","BiocManager"
)
ok <- vapply(required, function(p) requireNamespace(p, quietly = TRUE), logical(1))
cat("\n==== PACKAGE LOAD CHECK ====\n")
print(data.frame(package = required, available = ok, row.names = NULL))
missing <- required[!ok]
if (length(missing)) {
  stop("MISSING PACKAGES: ", paste(missing, collapse = ", "))
} else {
  cat("\nALL REQUIRED PACKAGES AVAILABLE\n")
}
