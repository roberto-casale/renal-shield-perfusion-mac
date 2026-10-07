# 01_register_kernel.R -- registers the Jupyter kernel "ir-renaliri", which runs the
# notebooks with the project-local R library (env/rlib).

project <- local({.r <- Sys.getenv("RENAL_IRI_ROOT", ""); if (nzchar(.r) && dir.exists(.r)) .r else { .d <- normalizePath(getwd(), winslash = "/"); while (!file.exists(file.path(.d, "run_all.ps1")) && dirname(.d) != .d) .d <- dirname(.d); .d }})
rlib    <- file.path(project, "env", "rlib")
venv    <- file.path(project, "env", "py")
.libPaths(c(rlib, .libPaths()))

old_path <- Sys.getenv("PATH")
Sys.setenv(PATH = paste(file.path(venv, "Scripts"), old_path, sep = .Platform$path.sep))

IRkernel::installspec(
  name        = "ir-renaliri",
  displayname = "R (renal IRI)",
  prefix      = venv
)

kdir <- file.path(venv, "share", "jupyter", "kernels", "ir-renaliri")
kj   <- file.path(kdir, "kernel.json")
spec <- jsonlite::fromJSON(kj, simplifyVector = FALSE)
spec$env <- list(R_LIBS_USER = rlib, R_LIBS_SITE = rlib)
writeLines(jsonlite::toJSON(spec, auto_unbox = TRUE, pretty = TRUE), kj)
cat("Kernel written to:", kj, "\n")
cat(readLines(kj), sep = "\n"); cat("\n")
