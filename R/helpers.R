# helpers.R -- shared settings for every notebook and module.
#
# Finds the project root (RENAL_IRI_ROOT, or the nearest folder upward that contains
# run_all.ps1), creates the working folders, and defines logging, CSV output and
# robust_fetch(): every network request goes through it, is retried on failure, and
# its answer is cached under data/cache/, so a re-run gives the same result even if
# the remote service has changed or is unavailable.

# Project root.
.renal_root <- function() {
  r <- Sys.getenv("RENAL_IRI_ROOT", unset = "")
  if (nzchar(r) && dir.exists(r)) return(normalizePath(r, winslash = "/"))
  d <- normalizePath(getwd(), winslash = "/")
  while (!file.exists(file.path(d, "run_all.ps1")) && dirname(d) != d) d <- dirname(d)
  d
}
PROJECT <- .renal_root()
RLIB    <- file.path(PROJECT, "env", "rlib")
if (dir.exists(RLIB)) .libPaths(c(RLIB, .libPaths()))

# Working folders, created if missing.
DIRS <- list(
  data    = file.path(PROJECT, "data"),
  cache   = file.path(PROJECT, "data", "cache"),
  geo     = file.path(PROJECT, "data", "geo"),
  holdout = file.path(PROJECT, "data", "holdout"),
  results = file.path(PROJECT, "results"),
  figures = file.path(PROJECT, "figures"),
  logs    = file.path(PROJECT, "logs")
)
for (d in DIRS) dir.create(d, recursive = TRUE, showWarnings = FALSE)

CACHE_DIR <- DIRS$cache
LOG_FILE  <- file.path(DIRS$logs, "run.log")
RES_LOG   <- file.path(DIRS$logs, "resource_versions.tsv")

# Fixed seed, so that any random step is reproducible.
set.seed(1234)

# Per-dataset differential-expression flag (P_ADJ_MAX, LFC_MIN; descriptive only),
# threshold of the preliminary consensus flag in build_consensus() (RRA_FDR_MAX),
# and retries of robust_fetch() (FETCH_*).
CONST <- list(
  P_ADJ_MAX       = 0.05,
  LFC_MIN         = 0.585,
  RRA_FDR_MAX     = 0.05,
  FETCH_MAX_TRIES = 4,
  FETCH_BACKOFF   = 3
)

# Timestamped message to the console and to logs/run.log.
log_msg <- function(msg, level = "INFO") {
  stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line  <- sprintf("[%s] %-5s %s", stamp, level, msg)
  cat(line, "\n", file = LOG_FILE, append = TRUE)
  message(line)
  invisible(line)
}

# Records each query to an external resource in logs/resource_versions.tsv.
log_resource <- function(resource, detail = "") {
  stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  if (!file.exists(RES_LOG)) {
    cat("timestamp\tresource\tdetail\n", file = RES_LOG)
  }
  cat(sprintf("%s\t%s\t%s\n", stamp, resource, detail),
      file = RES_LOG, append = TRUE)
  invisible(TRUE)
}

# File-system-safe cache key, shortened with a hash when too long.
.safe_key <- function(key) {
  k <- gsub("[^A-Za-z0-9._-]", "_", key)
  if (nchar(k) > 120) {
    k <- paste0(substr(k, 1, 90), "_", substr(digest_str(key), 1, 16))
  }
  k
}

# Small deterministic string hash used in cache keys.
digest_str <- function(x) {
  h <- 0
  for (b in utf8ToInt(x)) h <- (h * 31 + b) %% 2147483647
  sprintf("%08x", as.integer(h))
}

# Runs fn() (a network request), retrying with increasing waits. The answer is saved
# as data/cache/<key>.rds and read from there on every later run.
robust_fetch <- function(key, fn, parse = identity, refresh = FALSE,
                          max_tries = CONST$FETCH_MAX_TRIES,
                          backoff = CONST$FETCH_BACKOFF,
                          cache_dir = CACHE_DIR) {
  cache_file <- file.path(cache_dir, paste0(.safe_key(key), ".rds"))
  if (!refresh && file.exists(cache_file)) {
    obj <- tryCatch(readRDS(cache_file), error = function(e) NULL)
    if (!is.null(obj)) {
      log_msg(sprintf("CACHE HIT   %s (queried %s)", key, obj$date))
      return(list(ok = TRUE, value = obj$value, cached = TRUE,
                  date = obj$date, key = key))
    }
  }
  last_err <- NULL
  for (attempt in seq_len(max_tries)) {
    res <- tryCatch(parse(fn()),
                    error = function(e) structure(
                      list(msg = conditionMessage(e)), class = "fetch_error"))
    if (!inherits(res, "fetch_error")) {
      date <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      saveRDS(list(value = res, date = date, key = key), cache_file)
      log_msg(sprintf("LIVE OK     %s (attempt %d/%d)", key, attempt, max_tries))
      return(list(ok = TRUE, value = res, cached = FALSE, date = date, key = key))
    }
    last_err <- res$msg
    if (attempt < max_tries) {
      wait <- min(30, backoff * 2^(attempt - 1))
      log_msg(sprintf("RETRY %d/%d  %s : %s (wait %ds)",
                      attempt, max_tries, key, last_err, wait), "WARN")
      Sys.sleep(wait)
    }
  }
  log_msg(sprintf("SKIP+LOG    %s : FAILED after %d attempts : %s",
                  key, max_tries, last_err), "ERROR")
  list(ok = FALSE, value = NULL, cached = FALSE, date = NA, key = key,
       error = last_err)
}

# Default value for a NULL or empty argument.
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# Writes a data frame to results/<name>.
save_csv <- function(df, name) {
  path <- file.path(DIRS$results, name)
  utils::write.csv(df, path, row.names = FALSE)
  log_msg(sprintf("WROTE       results/%s (%d rows)", name, nrow(df)))
  invisible(path)
}

log_msg("helpers.R loaded; constants and robust_fetch() ready.")
