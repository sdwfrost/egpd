#!/usr/bin/env Rscript
# Precompute the Monte Carlo results displayed in
# vignettes/parameter-coverage.Rmd.
#
# The vignette shows its simulation chunks with `eval = FALSE` and loads the
# summaries saved by this script from vignettes/parameter-coverage-results.rds.
# The helper and simulation code is read from the vignette itself, so the
# stored numbers come from exactly the code the vignette displays.
#
# Each simulation chunk starts with its own set.seed(), so the six chunks are
# independent. They run in parallel (one chunk per process, via
# parallel::mclapply); within a chunk the replicates run serially, exactly as
# in the vignette, so the results do not depend on the number of cores.
#
# Usage, from the package root:
#   Rscript vignettes/precompute/parameter-coverage.R [--cores=4]
#       [--out=vignettes/parameter-coverage-results.rds] [--quick]
#
# --cores  number of processes (default 4, at most 6; forced to 1 on Windows).
# --out    output file.
# --quick  smoke test: 2 replicates per scenario (default output
#          parameter-coverage-quick.rds in the working directory).
#
# The package is loaded from the source tree with pkgload::load_all().
# The full run takes roughly an hour of CPU time.

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[[length(hit)]]) else default
}
quick <- "--quick" %in% args
cores <- as.integer(get_arg("cores", "4"))
out <- get_arg("out", if (quick) "parameter-coverage-quick.rds" else
  file.path("vignettes", "parameter-coverage-results.rds"))

rmd <- file.path("vignettes", "parameter-coverage.Rmd")
if (!file.exists(rmd) || !file.exists("DESCRIPTION"))
  stop("Run this script from the egpd package root.", call. = FALSE)
pkgload::load_all(".", quiet = TRUE)

read_chunks <- function(path) {
  lines <- readLines(path, warn = FALSE)
  starts <- grep("^```\\{r[ ,]", lines)
  ends <- grep("^```[[:space:]]*$", lines)
  chunks <- list()
  for (s in starts) {
    e <- ends[ends > s][1L]
    label <- sub("^```\\{r[ ,][[:space:]]*([^ ,}]+).*$", "\\1", lines[s])
    chunks[[label]] <- if (e > s + 1L) lines[(s + 1L):(e - 1L)] else character()
  }
  chunks
}

helper_labels <- c("sim-helper", "sim-both-helper")
sim_labels <- c("vary-sigma", "vary-xi", "vary-kappa",
                "both-sigma", "both-xi", "both-kappa")
chunks <- read_chunks(rmd)
missing <- setdiff(c(helper_labels, sim_labels), names(chunks))
if (length(missing))
  stop("Chunks not found in ", rmd, ": ", paste(missing, collapse = ", "),
       call. = FALSE)

env <- new.env(parent = globalenv())
for (lab in helper_labels)
  eval(parse(text = chunks[[lab]], keep.source = FALSE), envir = env)
helper_names <- ls(env)
if (quick) {
  formals(env$sim_param_coverage)$n_rep <- 2
  formals(env$sim_param_coverage_both)$n_rep <- 2
}

is_call_to <- function(e, fn) is.call(e) && identical(e[[1L]], as.name(fn))
is_seed <- function(e) is_call_to(e, "set.seed")
is_result <- function(e) {
  (is_call_to(e, "<-") || is_call_to(e, "=")) && is.name(e[[2L]]) &&
    startsWith(as.character(e[[2L]]), "res_")
}

# Grid constants are evaluated here, in chunk order (later chunks reuse
# constants from earlier ones). They must not consume random numbers, so the
# per-chunk set.seed() calls below reproduce the serial vignette exactly.
set.seed(1)
rng_before <- .Random.seed
jobs <- lapply(sim_labels, function(lab) {
  exprs <- as.list(parse(text = chunks[[lab]], keep.source = FALSE))
  seed <- Filter(is_seed, exprs)
  res <- Filter(is_result, exprs)
  if (length(seed) != 1L || length(res) != 1L)
    stop("Chunk '", lab, "' must contain one set.seed() call and one res_* ",
         "assignment.", call. = FALSE)
  for (e in Filter(function(e) !is_seed(e) && !is_result(e), exprs))
    eval(e, envir = env)
  list(label = lab, seed = seed[[1L]], res = res[[1L]],
       name = as.character(res[[1L]][[2L]]))
})
if (!identical(.Random.seed, rng_before))
  stop("Grid constants consumed random numbers.", call. = FALSE)

run_job <- function(job) {
  t0 <- proc.time()[["elapsed"]]
  jenv <- new.env(parent = env)
  eval(job$seed, envir = jenv)
  eval(job$res, envir = jenv)
  list(value = get(job$name, envir = jenv),
       seconds = proc.time()[["elapsed"]] - t0)
}

if (.Platform$OS.type == "windows") cores <- 1L
cores <- max(1L, min(cores, length(jobs)))
message(sprintf("Running %d simulation chunks on %d process(es)%s ...",
                length(jobs), cores, if (quick) " (quick)" else ""))
t_start <- proc.time()[["elapsed"]]
results <- if (cores > 1L) {
  parallel::mclapply(jobs, run_job, mc.cores = cores, mc.preschedule = FALSE)
} else {
  lapply(jobs, run_job)
}
failed <- vapply(results, function(r) inherits(r, "try-error") || is.null(r),
                 logical(1))
if (any(failed))
  stop("Simulation chunk(s) failed: ",
       paste(vapply(jobs[failed], `[[`, "", "label"), collapse = ", "), "\n",
       paste(unlist(lapply(results[failed], as.character)), collapse = "\n"),
       call. = FALSE)

# Keep only what the vignette displays: the per-replicate estimates are
# dropped, and the relative-bias matrices behind the boxplots are kept to
# 4 significant digits (which halves the file size).
trim <- function(res) lapply(res, function(r) {
  r$estimates <- NULL
  if (!is.null(r$rel_bias)) r$rel_bias <- signif(r$rel_bias, 4)
  r
})
values <- setNames(lapply(results, function(r) trim(r$value)),
                   vapply(jobs, `[[`, "", "name"))
constants <- setdiff(ls(env), helper_names)

git_head <- tryCatch(
  suppressWarnings(system2("git", c("describe", "--always", "--dirty"),
                           stdout = TRUE, stderr = FALSE)),
  error = function(e) NA_character_)
meta <- list(
  egpd_version = as.character(utils::packageVersion("egpd")),
  git_commit = if (length(git_head) == 1L) git_head else NA_character_,
  R_version = R.version.string,
  date = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  quick = quick,
  cores = cores,
  chunk_seconds = setNames(vapply(results, `[[`, numeric(1), "seconds"),
                           sim_labels),
  wall_seconds = proc.time()[["elapsed"]] - t_start
)

saveRDS(c(mget(constants, envir = env), values, list(meta = meta)), out,
        compress = "xz")
message(sprintf("Saved %s (%.1f KB) in %.1f min.", out,
                file.size(out) / 1024, meta$wall_seconds / 60))
