# odin_cache.R --------------------------------------------------------------
#
# Persistent compile cache for the odin2/dust2 models.
#
# odin2::odin() builds each model into a throw-away package under tempdir(),
# so every fresh R session pays the full C++ compile. For the six dust2
# models in this project that is ~3 minutes before any analysis starts.
#
# dust2 honours the DUST_WORKDIR_ROOT environment variable: with it set, a
# model is generated into <root>/dust_<hash> instead of a temporary
# directory. dust2 only rewrites the generated source when it has actually
# changed, so on the next session `make` sees an up-to-date .so and skips the
# compiler entirely. Cold build ~171s, warm start ~6s.
#
# The hash is taken over the *generated* C++, so editing a model file lands
# in a new cache entry and genuinely recompiles. The cache never needs to be
# invalidated by hand; odin_cache_clean() only exists to reclaim disk from
# entries belonging to older versions of a model.
#
# Sourced from .Rprofile (so plain `Rscript foo.R` from the project root gets
# the cache without any code change) and from model.R / det_model.R /
# stoch_model.R (so it still applies under `Rscript --vanilla`).
#
# Not cached: det_mod_cd.R and det_mod_ncd.R go through odin v1, which
# rewrites its generated C unconditionally and so always recompiles. That is
# ~5s total, against ~171s for the dust2 models.

odin_cache_models <- c("stoch_sir.R", "stoch_mod_adj.R", "stoch_mod_adj_sparse.R",
                       "stoch_linear.R", "stoch_ind.R", "det_mod_adj.R")


# Walk up from `start` to the directory holding the repo, so the cache lands
# in the same place whether a script was sourced with chdir=TRUE or not.
odin_project_root <- function(start = getwd()) {
  d <- normalizePath(start, mustWork = FALSE)
  repeat {
    if (file.exists(file.path(d, "CLAUDE.md")) || file.exists(file.path(d, ".git"))) {
      return(d)
    }
    parent <- dirname(d)
    if (identical(parent, d)) return(normalizePath(start, mustWork = FALSE))
    d <- parent
  }
}


odin_cache_root <- function() {
  root <- Sys.getenv("EATE_ODIN_CACHE", "")
  if (nzchar(root)) root else file.path(odin_project_root(), ".odin_cache")
}


# Point dust2 at the persistent cache. Safe to call repeatedly. The directory
# is not created here — dust2 and odin_cache_lock() create what they need, so
# merely starting R leaves no trace.
odin_cache_enable <- function(root = odin_cache_root()) {
  Sys.setenv(DUST_WORKDIR_ROOT = root)
  invisible(root)
}


odin_cache_disable <- function() {
  Sys.unsetenv("DUST_WORKDIR_ROOT")
  invisible(NULL)
}


# The slurm arrays start ~20 R processes at once. On a cold cache they would
# otherwise all generate and compile into the same <root>/dust_<hash>
# directory concurrently and clobber each other's object files. Hold a lock
# per model file: the first process compiles, the rest wait and then hit the
# up-to-date .so.
odin_cache_lock <- function(root, file) {
  if (!requireNamespace("filelock", quietly = TRUE)) {
    warning("Package 'filelock' is not installed: cannot serialise odin builds, ",
            "so this process will compile into a PRIVATE workdir (no cache ",
            "reuse). Install filelock to get the shared cache.", call. = FALSE)
    return(NULL)
  }
  dir <- file.path(root, "locks")
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  key <- sprintf("%s-%s.lock", tools::file_path_sans_ext(basename(file)),
                 substr(tools::md5sum(file), 1, 8))
  lock <- filelock::lock(file.path(dir, key), timeout = 15 * 60 * 1000)
  if (is.null(lock)) {
    warning(sprintf("Timed out waiting to build '%s'; compiling in a private workdir",
                    file), call. = FALSE)
  }
  lock
}


odin_cache_unlock <- function(lock) {
  if (!is.null(lock)) filelock::unlock(lock)
  invisible(NULL)
}


# A cache entry records nothing about where it came from: the directory is
# named after a hash and the generator is always called "odin_system". Keep a
# side index, one file per entry, so odin_cache_status() can say which model
# each entry belongs to. It lives outside the entry directory because dust2
# rejects a workdir containing files it did not write.
odin_cache_index <- function(root, entry) file.path(root, "index", entry)


odin_cache_note <- function(root, file, gen) {
  path <- attr(gen, "path", exact = TRUE)
  if (is.null(path)) return(invisible(NULL))
  dir.create(file.path(root, "index"), showWarnings = FALSE, recursive = TRUE)
  writeLines(basename(file), odin_cache_index(root, basename(path)))
  invisible(NULL)
}


# Drop-in replacement for odin2::odin(file, ...).
odin_cached <- function(file, ..., quiet = TRUE) {
  root <- odin_cache_enable()
  lock <- odin_cache_lock(root, file)
  if (is.null(lock)) {
    # No lock (filelock missing, or the wait timed out). Compiling into the
    # SHARED cache directory would then let concurrent processes clobber each
    # other's object files mid-build -- seen on the cluster as
    #   ld: cannot find cpp11.o: No such file or directory
    # after both cpp11.o and dust.o had apparently just been compiled. Fall
    # back to a per-process workdir: no cache reuse (slower), but correct.
    private <- file.path(tempdir(), "odin_cache_private")
    dir.create(private, showWarnings = FALSE, recursive = TRUE)
    Sys.setenv(DUST_WORKDIR_ROOT = private)
    on.exit(Sys.setenv(DUST_WORKDIR_ROOT = root), add = TRUE)
    return(odin2::odin(file, ..., quiet = quiet))
  }
  on.exit(odin_cache_unlock(lock), add = TRUE)
  gen <- odin2::odin(file, ..., quiet = quiet)
  odin_cache_note(root, file, gen)
  gen
}


# Compile every model into the cache. Run once after checking out or after
# editing a model file, so an interactive session or a slurm array starts
# warm: Rscript -e 'source("odin_cache.R"); odin_cache_warm()'
odin_cache_warm <- function(files = odin_cache_models) {
  for (f in files) {
    elapsed <- system.time(odin_cached(f))[["elapsed"]]
    message(sprintf("%-24s %6.2fs", f, elapsed))
  }
  invisible(odin_cache_root())
}


odin_cache_entries <- function(root = odin_cache_root()) {
  dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  dirs[grepl("^dust_", basename(dirs))]
}


odin_cache_status <- function(root = odin_cache_root()) {
  dirs <- odin_cache_entries(root)
  if (length(dirs) == 0) {
    return(data.frame(entry = character(), model = character(),
                      built = as.POSIXct(character()), mb = numeric()))
  }
  info <- lapply(dirs, function(d) {
    so <- list.files(file.path(d, "src"), pattern = "\\.so$", full.names = TRUE)
    idx <- odin_cache_index(root, basename(d))
    name <- if (file.exists(idx)) readLines(idx) else character()
    all <- list.files(d, recursive = TRUE, full.names = TRUE)
    data.frame(entry = basename(d),
               model = if (length(name)) name[[1]] else NA_character_,
               built = if (length(so)) file.mtime(so)[[1]] else as.POSIXct(NA),
               mb = round(sum(file.size(all), na.rm = TRUE) / 1e6, 1))
  })
  res <- do.call(rbind, info)
  res[order(res$built, decreasing = TRUE), ]
}


# Remove entries whose shared library was last built more than `days` ago.
# A cache hit does not touch the .so, so a model you still use will age out
# eventually; that only costs one recompile, after which it is warm again.
# odin_cache_clean(0) clears the cache completely.
odin_cache_clean <- function(days = 30, root = odin_cache_root()) {
  dirs <- odin_cache_entries(root)
  if (length(dirs) == 0) return(invisible(character()))
  age <- vapply(dirs, function(d) {
    so <- list.files(file.path(d, "src"), pattern = "\\.so$", full.names = TRUE)
    if (length(so) == 0) return(Inf)
    as.numeric(difftime(Sys.time(), file.mtime(so)[[1]], units = "days"))
  }, numeric(1))
  stale <- dirs[age > days]
  unlink(stale, recursive = TRUE)
  unlink(odin_cache_index(root, basename(stale)))
  message(sprintf("Removed %d of %d cache entries", length(stale), length(dirs)))
  invisible(basename(stale))
}


odin_cache_enable()
