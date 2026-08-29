## Project .Rprofile
##
## Turns on the persistent odin2/dust2 compile cache for any R session started
## from the project root (Rscript, R, RStudio). Without it a fresh session
## recompiles the six dust2 models from scratch (~3 min); with it, ~6 s.
## See odin_cache.R.
##
## Opt out with `R --vanilla`, or by setting EATE_ODIN_CACHE_DISABLE=1.
if (!nzchar(Sys.getenv("EATE_ODIN_CACHE_DISABLE")) && file.exists("odin_cache.R")) {
  try(source("odin_cache.R"), silent = TRUE)
}
