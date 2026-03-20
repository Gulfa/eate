# Run from the project root:
#   Rscript tests/run_tests.R

library(testthat)

test_dir(file.path(dirname(normalizePath("tests/run_tests.R")), "testthat"),
         reporter = "progress")
