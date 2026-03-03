# Run all DE-LIMP tests
# Usage: Rscript tests/testthat.R
#   or:  testthat::test_dir("tests/testthat") from the project root

library(testthat)
test_dir("tests/testthat")
