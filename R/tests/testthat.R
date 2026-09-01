# Prevent thread overuse (through data.table?) when running tests on CRAN.
Sys.setenv("OMP_THREAD_LIMIT" = 2L)

library(testthat)
library(power.transform)

testthat::test_check("power.transform")
