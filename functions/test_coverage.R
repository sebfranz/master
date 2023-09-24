# Install package from CRAN only if not installed, and load the library
if (!require(covr)) install.packages('covr')
if (!require(DT)) install.packages('DT')
if (!require(htmltools)) install.packages('htmltools')
if (!require(here)) install.packages('here')
library(covr)
library(here)
file_biclust <- here::here("functions", "biclust.R")
file_test_biclust <- here::here("functions", "test_biclust.R")

covr <- file_coverage(source_files=c(file_biclust),
                      test_files=c(file_test_biclust),
                      parent_env=parent.frame() )
covr
report(covr)
