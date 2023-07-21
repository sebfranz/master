# Install package from CRAN only if not installed, and load the library
if (!require(covr)) install.packages('covr')
library(covr)
source("./simulate_data/functions/generate_dummy_data_for_scregclust.R")

covr <- file_coverage(source_files=c("./simulate_data/functions/generate_dummy_data_for_scregclust.R"),
                      test_files=c("./simulate_data/functions/test-generate_dummy_data_for_scregclust.R"),
                      parent_env=parent.frame() )
covr
report(covr)
