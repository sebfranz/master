# Install package from CRAN only if not installed, and load the library
if (!require(covr)) install.packages('covr')
if (!require(DT)) install.packages('DT')
if (!require(htmltools)) install.packages('htmltools')
library(covr)

covr <- file_coverage(source_files=c("./simulate_data/functions/generate_dummy_data_for_scregclust.R"),
                      test_files=c("./simulate_data/functions/test_generate_dummy_data_for_scregclust.R"),
                      parent_env=parent.frame() )
covr
report(covr)
