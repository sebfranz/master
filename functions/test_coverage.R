# Install package from CRAN only if not installed, and load the library
if (!require(covr)) install.packages('covr')
if (!require(DT)) install.packages('DT')
if (!require(htmltools)) install.packages('htmltools')
library(covr)

covr <- file_coverage(source_files=c("./functions/biclust.R"),
                      test_files=c("./functions/test_biclust.R"),
                      parent_env=parent.frame() )
covr
report(covr)
