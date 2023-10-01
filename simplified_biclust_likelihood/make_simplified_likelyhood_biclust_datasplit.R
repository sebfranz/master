if (!require(here)) install.packages('here')
library(here)  # To work with paths
path_root <- here::here()

rmarkdown::render(file.path(path_root, "simplified_biclust_likelihood","simplified_likelyhood_biclust_datasplit.Rmd"))


