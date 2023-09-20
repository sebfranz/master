library(tidyverse)

execution_path <- dirname(rstudioapi::getActiveDocumentContext()$path)

rmarkdown::render(paste0(execution_path,"/simplified_likelyhood_biclust.Rmd"))


