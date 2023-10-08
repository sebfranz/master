library(tidyverse)

execution_path <- dirname(rstudioapi::getActiveDocumentContext()$path)

rmarkdown::render(paste0(execution_path,"/simplified_likelihood_biclust.Rmd"))


