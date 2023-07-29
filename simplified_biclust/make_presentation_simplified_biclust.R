library(tidyverse)

execution_path <- dirname(rstudioapi::getActiveDocumentContext()$path)

rmarkdown::render(paste0(execution_path,"/Simplified_biclust_presentation.Rmd"))


