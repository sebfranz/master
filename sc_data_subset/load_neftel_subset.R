rm(list = ls())
execution_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(execution_path)
# neftel_path <- paste0(execution_path, '/../../datasets_sctargettranslator/Neftel2019')
# setwd(neftel_path)

library(scregclust)
library('Seurat')
library(Matrix)

source(paste0(execution_path,"/../sc_data/setup.R"))

#cells are rows genes are columns
dim(mn_g1)
subsample <- mn_g1[sample(1:dim(mn_g1)[1], 1000), sample(1:dim(mn_g1)[2], 500)]

dim(subsample)
heatmap(subsample)


