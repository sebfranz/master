rm(list = ls())
execution_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(execution_path)
# neftel_path <- paste0(execution_path, '/../../datasets_sctargettranslator/Neftel2019')
# setwd(neftel_path)

library(scregclust)
library('Seurat')
library(Matrix)

#loads mn_g1
#cells are rows genes are columns
source(paste0(execution_path,"/../sc_data/setup.R"))
