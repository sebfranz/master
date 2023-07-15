rm(list = ls())
execution_path <- dirname(rstudioapi::getSourceEditorContext()$path)
# neftel_path <- paste0(execution_path, '/../../datasets_sctargettranslator/Neftel2019')
# setwd(neftel_path)

library(scregclust)
library('Seurat')
library(Matrix)

source(paste0(execution_path,"/setup.R"))


execution_path <- dirname(rstudioapi::getSourceEditorContext()$path)
source(paste0(execution_path,"/../functions/biclust.R"))
source(paste0(execution_path,"/../functions/plot_cluster_history.R"))

res <- biclust(max_iter=50,
               initial_cluster_history = cell_cluster_history,
               n_target_genes = sum(res1$is_regulator==0),
               n_regulator_genes = sum(res1$is_regulator==1),
               n_target_gene_clusters = c(1,2,4),
               n_cells = c(cells_per_cluster,cells_per_cluster,cells_per_cluster),
               train_dat = as.matrix(res))

plot_cluster_history(cell_cluster_history = res$cell_cluster_history)



