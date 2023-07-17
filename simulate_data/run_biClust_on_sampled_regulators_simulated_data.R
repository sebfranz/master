rm(list = ls())
library(scregclust)
library(plyr)
library(aricode)  # To calculate rand index
library(Seurat)

execution_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
sapply(list.files(paste0(execution_path,"/../functions/"),recursive = T),
       function(nm) source(paste0(execution_path,"/../functions/", nm)))


source(paste0(execution_path,"/functions/generate_scregclust_data_from_sampled_regulators.R"))
source(paste0(execution_path,"/functions/generate_biclust_data_from_sampled_regulators.R"))

set.seed(124)  # This seed crashes due to scregclust producing NULL in all target gene clusters in the only cell cluster left

# Set variables ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

n_target_genes = 100
n_cell_clusters = 2
n_target_gene_clusters = c(3,4)  # Number of target gene clusters in each cell cluster
regulator_expression_offset =  c(0,0)
n_cells = c(1000,1000)
coefficient_means = list(c(1,2,3), c(10,20,30,40))  # For generating dummy data, coefficient means in each cell clustertrue_cluster_allocation = rep(1:n_cell_clusters, times=n_cells)
# total_n_cells = sum(n_cells)

# Generate dummy data for each cell cluster that we want ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#load sampled data
execution_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
neftel_path <- paste0(execution_path, '/../../datasets_sctargettranslator/Neftel2019')
Neftel_g1 <- readRDS(file = paste0(neftel_path, '/r_files/', "neftel_seurat_group1"))

Neftel_g1<-SetIdent(Neftel_g1, value='malignant')
Neftel_g1_malignant<-subset(Neftel_g1, idents='yes')

z_g1<-GetAssayData(Neftel_g1_malignant, slot='scale.data')
dim(z_g1)
metaData<-Neftel_g1_malignant@meta.data$sample

out<-scregclust_format(z_g1)

genesymbols<-out[[1]]
sample_assignment<-out[[2]]
is_predictor<-out[[3]]

# z_g1 <- z_g1[order(is_predictor),]
# is_predictor <- is_predictor[order(is_predictor)]
regulator_expression <- t(z_g1[which(is_predictor == 1), ])
rm(z_g1)
rm(Neftel_g1)
rm(Neftel_g1_malignant)

res <- generate_biclust_data_from_sampled_regulators(
  n_target_genes = n_target_genes,
  n_cell_clusters = n_cell_clusters,
  n_target_gene_clusters = n_target_gene_clusters,  # Number of target gene clusters in each cell cluster
  regulator_expression_offset =  regulator_expression_offset,
  n_cells =n_cells,
  coefficient_means = coefficient_means,  # For generating dummy data, coefficient means in each cell cluster
  regulator_expression = regulator_expression
)


cell_cluster_history <- cbind(res$initial_cell_clust, res$disturbed_initial_cell_clust)
colnames(cell_cluster_history) <- c("True allocation", "Disturbed allocation")

res <- biclust(max_iter=50,
               initial_cluster_history = cell_cluster_history,
               is_regulator = c(rep(0, n_target_genes),rep(1,ncol(regulator_expression))),
               n_target_gene_clusters = n_target_gene_clusters,
               n_cells = n_cells,
               train_dat = res$train_dat)

plot_cluster_history(cell_cluster_history = res$cell_cluster_history)
