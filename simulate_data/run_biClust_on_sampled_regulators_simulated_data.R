rm(list = ls())
library(scregclust)
library(plyr)
library(aricode)  # To calculate rand index

execution_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
print(execution_path)
source(paste0(execution_path,"/functions/generate_scregclust_data_from_sampled_regulators.R"))
source(paste0(execution_path,"/functions/generate_biclust_data_from_sampled_regulators.R"))
source(paste0(execution_path,"/../functions/biclust.R"))
source(paste0(execution_path,"/../functions/plot_cluster_history.R"))

set.seed(1234)  # This seed crashes due to scregclust producing NULL in all target gene clusters in the only cell cluster left

# Set variables ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

n_target_genes = 100
n_cell_clusters = 3
n_target_gene_clusters = c(3,4,5)  # Number of target gene clusters in each cell cluster
regulator_expression_offset =  c(0,10,100)
n_cells = c(1000,5000,10000)
coefficient_means = list(c(1,20,30), c(1,2,3,4), c(1,2,3,4,5))  # For generating dummy data, coefficient means in each cell clustertrue_cluster_allocation = rep(1:n_cell_clusters, times=n_cells)
total_n_cells = sum(n_cells)

# Generate dummy data for each cell cluster that we want ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
res <- generate_biclust_data_from_sampled_regulators(
  n_target_genes = 100,
  n_cell_clusters = 3,
  n_target_gene_clusters = c(3,4,5),  # Number of target gene clusters in each cell cluster
  regulator_expression_offset =  c(0,10,100),
  n_cells = c(1000,5000,10000),
  coefficient_means = list(c(1,20,30), c(1,2,3,4), c(1,2,3,4,5))  # For generating dummy data, coefficient means in each cell cluster
)

cell_cluster_history <- cbind(res$initial_cell_clust, res$disturbed_initial_cell_clust)
colnames(cell_cluster_history) <- c("True allocation", "Disturbed allocation")

res <- biclust(max_iter=50,
               initial_cluster_history = cell_cluster_history,
               is_regulator = c(rep(0,n_target_genes),rep(1,n_regulator_genes)),
               n_target_gene_clusters = n_target_gene_clusters,
               n_cells = n_cells,
               train_dat = res$train_dat)

plot_cluster_history(cell_cluster_history = res$cell_cluster_history)
