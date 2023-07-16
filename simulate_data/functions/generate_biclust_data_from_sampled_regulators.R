execution_path <- dirname(rstudioapi::getActiveDocumentContext()$path)

library(plyr)

generate_biclust_data_from_sampled_regulators <- function(
    n_target_genes = 100,
    n_cell_clusters = 3,
    n_target_gene_clusters = c(3,4,5),  # Number of target gene clusters in each cell cluster
    regulator_expression_offset =  c(0,10,100),
    n_cells = c(1000,5000,10000),
    coefficient_means = list(c(1,20,30), c(1,2,3,4), c(1,2,3,4,5)),  # For generating dummy data, coefficient means in each cell cluster
    regulator_expression
){

  # Set variables ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  true_cluster_allocation = rep(1:n_cell_clusters, times=n_cells)
  total_n_cells = sum(n_cells)

  # Generate dummy data for each cell cluster that we want ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  dummy_data <- vector(mode = "list", length = n_cell_clusters)
  for(i_cluster in 1:n_cell_clusters){
    print(i_cluster)
    dummy_data[[i_cluster]] <- generate_scregclust_data_from_sampled_regulators(
      n_target_genes,
      n_cells = n_cells[i_cluster],
      regulator_expression_offset = regulator_expression_offset[i_cluster],
      n_target_gene_clusters = n_target_gene_clusters[i_cluster],
      coefficient_mean = coefficient_means[[i_cluster]],
      regulator_expression
      )
  }

  # Create Z_r and Z_t from dummy data
  Z_t <- dummy_data[[1]]$Z_t
  Z_r <- dummy_data[[1]]$Z_r
  if (n_cell_clusters>1){
    for(i_cluster in 2:n_cell_clusters){
      Z_t <- rbind(Z_t, dummy_data[[i_cluster]]$Z_t)
      Z_r <- rbind(Z_r, dummy_data[[i_cluster]]$Z_r)
    }
  }
  dat <- rbind(t(Z_t), t(Z_r)) # Columns are now cells

  # Split into train and test data for cell clustering
  # Skip for now
  # cell_data_split    <- sample(c(1,2), nrow(Z_t), replace = T)
  # train_indices      <- which(cell_data_split == 1)
  # train_dat          <- dat[,train_indices]

  train_dat <- dat

  # Get initial cell clustering
  initial_cell_clust <- kmeans(t(train_dat), n_cell_clusters)$cluster
  # initial_cell_clust <- sample(1:n_cell_clusters, n_cells, replace = T)


  # Kod för att flytta 1% av cellerna i varje kluster till ett annat kluster.
  disturbed_initial_cell_clust <- initial_cell_clust
  for(i_cluster in 1:n_cell_clusters){
    indexes_of_cluster <- which(initial_cell_clust == i_cluster)
    some_of_those_indexes <- sample(indexes_of_cluster, size=as.integer(length(indexes_of_cluster)*0.01), replace = F)
    disturbed_initial_cell_clust[some_of_those_indexes] <- sample(c(1:n_cell_clusters)[-i_cluster], size=length(some_of_those_indexes), replace=T)
  }

  return(list(disturbed_initial_cell_clust=disturbed_initial_cell_clust,
              initial_cell_clust=initial_cell_clust,
              train_dat=train_dat))
}
