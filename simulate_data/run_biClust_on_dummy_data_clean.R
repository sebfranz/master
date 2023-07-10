rm(list = ls())
execution_path <- dirname(rstudioapi::getSourceEditorContext()$path)
source(paste0(execution_path,"/dummy_data.R"))
library(scregclust)
library(plyr)
set.seed(1)  # This seed crashes due to scregclust producing NULL in all target gene clusters in the only cell cluster left

# Set variables ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

n_cell_clusters <- 3
n_target_gene_clusters <- c(3,4,5)  # Number of target gene clusters in each cell cluster
n_target_genes <- 20
n_regulator_genes <- 15
n_cells <- c(1000,5000,10000)
regulator_means = c(51,7435,1)  # For generating dummy data, regulator mean in each cell cluster
coefficient_means = list(c(10,20,30), c(1,200,300,400), c(1000,213,313,144,1245))  # For generating dummy data, coefficient means in each cell cluster
true_cluster_allocation = rep(1:n_cell_clusters, times=n_cells)
total_n_cells = sum(n_cells)

# Generate dummy data for each cell cluster that we want ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

dummy_data <- vector(mode = "list", length = n_cell_clusters)
for(i_cluster in 1:n_cell_clusters){
  print(i_cluster)
  dummy_data[[i_cluster]] <- generate_dummy_data(n_target_genes,
                                                 n_regulator_genes,
                                                 n_cells = n_cells[i_cluster],
                                                 n_target_gene_clusters[i_cluster],
                                                 regulator_mean = regulator_means[i_cluster],
                                                 coefficient_mean = coefficient_means[[i_cluster]])
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

previous_cell_clust <- initial_cell_clust


for (i_main in 1:50){
  print(paste("Iteration", i_main))
# Run scregclust for each cell cluster ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


  # Preallocate outputs
  out_list <- vector(mode = "list", length = n_cell_clusters)

  for(i_cluster in 1:n_cell_clusters){
    # Get data for cell cluster i_cluster specifically
    local_dat <- train_dat[,which(previous_cell_clust == i_cluster)]

    # Training data are represented by 1 and test data by 2
    cell_data_split    <- sample(c(1,2),ncol(local_dat), replace = T)

    # Initial target gene clustering
    gene_cluster_start <- kmeans(local_dat[1:n_target_genes,], n_target_gene_clusters[i_cluster])$cluster

    # Run scregclust
    scregclust(
      expression             = local_dat,               # scRegClust wants this form
      split_indices          = cell_data_split,         # Train data split
      genesymbols            = paste0('g', 1:(n_target_genes+n_regulator_genes)),  #gene row names
      is_regulator           = (1:(n_target_genes+n_regulator_genes) > n_target_genes) + 0,    #vector indicating which genes are regulators
      n_cl                   = n_target_gene_clusters[i_cluster],
      target_cluster_start   = gene_cluster_start,
      penalization           = 0.14,  # Maximal number of regulators for one cluster
      verbose                = FALSE
    ) -> out_list[[i_cluster]]
  }



# Calculate MSE -----------------------------------------------------------

  # Calculate MSE for per cells per cell cluster,
  # for all target gene clusters in those cell clusters
  # Missing clusters are rows with NA
  # matrix of cell cluster 1's target genes, then cell cluster 2's target genes - etc. In order. Columns are cells matching previous_cell_clust
  # Example output
  # > str(MSE)
  # num [1:12, 1:16000] 68.9 124.1 11.9 66.6 119.9 ...


  MSE <- matrix(nrow=sum(n_target_gene_clusters), ncol=total_n_cells)
  for(i_cell_cluster in 1:n_cell_clusters){
    # Expression of all cells
    # xvals expression of regulator cells
    # yvals expression of target cells
    xvals <- train_dat[(1:(n_target_genes+n_regulator_genes) > n_target_genes), which(previous_cell_clust == i_cell_cluster)]
    yvals <-  train_dat[(1:(n_target_genes+n_regulator_genes) <= n_target_genes), which(previous_cell_clust == i_cell_cluster)]
    i_total_target_geneclusters <- 0
    for(ii_cell_cluster in 1:n_cell_clusters){
      clustering <- out_list[[ii_cell_cluster]]$results[[1]]$output[[1]]$cluster[1:n_target_genes]
      for(i_target_gene_cluster in 1:n_target_gene_clusters[[ii_cell_cluster]]){
        i_total_target_geneclusters <- i_total_target_geneclusters + 1
        target_gene_ids_in_cluster_i <- which(clustering==i_target_gene_cluster)
        betas_for_gene_cluster_i <- out_list[[ii_cell_cluster]]$results[[1]]$output[[1]]$coeffs[[i_target_gene_cluster]]
        # If no target gene was assigned to this cluster we need to do something else
        if(!is.null(betas_for_gene_cluster_i)) {
          MSE[i_total_target_geneclusters, previous_cell_clust == i_cell_cluster] <- colMeans((yvals[target_gene_ids_in_cluster_i,] - t(betas_for_gene_cluster_i) %*% xvals)**2)
        }
      }
    }
  }

# Update cluster allocation to the appropriate cell cluster -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  # This can be done in some different ways,
  # One way would be to take mean squared error for the totality of the fitted model for that cell cluster
  # Another way would be to compare the average mse per gene cluster model
  # here we will compare the minimal mse per gene cluster model per cell cluster.
  # Could also use other metric than mse, e.g. r2

  if (all(is.na(MSE))){
    stop("scregclust put everything in noise cluster for all cellclusters. Exiting.")
  }
  updated_cell_clust <- rep(1:n_cell_clusters, n_target_gene_clusters)[apply(MSE, 2, which.min)]

  # Cross tabulation of clusters
  print("Table")
  print(data.frame(table(updated_cell_clust, previous_cell_clust)))
  flush.console()

  if (all(previous_cell_clust == updated_cell_clust)){
    print("Cell clustering same as last iteration. Exiting.")
    break
  }

  # Run next iteration with updated cell clustering
  previous_cell_clust <- updated_cell_clust

  # Handle case where cell cluster disappeared
  unique_cell_cluster_ids <- sort(unique(previous_cell_clust))
  if (length(unique_cell_cluster_ids)<n_cell_clusters){
    n_cell_clusters = length(unique_cell_cluster_ids)
    previous_cell_clust <- mapvalues(previous_cell_clust, from = unique_cell_cluster_ids, to = 1:n_cell_clusters)
    n_target_gene_clusters <- n_target_gene_clusters[unique_cell_cluster_ids]
  }

}
