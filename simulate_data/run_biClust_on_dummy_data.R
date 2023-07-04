rm(list = ls())
execution_path <- dirname(rstudioapi::getSourceEditorContext()$path)
source(paste0(execution_path,"/dummy_data.R"))
library(scregclust)
# set.seed(12411)

n_cell_clusters <- 3
n_target_gene_clusters <- c(3,4,5)  # Number of target gene clusters in each cell cluster
n_target_genes <- 20
n_regulator_genes <- 15
n_cells <- c(1000,5000,10000)
regulator_means = c(1,2,3)  # For generating dummy data, regulator mean in each cell cluster
coefficient_means = list(c(1,2,3), c(1,2,3,4), c(1,2,3,4,5))  # For generating dummy data, coefficient means in each cell cluster
true_cluster_allocation = rep(1:n_cell_clusters, times=n_cells)
total_n_cells = sum(n_cells)
############## preprocessing

# train/validation split of CELLS
# initialize with kmeans for all cells
# set train val split within cell clusters so cells start in same cell cluster

############### Outer loop

# check cell clustering fit
# adjusted r2?
# Check how well cells fit into any cell clusters, possibly update cell cluster membership

############### inner loop

#run scRegClust with proper initializations for each cell cluster.

#maybe to start

# 1 randomly split cells into train&val
# 2 for train cells, do a clustering thingy
# 3 for each cluster run scregclust
# 4 evaluate etc

# Generate dummy data for each cell cluster that we want
dummy_data <- vector(mode = "list", length = n_cell_clusters)
for(i_cluster in 1:n_cell_clusters){
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
    penalization           = n_regulator_genes  # Maximal number of regulators for one cluster
  ) -> out_list[[i_cluster]]
}

# Which target genes belong to which cluster>
str(out_list[[1]]$results[[1]])
str(out_list[[1]]$results[[1]]$output[[1]])



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
str(MSE)


# # Update the actual number of target gene clusters in each cell cluster
# n_target_gene_clusters_updated <- n_target_gene_clusters
# for(i_cell_cluster in 1:n_cell_clusters){
#   # This checks if each target gene clusters coeffs are numeric, if not, they are NULL.
#   # Then we just count those numeric ones.
#   n_target_gene_clusters_updated[i_cluster] <- sum(sapply(out_list[[2]]$results[[1]]$output[[1]]$coeffs, is.numeric))
# }




# Old MSE ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# MSE <- vector(mode = "list", length = n_cell_clusters)
# names(MSE) <- paste0('Cell cluster ', 1:n_cell_clusters)
# for(i_cell_cluster in 1:n_cell_clusters){
#
#   # Expression of all cells
#   # xvals expression of regulator cells
#   # yvals expression of target cells
#   xvals <- train_dat[(1:(n_target_genes+n_regulator_genes) > n_target_genes), which(previous_cell_clust == i_cell_cluster)]
#   yvals <-  train_dat[(1:(n_target_genes+n_regulator_genes) <= n_target_genes), which(previous_cell_clust == i_cell_cluster)]
#
#   MSE_in_cell_cluster_i <- vector(mode = "list", length = n_cell_clusters)
#   for(ii_cell_cluster in 1:n_cell_clusters){
#     clustering <- out_list[[ii_cell_cluster]]$results[[1]]$output[[1]]$cluster[1:n_target_genes]
#     MSE_in_cell_cluster_ii <- vector(mode = "list", length = n_target_gene_clusters[[ii_cell_cluster]])
#
#     for(i_target_gene_cluster in 1:n_target_gene_clusters[[ii_cell_cluster]]){
#       target_gene_ids_in_cluster_i <- which(clustering==i_target_gene_cluster)
#
#       betas_for_gene_cluster_i <- out_list[[ii_cell_cluster]]$results[[1]]$output[[1]]$coeffs[[i_target_gene_cluster]]
#       # If no target gene was assigned to this cluster we need to do something else
#       if(!is.null(betas_for_gene_cluster_i)) {
#         MSE_in_cell_cluster_ii[[i_target_gene_cluster]] <- colMeans((yvals[target_gene_ids_in_cluster_i,] - t(betas_for_gene_cluster_i) %*% xvals)**2)
#       }else{
#         print(paste0("Cluster doesn't exist ", i_target_gene_cluster, " in cell cluster ", ii_cell_cluster))
#       }
#     }
#     names(MSE_in_cell_cluster_ii) <- paste0('... MSE for target gene cluster ', 1:n_target_gene_clusters[[ii_cell_cluster]])
#     MSE_in_cell_cluster_i[[ii_cell_cluster]] <- MSE_in_cell_cluster_ii
#   }
#   names(MSE_in_cell_cluster_i) <- paste0('... MSE for cell cluster ', 1:n_cell_clusters)
#   # Store results in a list of size n_cell_clusters(cell clusters)
#   # of lists of n_target_gene_clusters(gene clusters)
#   MSE[[i_cell_cluster]] <- MSE_in_cell_cluster_i
# }
# str(MSE)



# Add each cell cluster after eachother

# Update cluster allocation to the appropriate cell cluster
# This can be done in some different ways,
# One way would be to take mean squared error for the totality of the fitted model for that cell cluster
# Another way would be to compare the average mse per gene cluster model
# here we will compare the minimal mse per gene cluster model per cell cluster.
# Could also use other metric than mse, e.g. r2

# updated_cell_clust <- previous_cell_clust + NA
# metrics <- matrix(data = NA,  nrow = n_cell_clusters, ncol = ncol(dat))
#
# for(cell in 1:ncol(dat)){
#   which_cell_cluster <- previous_cell_clust[cell]#should update to "current" cell cluster
#   which_within_cluster_index <- sum(previous_cell_clust[1:(cell)] == which_cell_cluster)
#
#   # for(cell_cluster in 1:n_cell_clusters){
#   #   #find cluster allocation metrics
#   #   metrics[cell_cluster,cell] <- min(do.call(rbind, MSE[[which_cell_cluster]][[cell_cluster]])[,which_within_cluster_index], na.rm = T)
#   # }
#
#   # updated_cell_clust[cell] <- which.min(metrics[,cell])
#   which.min(
#     sapply(
#       1:n_cell_clusters, function(j) min(sapply(1: n_cell_clusters, function(i) MSE[[which_cell_cluster]][[j]][[i]][which_within_cluster_index]), na.rm = T)
#     )
#   ) -> updated_cell_clust[cell]
#
# }
# updated_cell_clust <- sapply(1:ncol(dat), function(i) which.min(metrics[,i]))
# sapply(
#   1:ncol(dat),
#   function(cell) {
#     which_cell_cluster <- previous_cell_clust[cell]
#     which_within_cluster_index <- sum(previous_cell_clust[1:(cell)] == previous_cell_clust[cell])
#     which.min(
#       sapply(
#         1:n_cell_clusters, function(j) min(sapply(1: n_cell_clusters, function(i) MSE[[which_cell_cluster]][[j]][[i]][which_within_cluster_index]), na.rm = T)
#       )
#     )
#   }
# ) -> updated_cell_clust
#
# # cross tabulation of clusters
# data.frame( table(updated_cell_clust, previous_cell_clust))


for(i_cell_cluster in 1:n_cell_clusters){
  updated_cell_clust <- rep(1:n_cell_clusters, n_target_gene_clusters)[apply(MSE, 2, which.min)]
}

# sapply(
#   1:ncol(dat),
#   function(cell) {
#     which_cell_cluster <- previous_cell_clust[cell]
#     which_within_cluster_index <- sum(previous_cell_clust[1:(cell)] == previous_cell_clust[cell])
#     which.min(
#       sapply(
#         1:n_cell_clusters, function(j) min(sapply(1: n_cell_clusters, function(i) MSE[[which_cell_cluster]][[j]][[i]][which_within_cluster_index]), na.rm = T)
#       )
#     )
#   }
# ) -> updated_cell_clust

# Cross tabulation of clusters
data.frame(table(updated_cell_clust, previous_cell_clust))

# Notes -------------------------------------------------------------------





#now we arrive at the question, how do we update cell clusters given all this new information about the gene cluster structure within each previous cell cluster.

# Reb talked about using R2? what would that mean?
#if we can extract the regression model from each

# subproblem;
# for a cell, compute the MSE, and r2 of that one cell
# do this for ALL regression models in ALL cell clusters.
# use this to find which cell cluster a given cell fits best into
# sub-sub-problem;
# each cell cluster contains several regression models, one for each target
# gene cluster



min(do.call(rbind, MSE[[which_cell_cluster]][[cell_cluster]][,which_within_cluster_index]), na.rm = T)

