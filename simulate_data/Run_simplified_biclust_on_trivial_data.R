rm(list = ls())
library(scregclust)
library(plyr)
library(aricode)  # To calculate rand index
library(Seurat)

execution_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
sapply(list.files(paste0(execution_path,"/../functions/"),recursive = T),
       function(nm) source(paste0(execution_path,"/../functions/", nm)))

set.seed(1234)

# Set variables ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#
# n_target_genes = 100
# n_cell_clusters = 2
# n_target_gene_clusters = c(3,4)  # Number of target gene clusters in each cell cluster
# regulator_expression_offset =  c(0,0)
# coefficient_means = list(c(1,1,1), c(1,1,1,1))  # For generating dummy data, coefficient means in each cell clustertrue_cluster_allocation = rep(1:n_cell_clusters, times=n_cells)
# total_n_cells = sum(n_cells)

# Generate dummy data for each cell cluster that we want ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#generate some trivial data for K "cell" clusters

n_cell_clusters    <- 2
num_cells <- 100
intercepts <- rnorm(n_cell_clusters, mean = 10, sd = 5)
betas      <- rnorm(n_cell_clusters, mean = 1,  sd = 0.5)
regulator_expression <- rnorm(num_cells, mean = 1, sd = 0.1)
true_cell_clust <- c(rep(1,num_cells/2), rep(2,num_cells/2))

#build cell expression from linear model
cell_cluster_expression <- sapply(1:n_cell_clusters, function(cellCluster) c(intercepts[cellCluster], betas[cellCluster]) %*%
                                        t(cbind(rep(1, num_cells/2), regulator_expression[which(true_cell_clust == cellCluster)]))
      )

target_expression <- c(cell_cluster_expression[,1], cell_cluster_expression[,2])
dat <- cbind(target_expression, regulator_expression)

# Kod för att flytta 1% av cellerna i varje kluster till ett annat kluster.
disturbed_initial_cell_clust <- true_cell_clust
disturbed_fraction <- 0.20
for(i_cluster in 1:n_cell_clusters){
  indexes_of_cluster <- which(true_cell_clust == i_cluster)
  some_of_those_indexes <- sample(indexes_of_cluster, size=as.integer(length(indexes_of_cluster)*disturbed_fraction), replace = F)
  disturbed_initial_cell_clust[some_of_those_indexes] <-
    sample(c(1:n_cell_clusters)[-i_cluster],
           size=length(some_of_those_indexes), replace=T)
}


cell_cluster_history <- cbind(true_cell_clust, disturbed_initial_cell_clust)
colnames(cell_cluster_history) <- c("True allocation", "Disturbed allocation")

#first variant simplified_biclust()
max_iter <- 50

#find initial cluster labels
initial_clustering <- disturbed_initial_cell_clust
n_target_gene_clusters <- 1 #we are not clustering target genes for now and we only have one target gene


#set up some variables
n_cell_clusters = length(unique(initial_clustering))
n_target_genes = 1
n_regulator_genes = 1


#set up cluster history
initial_column_padding <- ncol(as.matrix(initial_clustering)) + 1  # +1 Because we have an index column that is not an index column it's an ID column
cell_cluster_history <- data.frame(matrix(NA, nrow = nrow(as.matrix(initial_clustering)), ncol = max_iter + initial_column_padding))
colnames(cell_cluster_history) <- c("Cell ID", 'Initial clustering', paste0("Iteration ", 1:max_iter))
cell_cluster_history[, 'Cell ID'] <-1:length(initial_clustering)  # Set cell names
cell_cluster_history[, 'Initial clustering'] <- initial_clustering

#preallocate all r2 matrices for later analysis if feasible
r2_all <- vector("list", length = max_iter)
#set flag for breaking out of the loop.
stop_iterating_flag = T

current_cell_cluster_allocation <- initial_clustering

# i_main  <- 1

for(i_main in 1:max_iter){
  # fit model to each cell cluster
  models <- vector("list", length = n_cell_clusters)
  for(cell_cluster in 1:n_cell_clusters){
    models[[cell_cluster]] <- lm(log(dat[,1]) ~ dat[,-1])
  }
  # plot(models[[1]])

  #for all cells, calculate the predicted r2 for all cell clusters.
  r2 <- matrix(0, nrow = nrow(dat), ncol = n_cell_clusters)

  #first calculate the mean target gene expression in these clusters
  #also the total sum of squares in that cluster
  target_gene_means <- vector('numeric', length = n_cell_clusters )
  SS_tot <- vector('numeric', length = n_cell_clusters )
  for(cell_cluster in 1:n_cell_clusters){
    target_gene_means[cell_cluster] <- mean(dat[which(current_cell_cluster_allocation == cell_cluster),1])
    SS_tot[cell_cluster]  <- sum((dat[,1] - target_gene_means[cell_cluster])^2)
  }

  #now to actually caltulate predicted or 'predicted' r2
  for(cell in 1:nrow(dat)){
    for(cell_cluster in 1:n_cell_clusters){
      #bug fix hack: remove NA coefficients
      if(any(is.na(models[[cell_cluster]]$coefficients))){
        NA_coeffs <-  unname(which(is.na(models[[cell_cluster]]$coefficients)))
        S_ERR <- (dat[cell,1] - as.vector(c(1,dat[cell,c(-1, -NA_coeffs)])) %*% models[[cell_cluster]]$coefficients[-NA_coeffs])^2
      }
      S_ERR <- (dat[cell,1] - as.vector(c(1,dat[cell,-1])) %*% models[[cell_cluster]]$coefficients)^2
      r2[cell,cell_cluster] <- 1-(S_ERR/SS_tot[[cell_cluster]])
    }
  }

  r2_all[[i_main]] <- r2


  #update cluster allocations
  updated_cell_clust <-  sapply(1:nrow(r2), function(row) which.min(r2[row,]))

  # If there is only 1 cell left in a cell cluster it doesn't work
  # Move that cell to the biggest cell cluster.
  updated_cell_clust_table = data.frame(table(updated_cell_clust))

  for(i_cell_cluster in 1:n_cell_clusters){
    if (length(which(updated_cell_clust==i_cell_cluster)) == 1){
      updated_cell_clust[which(updated_cell_clust==i_cell_cluster)] = which.max(updated_cell_clust_table$Freq)
    }
  }

  # Update data in cell_cluster_history
  cell_cluster_history[, i_main + initial_column_padding] <- updated_cell_clust

  # Compare clusters with with previous iterations so we can exit if we seen this allocation before
  for(prev_clustering in ((i_main-1):0) ){
    print(paste0('Comparing with iteration ', prev_clustering))
    if(RI(updated_cell_clust,cell_cluster_history[,prev_clustering + initial_column_padding]) == 1){
      print("Cell clustering from iteration same as some previous iteration. Exiting.")
      print(paste0("RI of ",
                   RI(updated_cell_clust, cell_cluster_history[,prev_clustering + initial_column_padding]),
                   " when comparing iteration ", i_main," to iteration ", prev_clustering))
      stop_iterating_flag = T
      break
    }
  }

  if(stop_iterating_flag == T){
    # Clean up cluster history
    cell_cluster_history <- cell_cluster_history[ , colSums(is.na(cell_cluster_history))==0, drop=FALSE]
    # Stop iterations/exit function
    break
  }

  current_cell_cluster_allocation <- updated_cell_clust

  # Handle case where cell cluster disappeared
  unique_cell_cluster_ids <- sort(unique(updated_cell_clust))
  if (length(unique_cell_cluster_ids)<n_cell_clusters){
    n_cell_clusters = length(unique_cell_cluster_ids)
    updated_cell_clust <- mapvalues(updated_cell_clust, from = unique_cell_cluster_ids, to = 1:n_cell_clusters)
    n_target_gene_clusters <- n_target_gene_clusters[unique_cell_cluster_ids]
  }
}


plot_cluster_history(cell_cluster_history = cell_cluster_history)


