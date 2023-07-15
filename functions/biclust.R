library(scregclust)
library(plyr)
library(aricode)  # To calculate rand index

biclust <- function(max_iter=50,
                    initial_cluster_history,
                    is_regulator, #input dataset has to have rows sorted so that targets are highest
                    n_target_gene_clusters = c(3,4,5),  #also necessary
                    n_cells = c(1000,5000,10000),
                    train_dat){

  #pre-setup
  n_cell_clusters <- length(n_cells)
  n_target_genes = sum(is_regulator==0)
  n_regulator_genes = sum(is_regulator==1)
  total_n_cells = sum(n_cells)

  # Preallocate memory
  if(!is.matrix(initial_cluster_history)){
    initial_cluster_history <- as.matrix(initial_cluster_history)
    colnames(initial_cluster_history) <- c("initial")
  }

  initial_column_padding <- ncol(initial_cluster_history) + 1  # +1 Because we have an index column that is not an index column it's an ID column
  cell_cluster_history <- data.frame(matrix(NA, nrow = nrow(initial_cluster_history), ncol = max_iter + initial_column_padding))
  colnames(cell_cluster_history) <- c("Cell ID", colnames(initial_cluster_history), paste0("Iteration ", 1:max_iter))
  cell_cluster_history[, 'Cell ID'] <-1:nrow(initial_cluster_history)  # Set cell names
  cell_cluster_history[, colnames(initial_cluster_history)] <- initial_cluster_history
  prev_cell_clust <- initial_cluster_history[, ncol(initial_cluster_history)]

  # Set exit flag
  stop_iterating_flag = F;

  for (i_main in 1:max_iter){
    print(paste("Iteration", i_main))
    # Run scregclust for each cell cluster ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    # Preallocate outputs
    out_list <- vector(mode = "list", length = n_cell_clusters)
    # List of last iterations gene clustering
    prev_target_gene_clusters <- vector(mode = "list", length = n_cell_clusters)

    for(i_cluster in 1:n_cell_clusters){
      # Get data for cell cluster i_cluster specifically
      local_dat <- train_dat[,which(prev_cell_clust == i_cluster)]

      # Training data are represented by 1 and test data by 2
      cell_data_split    <- sample(c(1,2),ncol(local_dat), replace = T)

      # The first iteration we initialize with kmeans, later iterations we start with
      # the same gene clustering as last time, possibly more efficient.
      if(i_main == 1){
        # Initial target gene clustering
        target_gene_cluster_start <- kmeans(local_dat[1:n_target_genes,], n_target_gene_clusters[i_cluster])$cluster
      } else {
        target_gene_cluster_start <- prev_target_gene_clusters[[i_cluster]]
      }

      # Run scregclust
      scregclust(
        expression             = local_dat,  # p rows of genes and n columns of cells of single cell expression data
        split_indices          = cell_data_split,  # Train/test data split indicated by 1s and 2s\
        genesymbols            = paste0('g', 1:(n_target_genes+n_regulator_genes)),  # Gene row names
        is_regulator           = (1:(n_target_genes+n_regulator_genes) > n_target_genes) + 0,  # Vectorindicating which genes are regulators
        n_cl                   = n_target_gene_clusters[i_cluster],
        target_cluster_start   = target_gene_cluster_start,
        penalization           = 0.14,
        verbose                = FALSE,
        max_optim_iter         = 100,
      ) -> out_list[[i_cluster]]

    # Store gene clustering for next iteration
    prev_target_gene_clusters[[i_cluster]] <- out_list[[i_cluster]]$results[[1]]$output[[1]]$cluster
    }



    # Calculate MSE -----------------------------------------------------------

    # Calculate MSE for per cells per cell cluster,
    # for all target gene clusters in those cell clusters
    # Missing clusters are rows with NA
    # matrix of cell cluster 1's target genes, then cell cluster 2's target genes - etc. In order. Columns are cells matching prev_cell_clust
    # Example output
    # > str(MSE)
    # num [1:12, 1:16000] 68.9 124.1 11.9 66.6 119.9 ...


    MSE <- matrix(nrow=sum(n_target_gene_clusters), ncol=total_n_cells)
    for(i_cell_cluster in 1:n_cell_clusters){
      # Expression of all cells
      # xvals expression of regulator cells
      # yvals expression of target cells
      xvals <- train_dat[(1:(n_target_genes+n_regulator_genes) > n_target_genes), which(prev_cell_clust == i_cell_cluster)]
      yvals <-  train_dat[(1:(n_target_genes+n_regulator_genes) <= n_target_genes), which(prev_cell_clust == i_cell_cluster)]
      i_total_target_geneclusters <- 0
      for(ii_cell_cluster in 1:n_cell_clusters){
        clustering <- out_list[[ii_cell_cluster]]$results[[1]]$output[[1]]$cluster[1:n_target_genes]
        for(i_target_gene_cluster in 1:n_target_gene_clusters[[ii_cell_cluster]]){
          i_total_target_geneclusters <- i_total_target_geneclusters + 1
          target_gene_ids_in_cluster_i <- which(clustering==i_target_gene_cluster)
          betas_for_gene_cluster_i <- out_list[[ii_cell_cluster]]$results[[1]]$output[[1]]$coeffs[[i_target_gene_cluster]]
          # If no target gene was assigned to this cluster we need to do something else
          if(!is.null(betas_for_gene_cluster_i)) {
            MSE[i_total_target_geneclusters, prev_cell_clust == i_cell_cluster] <- colMeans((yvals[target_gene_ids_in_cluster_i,] - t(betas_for_gene_cluster_i) %*% xvals)**2)
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

    # Update data in cell_cluster_history
    # skip
    cell_cluster_history[, i_main + initial_column_padding] <- updated_cell_clust

    # Cross tabulation of clusters
    print("Table")
    print(data.frame(table(updated_cell_clust, prev_cell_clust)))
    flush.console()

    # if (all(prev_cell_clust == updated_cell_clust)){
    #   print("Cell clustering same as last iteration. Exiting.")
    #   break
    # }

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
      cell_cluster_history <- cell_cluster_history[ , colSums(is.na(cell_cluster_history))==0]
      # Stop iterations/exit function
      break
    }

    # Run next iteration with updated cell clustering
    prev_cell_clust <- updated_cell_clust

    # Handle case where cell cluster disappeared
    unique_cell_cluster_ids <- sort(unique(prev_cell_clust))
    if (length(unique_cell_cluster_ids)<n_cell_clusters){
      n_cell_clusters = length(unique_cell_cluster_ids)
      prev_cell_clust <- mapvalues(prev_cell_clust, from = unique_cell_cluster_ids, to = 1:n_cell_clusters)
      n_target_gene_clusters <- n_target_gene_clusters[unique_cell_cluster_ids]
    }
  }

  return(list(cell_cluster_history=cell_cluster_history))
}
