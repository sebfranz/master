library(scregclust)
library(plyr)
library(aricode)  # To calculate rand index

biclust <- function(max_iter=50,
                    initial_cell_clust,
                    disturbed_initial_cell_clust,
                    train_dat){
  cell_cluster_history <- cbind(initial_cell_clust, disturbed_initial_cell_clust)
  previous_cell_clust <- disturbed_initial_cell_clust

  cell_cluster_history <- data.frame(matrix(NA, nrow = length(previous_cell_clust), ncol = max_iter + 2)) #preallocate memory
  colnames(cell_cluster_history) <- c("Cell_id", "Initial", paste0("Iteration_", 1:max_iter)) #set colnames
  cell_cluster_history$Cell_id <-1:length(previous_cell_clust) #set cell names
  cell_cluster_history$Initial <- previous_cell_clust

  # Set exit flag
  stop_iterating_flag = F;

  for (i_main in 1:max_iter){
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
        verbose                = FALSE,
        max_optim_iter         = 100000,
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

    # Update data in cell_cluster_history
    # skip
    cell_cluster_history[, i_main + 2] <- updated_cell_clust

    # Cross tabulation of clusters
    print("Table")
    print(data.frame(table(updated_cell_clust, previous_cell_clust)))
    flush.console()

    # if (all(previous_cell_clust == updated_cell_clust)){
    #   print("Cell clustering same as last iteration. Exiting.")
    #   break
    # }

    # Compare with previous iterations
    for(prev_clustering in ((i_main-1):0) ){
      print(paste0('comparing with iteration ', prev_clustering))
      if(RI(updated_cell_clust,cell_cluster_history[,prev_clustering + 2]) == 1){
        print("Cell clustering from iteration  same as some previous iteration. Exiting.")
        print(paste0("RI of ",
                     RI(updated_cell_clust,cell_cluster_history[,prev_clustering + 2]),
                     " when comparing iteration ", i_main," to iteration ", prev_clustering))
        stop_iterating_flag = T
        break
      }
    }

    if(stop_iterating_flag == T){
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

  return(list(cell_cluster_history=cell_cluster_history))
}
