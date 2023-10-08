#!/usr/bin/Rscript
library(scregclust)
library(plyr)
library(aricode)  # To calculate rand index

cluster_update <- function(metrics, n_cell_clusters, n_target_gene_clusters){
  # Minimises input matrix
  average_for_cell_cluster <- TRUE
  if (average_for_cell_cluster){
    d <- aggregate(metrics, list(rep(1:n_cell_clusters, n_target_gene_clusters)), FUN=function(x) mean(x, na.rm = TRUE))
    d <- as.matrix(d)
    d <- d[, 2:ncol(d), drop=FALSE]
    d <- as.vector(apply(d, 2, which.min))
  }else{
    d <- rep(1:n_cell_clusters, n_target_gene_clusters)[apply(metrics, 2, which.min)]
  }
  return(d)
}


calculate_R2_per_cluster <- function(n_target_gene_clusters,
                                     total_n_cells,
                                     n_cell_clusters,
                                     n_target_genes,
                                     n_regulator_genes,
                                     prev_cell_clust,
                                     out_list,
                                     train_dat,
                                     CELL_DATA_SPLIT){
  print("--- n_target_gene_clusters")
  print(str(n_target_gene_clusters))

  print("--- total_n_cells")
  print(str(total_n_cells))

  print("--- n_cell_clusters")
  print(str(n_cell_clusters))

  print("--- n_target_genes")
  print(str(n_target_genes))

  print("--- n_regulator_genes")
  print(str(n_regulator_genes))

  print("--- prev_cell_clust")
  print(str(prev_cell_clust))

  print("--- out_list")
  print(str(out_list))

  print("--- train_dat")
  print(str(train_dat))
  # Calculate R2 -----------------------------------------------------------

  # Calculate R2 for each cell for each target gene cluster
  # Output will be a matrix R2 with target gene clusters on rows and cells in columns.
  # Missing clusters are rows with NA.
  # which(rep(1:n_cell_clusters, n_target_gene_clusters) will be a vector telling you which rows will be in which cell cluster.
  # Example output R2:

  # Use a cell cluster as a 'from_cell_cluster' and then
  # calculate R2/predicitive R2 for that cluster applied to
  # all cell clusters' target gene clusters one by one.


  # R2 <- matrix(nrow=sum(n_target_gene_clusters), ncol=total_n_cells)
  adjusted_R2 <- matrix(nrow=n_cell_clusters, ncol=total_n_cells)
  # MSE <- matrix(nrow=sum(n_target_gene_clusters), ncol=total_n_cells)
  for(i_from_cell_cluster in 1:n_cell_clusters){
    # xvals expression of regulator cells in cell cluster i_from_cell_cluster
    # yvals expression of target cells in cell cluster i_from_cell_cluster
    regulator_gene_index <- (1:(n_target_genes+n_regulator_genes) > n_target_genes)
    target_gene_index <- (1:(n_target_genes+n_regulator_genes) <= n_target_genes)
    from_cell_cluster_index <- which(prev_cell_clust == i_from_cell_cluster)
    xvals <- train_dat[regulator_gene_index, from_cell_cluster_index, drop=FALSE]
    yvals <-  train_dat[target_gene_index, from_cell_cluster_index, drop=FALSE]
    # TODO: ska xvals och yvals alltid vara from_cell_cluster_index?
    center_all = FALSE
    if(center_all){
      yvals <- t(yvals)
      xvals <- t(xvals)
      yvals <- scale(yvals, scale=FALSE)
      xvals <- scale(xvals)
      yvals <- t(yvals)
      xvals <- t(xvals)
    }else{
      yvals <- t(yvals)
      xvals <- t(xvals)
      # print("cell data split")
      training_data_ind <- which(CELL_DATA_SPLIT[[i_from_cell_cluster]]==1)
      test_data_ind <- which(CELL_DATA_SPLIT[[i_from_cell_cluster]]==2)
      # if(length(training_data_ind)>length(test_data_ind)){
      #   test_data_ind <- c(test_data_ind, training_data_ind[1])
      # }
      # if(length(training_data_ind)<length(test_data_ind)){
      #   training_data_ind <- c(training_data_ind, test_data_ind[1])
      # }
      # print(str(training_data_ind))
      # print(str(test_data_ind))
      # print("yvals")
      # print(str(yvals[,training_data_ind]))
      org_yvals <- yvals
      org_xvals <- xvals
      # TODO: figure what to do when uneven train/test-split
      yvals[training_data_ind,] <- scale(yvals[training_data_ind,], scale = FALSE)
      xvals[training_data_ind,] <- scale(xvals[training_data_ind,])

      yvals[test_data_ind,] <- scale(yvals[test_data_ind,], colMeans(org_yvals[training_data_ind,]), scale = FALSE)
      xvals[test_data_ind,] <- scale(xvals[test_data_ind,], colMeans(org_xvals[training_data_ind,]), apply(org_xvals[training_data_ind,], 2, sd))

      yvals <- t(yvals)
      xvals <- t(xvals)
    }
    # i_total_target_geneclusters <- 0
    # Now we have xvals and yvals for a cell cluster we can find betas for each target gene cluster and calculate R2
    for(i_target_cell_cluster in 1:n_cell_clusters){
      indexes_of_not_deleted_target_gene_clusters <- which(sapply(out_list[[i_target_cell_cluster]]$results[[1]]$output[[1]]$coeffs, FUN=function(x) !is.null(x)))
      target_gene_cluster_index <- out_list[[i_target_cell_cluster]]$results[[1]]$output[[1]]$cluster[1:n_target_genes]  # n_target_genes, e.g. 1 1 1 1 3 3 3 3 3 3 3 3 1 3 1 3 3 3 3 3 1 3 1 1 3 1 3 3 1 3
      print("---- target_gene_cluster_index")
      print(target_gene_cluster_index)
      target_gene_cluster_index <- unlist(sapply(indexes_of_not_deleted_target_gene_clusters, FUN=function(x) which(target_gene_cluster_index==x)))
      # target_gene_cluster_index <- which(target_gene_cluster_index!=-1) # For target genes that wher put in noise cluster are -1
      # todo: sortera ---- target_gene_cluster_index så den är samma ordning för klustren som i betas_for_cells_cluster_i

      print("--- from to")
      print(paste(i_from_cell_cluster, i_target_cell_cluster))
      print("---- target_gene_cluster_index")
      print(target_gene_cluster_index)
      print("---- indexes_of_not_deleted_target_gene_clusters")
      print(indexes_of_not_deleted_target_gene_clusters)

      betas_for_cells_cluster_i <- out_list[[i_target_cell_cluster]]$results[[1]]$output[[1]]$coeffs[indexes_of_not_deleted_target_gene_clusters]
      print("---- betas_for_cells_cluster_i")
      print(str(betas_for_cells_cluster_i))

      betas_for_cells_cluster_i <- do.call(cbind, betas_for_cells_cluster_i)  # rows: n_regulator_genes, cols: n_target_genes

      print("---- betas_for_cells_cluster_i")
      print(str(betas_for_cells_cluster_i))

      # for(i_target_gene_cluster in 1:n_target_gene_clusters[[i_target_cell_cluster]]){
      #   i_total_target_geneclusters <- i_total_target_geneclusters + 1
      #   target_gene_ids_in_cluster_i <- which(target_gene_cluster_index==i_target_gene_cluster)
      #   betas_for_gene_cluster_i <- out_list[[i_target_cell_cluster]]$results[[1]]$output[[1]]$coeffs[[i_target_gene_cluster]]
      # If no target gene was assigned to this cluster we don't calculate anything.
      # This means rows of NA will represent target gene clusters that disappeared.
      if(!is.null(betas_for_cells_cluster_i)) {
        print(colSums(betas_for_cells_cluster_i))
        print(rowSums(betas_for_cells_cluster_i))
        target_gene_clusters_in_this_cell_cluster_i_yvals <- yvals[target_gene_cluster_index, , drop=FALSE]
        colmeaned_y <-  matrix(rep(colMeans(target_gene_clusters_in_this_cell_cluster_i_yvals), nrow(target_gene_clusters_in_this_cell_cluster_i_yvals)), nrow=nrow(target_gene_clusters_in_this_cell_cluster_i_yvals))
        SST <- (target_gene_clusters_in_this_cell_cluster_i_yvals - colmeaned_y)^2
        # SST <- (target_gene_clusters_in_this_cell_cluster_i_yvals - mean(target_gene_clusters_in_this_cell_cluster_i_yvals))^2
        SST <- as.matrix(SST, nrow=length(target_gene_ids_in_cluster_i), ncol=ncol(target_gene_clusters_in_this_cell_cluster_i_yvals))
        SST_sum_adjusted <- sum(SST)/(ncol(target_gene_clusters_in_this_cell_cluster_i_yvals) - nrow(xvals))

        sigmas <- sqrt(unlist(out_list[[i_target_cell_cluster]]$results[[1]]$output[[1]]$sigmas)^2)
        # print("Sigmas")
        # print(sigmas)
        # sigmas <- matrix(rep(sigmas,ncol(target_gene_clusters_in_this_cell_cluster_i_yvals)), nrow=length(sigmas), ncol=ncol(target_gene_clusters_in_this_cell_cluster_i_yvals))
        SSR <- (target_gene_clusters_in_this_cell_cluster_i_yvals - (t(betas_for_cells_cluster_i*0) %*% xvals))^2
        SSR <- as.matrix(SSR, nrow=length(target_gene_ids_in_cluster_i), ncol=ncol(target_gene_clusters_in_this_cell_cluster_i_yvals))/sigmas
        SSR_sum_adjusted <- colSums(SSR)/(ncol(target_gene_clusters_in_this_cell_cluster_i_yvals) - 1)

        # R2[i_total_target_geneclusters, prev_cell_clust == i_from_cell_cluster] <- 1 - colSums(SSR)/sum(SST)
        adjusted_R2[i_target_cell_cluster, prev_cell_clust == i_from_cell_cluster] <- 1 - SSR_sum_adjusted/SST_sum_adjusted
        # adjusted_R2[i_target_cell_cluster, prev_cell_clust == i_from_cell_cluster] <- 1 - SSR_sum_adjusted/SST_sum_adjusted
        # MSE[i_total_target_geneclusters, prev_cell_clust == i_from_cell_cluster] <- colMeans(SSR)
        # }
      }
    }
  }

  return(adjusted_R2)
}

calculate_R2 <- function(n_target_gene_clusters,
                         total_n_cells,
                         n_cell_clusters,
                         n_target_genes,
                         n_regulator_genes,
                         prev_cell_clust,
                         out_list,
                         train_dat){
  print("--- n_target_gene_clusters")
  print(str(n_target_gene_clusters))

  print("--- total_n_cells")
  print(str(total_n_cells))

  print("--- n_cell_clusters")
  print(str(n_cell_clusters))

  print("--- n_target_genes")
  print(str(n_target_genes))

  print("--- n_regulator_genes")
  print(str(n_regulator_genes))

  print("--- prev_cell_clust")
  print(str(prev_cell_clust))

  print("--- out_list")
  print(str(out_list))

  print("--- train_dat")
  print(str(train_dat))

  # save(list = c("train_dat", "out_list", "prev_cell_clust", "n_regulator_genes", "n_target_genes", "n_cell_clusters", "total_n_cells", "n_target_gene_clusters"), file = "test_data.Rdata")
  # Calculate R2 -----------------------------------------------------------

  # Calculate R2 for each cell for each target gene cluster
  # Output will be a matrix R2 with target gene clusters on rows and cells in columns.
  # Missing clusters are rows with NA.
  # which(rep(1:n_cell_clusters, n_target_gene_clusters) will be a vector telling you which rows will be in which cell cluster.
  # Example output R2:
  # num [1:6, 1:16000]
  #          X1        X2        X3        X4        X5        X6        X7
  # 1 0.9999561 0.9999349 0.9999627 0.9999567 0.9999449 0.9999573 0.9999398
  # 2        NA        NA        NA        NA        NA        NA        NA
  # 3 0.9999384 0.9999202 0.9999379 0.9999344 0.9999321 0.9999429 0.9999274
  # 4        NA        NA        NA        NA        NA        NA        NA
  # 5 0.9999168 0.9998943 0.9999142 0.9999105 0.9999101 0.9999234 0.9999042
  # 6        NA        NA        NA        NA        NA        NA        NA

  # Use a cell cluster as a 'from_cell_cluster' and then
  # calculate R2/predicitive R2 for that cluster applied to
  # all cell clusters' target gene clusters one by one.

  # R2 <- matrix(nrow=sum(n_target_gene_clusters), ncol=total_n_cells)
  adjusted_R2 <- matrix(nrow=sum(n_target_gene_clusters), ncol=total_n_cells)
  # MSE <- matrix(nrow=sum(n_target_gene_clusters), ncol=total_n_cells)
  for(i_from_cell_cluster in 1:n_cell_clusters){
    # xvals expression of regulator cells in cell cluster i_from_cell_cluster
    # yvals expression of target cells in cell cluster i_from_cell_cluster
    regulator_gene_index <- (1:(n_target_genes+n_regulator_genes) > n_target_genes)
    target_gene_index <- (1:(n_target_genes+n_regulator_genes) <= n_target_genes)
    from_cell_cluster_index <- which(prev_cell_clust == i_from_cell_cluster)
    xvals <- train_dat[regulator_gene_index, from_cell_cluster_index, drop=FALSE]
    yvals <-  train_dat[target_gene_index, from_cell_cluster_index, drop=FALSE]
    i_total_target_geneclusters <- 0
    # Now we have xvals and yvals for a cell cluster we can find betas for each target gene cluster and calculate R2
    for(i_target_cell_cluster in 1:n_cell_clusters){
      target_gene_cluster_index <- out_list[[i_target_cell_cluster]]$results[[1]]$output[[1]]$cluster[1:n_target_genes]
      for(i_target_gene_cluster in 1:n_target_gene_clusters[[i_target_cell_cluster]]){
        i_total_target_geneclusters <- i_total_target_geneclusters + 1
        target_gene_ids_in_cluster_i <- which(target_gene_cluster_index==i_target_gene_cluster)
        betas_for_gene_cluster_i <- out_list[[i_target_cell_cluster]]$results[[1]]$output[[1]]$coeffs[[i_target_gene_cluster]]
        # If no target gene was assigned to this cluster we don't calculate anything.
        # This means rows of NA will represent target gene clusters that disappeared.
        if(!is.null(betas_for_gene_cluster_i)) {
          target_gene_cluster_yvals <- yvals[target_gene_ids_in_cluster_i, , drop=FALSE]

          SST <- (target_gene_cluster_yvals - mean(target_gene_cluster_yvals))^2
          SST <- as.matrix(SST, nrow=length(target_gene_ids_in_cluster_i), ncol=ncol(target_gene_cluster_yvals))
          SST_sum_adjusted <- sum(SST)/(ncol(target_gene_cluster_yvals) - nrow(xvals))

          SSR <- (target_gene_cluster_yvals- t(betas_for_gene_cluster_i) %*% xvals)^2
          SSR <- as.matrix(SSR, nrow=length(target_gene_ids_in_cluster_i), ncol=ncol(target_gene_cluster_yvals))
          SSR_sum_adjusted <- colSums(SSR)/(ncol(target_gene_cluster_yvals) - 1)

          # R2[i_total_target_geneclusters, prev_cell_clust == i_from_cell_cluster] <- 1 - colSums(SSR)/sum(SST)
          adjusted_R2[i_total_target_geneclusters, prev_cell_clust == i_from_cell_cluster] <- 1 - SSR_sum_adjusted/SST_sum_adjusted
          # MSE[i_total_target_geneclusters, prev_cell_clust == i_from_cell_cluster] <- colMeans(SSR)
        }
      }
    }
  }
  return(adjusted_R2)
}

biclust <- function(max_iter=50,
                    initial_cluster_history,
                    is_regulator,  # Input dataset has to have rows sorted so that targets are highest
                    n_target_gene_clusters = c(3,4,5),
                    penalization_parameter = 0.14,
                    train_dat,
                    plot_r2 = TRUE,
                    ...){

  # Setup
  n_cell_clusters <- length(unique(initial_cluster_history[, 1]))
  n_target_genes <- sum(is_regulator==0)
  n_regulator_genes <- sum(is_regulator==1)
  total_n_cells <- ncol(train_dat)

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
  stop_iterating_flag <- F

  for (i_main in 1:max_iter){
    print(paste("Iteration", i_main))
    # Run scregclust for each cell cluster ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    # Preallocate outputs
    out_list <- vector(mode = "list", length = n_cell_clusters)
    # List of last iterations gene clustering
    prev_target_gene_clusters <- vector(mode = "list", length = n_cell_clusters)
    CELL_DATA_SPLIT <- vector(mode = "list", length = n_cell_clusters)
    for(i_cluster in 1:n_cell_clusters){
      # Get data for cell cluster i_cluster specifically
      local_dat <- train_dat[, which(prev_cell_clust == i_cluster), drop=FALSE]
      # local_dat <- as.matrix(local_dat, nrow=, ncol=length(which(prev_cell_clust == i_cluster)))
      # Training data are represented by 1 and test data by 2
      cell_data_split    <- sample(rep(1:2,length.out=ncol(local_dat)))
      CELL_DATA_SPLIT[[i_cluster]] <- cell_data_split
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
        is_regulator           = is_regulator,  # Vectorindicatsing which genes are regulators
        n_cl                   = n_target_gene_clusters[i_cluster],
        target_cluster_start   = target_gene_cluster_start,
        penalization           = penalization_parameter,
        verbose                = FALSE,
        max_optim_iter         = 10000,
        ...
      ) -> out_list[[i_cluster]]

      # Store gene clustering for next iteration
      prev_target_gene_clusters[[i_cluster]] <- out_list[[i_cluster]]$results[[1]]$output[[1]]$cluster
    }



    R2 <- calculate_R2_per_cluster(n_target_gene_clusters,
                                   total_n_cells,
                                   n_cell_clusters,
                                   n_target_genes,
                                   n_regulator_genes,
                                   prev_cell_clust,
                                   out_list,
                                   train_dat,
                                   CELL_DATA_SPLIT)

    print("R2----------")
    print(str(R2))
    print((data.frame(R2[, 1:7])))
    print("------------")

    # Plot histograms of r2
    if(plot_r2){
      par(mfrow=c(length(unique(prev_cell_clust)),n_cell_clusters))
      for(i_cells_from_cell_cluster in 1:length(unique(prev_cell_clust))){
        for(i_fits_into_cell_cluster in 1:length(unique(prev_cell_clust))){
          print(paste(i_cells_from_cell_cluster, i_fits_into_cell_cluster))
          ind_for_cell_cluster <- which(rep(1:n_cell_clusters, n_target_gene_clusters)==i_fits_into_cell_cluster)
          hist(R2[ind_for_cell_cluster, prev_cell_clust==i_cells_from_cell_cluster], breaks=10000, main=paste("Cells from cell cluster", i_cells_from_cell_cluster, "\nfits into cell cluster", i_fits_into_cell_cluster, "with r2:"))
        }
      }
      mtext(paste("Iteration", i_main), side = 3, line = -1, outer = TRUE)
    }
    # Update cluster allocation to the appropriate cell cluster -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    # This can be done in some different ways,
    # One way would be to take mean squared error for the totality of the fitted model for that cell cluster
    # Another way would be to compare the average mse per gene cluster model
    # here we will compare the minimal mse per gene cluster model per cell cluster.
    # Could also use other metric than mse, e.g. r2

    if (all(is.na(R2))){
      stop("scregclust put everything in noise cluster for all cellclusters. Exiting.")
    }
    # updated_cell_clust <-  cluster_update(-R2, n_cell_clusters, n_target_gene_clusters)
    updated_cell_clust <-  (1:n_cell_clusters)[apply(R2, 2, which.max)]

    # If there is only 1 cell left in a cell cluster it doesn't work
    # Move that cell to the biggest cell cluster.
    updated_cell_clust_table <- data.frame(table(updated_cell_clust))

    for(i_cell_cluster in 1:n_cell_clusters){
      if (length(which(updated_cell_clust==i_cell_cluster)) == 1){
        updated_cell_clust[which(updated_cell_clust==i_cell_cluster)] <- which.max(updated_cell_clust_table$Freq)
      }
    }

    # Update data in cell_cluster_history
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
        stop_iterating_flag <- T
        break
      }
    }

    if(stop_iterating_flag == T){
      # Clean up cluster history
      cell_cluster_history <- cell_cluster_history[ , colSums(is.na(cell_cluster_history))==0, drop=FALSE]
      # Stop iterations/exit function
      break
    }

    # Run next iteration with updated cell clustering
    prev_cell_clust <- updated_cell_clust

    # Handle case where cell cluster disappeared
    unique_cell_cluster_ids <- sort(unique(prev_cell_clust))
    if (length(unique_cell_cluster_ids)<n_cell_clusters){
      n_cell_clusters <- length(unique_cell_cluster_ids)
      prev_cell_clust <- mapvalues(prev_cell_clust, from = unique_cell_cluster_ids, to = 1:n_cell_clusters)
      n_target_gene_clusters <- n_target_gene_clusters[unique_cell_cluster_ids]
    }
  }

  return(list(cell_cluster_history=cell_cluster_history))
}


# runs only when script is run by itself
if (sys.nframe() == 0){
  # ... do main stuff
}
