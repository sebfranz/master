rm(list = ls())
library(scregclust)
library(plyr)
library(aricode)  # To calculate rand index
library(Seurat)
library(Matrix)

execution_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
master_folder <- dirname(execution_path)
function_folder <- file.path(master_folder, "functions")
setwd(master_folder)
all_function_files <- list.files(function_folder, recursive=T, full.names=T)
for(current_file in all_function_files){
    print(paste("Loading", current_file))
    source(current_file)
  }

setwd(execution_path)
set.seed(1234)

#load sampled data
neftel_path <- file.path(dirname(master_folder), "datasets_sctargettranslator", "Neftel2019", "r_files", "neftel_seurat_group1.rds")
print(paste("Loading Neftel data from: ", neftel_path))
Neftel_g1 <- readRDS(file = neftel_path)

Neftel_g1 <- SetIdent(Neftel_g1, value='malignant')

# rm(Neftel_g1)

#get malignant expression on proper form
z_g1<-GetAssayData(Neftel_g1, slot='scale.data')
is_malignant <- ("yes" == Neftel_g1@meta.data$malignant) + 0
metaData<-Neftel_g1@meta.data$sample


out<-scregclust_format(z_g1)

sample_assignment<-out[[2]]
is_predictor<-out[[3]]

#put some stuff in order
z_g1_ordered <- z_g1[order(is_predictor),order(is_malignant)]
metaData     <- metaData[order(is_malignant)]
is_predictor <- is_predictor[order(is_predictor)]
is_malignant <- is_malignant[order(is_malignant)]

rm(z_g1)
rm(Neftel_g1)

dim_(z_g1_ordered)



############### code to run scregclust if we want to ##################

# for (i in 1:length(unique(metaData))){
#   ix<-which(metaData==unique(metaData)[i])
#   sample_assignment[ix]<-i
# }
#
# if(!'neftel_scregfit_group1' %in% list.files(paste0(neftel_path,'/r_files/'))){
#   fit_g1 <- scregclust(
#     z_g1_ordered,
#     genesymbols = scregclust_format(z_g1_ordered)[[1]],
#     is_predictor,
#     n_cl = 10L,
#     penalization = 0.05,
#     total_proportion = 0.5,
#     n_cycles = 50L,
#     noise_threshold = 0.05,
#     min_cluster_size = 50,
#     sample_assignment = sample_assignment
#   )
#
#   saveRDS(fit_g1, file = paste0(path, '/r_files/', "neftel_scregfit_group1"))
# } else {
#   fit_g1 <- readRDS(file = paste0(path, '/r_files/', "neftel_scregfit_group1"))
# }
#



############################################

# True cell clustering here is the malignancy for example
true_cell_clust <- is_malignant + 1
n_cell_clusters<- length(unique(true_cell_clust))

cell_sample <- c(sample(which(is_malignant==0), 1500 ),sample(which(is_malignant==1), 1500 ) )
# cell_sample <- 1:length(is_malignant)

#optional line to make next few lines not take forever
cols2keep <- sample(cell_sample,100, replace = F)



# initial_cell_clust <- sample(1:n_cell_clusters, n_cells, replace = T)

#optional: slice out some rows for efficiency
dim(z_g1_ordered)

regulator_vars <- sapply(which(is_predictor==1), function(i) var(z_g1_ordered[i,cols2keep]))

#pick out top 2 variance regulator genes
high_variance_regulators <- which(regulator_vars >= sort( regulator_vars, decreasing = T)[2])

#shift to correct index in entire dataset
high_variance_regulators <- high_variance_regulators + sum(is_predictor == 0)

rm(regulator_vars)
# target_vars <- sapply(which(is_predictor==0), function(i) var(z_g1_ordered[i,]))

#pick out chich targets are highly correlated with either of the chosen regulators
target_vars <- sapply(which(is_predictor==0), function(i) max(abs(cor(t(z_g1_ordered[c(i,high_variance_regulators),cols2keep]))[1,c(2,3)])))

#pick out top 2 variance target genes
high_corr_targets <- which(target_vars >= sort( target_vars, decreasing = T)[2])

rm(target_vars)
#implement a simplified version of biclust here. We want to sidestep scregclust but keep the rest of the concept intact

#first variant simplified_biclust()
max_iter <- 50


#find initial cluster labels
true_cell_clust_sample <- true_cell_clust[cell_sample]
n_target_gene_clusters <- 1 #we are not clustering target genes for now

# Kod för att flytta % av cellerna i varje kluster till ett annat kluster.
disturbed_initial_cell_clust <- true_cell_clust_sample
disturbed_fraction <- 0.15
for(i_cluster in 1:n_cell_clusters){
  indexes_of_cluster <- which(true_cell_clust_sample == i_cluster)
  some_of_those_indexes <- sample(indexes_of_cluster, size=as.integer(length(indexes_of_cluster)*disturbed_fraction), replace = F)
  disturbed_initial_cell_clust[some_of_those_indexes] <-
    sample(c(1:n_cell_clusters)[-i_cluster],
           size=length(some_of_those_indexes), replace=T)
}


cell_cluster_history <- cbind(true_cell_clust_sample, disturbed_initial_cell_clust)
colnames(cell_cluster_history) <- c("True allocation", "Disturbed allocation")

# #find which target gene has highest difference in mean expression between noncancer and cancer
# means1  <- sapply(1:sum((is_predictor==0)), function(row) mean(z_g1_ordered[which(is_predictor==0)[row], cell_sample[true_cell_clust_sample ==1]]) )
# means2  <- sapply(1:sum((is_predictor==0)), function(row) mean(z_g1_ordered[which(is_predictor==0)[row], cell_sample[true_cell_clust_sample ==2]]) )
# max_var_target_gene <- which(is_predictor==0)[which.max(means1-means2)]

#gather data and transpose to put
# dat <- t(z_g1_ordered[c(max_var_target_gene,which(is_predictor==1)), cell_sample])
dat <- t(z_g1_ordered[c( high_corr_targets, high_variance_regulators), cell_sample])
ind_targetgenes <- c(1,2)
ind_reggenes    <- c(3,4)
dim_(dat)


initial_clustering <- disturbed_initial_cell_clust
dim_(initial_clustering)

#set up some variables
n_cell_clusters = length(unique(initial_clustering))
# n_target_genes = ncol(dat) - sum(is_predictor==1)
# n_regulator_genes = sum(is_predictor==1)
n_target_genes = 2
n_regulator_genes = 2


#set up cluster history
initial_column_padding <- ncol(as.matrix(initial_clustering)) + 1  # +1 Because we have an index column that is not an index column it's an ID column
cell_cluster_history <- data.frame(matrix(NA, nrow = nrow(as.matrix(initial_clustering)), ncol = max_iter + initial_column_padding))
colnames(cell_cluster_history) <- c("Cell ID", 'Initial clustering', paste0("Iteration ", 1:max_iter))
cell_cluster_history[, 'Cell ID'] <-1:length(initial_clustering)  # Set cell names
cell_cluster_history[, 'Initial clustering'] <- initial_clustering

#preallocate all like matrices for later analysis if feasible
likelyhood_all <- vector("list", length = max_iter)
#set flag for breaking out of the loop.
stop_iterating_flag = F

current_cell_cluster_allocation <- initial_clustering

# i_main  <- 1

for(i_main in 1:max_iter){
  # fit model to each cell cluster
  models <- vector("list", length = n_cell_clusters)

  #corresponds to running screg
  for(cell_cluster in 1:n_cell_clusters){
    current_rows <- which(current_cell_cluster_allocation == cell_cluster)
    models[[cell_cluster]] <- lm(dat[current_rows, ind_targetgenes] ~ dat[current_rows, ind_reggenes])
  }
  plot_lm_planes(dat, current_cell_cluster_allocation, ind_reggenes, ind_targetgenes, n_cell_clusters, title=paste("Iteration", i_main))

  #for all cells, calculate the predicted  likelythingy for all cell clusters.
  like <- matrix(0, nrow = nrow(dat), ncol = n_cell_clusters)

  # Calculate the residual target gene variance for each gene and cluster
  # (so just one gene).
  INV_TARGET_GENE_RESIDUAL_STD <-  matrix(0, nrow = n_target_genes, ncol = n_cell_clusters)

  lambda <- 0.1

  PENALTY <- matrix(0, nrow = 1, ncol = n_cell_clusters)

  for(cell_cluster in 1:n_cell_clusters){

    intercept_and_regulatorgenes <- as.matrix(cbind(rep(1,nrow(dat)),dat[,ind_reggenes]))  #
    current_rows <- which(current_cell_cluster_allocation == cell_cluster)

    std_of_residuals_of_target_genes_temp <- dat[current_rows,ind_targetgenes] -
                                               (intercept_and_regulatorgenes[current_rows,] %*%  models[[cell_cluster]]$coefficients)
    INV_TARGET_GENE_RESIDUAL_STD[,cell_cluster]  <- 1/diag(sqrt(var(std_of_residuals_of_target_genes_temp)))

    PENALTY[cell_cluster] <- lambda * sum(abs(models[[cell_cluster]]$coefficients)) # *  INV_TARGET_GENE_RESIDUAL_STD[,cell_cluster] )

  }

  for(cell in 1:nrow(dat)){
    for(cell_cluster in 1:n_cell_clusters){
      # #bug fix hack: remove NA coefficients
      # if(any(is.na(models[[cell_cluster]]$coefficients))){
      #   NA_coeffs <-  unname(which(is.na(models[[cell_cluster]]$coefficients)))
      #   S_ERR <- (dat[cell,1] - as.vector(c(1,dat[cell,c(-1, -NA_coeffs)])) %*% models[[cell_cluster]]$coefficients[-NA_coeffs])^2
      # }

      SQUARED_ERROR <- (dat[cell,ind_targetgenes] - as.vector(c(1,dat[cell,ind_reggenes])) %*% models[[cell_cluster]]$coefficients)^2


      like[cell,cell_cluster] <- SQUARED_ERROR %*% INV_TARGET_GENE_RESIDUAL_STD[,cell_cluster] + PENALTY[cell_cluster]

    }
  }

  likelyhood_all[[i_main]] <- like


  #update cluster allocations
  updated_cell_clust <-  sapply(1:nrow(like), function(row) which.min(like[row,]))

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

png("simplified biclust_simplified neftel.png")
plot_cluster_history(cell_cluster_history = cbind( true_cell_clust_sample, cell_cluster_history))
dev.off()
