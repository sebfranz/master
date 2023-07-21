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

#load sampled data
neftel_path <- paste0(getCurrentFileLocation(), '/../../datasets_sctargettranslator/Neftel2019')
Neftel_g1 <- readRDS(file = paste0(neftel_path, '/r_files/', "neftel_seurat_group1"))

Neftel_g1<-SetIdent(Neftel_g1, value='malignant')

# rm(Neftel_g1)

#get malignant expression on proper form
z_g1<-GetAssayData(Neftel_g1, slot='scale.data')
is_malignant <- ("yes" == Neftel_g1@meta.data$malignant) + 0
metaData<-Neftel_g1@meta.data$sample


out<-scregclust_format(z_g1)

sample_assignment<-out[[2]]
is_predictor<-out[[3]]

sum(is_predictor==1)
length(is_predictor)
length(is_malignant)
dim(z_g1)

#put some stuff in order
z_g1_ordered <- z_g1[order(is_predictor),order(is_malignant)]
metaData     <- metaData[order(is_malignant)]
is_predictor <- is_predictor[order(is_predictor)]
is_malignant <- is_malignant[order(is_malignant)]

rm(z_g1)
rm(Neftel_g1)

dim_(z_g1_ordered)
dim_(is_predictor)
dim_(is_malignant)


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

# initial_cell_clust <- sample(1:n_cell_clusters, n_cells, replace = T)


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

#optional: slice out some rows for efficiency
dim(z_g1_ordered)
cols2keep <- sample(ncol(z_g1_ordered), 4500)

vars <- sapply(which(is_predictor==0), function(i) var(z_g1_ordered[i,]))
#pick out top 500 variance target genes
vars_ <- which(vars>sort(vars,decreasing = T)[500])



#z_g1_ordered
#is_predictor
#is_malignant


#Make a simplified version of biclust here. We want to sidestep scregclust but keep the rest of the concept intact

#This essentially means more simply to fit a model for each cluster,
# get SOME metric for model fit for a new observation(predictive r2)
#and reassign clusters based on this


#first. the almost trivial case: a normal linear model between all regulatory genes and ONE target (response) gene

#second, a coop lasso using all target genes, but no clustering
#alternatively, a class of ordinary linear models
# install.packages('scoop')
# library(scoop)

#third step, would be scregclust, a coop lasso clustering method


#first variant simplified_biclust()
max_iter <- 50
# cell_sample <- c(sample(which(is_malignant==0), 1500 ),sample(which(is_malignant==1), 1500 ) )

cell_sample <- 1:length(is_malignant)

#find initial cluster labels
initial_clustering <- disturbed_initial_cell_clust[cell_sample]
true_cell_clust_sample <- true_cell_clust[cell_sample]
n_target_gene_clusters <- 1 #we are not clustering target genes for now

#find which target gene has highest difference in mean expression between noncancer and cancer
means1  <- sapply(1:sum((is_predictor==0)), function(row) mean(z_g1_ordered[which(is_predictor==0)[row], cell_sample[true_cell_clust_sample ==1]]) )
means2  <- sapply(1:sum((is_predictor==0)), function(row) mean(z_g1_ordered[which(is_predictor==0)[row], cell_sample[true_cell_clust_sample ==2]]) )
max_var_target_gene <- which(is_predictor==0)[which.max(means1-means2)]

#gather data and transpose to put
dat <- t(z_g1_ordered[c(max_var_target_gene,which(is_predictor==1)), cell_sample])
dim_(dat)

dim_(initial_clustering)

#set up some variables
n_cell_clusters = length(unique(initial_clustering))
n_target_genes = ncol(dat) - sum(is_predictor==1)
n_regulator_genes = sum(is_predictor==1)


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


