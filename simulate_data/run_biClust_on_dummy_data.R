rm(list = ls())
execution_path <- dirname(rstudioapi::getSourceEditorContext()$path)
source(paste0(execution_path,"/dummy_data.R"))
library(scregclust)
set.seed(12411)

Pt               = 10      #number of target genes
Pr               = 5       #number of regulator genes
n                = 10000   #number of cells
K                = 3       #Number of target gene clusters
regulator_mean   = 1
coefficient_mean = c(1,10,100)

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

K_cells <- 3 #number of cell clusters

#generate dummy data for each cell cluster that we want
dummy_res_1 <- generate_dummy_data(Pt, Pr, n, K,
                                   regulator_mean, coefficient_mean = 1)
dummy_res_2 <- generate_dummy_data(Pt, Pr, n, K,
                                   regulator_mean, coefficient_mean = 2)
dummy_res_3 <- generate_dummy_data(Pt, Pr, n, K,
                                   regulator_mean, coefficient_mean = 3)
# Append data
Z_t <- rbind( dummy_res_1$Z_t, dummy_res_2$Z_t, dummy_res_3$Z_t)
Z_r <- rbind( dummy_res_1$Z_r, dummy_res_2$Z_r, dummy_res_3$Z_r)
dat <- rbind(t(Z_t), t(Z_r)) #columns are now cells

#update parameters
n <- ncol(dat)

# Split into train and test data for cell clustering
# Skip for now
# cell_data_split    <- sample(c(1,2), nrow(Z_t), replace = T)
# train_indices      <- which(cell_data_split == 1)
# train_dat          <- dat[,train_indices]

train_dat <- dat

#get initial cell clustering
initial_cell_clust <- kmeans(t(train_dat), K_cells)$cluster
# initial_cell_clust <- sample(1:K_cells, n, replace = T)

previous_cell_clust <- initial_cell_clust

# preallocate outputs
out_list <- vector(mode = "list", length = K_cells)

for(inner_loop in 1:K_cells){
  # get local data
  local_dat <- train_dat[,which(previous_cell_clust == inner_loop)]

  cell_data_split    <- sample(c(1,2),ncol(local_dat), replace = T)    # train 1 val 2

  gene_cluster_start <- kmeans(local_dat[1:Pt,], K)$cluster        # initial target gene clustering

  #run scregclust
  scregclust(
    expression             = local_dat,               #scRegClust wants this form
    split_indices          = cell_data_split,         #train data split
    genesymbols            = paste0('g', 1:(Pt+Pr)),  #gene row names
    is_regulator           = (1:(Pt+Pr) > Pt) + 0,    #vector indicating which genes are regulators
    n_cl                   = K,
    target_cluster_start   = gene_cluster_start,
    penalization           = Pr
    #maximal number of regulators for one cluster
  ) -> out_list[[inner_loop]]
}

# which target genes belong to which cluster>
str(out_list[[1]]$results[[1]])
str(out_list[[1]]$results[[1]]$output[[1]])



# Calculate MSE -----------------------------------------------------------

# Calculate MSE for per cells per cell cluster,
# for all target gene clusters in those cell clusters

MSE <- vector(mode = "list", length = K_cells)
names(MSE) <- paste0('Cell cluster ', 1:K_cells)
for(i_cell_cluster in 1:K_cells){

  # Expression of all cells
  # xvals expression of regulator cells
  # yvals expression of target cells
  xvals <- train_dat[(1:(Pt+Pr) > Pt), which(previous_cell_clust == i_cell_cluster)]
  yvals <-  train_dat[(1:(Pt+Pr) <= Pt), which(previous_cell_clust == i_cell_cluster)]

  MSE_in_cell_cluster_i <- vector(mode = "list", length = K_cells)
  names(MSE_in_cell_cluster_i) <- paste0('... MSE for cell cluster ', 1:K_cells)
  for(ii_cell_cluster in 1:K_cells){
    clustering <- out_list[[ii_cell_cluster]]$results[[1]]$output[[1]]$cluster[1:Pt]
    MSE_in_cell_cluster_ii <- vector(mode = "list", length = K)
    names(MSE_in_cell_cluster_ii) <- paste0('... MSE for target gene cluster ', 1:K_cells)
    for(i_target_gene_cluster in 1:K){
      betas_for_gene_cluster_i <- out_list[[ii_cell_cluster]]$results[[1]]$output[[1]]$coeffs[[i_target_gene_cluster]]
      target_gene_ids_in_cluster_i <- which(clustering==i_target_gene_cluster)
      MSE_in_cell_cluster_ii[[i_target_gene_cluster]] <- colMeans((yvals[target_gene_ids_in_cluster_i,] - t(betas_for_gene_cluster_i) %*% xvals)**2)
    }
    MSE_in_cell_cluster_i[[ii_cell_cluster]] <- MSE_in_cell_cluster_ii
  }
  #store results in a list of size K_cells(cell clusters)
    #of lists of K(gene clusters)
  MSE[[i_cell_cluster]] <- MSE_in_cell_cluster_i
}
str(MSE)

#update cluster allocation to the appropriate cell cluster
  # This can be done in some different ways,
  # One way would be to take mean squared error for the totality of the fitted model for that cell cluster
  # Another way would be to compare the average mse per gene cluster model
  # here we will compare the minimal mse per gene cluster model per cell cluster.
  # Could also use other metric than mse, e.g. r2

updated_cell_clust <- previous_cell_clust + NA

for(cell in 1:ncol(dat)){
  which_cell_cluster <- previous_cell_clust[cell]#should update to "current" cell cluster
  which_within_cluster_index <- sum(previous_cell_clust[1:(cell)] == previous_cell_clust[cell])
  metrics <- matrix(data = NA,  nrow = K_cells, ncol = 1)

  for(cell_cluster in 1:K_cells){
    #find cluster allocation metrics
    metrics[cell_cluster] <- min(do.call(rbind, MSE[[which_cell_cluster]][[cell_cluster]])[,which_within_cluster_index])
  }
  updated_cell_clust[cell] <- which.min(metrics)
}

# cross tabulation of clusters
data.frame( table(updated_cell_clust, previous_cell_clust))

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
