execution_path <- dirname(rstudioapi::getSourceEditorContext()$path)
source(paste0(execution_path,"/dummy_data.R"))


Pt               = 10      #number of target genes
Pr               = 5       #number of regulator genes
n                = 10000   #number of cells
K                = 3       #Number of target gene clusters
regulator_mean   = 1
coefficient_mean = 1

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
#append data
Z_t <- rbind( dummy_res_1$Z_t, dummy_res_2$Z_t, dummy_res_3$Z_t)
Z_r <- rbind( dummy_res_1$Z_r, dummy_res_2$Z_r, dummy_res_3$Z_r)
dat <- rbind(t(Z_t), t(Z_r)) #columns are now cells

cell_data_split    <- sample(c(1,2), nrow(Z_t), replace = T)
train_indices      <- which(cell_data_split == 1)
train_dat          <- dat[,train_indices]
initial_cell_clust <- kmeans(t(train_dat), K_cells)$cluster #get initial cell clustering

#preallocate outputs
out_list <- vector(mode = "list", length = K_cells)

for(inner_loop in 1:K_cells){
  #get local data
  local_dat <- train_dat[,which(initial_cell_clust == inner_loop)]

  cell_data_split    <- sample(c(1,2),ncol(local_dat), replace = T)    # train 1 val 2

  gene_cluster_start <- kmeans(local_dat[1:Pt,], K)$cluster        # initial target gene clustering

  #run scregclust
  scregclust(
    expression             = local_dat,               #scRegClust wants this form
    split_indices          = cell_data_split,         #train data split
    genesymbols            = 1:(Pt+Pr),               #gene row numbers
    is_regulator           = (1:(Pt+Pr) > Pt) + 0,    #vector indicating which genes are regulators
    n_cl                   = K,
    target_cluster_start   = gene_cluster_start,
    penalization           = Pr
    #maximal number of regulators for one cluster
  ) -> out_list[[inner_loop]]
}


out_list[[1]]$results

str(out_list[[1]], vec.len=Inf)

out_list[[2]]$results

out_list[[3]]$results

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






