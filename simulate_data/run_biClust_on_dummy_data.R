execution_path <- dirname(rstudioapi::getSourceEditorContext()$path)
source(paste0(execution_path,"/dummy_data.R"))


Pt = 10      #number of target genes
Pr = 5       #number of regulator genes
n  = 10000   #number of cells
K  = 3       #Number of target gene clusters


#generate dummy data for each cell cluster that we want
dummy_res <- generate_dummy_data(Pt, Pr, n,K)

#hack to get the list dummy_res into the global environment
for(iter in 1:length(names(dummy_res))) {
  assign(eval(names(dummy_res)[iter]),dummy_res[[iter]])
}

# scRegOut$results
# plot(scRegOut)$data


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

K_cells <- 3

cell_data_split    <- sample(c(1,2),n, replace = T)
train_indices      <- which(cell_data_split == 1)
train_dat          <- dat[,train_indices]
initial_cell_clust <- kmeans(t(train_dat), K_cells)$cluster #get initial cell clustering

#preallocate outputs
out_list <- vector(mode = "list", length = K_cells)

for(inner_loop in 1:K_cells){
  #get local data
  local_dat <- train_dat[,which(initial_cell_clust == inner_loop)]

  cell_data_split    <- sample(c(1,2),ncol(local_dat), replace = T)    # train 1 val 2

  gene_cluster_start <- kmeans(local_dat[1:Pt,], K)$cluster        # initial clustering

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

out_list[[2]]$results

out_list[[3]]$results
