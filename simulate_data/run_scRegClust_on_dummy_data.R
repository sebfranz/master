#!/usr/bin/Rscript
# This script check how many times out of n_loops we have rand index 1
# It generates data and then runs scregclust and then calculate rand index, e.g. 1000 times.
if (!require(here)) install.packages('here')
if (!require(here)) install.packages('aricode')
library(scregclust)
library(here)
library(aricode)  # To calculate rand index

root_path <- here::here()
execution_path <- file.path(root_path, "simulate_data")
generate_dummy_data_function_path <- file.path(execution_path,"functions", "generate_dummy_data_for_scregclust.R")
source(generate_dummy_data_function_path)

# set.seed(14)
Pt = 10      #number of target genes
Pr = 5       #number of regulator genes
n  = 10000   #number of cells
K  = 3       #Number of target gene clusters
regulator_mean   = 1
coefficient_mean = c(1,10,100)
n_loops = 1000


rands <- vector(length = n_loops)
for(i_loop in 1:n_loops){

  #generate dummy data for each cell cluster that we want
  dummy_res <- generate_dummy_data_for_scregclust(n_target_genes = Pt,  # Number of target genes
                                                      n_regulator_genes = Pr,  # Number of regulator genes
                                                      n_cells  = n,  # Number of cells
                                                      n_target_gene_clusters  = K,  # Number of target gene clusters
                                                      regulator_mean   = regulator_mean,  # Mean expression of regulator genes
                                                      coefficient_mean = coefficient_mean)  # Mean coefficients in true model, length n_target_gene_clusters)

  #hack to get the list dummy_res into the global environment
  for(iter in 1:length(names(dummy_res))) {
    assign(eval(names(dummy_res)[iter]),dummy_res[[iter]])
  }

  # apply simulated data to scregclust ---------------------------------------

  scregclust(
    expression = rbind(t(Z_t), t(Z_r)),    #scRegClust wants this form
    genesymbols = 1:(Pt+Pr),               #gene row numbers
    is_regulator = (1:(Pt+Pr) > Pt) + 0, #vector indicating which genes are regulators
    n_cl        = K,
    # target_cluster_start = K,
    penalization = max(sapply(seq_along(R[,1]), function(i) sum(R[i,]))) + 1,
    #maximal number of regulators for one cluster
    verbose = FALSE
  )-> scRegOut

  # scRegOut$results
  # rowSums(Pi)
  # rowSums(R)
  # plot(scRegOut)$data

  true_clust_allocation <- apply(X=Pi, MARGIN=2, FUN=function(x) which(x==1))
  predicted_cluster_allocation <- scRegOut$results[[1]]$output[[1]]$cluster[1:Pt]
  rand_index <- aricode::RI(true_clust_allocation, predicted_cluster_allocation)
  print(paste(i_loop, "Rand Index:", rand_index))
  rands[i_loop] <- rand_index
  if(rand_index!=1){
    break()
  }
  # Create an sB matrix that looks like generated dummy data's B
  # A regulator gene can affect all or none target genes. A target gene can only be in one target gene cluster.
  # The coeffs from scregclust are probably normalised so they can probably not be compared to the generated datas B.
  # sB_temp <- scRegOut$results[[1]]$output[[1]]$coeffs
  # cluster <- scRegOut$results[[1]]$output[[1]]$cluster[1:Pt]
  # sB <-  vector("list", length = K)
  # for(i_cell_cluster in 1:K){
  #   sB[[i_cell_cluster]] <- matrix(data = 0, nrow = Pr, ncol=Pt)
  #   inds <- which(cluster==i_cell_cluster)
  #   sB[[i_cell_cluster]][,inds] <- sB_temp[[i_cell_cluster]]
  # }
}
hist(rands, bins=100)

print("Percent rand index=1:", sum(rands==1)/n_loops*100)
