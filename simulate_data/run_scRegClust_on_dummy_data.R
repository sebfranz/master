#!/usr/bin/Rscript
rm(list = ls())
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

set.seed(14)
Pt = 50      #number of target genes
Pr = 20       #number of regulator genes
n  = 10000   #number of cells
K  = 3       #Number of target gene clusters
regulator_mean   = 1
coefficient_mean = c(1,10,100)
n_loops = 1000


rands <- vector(length = n_loops)
for(i_loop in 1:10){

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
    verbose = TRUE
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
  # if(rand_index!=1){
  #   break()
  # }
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

  b <- scRegOut$results[[1]]$output[[1]]$coeffs
  total_targets <- 0
  total_regulator_vector <- vector(length=Pr)
  total_used_target_gene_clusters <- 0
  for(i in 1:length(b)){
    # print(sum(apply(X = !is.na(b[[i]]), MARGIN = 1, FUN = all)))
    if(is.null(b[[i]])){
      targets <- -
      regulators <- -
    }else{
      targets <- (ncol(b[[i]]))
      regulator_vector <- apply(X = b[[i]]!=0, MARGIN = 1, FUN = all)
      regulators <- sum(regulator_vector)
      total_regulator_vector <- total_regulator_vector + regulator_vector
      total_used_target_gene_clusters <- total_used_target_gene_clusters + 1
    }
    total_targets <- total_targets + targets

    print(paste("Cluster", i, "T", targets, "R", regulators))

  }
  total_used_regulators <- sum(regulator_vector!=0)
  print(paste0("Total used target genes ", total_targets, "/", Pt, ". Total used regulators ", total_used_regulators, "/", Pr, ". Total used target gene clusters ", total_used_target_gene_clusters, "/", length(b)))
}



# hist(rands, bins=100)

# print("Percent rand index=1:", sum(rands==1)/n_loops*100)
