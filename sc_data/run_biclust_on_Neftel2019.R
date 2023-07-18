rm(list = ls())
library(scregclust)
library(plyr)
library(aricode)  # To calculate rand index
library(Seurat)

execution_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
sapply(list.files(paste0(execution_path,"/../functions/"),recursive = T),
       function(nm) source(paste0(execution_path,"/../functions/", nm)))

set.seed(124)  # This seed crashes due to scregclust producing NULL in all target gene clusters in the only cell cluster left

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
initial_cell_clust <- is_malignant + 1
n_cell_clusters<- length(unique(initial_cell_clust))

# initial_cell_clust <- sample(1:n_cell_clusters, n_cells, replace = T)


# Kod för att flytta 1% av cellerna i varje kluster till ett annat kluster.
disturbed_initial_cell_clust <- initial_cell_clust
disturbed_fraction <- 0.20
for(i_cluster in 1:n_cell_clusters){
  indexes_of_cluster <- which(initial_cell_clust == i_cluster)
  some_of_those_indexes <- sample(indexes_of_cluster, size=as.integer(length(indexes_of_cluster)*disturbed_fraction), replace = F)
  disturbed_initial_cell_clust[some_of_those_indexes] <-
    sample(c(1:n_cell_clusters)[-i_cluster],
           size=length(some_of_those_indexes), replace=T)
}


cell_cluster_history <- cbind(initial_cell_clust, disturbed_initial_cell_clust)
colnames(cell_cluster_history) <- c("True allocation", "Disturbed allocation")

set.seed(1234)
#optionally create a vector of cells to include, if we want TIME to PASS
cols2keep <- sample(ncol(z_g1_ordered), 4500)


# if(!'neftel_biclustfit_group1' %in% list.files(paste0(execution_path,'/r_files/'))){

#thiw might be optional not sure, but ida does it and we keep getting stuck in loops
for (i in 1:length(unique(metaData))){
  ix<-which(metaData==unique(metaData)[i])
  sample_assignment[ix]<-i
}

res <- biclust(max_iter=50,
               initial_cluster_history = cell_cluster_history[cols2keep,],
               is_regulator =  is_predictor,
               n_target_gene_clusters = rep(10,n_cell_clusters),
               penalization_parameter = 0.005,
               train_dat = z_g1_ordered[,cols2keep]
               # sample_assignment = sample_assignment[cols2keep]  can't just send in requires special care
)

  saveRDS(res, file = paste0(execution_path, '/r_files/', "neftel_biclustfit_group1"))
# } else {
#   res <- readRDS(file = paste0(execution_path, '/r_files/', "neftel_biclustfit_group1"))
# }


# res <- biclust(max_iter=50,
#                initial_cluster_history = cell_cluster_history,
#                is_regulator =  is_predictor,
#                n_target_gene_clusters = n_target_gene_clusters,
#                train_dat = z_g1_ordered
# )

plot_cluster_history(cell_cluster_history = res$cell_cluster_history)

