rm(list = ls())
library(scregclust)
library(plyr)
library(aricode)  # To calculate rand index
library(Seurat)

execution_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
sapply(list.files(paste0(execution_path,"/../functions/"),recursive = T),
       function(nm) source(paste0(execution_path,"/../functions/", nm)))


source(paste0(execution_path,"/functions/generate_scregclust_data_from_sampled_regulators.R"))
source(paste0(execution_path,"/functions/generate_biclust_data_from_sampled_regulators.R"))

set.seed(124)  # This seed crashes due to scregclust producing NULL in all target gene clusters in the only cell cluster left

# Set variables ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

n_target_genes = 100
n_cell_clusters = 2
n_target_gene_clusters = c(3,4)  # Number of target gene clusters in each cell cluster
regulator_expression_offset =  c(0,0)
coefficient_means = list(c(1,1,1), c(1,1,1,1))  # For generating dummy data, coefficient means in each cell clustertrue_cluster_allocation = rep(1:n_cell_clusters, times=n_cells)
# total_n_cells = sum(n_cells)

# Generate dummy data for each cell cluster that we want ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#load sampled data
execution_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
neftel_path <- paste0(execution_path, '/../../datasets_sctargettranslator/Neftel2019')
Neftel_g1 <- readRDS(file = paste0(neftel_path, '/r_files/', "neftel_seurat_group1"))

Neftel_g1<-SetIdent(Neftel_g1, value='malignant')
Neftel_g1_malignant<-subset(Neftel_g1, idents='yes')
Neftel_g1_non_malignant<-subset(Neftel_g1, idents='non_cancer')
rm(Neftel_g1)

#get malignant expression on proper form
z_g1<-GetAssayData(Neftel_g1_malignant, slot='scale.data')
# str(z_g1)
metaData<-Neftel_g1_malignant@meta.data$sample

out<-scregclust_format(z_g1)

genesymbols<-out[[1]]
sample_assignment<-out[[2]]
is_predictor<-out[[3]]
sum(is_predictor==1)

# z_g1 <- z_g1[order(is_predictor),]
# is_predictor <- is_predictor[order(is_predictor)]
malignant_regulator_expression <- t(z_g1[which(is_predictor == 1), ])
rm(z_g1)
rm(Neftel_g1_malignant)


#get nonmalignant expression
z_g1<-GetAssayData(Neftel_g1_non_malignant, slot='scale.data')
metaData<-Neftel_g1_non_malignant@meta.data$sample

out<-scregclust_format(z_g1)

genesymbols<-out[[1]]
sample_assignment<-out[[2]]
is_predictor<-out[[3]]
sum(is_predictor == 1)
# z_g1 <- z_g1[order(is_predictor),]
# is_predictor <- is_predictor[order(is_predictor)]
non_malignant_regulator_expression <- t(z_g1[which(is_predictor == 1), ])
rm(z_g1)
rm(Neftel_g1_non_malignant)


dim_(malignant_regulator_expression)

dim_(non_malignant_regulator_expression)

# res <- generate_biclust_data_from_sampled_regulators(
#   n_target_genes = n_target_genes,
#   n_cell_clusters = n_cell_clusters,
#   n_target_gene_clusters = n_target_gene_clusters,  # Number of target gene clusters in each cell cluster
#   regulator_expression_offset =  regulator_expression_offset,
#   n_cells =n_cells,
#   coefficient_means = coefficient_means,  # For generating dummy data, coefficient means in each cell cluster
#   regulator_expression = regulator_expression
# )

res1 <- generate_scregclust_data_from_sampled_regulators(
    n_target_genes = n_target_genes,              #number of target genes
    regulator_expression_offset = 0,
    n_target_gene_clusters  = 2,               #Number of target gene clusters
    # coefficient_mean = c(1,2), #mean coefficients in true  model, length n_target_gene_clusters
    coefficient_min = c(0.01,0.01),
    coefficient_max = c(0.1,0.1),
    regulator_expression = malignant_regulator_expression
)

res2 <- generate_scregclust_data_from_sampled_regulators(
  n_target_genes = n_target_genes,              #number of target genes
  regulator_expression_offset = 0,
  n_target_gene_clusters  = 2,               #Number of target gene clusters
  # coefficient_mean = c(1,2), #mean coefficients in true  model, length n_target_gene_clusters
  coefficient_min = c(0.01,0.01),
  coefficient_max = c(0.1,0.1),
  regulator_expression = non_malignant_regulator_expression
)

# is_regulator = c(rep(0, n_target_genes),rep(1,ncol(malignant_regulator_expression)))
is_regulator =  c(rep(0, ncol(res1$Z_t)),rep(1,ncol(res1$Z_r)))

# up to c(nrow(res1$Z_t), nrow(res2$Z_t))
cell_subsample <- c(1000, 1000)

train_dat <- rbind(cbind(res1$Z_t,res1$Z_r)[1:cell_subsample[1],],
                   cbind(res2$Z_t,res2$Z_r)[1:cell_subsample[2],])


# Get initial cell clustering
initial_cell_clust <- c(rep(1, cell_subsample[1]),rep(2,cell_subsample[2]))
n_cell_clusters<- 2
# initial_cell_clust <- sample(1:n_cell_clusters, n_cells, replace = T)


# Kod för att flytta 1% av cellerna i varje kluster till ett annat kluster.
disturbed_initial_cell_clust <- initial_cell_clust
disturbed_fraction <- 0.05
for(i_cluster in 1:n_cell_clusters){
  indexes_of_cluster <- which(initial_cell_clust == i_cluster)
  some_of_those_indexes <- sample(indexes_of_cluster, size=as.integer(length(indexes_of_cluster)*disturbed_fraction), replace = F)
  disturbed_initial_cell_clust[some_of_those_indexes] <-
    sample(c(1:n_cell_clusters)[-i_cluster],
           size=length(some_of_those_indexes), replace=T)
}


cell_cluster_history <- cbind(initial_cell_clust, disturbed_initial_cell_clust)
colnames(cell_cluster_history) <- c("True allocation", "Disturbed allocation")

res <- biclust(max_iter=50,
               initial_cluster_history = cell_cluster_history,
               is_regulator =  is_regulator,
               n_target_gene_clusters = n_target_gene_clusters,
               n_cells = nrow(train_dat),
               train_dat = t(train_dat)
               )

plot_cluster_history(cell_cluster_history = res$cell_cluster_history)
