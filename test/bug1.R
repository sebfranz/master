library(scregclust)
set.seed(1)
load("test_data.RData")
# If scregclust fails save all data with
# n_target_gene_clusters <- n_target_gene_clusters[i_cluster]
# save(list=c("local_dat", "cell_data_split", "n_target_genes", "n_regulator_genes", "n_target_gene_clusters", "gene_cluster_start"), file="test_data.RData")
# and send a bug report

scregclust(
  expression             = local_dat,  # "...p rows of genes and n columns of cells."
  split_indices          = cell_data_split,  # Training/test data split indicating by 1 and 2
  genesymbols            = paste0('g', 1:(n_target_genes+n_regulator_genes)),  # Gene row names
  is_regulator           = (1:(n_target_genes+n_regulator_genes) > n_target_genes) + 0,  # Vector indicating which genes are regulators
  n_cl                   = n_target_gene_clusters,
  target_cluster_start   = gene_cluster_start,
  penalization           = 0.14,
  verbose                = TRUE
) -> result

