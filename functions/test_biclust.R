#!/usr/bin/Rscript
# This file is to showcase how to write unit tests.
if (!require(testthat)) install.packages('testthat')
library(testthat)

# This is to setup testthat version 3e. It's opt in.
# There is an automatic command to do it, but it is expected to run
# inside of a package. So it throws an error but works anyway.
# Because it throws an error I have it surrounded by a tryCatch.
tryCatch(
  expr = {
    usethis::use_testthat(3)  # Use testthat version 3e
  },
  error = function(e){
    # (Optional)
    # Do this if an error is caught...
  },
  warning = function(w){
    # (Optional)
    # Do this if an warning is caught...
  },
  finally = {
    # (Optional)
    # Do this at the end before quitting the tryCatch structure...
  }
)


# The tests -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


filename <- "./functions/biclust.R"
if(!file.exists(filename)){
  filename <- "biclust.R"
}
source(filename)

test_that("Check if list is returned", {

  expect_length(1,1)
#
#   n_target_gene_clusters <- c(1)  # c(2, 3)  # length is n_cell_clusters
#   total_n_cells <- 2*10  # 16000
#   n_cell_clusters <- 2  # 2
#   n_target_genes <- 10  # 30
#   n_regulator_genes <- 5  # 20
#   # Assign half of the cells to be in either cluster
#   prev_cell_clust <- c(rep(1,total_n_cells/2), rep(2,total_n_cells/2))
#   # train_dat <- num [1:50, 1:16000]
#
#   # Generate data
#   set.seed(1234)
#   regulator_expression <- rnorm(total_n_cells, mean = 1, sd = 0.1)
#   # Sample intercepts
#   intercepts <- rnorm(n_cell_clusters, mean = 10, sd = 5)
#   # Sample slopes
#   betas      <- rnorm(n_cell_clusters, mean = 1,  sd = 0.5)
#   # Build cell expression from linear model
#   cell_cluster_expression <- sapply(1:n_cell_clusters, function(cellCluster) c(intercepts[cellCluster], betas[cellCluster]) %*%
#                                       t(cbind(rep(1, total_n_cells/2), regulator_expression[which(prev_cell_clust == cellCluster)]))
#   )
#
#   target_expression <- c(cell_cluster_expression[,1], cell_cluster_expression[,2])
#   dat <- cbind(target_expression, regulator_expression)
#
#   models <- vector("list", length = n_cell_clusters)
#   dat <- cbind(dat, NA)  # Add columns so we can have colors for clusters in the plot
#   colnames(dat)[3] <- 'cluster'  # Name that column
#   for(i_cell_cluster in 1:n_cell_clusters){
#     current_rows <- which(prev_cell_clust == i_cell_cluster)
#     dat[current_rows, 'cluster'] <- i_cell_cluster  # Again, so we can plot the colors for each cluster
#     models[[i_cell_cluster]] <- lm(dat[current_rows, 'target_expression'] ~ dat[current_rows,'regulator_expression'])
#   }
#
#   # coeffs is a list of n_target_gene_clusters
#   # coeffs is n_regulator_genes rows and n_target_genes cols (n_target_genes cols is split across different clusters/coeffs-elements)
#   coeffs_a <- matrix(rnorm(n_regulator_genes*n_target_genes, mean = 1, sd = 0.1), nrow=n_regulator_genes, ncol=n_target_genes)
#   coeffs_b <- matrix(rnorm(n_regulator_genes*n_target_genes, mean = 1, sd = 0.1), nrow=n_regulator_genes, ncol=n_target_genes)
#   i_target_gene_cluster <- 1
#   i_target_cell_cluster <- 1
#   out_list[[i_target_cell_cluster]]$results[[1]]$output[[1]]$coeffs[[i_target_gene_cluster]] <- coeffs_a
#   i_target_gene_cluster <- 1
#   i_target_cell_cluster <- 2
#   out_list[[i_target_cell_cluster]]$results[[1]]$output[[1]]$coeffs[[i_target_gene_cluster]] <- coeffs_b
#
#   R2 <- calculate_R2(n_target_gene_clusters,
#                      total_n_cells,
#                      n_cell_clusters,
#                      n_target_genes,
#                      n_regulator_genes,
#                      prev_cell_clust,
#                      out_list,
#                      train_dat)
#
#   # Test that the result is the correct value
#   expect_length(R2, 5)
})



