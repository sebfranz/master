#!/usr/bin/Rscript
rm(list = ls())

library(tidyverse)
library(aricode)  # To calculate rand index
library(ggplot2)  # To plot things #TODO: What things?
library(ggalluvial)  # To plot things #TODO: What things?
library(reshape2)
library(here)  # To work with paths
library(ggfortify)  # For pca-plot

# Get absolute path where script is located, by using relative paths.
path_root <- here::here()
execution_path <- here::here("simplified_biclust_likelihood",
                             "simplified_iterating_likelihood_biclust")
output_path <- execution_path
function_folder <- here::here("functions")
all_function_files <- list.files(function_folder, recursive = T, full.names = T)
for (current_file in all_function_files) {
  print(paste("Loading", current_file))
  source(current_file)
}
source(file.path(execution_path, "function_generate_data_lm.R"))

# Set seed for example
set.seed(1234)


dat <- generate_data_lm(n_cell_clusters = 3,
                        n_target_gene_type = 4,  # We have x named target genes that have one expression per cell
                        n_regulator_gene_type = 5,  # We have x named regulator genes that have one expression per cell
                        n_cells = c(1000, 5000, 10000),
                        regulator_means = c(1, 2, 5),  # Regulator mean expression in each cell cluster.
                        regulator_standard_deviations = c(0.1, 0.2, 0.3),  # Regulator sd for expression in each cell cluster.
                        coefficients_standard_deviation = 100, # 'The betas/slopes'. One per target gene. Instead of providing mean and standard deviation for each target gene, we provide the standard deviation from which these will be generated. Mean will be 0.
                        target_gene_type_standard_deviation = 3
)

# We assume target genes start with t then a number. And likewise for regulator genes.
n_total_cells <- nrow(dat)
ind_targetgenes <- which(str_detect(colnames(dat), "t\\d"))
ind_reggenes <- which(str_detect(colnames(dat), "r\\d"))


# Plot original data ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(reshape2)
d <- reshape2::melt(dat[, 2:ncol(dat)], id.vars = "true_cell_cluster_allocation")

if (length(ind_targetgenes) == 1 && length(ind_reggenes) == 1) {
  p <- ggplot(dat,
              aes(x = r1,
                  y = t1,
                  color = true_cell_cluster_allocation)
  )
  png(file.path(execution_path, "scatterplot_data.png"))
  p + geom_point()
  dev.off()
}


# Randomise cluster labels ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
randomise_cluster_labels <- function(cluster_labels = dat$true_cell_cluster_allocation,
                                     fraction_randomised = 0.05) {

  disturbed_initial_cell_clust <- cluster_labels
  n_cell_clusters <- length(unique(cluster_labels))

  for (i_cluster in 1:n_cell_clusters) {
    indexes_of_cluster <- which(cluster_labels == i_cluster)
    some_of_those_indexes <- sample(indexes_of_cluster,
                                    size = as.integer(length(indexes_of_cluster) * fraction_randomised),
                                    replace = F)
    disturbed_initial_cell_clust[some_of_those_indexes] <-
      sample((1:n_cell_clusters)[-i_cluster],
             size = length(some_of_those_indexes), replace = T)
  }
  return(disturbed_initial_cell_clust)
}

disturbed_initial_cell_clust <- randomise_cluster_labels(cluster_labels = dat$true_cell_cluster_allocation)

cell_cluster_history <- tibble::tibble(dat$cell_id, dat$true_cell_cluster_allocation, disturbed_initial_cell_clust)
colnames(cell_cluster_history) <- c("cell_id", "True allocation", "Disturbed allocation")


# Begin iteration etc
stop_iterating_flag <- 0  # Flag if we have converged
# i_main <- 1  # Main iteration variable, just set it to one for now.
max_iter <- 40  # This one would usually cap number of iterations,

# Find initial cluster labels
initial_clustering <- disturbed_initial_cell_clust
n_target_gene_clusters <- 1  # We are not clustering target genes for now and we only have one target gene


# Set up some variables
n_cell_clusters <- length(unique(initial_clustering))
n_target_genes <- length(ind_targetgenes)
n_regulator_genes <- length(ind_reggenes)


# Preallocate cluster history
initial_column_padding <- ncol(as.matrix(initial_clustering)) + 1  # +1 Because we have an index column that is not an index column it's an ID column
cell_cluster_history <- data.frame(matrix(NA, nrow = nrow(as.matrix(initial_clustering)), ncol = max_iter + initial_column_padding))
colnames(cell_cluster_history) <- c("Cell ID", 'Initial clustering', paste0("Iteration ", 1:max_iter))
cell_cluster_history[, 'Cell ID'] <- seq_along(initial_clustering)  # Set cell names
cell_cluster_history[, 'Initial clustering'] <- initial_clustering
cell_cluster_history <- tibble::as_tibble(cell_cluster_history)

# Pre-allocate all r2 matrices for later analysis if feasible
likelihood_all <- vector(mode = "list", length = max_iter)

# Set the current cell clustering
current_cell_cluster_allocation <- initial_clustering

for (i_main in 1:max_iter) {

  # Fit model to each cell cluster
  models <- vector(mode = "list", length = n_cell_clusters)

  for (cell_cluster in 1:n_cell_clusters) {
    cell_cluster_rows <- which(current_cell_cluster_allocation == cell_cluster)
    cell_cluster_target_genes <- as.matrix(dat[cell_cluster_rows, ind_targetgenes])
    cell_cluster_regulator_genes <- as.matrix(dat[cell_cluster_rows, ind_reggenes])
    models[[cell_cluster]] <- lm(formula = 'cell_cluster_target_genes ~ 0 + cell_cluster_regulator_genes')
  }

  # For all cells, calculate the likelihood of coming from the model corresponding to each
  likelihood <- matrix(data = 0, nrow = nrow(dat), ncol = n_cell_clusters)

  penalization_lambda <- 0.5

  # Calculate the residual target gene variance for each gene and cluster
  # (so just one gene).

  # Pre-allocate residual variance estimates
  target_genes_residual_var <- matrix(data = 0, nrow = n_cell_clusters, ncol = n_target_genes)
  # dat is dat <- cbind(target_expression, regulator_expression), e.g. a 2x100, with e.g. the first 50 rows being true cell cluster 1
  # 100x2 * 2x1

  # TODO: This calculates one variance value for each and every target gene type for every cell cluster. Is that correct?
  # Output: n_cell_clusters x n_target_genes
  for (cell_cluster in 1:n_cell_clusters) {
    current_rows <- which(current_cell_cluster_allocation == cell_cluster)
    current_regulator_genes <- as.matrix(dat[current_rows, ind_reggenes])
    current_target_genes <- as.matrix(dat[current_rows, ind_targetgenes])
    cell_cluster_betas <- models[[cell_cluster]]$coefficients

    predicted_values <- current_regulator_genes %*% cell_cluster_betas

    residuals <- current_target_genes - predicted_values

    target_genes_residual_var[cell_cluster,] <- diag(var(residuals))
  }

  # Now to actually calculate predicted or 'predicted' r2
  for (cell in seq_len(nrow(dat))) {
    for (cell_cluster in seq_len(n_cell_clusters)) {
      cell_regulator_genes <- as.matrix(dat[cell, ind_reggenes])
      cell_cluster_betas <- models[[cell_cluster]]$coefficients
      observed_value <- as.matrix(dat[cell, ind_targetgenes])
      cell_cluster_target_genes_residual_var <- target_genes_residual_var[cell_cluster, , drop = FALSE]
      # Bug fix hack: remove NA coefficients
      # if (any(is.na(current_betas))) {
      #   NA_coeffs <- unname(which(is.na(current_betas)))
      #   S_ERR <- (dat[cell, ind_targetgenes] - as.vector(c(1, dat[cell, c(-1, -NA_coeffs)])) %*% current_betas[-NA_coeffs])^2
      # }

      predicted_value <- cell_regulator_genes %*% cell_cluster_betas

      cell_cluster_betas_vector_1norm <- sum(abs(cell_cluster_betas))

      cell_squared_error <- (observed_value - predicted_value)^2

      penalization <- penalization_lambda * cell_cluster_betas_vector_1norm / cell_cluster_target_genes_residual_var

      # TODO: Figure out what the formula should be, and make sure you have which.min or which.max that is correct later around line 220.
      # likelihood[cell,cell_cluster] <- squared_error / 2 / target_genes_residual_var[cell_cluster] - penalization #negative penalty as higher likelihood is better

      # Here we are optimizing the penalized NEGATIVE likelyhood, so penalty is positive
      temp_likelihood <- as.numeric(log(cell_cluster_target_genes_residual_var) / 2 +
                                      cell_squared_error / (2 * target_genes_residual_var[cell_cluster]) +
                                      penalization)
      likelihood[cell, cell_cluster] <- sum(temp_likelihood)

    }
  }

  likelihood_all[[i_main]] <- likelihood

  # r2plot(iteration = i_main,
  #      r2 = likelihood,
  #      prev_cell_clust = dat$true_cell_cluster_allocation)

  # Scatter-plot the log-liklihood on each axis, color with true allocation
  # If more than 2 dim make pca plot
  true_cell_cluster_allocation_vector <- paste("Cluster", pull(dat, var = 'true_cell_cluster_allocation'))  # These needs to be strings for discrete labels in pca plot
  colnames(likelihood) <- paste("Likelihood cell cluster", seq_len(n_cell_clusters))
  likelihood_tibble <- tibble::as_tibble(likelihood)
  data_for_plotting <- tibble::tibble(likelihood_tibble, true_cell_cluster_allocation = true_cell_cluster_allocation_vector)

  filename_plot <- paste0("Decision_line_lambda_", round(penalization_lambda, digits = 2), "_iteration_", i_main, ".png")
  if (ncol(likelihood) == 2) {
    p <- ggplot2::ggplot(data = data_for_plotting, ggplot2::aes(x = "Likelihood cell cluster 1", y = "Likelihood cell cluster 2", color = true_cell_cluster_allocation)) +
      ggplot2::geom_point() +
      ggplot2::geom_abline(intercept = 0, slope = 1) +
      ggplot2::labs(x = "Log-likelihood for fitting into cluster 1", y = "Log-likelihood for fitting into cluster 2")
    png(file.path(execution_path, filename_plot))
    p + ggplot2::labs(color = "True cell cluster")
    dev.off()
  } else {
    pca_res <- prcomp(data_for_plotting[, seq_len(ncol(data_for_plotting) - 1)], scale. = TRUE)
    p <- ggplot2::autoplot(pca_res, data = data_for_plotting, colour = 'true_cell_cluster_allocation')
    png(file.path(execution_path, filename_plot))
    plot(p)
    dev.off()
  }

  # true_cell_cluster_allocation_vector <- pull(dat, var = 'true_cell_cluster_allocation')
  # likelihood_tibble['cell_id'] <- seq_len(nrow(likelihood_tibble))
  # for (cell_cluster in seq_len(n_cell_clusters)) {
  #   cell_cluster_rows <- which(true_cell_cluster_allocation_vector == cell_cluster)
  #
  #   cell_cluster_likelihood <- likelihood_tibble[cell_cluster_rows, ]
  #   # data_for_plotting <- tibble::tibble(cell_cluster_likelihood, true_cell_cluster_allocation = true_cell_cluster_allocation_vector)
  #   cell_cluster_likelihood <- reshape2::melt(cell_cluster_likelihood, id.vars = "cell_id")
  #   # Change histogram plot line colors by groups
  #   ggplot(cell_cluster_likelihood, aes(x = cell_id, color = variable)) +
  #     geom_histogram(fill = "white") +
  #     stat_bin(bins=100)
  #   # Overlaid histograms
  #   ggplot(df, aes(x = cell_id, color = variable)) +
  #     geom_histogram(fill = "white", alpha = 0.5, position = "identity")
  # }


  # Update cluster allocations
  updated_cell_clust <- sapply(seq_len(nrow(likelihood)), function(row) which.min(likelihood[row,]))

  # Update data in cell_cluster_history
  cell_cluster_history[, i_main + initial_column_padding] <- updated_cell_clust
  aricode::RI(dat$true_cell_cluster_allocation, updated_cell_clust)

  # Check convergence of cluster labels
  # Compare clusters with with previous iterations so we can exit if we seen this allocation before
  for (prev_clustering in ((i_main - 1):0)) {
    print(paste0('Comparing with iteration ', prev_clustering))
    rand_index <- aricode::RI(updated_cell_clust,
                              as.matrix(cell_cluster_history)[, prev_clustering + initial_column_padding]
    )
    if (rand_index == 1) {
      print("Cell clustering from iteration same as some previous iteration. Exiting.")
      print(paste0("RI of ", rand_index,
                   " when comparing iteration ", i_main,
                   " to iteration ", prev_clustering))
      stop_iterating_flag <- T
      break
    }
  }


  if (stop_iterating_flag) {
    # Clean up cluster history
    cell_cluster_history <- cell_cluster_history[, colSums(is.na(cell_cluster_history)) == 0, drop = FALSE]
    # Stop iterations/exit function
    break
  }

}
cell_cluster_history_plotting <- cbind('Cell ID' = cell_cluster_history[, 1],
                                       dat$true_cell_cluster_allocation,
                                       cell_cluster_history[, c(2, 3, 4)])
png(file.path(execution_path, paste0("Alluvial_diag_", "_lambda_",
                                     round(penalization_lambda, 0.1), ".png")))
plot_cluster_history(cell_cluster_history = cell_cluster_history_plotting)
dev.off()
