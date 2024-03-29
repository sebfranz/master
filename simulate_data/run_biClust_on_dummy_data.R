rm(list = ls())
if (!require(plyr)) install.packages('plyr')
if (!require(here)) install.packages('aricode')
if (!require(here)) install.packages('here')
library(scregclust)
library(plyr)
library(here)
library(aricode)  # To calculate rand index

root_path <- here::here()
execution_path <- file.path(root_path, "simulate_data")
source(file.path(execution_path,"functions", "generate_dummy_data_for_cell_clustering.R"))
source(file.path(root_path,"functions", "biclust.R"))
source(file.path(root_path,"functions", "plot_cluster_history.R"))

set.seed(1)

# Set variables ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

n_cell_clusters <- 3
n_target_gene_clusters <- c(2,3,4)  # Number of target gene clusters in each cell cluster
n_target_genes <- 50
n_regulator_genes <- 30
n_cells <- c(1000,5000,7000)
regulator_means <- c(1, 5, 10)  # For generating dummy data, regulator mean in each cell cluster
coefficient_means <- list(c(1, 20), c(4, 5, 6), c(10, 12, 14, 16))  # For generating dummy data, coefficient means in each cell cluster
true_cluster_allocation <- rep(1:n_cell_clusters, times=n_cells)
total_n_cells <- sum(n_cells)

# Generate dummy data for each cell cluster that we want ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
res <- generate_dummy_data_for_cell_clustering(
    n_cell_clusters = n_cell_clusters,
    n_target_gene_clusters = n_target_gene_clusters,  # Number of target gene clusters in each cell cluster
    n_target_genes = n_target_genes,
    n_regulator_genes = n_regulator_genes,
    n_cells = n_cells,
    regulator_means = regulator_means,  # For generating dummy data, regulator mean in each cell cluster
    coefficient_means = coefficient_means,  # For generating dummy data, coefficient means in each cell cluster
    disturbed_fraction = 0
)

cell_cluster_history <- cbind(res$initial_cell_clust, res$disturbed_initial_cell_clust)
colnames(cell_cluster_history) <- c("True allocation", "Disturbed allocation")

# install.packages("gplots")
# library("gplots")
# heatmap.2(res$train_dat/max(res$train_dat), scale = "none", col = bluered(100),
#           trace = "none", density.info = "none")

res <- biclust(max_iter=50,
               initial_cluster_history = cell_cluster_history,
               is_regulator = c(rep(0,n_target_genes),rep(1,n_regulator_genes)),
               n_target_gene_clusters = n_target_gene_clusters,
               train_dat = res$train_dat,
               penalization_parameter = 1,
               plot_r2 = FALSE)
  print(res$stats_scregclust)
  print(str(res$stats_scregclust))
  for(i in 1:length(res$stats_scregclust)){
    b <- res$stats_scregclust[[i]]
    if(!is.null(b)){
      cat(paste0("Iteration ", i, "\n"))
      for(ii in 1:length(b)){
        cat(paste0(" Cell cluster ", ii, "\n"))
        a <- b[[ii]]
        if(!is.null(a)){
          cat(paste0("  Total used target genes ", as.character(a["total_target_genes_used"]), "/", as.character(a["total_target_genes"])))
          cat(paste0(". Total used regulators ", as.character(a["total_used_regulators"]), "/", as.character(a["total_regulator_genes"])))
          cat(paste0(". Total used target gene clusters ", as.character(a["total_used_target_gene_clusters"]), "/", as.character(a["total_target_gene_clusters"])))
          cat("\n")
        }
      }
    }
  }

  plot_cluster_history(cell_cluster_history = res$cell_cluster_history, correct_plot=FALSE)
