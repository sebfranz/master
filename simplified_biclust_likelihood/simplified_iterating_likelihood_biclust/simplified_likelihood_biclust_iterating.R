library(tidyverse)
library(aricode)  # To calculate rand index
library(ggplot2)
library(ggalluvial)
library(reshape2)
library(here)  # To work with paths

path_root <- here::here()
execution_path <- here::here("simplified_biclust_likelihood",
                             "simplified_iterating_likelihood_biclust")
output_path <- execution_path
function_folder <- here::here("functions")
all_function_files <- list.files(function_folder, recursive=T, full.names=T)
for(current_file in all_function_files){
    print(paste("Loading", current_file))
    source(current_file)
  }

# Simulate appropriate data from linear model

set.seed(1234)

generate_data <- function(){
  num_cells <- 100
  regulator_expression <- rnorm(num_cells, mean = 1, sd = 0.1)
  # hist(regulator_expression)

  n_cell_clusters    <- 2

  # Sample intercepts
  intercepts <- rnorm(n_cell_clusters, mean = 10, sd = 5)  * 0 #can just set to zero without loss of generality as data will be

  # Sample slopes
  betas      <- rnorm(n_cell_clusters, mean = 0,  sd =  5)

  # target-gene variance per target gene cluster
  sigmas     <- rnorm(n_cell_clusters, mean = 1,  sd = 0.1)


  # Assign half of the cells to be in either cluster
  true_cell_clust <- c(rep(1,num_cells/2), rep(2,num_cells/2))

  # Build cell expression from linear model
  # | intercept, beta | x | r11, r12 |
  #                       | r21, r22 |
  #                       |    ...   |
  # For every cell cluster, 1x2 * 2x50

  cell_cluster_target_gene_expression <- sapply(1:n_cell_clusters, function(cellCluster){
    intercept_and_coefficient <- c(intercepts[cellCluster], betas[cellCluster])
    design_matrix <- t(cbind(rep(1, num_cells/2),  regulator_expression[which(true_cell_clust == cellCluster)]))
    residuals     <- rnorm(num_cells/2, 0, sd = sqrt(sigmas[cellCluster]))

    intercept_and_coefficient %*% design_matrix + residuals
  }
  )

  target_expression <- c(cell_cluster_target_gene_expression[,1],
                         cell_cluster_target_gene_expression[,2])
  dat <- tibble(cell_id = (1:num_cells), true_cell_clust,
               target_expression, regulator_expression)
  dat
}

dat <- generate_data()

ind_targetgenes <- c(3)

ind_reggenes <- c(4)

#plot original data

p <- ggplot(dat,
            aes(x = regulator_expression,
                y = target_expression,
                color=factor(true_cell_clust)))
png(file.path(execution_path,"scatterplot_data.png"))
p + geom_point()
dev.off()

#randomise cluster labels
randomise_cluster_labels <- function(cluster_labels = dat$true_cell_clust,
                                     fraction_randomised = 0.20){

  disturbed_initial_cell_clust <- cluster_labels
  n_cell_clusters <- length(unique(cluster_labels))

  for(i_cluster in 1:n_cell_clusters){
    indexes_of_cluster <- which(cluster_labels == i_cluster)
    some_of_those_indexes <- sample(indexes_of_cluster,
                                    size=as.integer(length(indexes_of_cluster)*fraction_randomised),
                                    replace = F)
    disturbed_initial_cell_clust[some_of_those_indexes] <-
      sample(c(1:n_cell_clusters)[-i_cluster],
             size=length(some_of_those_indexes), replace=T)
  }
  return(disturbed_initial_cell_clust)
}

disturbed_initial_cell_clust <- randomise_cluster_labels(cluster_labels = dat$true_cell_clust)

cell_cluster_history <- tibble(dat$cell_id,dat$true_cell_clust, disturbed_initial_cell_clust)
colnames(cell_cluster_history) <- c("cell_id", "True allocation", "Disturbed allocation")

#Begin iteration etc

stop_iterating_flag <- 0 #flag if we have converged
# i_main <- 1  # Main iteration variable, just set it to one for now.
max_iter <- 40  # This one would usually cap number of iterations,

# Find initial cluster labels
initial_clustering <- disturbed_initial_cell_clust
n_target_gene_clusters <- 1  # We are not clustering target genes for now and we only have one target gene


# Set up some variables
n_cell_clusters = length(unique(initial_clustering))
n_target_genes = 1
n_regulator_genes = 1


# Preallocate cluster history
initial_column_padding <- ncol(as.matrix(initial_clustering)) + 1  # +1 Because we have an index column that is not an index column it's an ID column
cell_cluster_history <- data.frame(matrix(NA, nrow = nrow(as.matrix(initial_clustering)), ncol = max_iter + initial_column_padding))
colnames(cell_cluster_history) <- c("Cell ID", 'Initial clustering', paste0("Iteration ", 1:max_iter))
cell_cluster_history[, 'Cell ID'] <-1:length(initial_clustering)  # Set cell names
cell_cluster_history[, 'Initial clustering'] <- initial_clustering
cell_cluster_history <- as_tibble(cell_cluster_history)

# Pre-allocate all r2 matrices for later analysis if feasible
likelihood_all <- vector("list", length = max_iter)

# Set the current cell clustering
current_cell_cluster_allocation <- initial_clustering

for(i_main in (1:max_iter)){

# Fit model to each cell cluster
models <- vector("list", length = n_cell_clusters)

for(cell_cluster in 1:n_cell_clusters){
  current_rows <- which(current_cell_cluster_allocation == cell_cluster)
  models[[cell_cluster]] <- lm(formula = 'target_expression ~ 0 + regulator_expression',
                               data = dat,
                               subset = current_rows)
}

# For all cells, calculate the likelihood of coming from the model corresponding to each
like <- matrix(0, nrow = nrow(dat), ncol = n_cell_clusters)

penalization_LAMBDA <-  0

# Calculate the residual target gene variance for each gene and cluster
# (so just one gene).

#preallocate residual variance estimates
TARGET_GENE_RESIDUAL_var <-  matrix(0, nrow = n_target_genes, ncol = n_cell_clusters)
# dat is dat <- cbind(target_expression, regulator_expression), e.g. a 2x100, with e.g. the first 50 rows being true cell cluster 1
# 100x2 * 2x1

  for(cell_cluster in 1:n_cell_clusters){


    intercept_and_regulatorgenes <- as.matrix(cbind(rep(1,nrow(dat)),dat[,ind_reggenes]))  # 100x2
    current_rows <- which(current_cell_cluster_allocation == cell_cluster)

    if(length(models[[cell_cluster]]$coefficients) == 1){
      PREDICTED_VALUES <- (as.matrix(intercept_and_regulatorgenes[current_rows, -1]) %*%
                             models[[cell_cluster]]$coefficients)
    }else{
      PREDICTED_VALUES <- (intercept_and_regulatorgenes[current_rows,] %*%
                             models[[cell_cluster]]$coefficients)
    }

    residuals <- dat[current_rows,ind_targetgenes] - PREDICTED_VALUES

    TARGET_GENE_RESIDUAL_var[cell_cluster]  <- var(residuals)
  }

  # Now to actually calculate predicted or 'predicted' r2
  for(cell in 1:nrow(dat)){
    for(cell_cluster in 1:n_cell_clusters){
      # Bug fix hack: remove NA coefficients
      if(any(is.na(models[[cell_cluster]]$coefficients))){
        NA_coeffs <-  unname(which(is.na(models[[cell_cluster]]$coefficients)))
        S_ERR <- (dat[cell,ind_targetgenes] - as.vector(c(1,dat[cell,c(-1, -NA_coeffs)])) %*% models[[cell_cluster]]$coefficients[-NA_coeffs])^2
      }

      if(length(models[[cell_cluster]]$coefficients) == 1){
        PREDICTED_VALUE <- as.matrix(dat[cell,ind_reggenes]) %*% models[[cell_cluster]]$coefficients
      }else{
        PREDICTED_VALUE <- as.vector(c(1,dat[cell,ind_reggenes])) %*% models[[cell_cluster]]$coefficients
      }

      if(length(models[[cell_cluster]]$coefficients) == 1){
        coefficient_vector_1norm <- sum(abs(models[[cell_cluster]]$coefficients))
        }else{
         coefficient_vector_1norm <- sum(abs(models[[cell_cluster]]$coefficients[-1]))
      }

      OBSERVED_VALUE <- dat[cell,ind_targetgenes]

      SQUARED_ERROR <- (OBSERVED_VALUE - PREDICTED_VALUE)^2

      penalization <- penalization_LAMBDA * coefficient_vector_1norm/TARGET_GENE_RESIDUAL_var[cell_cluster]

      # like[cell,cell_cluster] <- SQUARED_ERROR / 2 / TARGET_GENE_RESIDUAL_var[cell_cluster] - penalization #negative penalty as higher likelihood is better

      #here we are optimizing the penalized NEGATIVE likelyhood,
      #so penalty is positive

      like[cell,cell_cluster] <- as.numeric(log(TARGET_GENE_RESIDUAL_var[cell_cluster])/2 +
                                  SQUARED_ERROR/2/TARGET_GENE_RESIDUAL_var[cell_cluster] +
                                  penalization)

    }
  }

  likelihood_all[[i_main]] <- like

  # r2plot(iteration = i_main,
  #      r2 = like,
  #      prev_cell_clust = dat$true_cell_clust)

  p <- ggplot(data = as.tibble(cbind(like, dat$true_cell_clust)),
       aes(x = V1, y = V2, color = as.factor(dat$true_cell_clust))
      ) +
      geom_point() +
        geom_abline(intercept = 0, slope = 1)+
        labs(x = "Log-likelihood for fitting into cluster 1", y = "Log-likelihood for fitting into cluster 2")

  png(file.path(execution_path,paste0("Decision_line","_lambda_",
                                      round(penalization_LAMBDA, 0.1),".png")))
  p + labs(color = "True cell cluster")
  dev.off()

# Update cluster allocations
updated_cell_clust <-  sapply(1:nrow(like), function(row) which.min(like[row,]))

# Update data in cell_cluster_history
cell_cluster_history[, i_main + initial_column_padding] <- updated_cell_clust
RI(dat$true_cell_clust,updated_cell_clust)

#check convergence of cluster labels
# Compare clusters with with previous iterations so we can exit if we seen this allocation before
for(prev_clustering in ((i_main-1):0) ){
  print(paste0('Comparing with iteration ', prev_clustering))
  RAND_INDEX <- RI(updated_cell_clust,
                   as.matrix(cell_cluster_history)[,prev_clustering + initial_column_padding]
  )
    if(RAND_INDEX == 1){
    print("Cell clustering from iteration same as some previous iteration. Exiting.")
    print(paste0("RI of ", RAND_INDEX,
                 " when comparing iteration ", i_main,
                 " to iteration ", prev_clustering))
    stop_iterating_flag <- T
    break
  }
}


if(stop_iterating_flag == T){
  # Clean up cluster history
  cell_cluster_history <- cell_cluster_history[ , colSums(is.na(cell_cluster_history))==0, drop=FALSE]
  # Stop iterations/exit function
  break
}

}
cell_cluster_history_plotting <-cbind('Cell ID' = cell_cluster_history[,1],
                             dat$true_cell_clust,
                             cell_cluster_history[,c(2,3,4)])
png(file.path(execution_path,paste0("Alluvial_diag_","_lambda_",
                                    round(penalization_LAMBDA, 0.1),".png")))
plot_cluster_history(cell_cluster_history = cell_cluster_history_plotting)
dev.off()
