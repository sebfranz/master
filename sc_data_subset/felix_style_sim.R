rm(list = ls())
execution_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(execution_path)
# neftel_path <- paste0(execution_path, '/../../datasets_sctargettranslator/Neftel2019')
# setwd(neftel_path)
set.seed(1)

library(scregclust)
library('Seurat')
library(Matrix)



felix_style_sim <- function(mn_g1,
                            number_regulators = NA, #only applicable if NOT including sel_corr_idx
                            num_targets = 200,
                            num_clusters,
                            group_mean_min = 0.01,
                            group_mean_max = 0.1,
                            sigma_of_betas = 0.1,
                            sel_corr_idx = NA#optional selection of regulator genes, important if calling several times and planning to append (union)
                            ){
  # Assume there are 5 true clusters in the data
  n_cl <- num_clusters

  if(all(!is.na(sel_corr_idx))){
    mn_g1 <- t(mn_g1[sel_corr_idx, ]) #picks out genes from rows and transposes
  }else{

    #number of regulators to keep
    n_reg <- number_regulators

    sel_corr_idx <- sample(nrow(mn_g1), n_reg)

    mn_g1 <- t(mn_g1[sel_corr_idx, ]) #picks out genes from rows and transposes
  }
  # Each cluster is connected to a set of regulators and respective signs

  # Select 5-15% of active regulators in each cluster to keep
  # the models sparse
  reg_idx <- lapply(seq_len(n_cl), function(i) {
    sort(sample(n_reg, sample(seq(0.05 * n_reg, 0.15 * n_reg), 1)))
  })

  # Choose signs for each cluster
  signs_og <- lapply(
    reg_idx, function(idx) sample(c(-1, 1), length(idx), replace = TRUE)
  )

  # Simulate group means on a uniform scale
  grp_means <- lapply(reg_idx, function(idx) runif(length(idx),
                                                   group_mean_min, group_mean_max))

  # Simulate cluster sizes (make them all equal size for now)
  n_target <- num_targets
  k_true <- rep(seq_len(n_cl), each = n_target / n_cl)

  # Generate coefficients for each target gene
  sigma_beta <- sigma_of_betas
  beta_og <- mapply(function(i, m_vec, s_vec) {
    # beta_og_pre <- mapply(function(i, m_vec, s_vec) {
    t(mapply(
      function(m, s) s * rnorm(sum(k_true == i), m, sigma_beta), m_vec, s_vec
    ))
  }, seq_len(n_cl), grp_means, signs_og, SIMPLIFY = FALSE)

  mn_target_nf_list <- mapply(function(beta, idx) {
    mn_g1[, idx] %*% beta
  }, beta_og, reg_idx, SIMPLIFY = FALSE)

  # Simple uncorrelated white noise for now with control
  # over the Signal-to-Noise ratio
  avg_signal_strength <- mean(sapply(mn_target_nf_list, function(z_) {
    mean(apply(z_, 2, sd))
  }))

  signal_to_noise_ratio <- 0.8

  sigma_noise <- avg_signal_strength / signal_to_noise_ratio

  noise_target_list <- lapply(mn_target_nf_list, function(z) {
    matrix(
      rnorm(prod(dim(z)), 0, sigma_noise),
      nrow = nrow(z),
      ncol = ncol(z)
    )
  })

  mn_target <- (
    do.call(cbind, mn_target_nf_list) + do.call(cbind, noise_target_list)
  )
  #
  # mn_target_scaled <- scale(mn_target, center = FALSE)
  # apply(mn_target_scaled, 2, sd)

  expression <- t(cbind(mn_target, mn_g1))
  is_regulator <- c(rep(0,ncol(mn_target)),rep(1,ncol(mn_g1)))
  # Without prior information gene symbols do not matter.
  # However, we have to supply them.
  genesymbols <- c(
    sprintf("T%d", seq_len(n_target)),
    sprintf("R%d", seq_len(n_reg))
  )
  # is_regulator <- c(rep.int(0, n_target), rep.int(1, n_reg))
  return(list(expression=as.data.frame(expression),
              genesymbols = genesymbols,
              is_regulator = is_regulator
              ))
}

#loads mn_g1
#cells are rows genes are columns
source(paste0(execution_path,"/../functions/biclust.R"))
source(paste0(execution_path,"/../functions/plot_cluster_history.R"))

neftel_path <- paste0(execution_path, '/../../datasets_sctargettranslator/Neftel2019')
Neftel_g1 <- readRDS(file = paste0(neftel_path, '/r_files/', "neftel_seurat_group1"))


# p1 <- DimPlot(Neftel_g1,reduction = "umap", group.by='malignant')
# p1

Neftel_g1<-SetIdent(Neftel_g1, value='malignant')
Neftel_g1_malignant<-subset(Neftel_g1, idents='yes')

z_g1<-GetAssayData(Neftel_g1_malignant, slot='scale.data')
dim(z_g1)
metaData<-Neftel_g1_malignant@meta.data$sample

out<-scregclust_format(z_g1)

genesymbols<-out[[1]]
sample_assignment<-out[[2]]
is_predictor<-out[[3]]

z_g1 <- z_g1[order(is_predictor),]
is_predictor <- is_predictor[order(is_predictor)]

#now cells are cols genes are row
dim(z_g1)

sel_corr_idx <- which(is_predictor == 1)

set.seed(10)

res1 <- felix_style_sim(mn_g1 = z_g1 ,
                number_regulators = NA, #only applicable if NOT including sel_corr_idx
                num_targets = 200,
                num_clusters = 1,
                group_mean_min = 0.01,
                group_mean_max = 0.033,
                sigma_of_betas = 0.1,
                sel_corr_idx = sel_corr_idx#optional selection of regulator genes, important if calling several times and planning to append (union)
)
res2 <- felix_style_sim(mn_g1 = z_g1 ,
                        number_regulators = NA, #only applicable if NOT including sel_corr_idx
                        num_targets = 200,
                        num_clusters = 2,
                        group_mean_min = 0.033,
                        group_mean_max = 0.066,
                        sigma_of_betas = 1,
                        sel_corr_idx = sel_corr_idx#optional selection of regulator genes, important if calling several times and planning to append (union)
)
res3 <- felix_style_sim(mn_g1 = z_g1 ,
                        number_regulators = NA, #only applicable if NOT including sel_corr_idx
                        num_targets = 200,
                        num_clusters = 4,
                        group_mean_min = 0.066,
                        group_mean_max = 0.1,
                        sigma_of_betas = 10,
                        sel_corr_idx = sel_corr_idx#optional selection of regulator genes, important if calling several times and planning to append (union)
)

dim(res1$expression)

res1$is_regulator

#set things up
n_cell_clusters <- 3

#toss out a bunch of cells
cells_per_cluster <- 3000
res1$expression <-res1$expression[,1:cells_per_cluster]
res2$expression <-res2$expression[,1:cells_per_cluster]
res3$expression <-res3$expression[,1:cells_per_cluster]

res<-merge(res1$expression, res2$expression, by = 0, all = TRUE)
rownames(res)<- rownames(res1$expression) #merge is a stupid function that merges in an id and then tosses the id for some reason
res<-merge(res[,-1], res3$expression, by = 0, all = TRUE)
res<- res[,-1] #merge is frikkin dumb
rownames(res)<- rownames(res1$expression)

#save true cell clustering
initial_cell_clust <- rep(1:n_cell_clusters, each = cells_per_cluster)

# Kod för att flytta 1% av cellerna i varje kluster till ett annat kluster.
disturbed_initial_cell_clust <- initial_cell_clust
for(i_cluster in 1:n_cell_clusters){
  indexes_of_cluster <- which(initial_cell_clust == i_cluster)
  some_of_those_indexes <- sample(indexes_of_cluster, size=as.integer(length(indexes_of_cluster)*0.01), replace = F)
  disturbed_initial_cell_clust[some_of_those_indexes] <- sample(c(1:n_cell_clusters)[-i_cluster], size=length(some_of_those_indexes), replace=T)
}

cell_cluster_history <- cbind(initial_cell_clust, disturbed_initial_cell_clust)
colnames(cell_cluster_history) <- c("True allocation", "Disturbed allocation")

#run biclust

execution_path <- dirname(rstudioapi::getSourceEditorContext()$path)
source(paste0(execution_path,"/../functions/biclust.R"))
source(paste0(execution_path,"/../functions/plot_cluster_history.R"))

res <- biclust(max_iter=50,
               initial_cluster_history = cell_cluster_history,
               is_regulator = is_regulator,
               n_target_gene_clusters = c(1,2,4),
               n_cells = c(cells_per_cluster,cells_per_cluster,cells_per_cluster),
               train_dat = as.matrix(res))

plot_cluster_history(cell_cluster_history = res$cell_cluster_history)


