rm(list = ls())
execution_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(execution_path)
# neftel_path <- paste0(execution_path, '/../../datasets_sctargettranslator/Neftel2019')
# setwd(neftel_path)

library(scregclust)
library('Seurat')
library(Matrix)

#loads mn_g1
#cells are rows genes are columns
source(paste0(execution_path,"/../sc_data/setup.R"))



#number of regulators to keep
n_reg <- 100
# Assume there are 5 true clusters in the data
n_cl <- 5

# sel_corr_idx <- sample(corr_idx, n_reg)
sel_corr_idx <- sample(nrow(mn_g1), n_reg)

mn_g1 <- t(mn_g1[sel_corr_idx, ])
dim(mn_g1)


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
grp_means <- lapply(reg_idx, function(idx) runif(length(idx), 0.01, 0.1))

# Simulate cluster sizes (make them all equal size for now)
n_target <- 500
k_true <- rep(seq_len(n_cl), each = n_target / n_cl)

# Generate coefficients for each target gene
sigma_beta <- 0.1
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

# Without prior information gene symbols do not matter.
# However, we have to supply them.
genesymbols <- c(
  sprintf("T%d", seq_len(n_target)),
  sprintf("R%d", seq_len(n_reg))
)
is_regulator <- c(rep.int(0, n_target), rep.int(1, n_reg))
