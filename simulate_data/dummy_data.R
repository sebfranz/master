# Lite R-kod för att testköra scregclust

# n_target_genes <- 10   #number of target genes
# n_regulator_genes <- 5    #number of regulator genes
# n_cells <- 100   #number of cells
# n_target_gene_clusters <- 3     #Number of target gene clusters
# regulator_mean   = 1
# coefficient_mean = c(1,10,100)

# Binary matrix Pi --------------------------------------------------------
# Which target gene is allocated to which cluster.
# Here it's randomly generated, for real data it would be smartly guessed.
# Rows are cluster index.
# Cols are target gene index.
generate_dummy_data <- function(
    n_target_genes = 10,              #number of target genes
    n_regulator_genes = 5,               #number of regulator genes
    n_cells  = 100,             #number of cells
    n_target_gene_clusters  = 3,               #Number of target gene clusters
    regulator_mean   = 1, #mean expression of regulator genes
    coefficient_mean = c(1,10,100) #mean coefficients in true  model, length n_target_gene_clusters
)
{
  # Check arguments
  check_positive_scalar <- function(x){
    if (!(length(x) == 1L &&
        is.atomic(x) &&
        is.numeric(x) &&
        x == round(x) &&
        x >= 1)) {
      stop(paste0(deparse(substitute(a)), " must be positive scalar."))
    }
  }

  check_positive_scalar(n_target_genes)
  check_positive_scalar(n_regulator_genes)
  check_positive_scalar(n_cells)
  check_positive_scalar(n_target_gene_clusters)
  check_positive_scalar(regulator_mean)
  if (length(coefficient_mean) != n_target_gene_clusters) {
    stop("coefficient_mean must be a vector of positive numbers, and have length n_target_gene_clusters.")
  }


  diag_ <- function(vec)diag(x = vec, nrow = length(vec))
  # set.seed(3333) #To get a nice matrix
  Pi <- matrix(0, n_target_gene_clusters, n_target_genes)

  # To guarantee at least one target gene per cluster, assign the first n_target_gene_clusters target
  # gene to different clusters
  for(index in 1:n_target_gene_clusters){
    Pi[index, index] <- 1
  }
  # Set rest of matrix randomly
  for(index in (n_target_gene_clusters + 1):n_target_genes){
    Pi[sample.int(n = n_target_gene_clusters, size = 1), index] <- 1
  }

  Pi

  # Binary matrix R ---------------------------------------------------------
  # The regulators (columns) that affect each cluster (rows, i in the manuscript).
  # Note that in the paper this is a vector with set of indexes,
  # not a binary matrix, see code for an example.
  # Note that regulators can affect any number of clusters.

  # set.seed(12345)  # To get a nice matrix
  R <- matrix(rbinom(n_target_gene_clusters * n_regulator_genes, 1, 1/n_target_gene_clusters), n_target_gene_clusters, n_regulator_genes)  # Row i is an indicator version of R_i in manuscript

  #check if any cluster (row) has no regulator, and if so assign one
  for(row in 1:n_target_gene_clusters){
    if(sum(R[row,]) == 0){
      R[row, sample.int(n_regulator_genes, 1)] <- 1
    }
  }

  R

  R[1,]  # Cluster 1 is affected by these regulators
  sum(R[1,])  # R_1 in the manuscript is which regulators affect cluster 1

  R2R_i <- function(i) {
    which(R[i,]!= 0)
  }
  R2R_i(1) # R_1 in the manuscript is which regulators affect cluster 1

  # Matrix S ----------------------------------------------------------------
  # A n_target_gene_clusters x n_regulator_genes matrix with 1 or -1 if the regulator (columns)
  # of how (1 stimulating, -1 repressing) a regulator affects a cluster (rows),
  # 0 if it doesn't affect it.
  # This has the same information as the manuscript's s_i
  # set.seed(10) # to get a nice matrix
  S <- R * matrix(rbinom(n = n_target_gene_clusters * n_regulator_genes, 1, 0.8)*2-1 , n_target_gene_clusters, n_regulator_genes)  # Just randomize signs

  S2S_i <- function(i, mat=S) {
    mat[i, which(mat[i,]!= 0) ]
  }
  # S2S_i(2)  # Non-zero entries of this is s_i in manuscript

  # For each regulator j in module i, a mean value was
  # chosen uniformly at random between 0.01 and 0.1
  # Regulator_means<- runif(n_regulator_genes, 0.01, 0.1)
  # Regulator_means

  # Matrix Zr ---------------------------------------------------------------
  # n_cells x n_regulator_genes, cells are rows, regulator genes are columns
  # Just get some random expression for regulator genes for now
  Z_r <- matrix(data = rnorm(n_cells * n_regulator_genes, mean = regulator_mean, sd = 0.1),
                nrow = n_cells, ncol=n_regulator_genes)
  # Z_r <- sapply(1:n_regulator_genes, function(i)  rnorm(n_cells, mean = runif(1, 0.1,1), sd = 0.1) )

  Z_r
  dim(Z_r)

  # Array 𝚩 ---------------------------------------------------------------
  # Now we want to build 𝚩 and use 𝚩 to build Z_t.
  # For that we need some coefficients from our regression models.
  # These coefficients are stored in the array Beta, which has
  # dimension  n_regulator_genes x n_target_genes x n_target_gene_clusters with nonegative entries
  # in this we store the (in the manuscript only the non-zero) coefficients
  # describing how the regulator genes affect the target genes.

  # Beta <- array(data = abs(rnorm(n_regulator_genes * n_target_genes * n_target_gene_clusters,
  #                                 mean = coefficient_mean, sd = 0.1)),
  #               c(n_regulator_genes,n_target_genes,n_target_gene_clusters))

  Beta <- array(
    data = sapply(1:n_target_gene_clusters, function(i) rnorm(n_regulator_genes*n_target_genes, mean = coefficient_mean[i], sd = 0.1)),
    dim = c(n_regulator_genes,n_target_genes,n_target_gene_clusters)
  )

  # Beta <- array(data = unlist(
  #   sapply(1:n_target_gene_clusters, function(i) rnorm(n_regulator_genes*n_target_genes, mean = runif(1,min = 1, max = 1), sd = 0.1), simplify = F)
  # ), dim = c(n_regulator_genes,n_target_genes,n_target_gene_clusters) )

  # Make 𝚩 zero in appropriate spots
  for (clust in 1:n_target_gene_clusters){
    Beta[,,clust] <-  diag_(R[clust,]) %*% Beta[,,clust]
  }
  Beta
  dim(Beta)

  # In the manuscript the zero rows are just dropped

  Beta2Beta_i <- function(i){
    matrix(data = Beta[R2R_i(i),,i], nrow = length(R2R_i(i)), ncol = n_target_genes)
  }

  Beta2Beta_i(3) #Beta_i as in the manuscript, has dimension |R_i| x n_target_genes
  Beta[,,3]

  # Matrix Z_t --------------------------------------------------------------
  # If j is one target cell, that cells expression should then be,
  # according to (1) in the manuscript.

  # Z_t could below be initialized to something nonzero, and in the next step the
  # right hand side could be added instead of just inserted, this would make the
  # initialisation similar to some baseline exposure, or intercept in the model.

  Z_t <- matrix(data = 0, nrow = n_cells, ncol=n_target_genes)

  #todo:  vectorize this
  for(i in 1:n_target_gene_clusters){
    for(j in 1 : n_target_genes){
      Z_t[,j] <-
        Z_t[,j] +
        Pi[i,j] *               #  True cluster allocation, zero if Z_t[,j] is
        (                       # not in cluster i
          Z_r[,R2R_i(i)] %*%    #  Gene expression of regulators of cluster i
            diag_(S2S_i(i)) %*%  #  signs for wether regulators are stim or repress
            Beta2Beta_i(i)[,j]  #  how much reg of cluster i affects target j
        )
    }
    cat(paste0("building cluster ", i,"\n_cells"))
  }
  # This can probably be vectorized
  # For this we are omitting the variance terms.
  # Z_t, Z_r, S_i, and B_i as here will minimize (1) in the manuscript
  list(
    Z_t = Z_t,
    Z_r = Z_r,
    Pi  = Pi,
    R   = R,
    S   = S
  )
}

# Z_t
# dim(Z_t)

# # apply simulated data to scregclust ---------------------------------------
# library(scregclust)
# ?scregclust
#
#
# scregclust(
#   expression = rbind(t(Z_t), t(Z_r)),    #scRegClust wants this form
#   genesymbols = 1:(n_target_genes+n_regulator_genes),               #gene row numbers
#   is_regulator = (1:(n_target_genes+n_regulator_genes) > n_target_genes) + 0, #vector indicating which genes are regulators
#   target_cluster_start = n_target_gene_clusters,
#   penalization = max(sapply(seq_along(R[,1]), function(i) sum(R[i,]))) + 1
#   #maximal number of regulators for one cluster
# )-> scRegOut
#
# scRegOut$results
#
# R
#
# Pi
