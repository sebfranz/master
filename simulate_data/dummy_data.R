# Lite R-kod för att testköra scregclust

Pt <- 10   #number of target genes
Pr <- 5    #number of regulator genes
n <- 100   #number of cells
K <- 3     #Number of target gene clusters

# Binary matrix Pi --------------------------------------------------------
# Which target gene is allocated to which cluster.
# Here it's randomly generated, for real data it would be smartly guessed.
# Rows are cluster index.
# Cols are target gene index.

set.seed(3333) #To get a nice matrix
Pi <- matrix(0, K, Pt)

for(index in 1:Pt){
  Pi[ sample.int( n = K, size = 1), index] <- 1
}
Pi

# Binary matrix R ---------------------------------------------------------
# The regulators (columns) that affect each cluster (rows, i in the manuscript).
# Note that in the paper this is a vector with set of indexes,
# not a binary matrix, see code for an example.
# Note that regulators can affect any number of clusters.

set.seed(12345)  # To get a nice matrix
R <- matrix(rbinom( K * Pr, 1, 1/K), K, Pr)  # Row i is an indicator version of R_i in manuscript
R

R[1,]  # Cluster 1 is affected by these regulators
sum(R[1,])  # R_1 in the manuscript is which regulators affect cluster 1

R2R_i <- function(i) {
  which(R[i,]!= 0)
}
R2R_i(1) # R_1 in the manuscript is which regulators affect cluster 1


# Matrix S ----------------------------------------------------------------
# A K x Pr matrix with 1 or -1 if the regulator (columns)
# of how (1 stimulating, -1 repressing) a regulator affects a cluster (rows),
# 0 if it doesn't affect it.
# This has the same information as the manuscript's s_i
set.seed(10) # to get a nice matrix
S <- R * matrix(rbinom(n = K * Pr, 1, 0.8)*2-1 , K, Pr)  # Just randomize signs
S

S2S_i <- function(i) {
  S[i, which(S[i,]!= 0) ]
}
S2S_i(2)  # Non-zero entries of this is s_i in manuscript

# For each regulator j in module i, a mean value was
# chosen uniformly at random between 0.01 and 0.1
# Regulator_means<- runif(Pr, 0.01, 0.1)
# Regulator_means

# Matrix Zr ---------------------------------------------------------------
# n x Pr, cells are rows, regulator genes are columns
# Just get some random expression for regulator genes for now
Z_r <- matrix( data = rnorm(n * Pr, mean = 1, sd = 0.1),
               nrow = n, ncol=Pr)
Z_r
dim(Z_r)

# Array 𝚩 ---------------------------------------------------------------
# Now we want to build 𝚩 and use 𝚩 to build Z_t.
# For that we need some coefficients from our regression models.
# These coefficients are stored in the array Beta, which has
# dimension Pt x Pr x K
# in this we store the (in the manuscript only the non-zero) coefficients
# describing how the regulator genes affect the target genes.

# For now its just one distr could be made more sophisticated

Beta <- array(data = rnorm( Pt * Pr * K, mean = 1, sd = 0.1), c(Pr,Pt,K))

# Make 𝚩 zero in appropriate spots
for (clust in 1:K){
  Beta[,,clust] <-  diag(R[clust,]) %*% Beta[,,clust]
}
Beta
dim(Beta)

# In the manuscript the zero rows are just dropped

Beta2Beta_i <- function(i){
  matrix(data = Beta[R2R_i(i),,i], nrow = length(R2R_i(i)), ncol = Pt)
}

Beta2Beta_i(3) #Beta_i as in the manuscript, has dimension |R_i| x Pt
Beta[,,3]

# Matrix Z_t --------------------------------------------------------------
# If j is one target cell, that cells expression should then be,
# according to (1) in the manuscript.

# Z_t could below be initialized to something nonzero, and in the next step the
# right hand side could be added instead of just inserted, this would make the
# initialisation similar to some baseline exposure, or intercept in the model.

Z_t <- matrix( data = 0, nrow = n, ncol=Pt)

#todo:  vectorize this
for(i in 1:K){
  for(j in 1 : Pt){
    Z_t[,j] <-
      Z_t[,j] +
      Pi[i,j] *
      (
        Z_r[,R2R_i(i)] %*%
          diag(S2S_i(i)) %*%
          Beta2Beta_i(i)[,j]
      )
  }
  cat(paste0("building cluster ", i,"\n"))
}
# This can probably be vectorized
# For this we are omitting the variance terms.
# Z_t, Z_r, S_i, and B_i as here will minimize (1) in the manuscript


# Z_t
dim(Z_t)

# apply simulated data to scregclust ---------------------------------------
library(scregclust)
?scregclust


scregclust(
  expression = rbind(t(Z_t), t(Z_r)),    #scRegClust wants this form
  genesymbols = 1:(Pt+Pr),               #gene row numbers
  is_regulator = (genesymbols > Pt) + 0, #vector indicating which genes are regulators
  target_cluster_start = K,
  penalization = max(sapply(seq_along(R[,1]), function(i) sum(R[i,])))
  #maximal number of regulators for one cluster
)-> scRegOut

names(scRegOut)

scRegOut$results

R

Pi
