#lite rkod för att testköra scregclust

library(scregclust)
#?scregclust


#TODO make below a callable function

set.seed(3333)

Pt <- 10   #number of target genes
Pr <- 5    #number of regulator genes
n <- 100   #number of cells
K <- 3     #Number of target gene clusters

#matrix of true cluster allocations
Pi <- matrix(0, K, Pt)

for(index in 1:Pt){
  Pi[ sample.int( n = K, size = 1), index] <- 1
}
Pi

# for each cluster, the subset of
# regulators that are affecting that cluster R_i
# note that regulators can affect any number of clusters

set.seed(12345)
R <- matrix(rbinom( K * Pr, 1, 1/K), K, Pr) #row i is an indicator version of R_i in manuscript
R

R[1,] #cluster 1 is affected by these regulators
sum(R[1,]) #|R_1| is the amount of regulator genes affecting cluster 1

R2R_i <- function(i) {
  which(R[i,]!= 0)
}
R2R_i(1) # R_1 in the manuscript is which regulators affect cluster 1


#we start by making S, a K x Pr matrix with ones/minus ones if the columns
#regulator affects that rows cluster, zeroes otherwise.
# This has the same information as the publications' s_i
set.seed(10)
S <- R * matrix(rbinom(n = K * Pr, 1, 0.8)*2-1 , K, Pr) #just randomize signs
S

S2S_i <- function(i) {
  S[i, which(S[i,]!= 0) ]
}
S2S_i(2) #nonzero entries of this is s_i in manuscript

#For each regulator j in module i, a mean value was
#chosen uniformly at random between 0.01 and 0.1
# Regulator_means<- runif(Pr, 0.01, 0.1)
# Regulator_means

#now generate Zr

#just get some random expression for regulator genes for now
Z_r <- matrix( data = rnorm(n * Pr, mean = 1, sd = 0.1),
               nrow = n, ncol=Pr)
Z_r
dim(Z_r)

#now we want to build beta and use beta to build Z_t

#For that we need some coefficients from our regression models
# These coefficients are stored in the array Beta, which has
#dimension Pr x Pt x K
#in this we store the (in the manuscript only the nonzero) coefficients
#describing how the regulator genes affect the target genes.

#for now its just one distr could be made more sophisticated

Beta <- array(data = rnorm(K * Pt * Pr, mean = 1, sd = 0.1), c(Pr,Pt,K))
#make beta zero in appropriate spots
for (clust in 1:K){
  Beta[,,clust] <-  diag(R[clust,]) %*% Beta[,,clust]
}
Beta
dim(Beta)


# Beta_alt <- x <- vector(mode = "list", length = K)
# for (clust in 1:K){
#   Beta_alt[[clust]] <-  as.matrix( Beta[which(Beta[,1,clust] != 0),,clust])
# } # this can probably be vectorized
# Beta_alt


#in the manuscript the zero rows are just dropped


Beta2Beta_i <- function(i){
  matrix(data = Beta[R2R_i(i),,i], nrow = length(R2R_i(i)), ncol = Pt)
}

Beta2Beta_i(3)
Beta[,,3]

# if j is one target cell, that cells expression should then be,
# according to (1) in the manuscript

  # Z_t could below be initialized to something nonzero, and in the next step the
  # right hand side could be added instead of just inserted, this would make the
    # initialisation similar to some baseline exposure, or intercept in the model.

Z_t <- matrix( data = 0,
               nrow = n, ncol=Pt)

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
} #this can probably be vectorized

# Z_t
dim(Z_t)





