# lite rkod för att testköra scregclust

library(scregclust)
library(tidyverse)
#?scregclust

set.seed(3333)

Pt <- 10   # number of target genes
Pr <- 5    # number of regulator genes
n <- 100   # number of cells
K <- 3     # Number of target gene clusters


# Binary matrix Pi --------------------------------------------------------
# Which target gene is allocated to which cluster.
# Here it's randomly generated, for real data it would be smartly guessed.
# Rows are cluster index.
# Cols are target gene index.
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
R_init <- rbinom( K * Pr, 1, 1/K)  # needed later

R <- matrix(R_init, K, Pr)  # Row i is an indicator version of R_i in manuscript
R
R[1,]  # Cluster 1 is affected by these regulators
which(R[1,]!= 0)  # R_1 in the manuscript is which regulators affect cluster 1
sum(R[1,])  # |R_1| is the amount of regulator genes affecting cluster 1


# Matrix S ----------------------------------------------------------------
# A K x Pr matrix with 1 or -1 if the regulator (columns)
# of how (1 stimulating, -1 repressing) a regulator affects a cluster (rows),
# 0 if it doesn't affect it.
# This has the same information as the manuscript's s_i

# Create S by randomizing the signs of R_init.
S <- matrix(R_init * (rbinom(K * Pr, 1, 1/K)*2-1) , K, Pr)
S
# Non-zero entries of this is s_i in manuscript
S[1, which(S[1,]!= 0) ]
S[2, which(S[2,]!= 0) ]


# No idea -----------------------------------------------------------------

# For each regulator j in module i, a mean value was
# chosen uniformly at random between 0.01 and 0.1
# Regulator_means<- runif(Pr, 0.01, 0.1)
# Regulator_means

# Matrix Zr ---------------------------------------------------------------

# n x Pr, cells are rows, regulator genes are columns
# Just get some random expression for regulator genes for now
Z_r <- matrix( data = rnorm(n * Pr, mean = 100, sd = 10),
               nrow = n, ncol=Pr)
Z_r
dim(Z_r)

# Now we want to build 𝚩 and use 𝚩 to build Z_t

# Similar notation as publication
s_1 <- S[1,which(S[1,]!= 0) ]
R_1 <- which(R[1,]!= 0)
Z_r[,R_1] %*% diag(s_1)

# Needs work from here

# So the restricted 𝚩 in (1) has dimension |R_1| x Pt
# Complete matrix 𝚩 then has dimension Pr x Pt!
# this might actually need to be an array of dim K x Pr x Pt
# with separate depth representing different model
# todo: instead use R or S, create a 𝚩 matrix for each row
# Beta_i then should only have barameters for the non-zero
# entries of R
Beta <- matrix( data = rnorm(Pt * Pr, mean = 100, sd = 10),
               nrow = Pr, ncol=Pt)

dim(Beta)

Beta[R_i,]

(Z_r[,R_1] %*% diag(s_1 ) %*%  Beta[R_1,])[,which(Pi[1,] == 1)]

# if j is one target cell, that cells expression should then be,
# according to (1) in the manuscript

# Z_t could below be initialized to something non-zero, and in the next step the
# right hand side could be added instead of just inserted, this would make the
# initialisation similar to some baseline exposure, or intercept in the model.

Z_t <- matrix( data = 0,
               nrow = n, ncol=Pt)
dim(Z_t)

i <- 1

for(j in 1 : Pt){
  Z_t[,j] <-
    (Z_r[,which(R[i,]!= 0)] %*%
     diag(S[i,which(S[i,]!= 0) ] ) %*%
     Beta[which(R[i,]!= 0),j])
}  #todo:  vectorize this
Z_t
dim(Z_t)
