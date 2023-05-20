#lite rkod för att testköra scregclust

library(scregclust)
library(tidyverse)
#?scregclust

set.seed(3333)

Pt <- 10   #number of target genes
Pr <- 5    #number of regulator genes
n <- 100   #number of cells
K <- 3     #Number of target gene clusters

#matrix of true cluster allocations
Pi <- matrix(0, K, Pt)
  #fan va fult refrakturera

for(index in 1:Pt){
  Pi[ sample.int( n = K, size = 1), index] <- 1
}
Pi

# for each cluster, the subset of
# regulators that are affecting that cluster R_i
# note that regulators can affect any number of clusters

R_init <- rbinom( K * Pr, 1, 1/K) #needed later

R <- matrix(R_init, K, Pr) #row i is an indicator version of R_i in manuscript
R
R[1,] #cluster 1 is affected by these regulators
which(R[1,]!= 0) # R_1 in the manuscript is which regulators affect cluster 1
sum(R[1,]) #|R_1| is the amount of regulator genes affecting cluster 1


#we start by making S, a K x Pr matrix with ones/minus ones if the columns
#regulator affects that rows cluster, zeroes otherwise.
# This has the same information as the publications' s_i
S <- matrix(R_init * (rbinom(K * Pr, 1, 1/K)*2-1) , K, Pr)#just randomize signs
S
S[1, which(S[1,]!= 0) ] #nonzero entries of this is s_i in manuscript
S[2, which(S[2,]!= 0) ]

#For each regulator j in module i, a mean value was
#chosen uniformly at random between 0.01 and 0.1
# Regulator_means<- runif(Pr, 0.01, 0.1)
# Regulator_means

#now generate Zr

#just get some random expression for regulator genes for now
Z_r <- matrix( data = rnorm(n * Pr, mean = 100, sd = 10),
               nrow = n, ncol=Pr)
Z_r
dim(Z_r)

#now we want to build beta and use beta to build Z_t

# similar notation as publication
s_1 <- S[1,which(S[1,]!= 0) ]
R_1 <- which(R[1,]!= 0)


Z_r[,R_1] %*% diag(s_1 )

#needs work from here

#so the restricted beta in (1) has dimention |R_1| x Pt
# complete matrix beta then has dimension Pr x Pt!
    # this might actually need to be an array of dim K x Pr x Pt
    # with separate depth representing different model
Beta <- matrix( data = rnorm(Pt * Pr, mean = 100, sd = 10),
               nrow = Pr, ncol=Pt)

dim(Beta)

Beta[R_i,]

(Z_r[,R_1] %*% diag(s_1 ) %*%  Beta[R_1,])[,which(Pi[1,] == 1)]

# if j is one target cell, that cells expression should then be,
# according to (1) in the manuscript

# Z_t could below be initialized to something nonzero, and in the next step the
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