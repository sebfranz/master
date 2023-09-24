
execution_path <- dirname(rstudioapi::getSourceEditorContext()$path)
source(paste0(execution_path,"/functions/generate_dummy_data_for_scregclust.R"))


Pt = 10      #number of target genes
Pr = 5       #number of regulator genes
n  = 10000   #number of cells
K  = 3       #Number of target gene clusters
regulator_mean   = 1
coefficient_mean = c(1,10,100)

#generate dummy data for each cell cluster that we want
dummy_res <- generate_dummy_data_for_scregclust(n_target_genes = Pt,  # Number of target genes
                                                    n_regulator_genes = Pr,  # Number of regulator genes
                                                    n_cells  = n,  # Number of cells
                                                    n_target_gene_clusters  = K,  # Number of target gene clusters
                                                    regulator_mean   = regulator_mean,  # Mean expression of regulator genes
                                                    coefficient_mean = coefficient_mean)  # Mean coefficients in true model, length n_target_gene_clusters)

#hack to get the list dummy_res into the global environment
for(iter in 1:length(names(dummy_res))) {
  assign(eval(names(dummy_res)[iter]),dummy_res[[iter]])
}

# apply simulated data to scregclust ---------------------------------------
library(scregclust)
?scregclust

scregclust(
  expression = rbind(t(Z_t), t(Z_r)),    #scRegClust wants this form
  genesymbols = 1:(Pt+Pr),               #gene row numbers
  is_regulator = (1:(Pt+Pr) > Pt) + 0, #vector indicating which genes are regulators
  n_cl        = K,
  # target_cluster_start = K,
  penalization = max(sapply(seq_along(R[,1]), function(i) sum(R[i,]))) + 1
  #maximal number of regulators for one cluster
)-> scRegOut

scRegOut$results
rowSums(Pi)
rowSums(R)
plot(scRegOut)$data


