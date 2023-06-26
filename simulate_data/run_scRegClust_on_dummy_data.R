
execution_path <- dirname(rstudioapi::getSourceEditorContext()$path)
source(paste0(execution_path,"/dummy_data.R"))
generate_dummy_data()

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
plot(scRegOut)$data
R

Pi
