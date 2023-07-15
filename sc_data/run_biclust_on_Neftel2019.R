rm(list = ls())
execution_path <- dirname(rstudioapi::getSourceEditorContext()$path)
# neftel_path <- paste0(execution_path, '/../../datasets_sctargettranslator/Neftel2019')
# setwd(neftel_path)

library(scregclust)
library('Seurat')
library(Matrix)

#cells are rows genes are columns
source(paste0(execution_path,"/setup.R"))
dim(mn_g1)

execution_path <- dirname(rstudioapi::getSourceEditorContext()$path)
source(paste0(execution_path,"/../functions/biclust.R"))
source(paste0(execution_path,"/../functions/plot_cluster_history.R"))

neftel_path <- paste0(execution_path, '/../../datasets_sctargettranslator/Neftel2019')
Neftel_g1 <- readRDS(file = paste0(neftel_path, '/r_files/', "neftel_seurat_group1"))


p1 <- DimPlot(Neftel_g1,reduction = "umap", group.by='malignant')
p1

Neftel_g1<-SetIdent(Neftel_g1, value='malignant')
Neftel_g1_malignant<-subset(Neftel_g1, idents='yes')

z_g1<-GetAssayData(Neftel_g1_malignant, slot='scale.data')
dim(z_g1)
metaData<-Neftel_g1_malignant@meta.data$sample

out<-scregclust_format(z_g1)

genesymbols<-out[[1]]
sample_assignment<-out[[2]]
is_predictor<-out[[3]]

rm(mn_g1)
z_g1 <- z_g1[order(is_predictor),]
is_predictor <- is_predictor[order(is_predictor)]

#get initial cell clustering
cell_cluster_history <- kmeans(t(z_g1), 3)$cluster


res <- biclust(max_iter=50,
               initial_cluster_history = cell_cluster_history,
               is_regulator = is_predictor,
               n_target_gene_clusters = c(3,3,3),
               n_cells = c(sum(cell_cluster_history == 1),
                           sum(cell_cluster_history == 2),
                           sum(cell_cluster_history == 3)),
               train_dat = z_g1)

plot_cluster_history(cell_cluster_history = res$cell_cluster_history)



