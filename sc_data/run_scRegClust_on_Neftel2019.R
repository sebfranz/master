rm(list = ls())
execution_path <- dirname(rstudioapi::getSourceEditorContext()$path)

path <- paste0(execution_path, '/../../datasets_sctargettranslator/Neftel2019')
setwd(path)

library(scregclust)
library('Seurat')
library(Matrix)

source(paste0(execution_path,"/setup.R"))

if(!'neftel_seurat_group1' %in% list.files(paste0(path,'/r_files/'))){
  Neftel_g1 <- CreateSeuratObject(mn_g1,min.cells = 3,min.features =500,meta.data = cells)
  Neftel_g1 <- SCTransform(Neftel_g1)

  Neftel_g1 <- RunPCA(Neftel_g1, verbose = FALSE)
  Neftel_g1 <- RunUMAP(Neftel_g1, dims = 1:30, verbose = FALSE)

  saveRDS(Neftel_g1, file = paste0(path, '/r_files/', "neftel_seurat_group1"))
}else{
  Neftel_g1 <- readRDS(file = paste0(path, '/r_files/', "neftel_seurat_group1"))
}


p1 <- DimPlot(Neftel_g1,reduction = "umap", group.by='malignant')
p1

Neftel_g1<-SetIdent(Neftel_g1, value='malignant')
Neftel_g1_malignant<-subset(Neftel_g1, idents='yes')

z_g1<-GetAssayData(Neftel_g1_malignant, slot='scale.data')
metaData<-Neftel_g1_malignant@meta.data$sample

out<-scregclust_format(z_g1)

genesymbols<-out[[1]]
sample_assignment<-out[[2]]
is_predictor<-out[[3]]

for (i in 1:length(unique(metaData))){
  ix<-which(metaData==unique(metaData)[i])
  sample_assignment[ix]<-i
}

if(!'neftel_scregfit_group1' %in% list.files(paste0(path,'/r_files/'))){
  fit_g1 <- scregclust(
    z_g1,
    genesymbols,
    is_predictor,
    n_cl = 10L,
    penalization = 0.05,
    total_proportion = 0.5,
    n_cycles = 50L,
    noise_threshold = 0.05,
    min_cluster_size = 50,
    sample_assignment = sample_assignment
  )

  saveRDS(fit_g1, file = paste0(path, '/r_files/', "neftel_scregfit_group1"))
} else {
  fit_g1 <- readRDS(file = paste0(path, '/r_files/', "neftel_scregfit_group1"))
}
# library(regnet)
# plotRegNet(fit)

fit_g1$results
# rowSums(Pi)
# rowSums(R)
plot(fit_g1)$data

#could do some smarter visualisation of output


group_2_flag = F
if(group_2_flag){

  if(!'neftel_seurat_group2' %in% list.files(paste0(path,'/r_files/'))){
    Neftel_g2<-CreateSeuratObject(mn_g2,min.cells = 3,min.features =500,meta.data=cells)

    Neftel_g2<-ScaleData(Neftel_g2, do.scale=FALSE, do.center=FALSE)
    Neftel_g2<-FindVariableFeatures(Neftel_g2)

    Neftel_g2 <- RunPCA(Neftel_g2, npcs = 30, verbose = FALSE)
    Neftel_g2 <- RunUMAP(Neftel_g2, reduction = "pca", dims = 1:30)

    saveRDS(Neftel_g2, file = paste0(path, '/r_files/', "neftel_seurat_group2"))
  }else{
    Neftel_g2 <- readRDS(file = paste0(path, '/r_files/', "neftel_seurat_group2"))
  }

  p1 <- DimPlot(Neftel_g2, reduction = "umap", group.by = "malignant")
  p1

  Neftel_g2<-SetIdent(Neftel_g2,value='malignant')
  Neftel_g2_malignant<-subset(Neftel_g2, idents='yes')

  ix<-which(is.na(Neftel_g2_malignant@meta.data$MESlike2))
  Neftel_g2_malignant <- Neftel_g2_malignant[,!colnames(Neftel_g2_malignant) %in% colnames(Neftel_g2_malignant)[ix]]

  z_g2<-GetAssayData(Neftel_g2_malignant,slot='scale.data')
  metaData<-Neftel_g2_malignant@meta.data$sample

  out<-scregclust_format(z_g2)

  genesymbols<-out[[1]]
  sample_assignment<-out[[2]]
  is_predictor<-out[[3]]
  # skip this for now, it gives errors and is probably notn ecessary for our purposes
  # # this gives:
  # # Error in solve.default(xtx, xty) : ──────────────| 1/10 clusters
  # # Lapack routine dgesv: system is exactly singular: U[1596,1596] = 0
  # # In addition: Warning messages:
  # #   1: Quick-TRANSfer stage steps exceeded maximum (= 990350)
  # # 2: Quick-TRANSfer stage steps exceeded maximum (= 990350)
  #
  # if(!'neftel_scregfit_group2' %in% list.files(paste0(path,'/r_files/'))){
  #   fit_g2 <- scregclust(
  #     z_g2,
  #     genesymbols,
  #     is_predictor,
  #     n_cl = 10L,
  #     penalization = 0.03,
  #     total_proportion = 0.5,
  #     n_cycles = 50L,
  #     noise_threshold = 0.05,
  #     min_cluster_size = 20,
  #     sample_assignment = sample_assignment
  #   )
  #
  #   saveRDS(fit_g2, file = paste0(path, '/r_files/', "neftel_scregfit_group2"))
  # } else {
  #   fit_g2 <- readRDS(file = paste0(path, '/r_files/', "neftel_scregfit_group2"))
  # }
  # library(regnet)
  # plotRegNet(fit_g2)
}
