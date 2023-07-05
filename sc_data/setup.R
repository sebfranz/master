
execution_path <- dirname(rstudioapi::getSourceEditorContext()$path)

if( !'datasets_sctargettranslator' %in% list.files(paste0(execution_path, '/../../'))){
  #create folder and put data there
  dir.create((paste0(execution_path, '/../../datasets_sctargettranslator')) )
}

setwd(paste0(execution_path, '/../../datasets_sctargettranslator'))
# download data
# download.file('https://uppsala.box.com/s/aqf3fyeb43i5jy8v4s527u0l262virue',
#               'C:/Users/Filip/repos/master/sc_data/test' )




#make things prettier

## Neftel 10X ##
setwd(paste0(execution_path, '/../../datasets_sctargettranslator/Neftel2019'))

library(Matrix)
mn<-readMM(file='./Group1/Exp_data_UMIcounts_10X.mtx')
mn<-as.matrix(mn)

cells<-read.table(file='./Group1/Cells_10X.txt',sep=' ',header=TRUE,stringsAsFactors = FALSE)
genes<-read.table(file='./Group1/Genes_10X.txt',sep='\t', header=FALSE,stringsAsFactors = FALSE)
genes<-genes[,1]

rownames(cells)<-cells[,1]

rownames(mn)<-genes
colnames(mn)<-cells$cell_name

library('Seurat')
Neftel <- CreateSeuratObject(mn,min.cells = 3,min.features =500,meta.data = cells)
Neftel <- SCTransform(Neftel)

Neftel <- RunPCA(Neftel, verbose = FALSE)
Neftel <- RunUMAP(Neftel, dims = 1:30, verbose = FALSE)

p1 <- DimPlot(Neftel,reduction = "umap",group.by='malignant')
p1

Neftel<-SetIdent(Neftel,value='malignant')
Neftel_malignant<-subset(Neftel, idents='yes')

z<-GetAssayData(Neftel_malignant,slot='scale.data')
metaData<-Neftel_malignant@meta.data$sample

setwd('/Users/idala384/Desktop/scEM_package/scripts/')
devtools::load_all()

out<-scregclust_format(z)

genesymbols<-out[[1]]
sample_assignment<-out[[2]]
is_predictor<-out[[3]]

for (i in 1:length(unique(metaData))){
  ix<-which(metaData==unique(metaData)[i])
  sample_assignment[ix]<-i
}

fit <- scregclust(
  z, genesymbols, is_predictor,
  k_start = 10L, lambda = 0.05, total_proportion = 0.5,
  n_cycles = 50L, noise_threshold = 0.05, min_cluster_size = 50,
  sample_assignment = sample_assignment
)

plotRegNet(fit)

## Neftel SmartSeq2 ##
setwd('/Users/idala384/Desktop/Neftel2019')

cells<-read.table(file='./Group2/cells.txt',sep=' ',header=TRUE,stringsAsFactors = FALSE)
rownames(cells)<-cells[,1]

mn<-readMM(file='./exp_data_TPM.mtx')
mn<-as.matrix(mn)

Neftel<-CreateSeuratObject(mn,min.cells = 3,min.features =500,meta.data=cells)

Neftel<-ScaleData(Neftel, do.scale=FALSE, do.center=FALSE)
Neftel<-FindVariableFeatures(Neftel)

Neftel <- RunPCA(Neftel, npcs = 30, verbose = FALSE)
Neftel <- RunUMAP(Neftel, reduction = "pca", dims = 1:30)

p1 <- DimPlot(Neftel, reduction = "umap", group.by = "malignant")
p1

Neftel<-SetIdent(Neftel,value='malignant')
Neftel_malignant<-subset(Neftel, idents='yes')

ix<-which(is.na(Neftel_malignant@meta.data$MESlike2))
Neftel_malignant <- Neftel_malignant[,!colnames(Neftel_malignant) %in% colnames(Neftel_malignant)[ix]]

z<-GetAssayData(Neftel_malignant,slot='scale.data')
metaData<-Neftel_malignant@meta.data$sample

setwd('/Users/idala384/Desktop/scEM_package/scripts/')
devtools::load_all()

out<-scregclust_format(z)

genesymbols<-out[[1]]
sample_assignment<-out[[2]]
is_predictor<-out[[3]]

for (i in 1:length(unique(metaData))){
  ix<-which(metaData==unique(metaData)[i])
  sample_assignment[ix]<-i
}

fit <- scregclust(
  z, genesymbols, is_predictor,
  k_start = 10L, lambda = 0.03, total_proportion = 0.5,
  n_cycles = 50L, noise_threshold = 0.05, min_cluster_size = 20,
  sample_assignment = sample_assignment
)

plotRegNet(fit)
