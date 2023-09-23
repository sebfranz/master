library(Matrix)
library(Seurat)
# This script requires that you have downloaded the neftel2019 dataset and
# added it to ~\repos\datasets_sctargettranslator\Neftel2019
# Download Group1 from https://uppsala.box.com/s/aqf3fyeb43i5jy8v4s527u0l262virue

# All the folders and filenames
execution_path <- dirname(rstudioapi::getSourceEditorContext()$path)
subfolder_to_master <- dirname(dirname(execution_path))
datasets_sctargettranslator <- file.path(subfolder_to_master, "datasets_sctargettranslator")
Neftel2019 <- file.path(datasets_sctargettranslator, "Neftel2019")
r_files <- file.path(Neftel2019, "r_files")
Group1 <- file.path(r_files, "Group1")
neftel_mn_group1 <- file.path(r_files, "neftel_mn_group1.rds")
mtx_file <- file.path(Group1, 'Exp_data_UMIcounts_10X.mtx')
cell_file <- file.path(Group1, 'Cells_10X.txt')
gene_file <- file.path(Group1, 'Genes_10X.txt')
neftel_seurat_group1 <- file.path(r_files, "neftel_seurat_group1.rds")

# Create folder tree and check if user downloaded the Group 1 Neftel data
if(!file.exists(datasets_sctargettranslator)){
  dir.create(datasets_sctargettranslator)
}
if(!file.exists(Neftel2019)){
  dir.create(Neftel2019)
}
if(!file.exists(r_files)){
  dir.create(r_files)
  stop(paste("Please download Group1 from https://uppsala.box.com/s/aqf3fyeb43i5jy8v4s527u0l262virue and put them in", r_files))
}


if(!file.exists(neftel_mn_group1)){
  mn_g1<-readMM(file=mtx_file)
  mn_g1<-as.matrix(mn_g1)

  cells<-read.table(file=cell_file, sep=' ', header=TRUE, stringsAsFactors = FALSE)
  genes<-read.table(file=gene_file, sep='\t', header=FALSE, stringsAsFactors = FALSE)
  genes<-genes[,1]

  rownames(cells)<-cells[,1]

  rownames(mn_g1)<-genes
  colnames(mn_g1)<-cells$cell_name

  saveRDS(mn_g1, file = file.path(neftel_mn_group1))

}else{
  mn_g1 <- readRDS(file = neftel_mn_group1)

  cells<-read.table(file=cell_file, sep=' ', header=TRUE, stringsAsFactors = FALSE)
  genes<-read.table(file=gene_file, sep='\t', header=FALSE, stringsAsFactors = FALSE)
  genes<-genes[,1]

  rownames(cells)<-cells[,1]
}

if(!file.exists(neftel_seurat_group1)){
  Neftel_g1 <- CreateSeuratObject(mn_g1, min.cells=3, min.features=500, meta.data=cells)
  Neftel_g1 <- SCTransform(Neftel_g1)

  Neftel_g1 <- RunPCA(Neftel_g1, verbose = FALSE)
  Neftel_g1 <- RunUMAP(Neftel_g1, dims = 1:30, verbose = FALSE)

  saveRDS(Neftel_g1, file=neftel_seurat_group1)
}else{
  Neftel_g1 <- readRDS(file=neftel_seurat_group1)
}




# Group 2 ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

group_2_flag = F
if(group_2_flag){

  library(Matrix)
  if(!'neftel_mn_group2' %in% list.files(paste0(path,'/r_files/'))){
    mn_g2<-readMM(file='./Group2/exp_data_TPM.mtx')
    mn_g2<-as.matrix(mn_g2)

    cells<-read.table(file='./Group2/Cells.txt',sep=' ',header=TRUE,stringsAsFactors = FALSE)
    genes<-read.table(file='./Group2/Genes.txt',sep='\t', header=FALSE,stringsAsFactors = FALSE)
    genes<-genes[,1]

    rownames(cells)<-cells[,1]

    rownames(mn_g2)<-genes
    colnames(mn_g2)<-cells$cell_name

    saveRDS(mn_g2, file = paste0(path, '/r_files/', "neftel_mn_group2"))
    rm(cells, genes,group_2_flag)
  }else{
    mn_g2 <- readRDS(file = paste0(path, '/r_files/', "neftel_mn_group2"))
  }
}


