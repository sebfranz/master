if (!require(Matrix)) install.packages('Matrix')
if (!require(Seurat)) install.packages('Seurat')
if (!require(here)) install.packages('here')
library(Matrix)
library(Seurat)  # To work with Neftel data
library(here  # To handle paths
# This script requires that you have downloaded the neftel2019 dataset and
# added it to ~\repos\datasets_sctargettranslator\Neftel2019
# Download Group1 and/or Group2 from https://uppsala.box.com/s/aqf3fyeb43i5jy8v4s527u0l262virue
# This creates and/or loads the datasets

# Which dataset to load
group_1_flag <- T
group_2_flag <- F

# Folders and filenames
root_path <- here::here()
execution_path <- here::here("sc_data")
subfolder_to_root <- dirname(root_path)
path_datasets_sctargettranslator <- file.path(subfolder_to_root, "datasets_sctargettranslator")
path_Neftel2019 <- file.path(path_datasets_sctargettranslator, "Neftel2019")
path_r_files <- file.path(path_Neftel2019, "r_files")

# Create folder tree if it doesn't exist
if(!file.exists(path_datasets_sctargettranslator)){
  dir.create(path_datasets_sctargettranslator)
}
if(!file.exists(path_Neftel2019)){
  dir.create(path_Neftel2019)
}
if(!file.exists(path_r_files)){
  dir.create(path_r_files)
}

# Group 1 ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if(group_1_flag){
  path_Group1 <- file.path(path_r_files, "Group1")
  path_mtx_file <- file.path(path_Group1, 'Exp_data_UMIcounts_10X.mtx')
  path_cell_file <- file.path(path_Group1, 'Cells_10X.txt')
  path_gene_file <- file.path(path_Group1, 'Genes_10X.txt')
  path_neftel_mn_group1 <- file.path(path_r_files, "neftel_mn_group1.rds")
  path_neftel_seurat_group1 <- file.path(path_r_files, "neftel_seurat_group1.rds")

  if(!file.exists(path_Group1)){
    dir.create(path_Group1)
    stop(paste("Please download Group1 from https://uppsala.box.com/s/aqf3fyeb43i5jy8v4s527u0l262virue and put them in", path_r_files))
  }

  if(!file.exists(path_neftel_mn_group1)){
    mn_g1<-readMM(file=path_mtx_file)
    mn_g1<-as.matrix(mn_g1)

    cells<-read.table(file=path_cell_file, sep=' ', header=TRUE, stringsAsFactors = FALSE)
    genes<-read.table(file=path_gene_file, sep='\t', header=FALSE, stringsAsFactors = FALSE)
    genes<-genes[,1]

    rownames(cells)<-cells[,1]

    rownames(mn_g1)<-genes
    colnames(mn_g1)<-cells$cell_name

    saveRDS(mn_g1, file = file.path(path_neftel_mn_group1))

  }else{
    mn_g1 <- readRDS(file = path_neftel_mn_group1)

    cells<-read.table(file=path_cell_file, sep=' ', header=TRUE, stringsAsFactors = FALSE)
    genes<-read.table(file=path_gene_file, sep='\t', header=FALSE, stringsAsFactors = FALSE)
    genes<-genes[,1]

    rownames(cells)<-cells[,1]
  }

  if(!file.exists(path_neftel_seurat_group1)){
    Neftel_g1 <- CreateSeuratObject(mn_g1, min.cells=3, min.features=500, meta.data=cells)
    Neftel_g1 <- SCTransform(Neftel_g1)

    Neftel_g1 <- RunPCA(Neftel_g1, verbose = FALSE)
    Neftel_g1 <- RunUMAP(Neftel_g1, dims = 1:30, verbose = FALSE)

    saveRDS(Neftel_g1, file=path_neftel_seurat_group1)
  }else{
    Neftel_g1 <- readRDS(file=path_neftel_seurat_group1)
  }
}




# Group 2 ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


path_neftel_mn_group2 <- file.path(path_r_files, "neftel_mn_group2.rds")
path_Group2 <- file.path(path_r_files, "Group2")
if(group_2_flag){
  if(!file.exists(path_neftel_mn_group2)){
    if(!file.exists(path_datasets_sctargettranslator)){
     dir.create(path_datasets_sctargettranslator)
    }
    if(!file.exists(path_Neftel2019)){
     dir.create(path_Neftel2019)
    }
    if(!file.exists(path_r_files)){
     dir.create(path_r_files)
    }
    if(!file.exists(path_Group2)){
      dir.create(path_Group2)
      stop(paste("Please download Group2 from https://uppsala.box.com/s/aqf3fyeb43i5jy8v4s527u0l262virue and put them in", path_r_files))
    }

    mn_g2<-readMM(file=file.path(path_Group2, 'exp_data_TPM.mtx'))
    mn_g2<-as.matrix(mn_g2)

    cells<-read.table(file=file.path(path_Group2, 'Cells.txt'),sep=' ',header=TRUE,stringsAsFactors = FALSE)
    genes<-read.table(file=file.path(path_Group2, 'Genes.txt'),sep='\t', header=FALSE,stringsAsFactors = FALSE)
    genes<-genes[,1]

    rownames(cells)<-cells[,1]

    rownames(mn_g2)<-genes
    colnames(mn_g2)<-cells$cell_name

    saveRDS(mn_g2, file = file.path(path_r_files, "neftel_mn_group2.rds"))
    rm(cells, genes,group_2_flag)
  }else{
    mn_g2 <- readRDS(file = file.path(path_r_files, "neftel_mn_group2.rds"))
  }

  path_neftel_seurat_group2 <- file.path(path_r_files, "neftel_seurat_group2.rds")
  if(!file.exists()){
    Neftel_g2 <- CreateSeuratObject(mn_g2, min.cells=3, min.features=500, meta.data=cells)
    Neftel_g2 <- SCTransform(Neftel_g2)

    Neftel_g2 <- RunPCA(Neftel_g2, verbose = FALSE)
    Neftel_g2 <- RunUMAP(Neftel_g2, dims = 1:30, verbose = FALSE)

    saveRDS(Neftel_g2, file=path_neftel_seurat_group2)
  }else{
    Neftel_g2 <- readRDS(file=path_neftel_seurat_group2)
  }
}


