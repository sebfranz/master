library(Matrix)
# This script requires that you have downloaded the neftel2019 dataset and
# added it to ~\repos\datasets_sctargettranslator\Neftel2019
# Download Group1 from https://uppsala.box.com/s/aqf3fyeb43i5jy8v4s527u0l262virue

execution_path <- dirname(rstudioapi::getSourceEditorContext()$path)
subfolder_to_master <- dirname(dirname(execution_path))
datasets_sctargettranslator <- file.path(subfolder_to_master, "datasets_sctargettranslator")
Neftel2019 <- file.path(datasets_sctargettranslator, "Neftel2019")
r_files <- file.path(Neftel2019, "r_files")
Group1 <- file.path(r_files, "Group1")
neftel_mn_group1 <- file.path(r_files, "neftel_mn_group1")
if(!file.exists(datasets_sctargettranslator)){
  # Create folder and put data there
  dir.create(datasets_sctargettranslator)
}
if(!file.exists(Neftel2019)){
  # Create folder and put data there
  dir.create(Neftel2019)
}
if(!file.exists(r_files)){
  # Create folder and put data there
  dir.create(r_files)
  stop(paste("Please download Group1 from https://uppsala.box.com/s/aqf3fyeb43i5jy8v4s527u0l262virue and put them in", r_files))
}

mtx_file <- file.path(Group1, 'Exp_data_UMIcounts_10X.mtx')
cell_file <- file.path(Group1, 'Cells_10X.txt')
gene_file <- file.path(Group1, 'Genes_10X.txt')
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


## Neftel SmartSeq2 ##

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


