
execution_path <- dirname(rstudioapi::getSourceEditorContext()$path)

#This script requires that you have downloaded the neftel2019 dataset and
# added it to ~\repos\datasets_sctargettranslator\Neftel2019

if( !'datasets_sctargettranslator' %in% list.files(paste0(execution_path, '/../../'))){
  #create folder and put data there
  # dir.create((paste0(execution_path, '/../../datasets_sctargettranslator')) )

  stop(paste0('please check that you have a folder called "datasets_sctargettranslator" with subfolder "Neftel2019" in ', execution_path))

}

setwd(paste0(execution_path, '/../../datasets_sctargettranslator'))
# todo: download data automatically from the box
# download.file('https://uppsala.box.com/s/aqf3fyeb43i5jy8v4s527u0l262virue',
#               'C:/Users/Filip/repos/master/sc_data/test' )

#make things prettier

## Neftel 10X ##
path <- paste0(execution_path, '/../../datasets_sctargettranslator/Neftel2019')
setwd(path)

if( !'r_files' %in% dir(path)){
  #create folder and put data there
  stop(paste0('please check that you have a folder called r_files in ', path))
}

library(Matrix)
if(!'neftel_mn_group1' %in% list.files(paste0(path,'/r_files/'))){
  mn_g1<-readMM(file='./Group1/Exp_data_UMIcounts_10X.mtx')
  mn_g1<-as.matrix(mn_g1)

  cells<-read.table(file='./Group1/Cells_10X.txt',sep=' ',header=TRUE,stringsAsFactors = FALSE)
  genes<-read.table(file='./Group1/Genes_10X.txt',sep='\t', header=FALSE,stringsAsFactors = FALSE)
  genes<-genes[,1]

  rownames(cells)<-cells[,1]

  rownames(mn_g1)<-genes
  colnames(mn_g1)<-cells$cell_name

  saveRDS(mn_g1, file = paste0(path, '/r_files/', "neftel_mn_group1"))
  rm(cells, genes)

}else{
  mn_g1 <- readRDS(file = paste0(path, '/r_files/', "neftel_mn_group1"))
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


