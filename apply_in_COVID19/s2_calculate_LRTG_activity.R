#############
## library ##
#############

library(dplyr)
library(Seurat)
library(SeuratWrappers)

rm(list=ls())
gc()
setwd("E:/stMLnet/apply_in_COVID19/")

source('../code/code.R')

##########
## main ##
##########

## load data

st_rctd <- readRDS('./input/st_rctd.rds')

annoMat <- data.frame(colnames(st_covid_rctd),st_covid_rctd$Cluster)
colnames(annoMat) <- c('Barcode','Cluster')

seed <- 4321
exprMat = st_covid_rctd@assays$SCT@data
exprMat.Impute <- run_Imputation(exprMat, use.seed = T,seed = seed)

locaMat <- data.frame(st_covid_rctd@meta.data[,c('xdim','ydim')])
rownames(locaMat) <- colnames(st_covid_rctd)
distMat <- as.matrix(dist(locaMat))
distMat[distMat==0] <- 1

## load network

wd <- "./runscMLnet/"
files = list.files(wd)
files_tumor = files[grep(".csv",files,invert = T)]
mulNetList = lapply(files_tumor, function(file_tumor){
  
  readRDS(paste0(wd,file_tumor,"/scMLnet.rds"))
  
})
names(mulNetList) = files_tumor
mulNetList = mulNetList[!unlist(lapply(mulNetList, function(mulnet){nrow(mulnet$LigRec)==0}))]
str(mulNetList,max.level = 2)

## Sender-Receiver

for (cellpair in names(mulNetList)) {
  
  Receiver <- gsub('.*_','',cellpair)
  Sender <- gsub('_.*','',cellpair)
  
  LRTG_allscore = calculate_LRTG_allscore_V2(exprMat = exprMat.Impute, distMat = distMat, annoMat = annoMat,
                                             mulNetList = mulNetList, Receiver = Receiver, Sender = Sender)
  
  if(length(LRTG_allscore$LRs_score)!=0){
    
    saveRDS(LRTG_allscore, paste0("./runModel/LRTG_allscore_",Sender,'-',Receiver,'.rds'))
    
  }
  
}

