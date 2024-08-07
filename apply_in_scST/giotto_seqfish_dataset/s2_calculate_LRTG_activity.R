#############
#  library  #
#############

library(dplyr)
library(Seurat)
library(SeuratWrappers)

rm(list=ls())
gc()
setwd("/home/cjy/project/giotto_seqfish_dataset/")

source('../code/code.R')

##########
## main ##
##########

## load data

load('./giotto_seqfish_output.rda')

annoMat <- df_anno
head(annoMat)

seed <- 4321
exprMat <- df_norm
exprMat.Impute <- run_Imputation(exprMat, use.seed = T,seed = seed)

locaMat <- df_loca
head(locaMat)

distMat <- as.matrix(dist(locaMat))
distMat[distMat==0] <- 1
distMat[1:4,1:4]

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

wd <- paste0("./runModel/")
dir.create(wd,recursive = T)

for (cellpair in names(mulNetList)) {
  
  Receiver <- gsub('.*_','',cellpair)
  Sender <- gsub('_.*','',cellpair)
  
  LRTG_allscore = calculate_LRTG_allscore_V2(exprMat = exprMat.Impute, distMat = distMat, annoMat = annoMat,
                                             mulNetList = mulNetList, Receiver = Receiver, Sender = Sender)
  
  if(length(LRTG_allscore$LRs_score)!=0){
    
    saveRDS(LRTG_allscore, paste0("./runModel/LRTG_allscore_",Sender,'_',Receiver,'.rds'))
    
  }
  
}

