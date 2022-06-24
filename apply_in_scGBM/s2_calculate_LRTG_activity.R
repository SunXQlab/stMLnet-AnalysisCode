#############
## library ##
#############

library(dplyr)
library(Seurat)
library(SeuratWrappers)

rm(list=ls())
gc()
setwd("E:/stMLnet/apply_in_scGBM/")

source('../code/code.R')

##########
## load ##
##########

seur <- readRDS("F:/finalVersion/vaild_scRNAseq/2019_cell_scRNAseq/result/seur.rds")
annoMat <- data.frame(Barcode = rownames(seur@meta.data),
                      Cluster = seur$orig.annotation)
distMat <- readRDS("F:/finalVersion/vaild_scRNAseq/input/euclidean_distMat_hvg_alra.rds")
exprMat.Impu <- seur@assays$alra@data

##################
## TEM-Receiver ##
##################

# load

wd = './runscMLnet/'
folders = list.dirs(wd)
mulNetAllList <- lapply(folders[-1], function(fd){
  
  readRDS(paste0(fd,"/scMLnet.rds"))
  
})
names(mulNetAllList) <- stringr::str_split(folders[-1],pattern = "/",simplify = T)[,9] %>% as.vector()
mulNetAllList = mulNetAllList[!unlist(lapply(mulNetAllList, function(mulnet){nrow(mulnet$LigRec)==0}))]

# main

for (cellpair in names(mulNetAllList)) {
  
  message(cellpair)
  Sender = strsplit(cellpair,"_")[[1]][1]
  Receiver = strsplit(cellpair,"_")[[1]][2]
  
  LRTG_allscore_tumor_merge = calculate_LRTG_allscore_V2(
    exprMat = exprMat.Impu, distMat = distMat, annoMat = annoMat,
    mulNetList = mulNetAllList, Receiver = Receiver, Sender = Sender)
  
  saveRDS(LRTG_allscore_tumor_merge,
          paste0("./runModel/LRTG_allscore_",cellpair,".rds"))
  
}

#####################
## TME-macrophages ##
#####################

Receiver = 'macrophages'
Sender = NULL

wd = './runscMLnet/'
folders = list.dirs(wd)
folders = folders[grep('_macrophages',folders)] 
mulNetAllList <- lapply(folders, function(fd){
  
  readRDS(paste0(fd,"/scMLnet.rds"))
  
})
names(mulNetAllList) <- stringr::str_split(folders,pattern = "/",simplify = T)[,4] %>% as.vector()
mulNetAllList = mulNetAllList[!unlist(lapply(mulNetAllList, function(mulnet){nrow(mulnet$LigRec)==0}))]

LRTGscore_TC_TAM <- calculate_LRTG_allscore_V2(
  exprMat = exprMat.Impu, distMat = distMat, annoMat = annoMat,
  mulNetList = mulNetAllList, Receiver = Receiver, Sender = Sender)

saveRDS(LRTGscore_TC_TAM, paste0("./runModel/LRTG_allscore_TME_TAM.rds"))

