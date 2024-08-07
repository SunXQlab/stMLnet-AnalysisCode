#############
## library ##
#############

library(Matrix)
library(dplyr)
library(parallel)
library(Seurat)

rm(list=ls())
gc()

setwd("./stMLnet/apply_in_COVID19/")

source('../code/code.R')

###############
## get MLnet ##
###############

## load

st_rctd <- readRDS('./input/st_rctd.rds')
Ligs_up_list <- readRDS("./input/seur_Ligs_up.rds")
Recs_expr_list <- readRDS("./input/seur_Recs_expr.rds")
ICGs_list <- readRDS("./input/seur_TGs_diff.rds")
Databases <- readRDS('../prior_knowledge/output/Databases.rds')

## data

GCMat <- st_rctd@assays$SCT@data
BarCluTable <- data.frame(colnames(st_rctd),st_rctd$Cluster)
colnames(BarCluTable) <- c('Barcode','Cluster')
clusters <- BarCluTable$Cluster %>% as.character() %>% unique()

## parameters

wd <- paste0("./runscMLnet/")
dir.create(wd,recursive = T)

## prior database

quan.cutoff = 0.98

RecTF.DB <- Databases$RecTF.DB %>% 
  .[.$score > quantile(.$score, quan.cutoff),] %>%
  dplyr::distinct(source, target)

LigRec.DB <- Databases$LigRec.DB %>%
  dplyr::distinct(source, target) %>%
  dplyr::filter(target %in% RecTF.DB$source)

TFTG.DB <- Databases$TFTG.DB %>%
  dplyr::distinct(source, target) %>%
  dplyr::filter(source %in% RecTF.DB$target)

## get multi-layer

for(RecClu in clusters){
  
  LigClus = clusters[clusters != RecClu]
  Output <- matrix(ncol = length(LigClus), nrow = 10) %>% as.data.frame()
  MLnet_list <- list()
  
  for(i in 1:length(LigClus)){
    
    LigClu <- LigClus[i]
    message(paste(LigClu,RecClu,sep = '-'))
    
    MLnet <- mainfunc(LigClu, RecClu, wd)
    MLnet_list[[i]] <- MLnet
    Output[,i] <- c(Ligs_up_list[[LigClu]] %>% length(),
                    Recs_expr_list[[RecClu]] %>% length(),
                    ICGs_list[[RecClu]] %>% length(),
                    nrow(MLnet$LigRec),nrow(MLnet$RecTF),nrow(MLnet$TFTar),
                    ifelse(nrow(MLnet$LigRec)==0,0,MLnet$LigRec$source %>% unique() %>% length()),
                    ifelse(nrow(MLnet$LigRec)==0,0,MLnet$LigRec$target %>% unique() %>% length()),
                    ifelse(nrow(MLnet$TFTar)==0,0,MLnet$TFTar$source %>% unique() %>% length()),
                    ifelse(nrow(MLnet$TFTar)==0,0,MLnet$TFTar$target %>% unique() %>% length()))
    
    
  }
  names(MLnet_list) <- paste(LigClus,RecClu,sep = "_")
  colnames(Output) <- paste(LigClus,RecClu,sep = "_")
  rownames(Output) <- c('Lig_bk','Rec_bk','ICG_bk',
                        "LRpair","RecTFpair","TFTGpair",
                        "Ligand", "Receptor", "TF", "TG")
  
  write.csv(Output, file = paste0(wd,"/TME_",RecClu,".csv"))
  
}
