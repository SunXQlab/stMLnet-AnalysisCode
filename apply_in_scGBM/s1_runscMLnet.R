#############
## library ##
#############

library(dplyr)
library(Seurat)

rm(list=ls())
gc()
setwd("E:/stMLnet/apply_in_scGBM/")

source('../code/code.R')

###############
## get MLnet ##
###############

## load

load("./input/input.rda")
GCMat <- exprMat.Impu
BarCluTable <- annoMat
clusters <- BarCluTable$Cluster %>% as.character() %>% unique()

Ligs_up_list <- readRDS("./input/Ligs_up_list.rds")
Recs_expr_list <- readRDS("./input/Recs_expr_list.rds")

## load DEGs

TGs_list <- readRDS('./data/2016Science/TGs_list.rds')
ICGs_list <- TGs_list[grep("RebEP",names(TGs_list))]
names(ICGs_list) <- gsub('RebEP_TAM','macrophages',names(ICGs_list))
names(ICGs_list) <- gsub('RebEP_TC','Malignant',names(ICGs_list))

DEGs_list <- readRDS("./input/DEGs_lists.rds")
ICGs_list[['oligodendrocytes']] <- DEGs_list$oligodendrocytes
ICGs_list[['Tcell']] <- DEGs_list$Tcell

## database

quan.cutoff = 0.98
Databases <- readRDS('../prior_knowledge/output/Databases.rds')

RecTF.DB <- Databases$RecTF.DB %>% 
  .[.$score > quantile(.$score, quan.cutoff),] %>%
  dplyr::distinct(source, target) %>%
  as.data.frame()

LigRec.DB <- Databases$LigRec.DB %>%
  dplyr::distinct(source, target) %>%
  dplyr::filter(target %in% RecTF.DB$source) %>%
  as.data.frame()

TFTG.DB <- Databases$TFTG.DB %>%
  dplyr::distinct(source, target) %>%
  dplyr::filter(source %in% RecTF.DB$target) %>%
  as.data.frame()

## get multi-layer

wd <- './runscMLnet/'
for(RecClu in clusters){
  
  LigClus = clusters[clusters != RecClu]
  Output <- matrix(ncol = length(LigClus), nrow = 10) %>% as.data.frame()
  MLnet_list <- list()
  
  for(i in 1:length(LigClus)){
    
    LigClu <- LigClus[i]
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