#############
## prepare ##
#############

library(Matrix)
library(dplyr)
library(ggsci)
library(ggplot2)
library(igraph)
library(plotrix)
library(ggraph)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggalluvial)

rm(list = ls())
gc()

setwd("~/cell_cell_interaction/apply_in_stGBM/UKF_304_T/")

source('../code/code.R')

###########
## color ##
###########

scales::show_col(pal_nejm(palette = "default", alpha = 0.7)(8))
mycolors_nejm <- pal_nejm(palette = "default", alpha = 0.7)(8)

celltype <- c("Malignant","macrophages","Astrocytes","oligodendrocytes","endothelial")
mycolor_ct <- mycolors_nejm[1:length(celltype)]
names(mycolor_ct) <- celltype
scales::show_col(mycolor_ct)

scales::show_col(pal_locuszoom(palette = "default", alpha = 0.8)(7))
mycolors_locus <- pal_locuszoom(palette = "default", alpha = 0.8)(7)

nodekey <- c("Ligand","Receptor","TF","Target")
mycolor_key <- mycolors_locus[1:length(nodekey)]
names(mycolor_key) <- nodekey
scales::show_col(mycolor_key)

nodetype <- c("cell","Sender","Receiver")
mycolor_nt <- mycolors_locus[1:3]
names(mycolor_nt) <- nodetype
scales::show_col(mycolor_nt)

#############
## workdir ##
#############

plotdir = './visualize_CCI/'

#################
## NetworkPlot ##
#################

inputdir <- './getPIM/'
files <- list.files(inputdir)[grep('_im_',list.files(inputdir))]
files <- files[!grepl('TME',files)]
LRTG_detail <- lapply(files, function(f){
  
  # f = files[1]
  LRTG_im <- readRDS(paste0(inputdir,f))
  c(length(unique(LRTG_im$LRpair)),length(unique(LRTG_im$Target)),
    sum(LRTG_im$IM),sum(LRTG_im$im_norm),
    sum(LRTG_im$IM)/nrow(LRTG_im),sum(LRTG_im$im_norm)/nrow(LRTG_im))
  
}) %>% do.call('rbind',.) %>% as.data.frame()
df_cellpair <- gsub('LRTG_im_clean_|.rds','',files) %>% strsplit(.,"_") %>% do.call('rbind',.) %>% as.data.frame()
LRTG_detail <- cbind(df_cellpair,LRTG_detail)
colnames(LRTG_detail) <- c('cell_from','cell_to','n_LRs','n_TGs','IM','IM_norm','mean_IM','mean_IM_norm')

for (key in colnames(LRTG_detail)[3:8]) {
  
  # key <- 'n_LRs'
  tmeTab <- LRTG_detail[,c('cell_from','cell_to',key)]
  colnames(tmeTab) <- c('cell_from','cell_to','n')
  
  png(paste0("./visualize_CCI/networkPlot_",key,".png"),
      height = 6,width = 8, units = 'in', res = 300)
  DrawCellComm(tmeTab,mycolor_ct,gtitle = key)
  dev.off()
  
  pdf(paste0("./visualize_CCI/networkPlot_",key,".pdf"),height = 6,width = 8)
  DrawCellComm(tmeTab,mycolor_ct,gtitle = key)
  dev.off()
  
}

############################
## MLnetPlot: directedCCI ##
############################

## TAM_TC

wd <- "./runscMLnet/"
MLnet <- readRDS(paste0(wd,"macrophages_Malignant/scMLnet.rds"))
MLnet$LigRec

wd <- "./getPIM/"
LRTG_im <- readRDS(paste0(wd,"LRTG_im_clean_macrophages_Malignant.rds"))

Key <- c('IGF1')
Type <- 'Ligand'
#MLnet_key <- prepareMLnetworkPlotData_V3(mlnet = MLnet,lrtg_im = LRTG_im,Key=Key,Type=Type,do.check = T)
MLnet_key <- prepareMLnetworkPlotData(mlnet = MLnet,lrtg_im = LRTG_im,Key=Key,Type=Type,do.check = T)
str(MLnet_key)

# Key2 <- c('ITGAV')
# Type2 <- 'Receptor'
# #MLnet_key <- prepareMLnetworkPlotData_V3(mlnet = MLnet,lrtg_im = LRTG_im,Key=Key,Type=Type,do.check = T)
# MLnet_key2 <- prepareMLnetworkPlotData(mlnet = MLnet_key1,lrtg_im = LRTG_im,Key=Key2,Type=Type2,do.check = T)
# str(MLnet_key2)

colodb = pal_locuszoom(palette = "default", alpha = 0.5)(4)
names(colodb) <- nodekey
scales::show_col(colodb)

downstream <- 'Target'
gtitle <- 'TAM_TC_IGF1_ITGAV'
wd <- './visualize_CCI/'
# drawMLnetworkPlot_V4(mlnet=MLnet_key,colodb=colodb,downstream = downstream,
#                      gtitle=gtitle,wd=wd,p_height = 4,p_width = 7)
drawMLnetworkPlot(mlnet=MLnet_key,colodb=colodb,downstream = downstream,
                     gtitle=gtitle,wd=wd,p_height = 3.7,p_width = 6)


## TC_TAM

wd <- "./runscMLnet/"
MLnet <- readRDS(paste0(wd,"Malignant_macrophages/scMLnet.rds"))
MLnet$LigRec

wd <- "./getPIM/"
LRTG_im <- readRDS(paste0(wd,"LRTG_im_clean_Malignant_macrophages.rds"))

Key <- c('IL4R','CSF1R')
Type <- 'Receptor'
#MLnet_key <- prepareMLnetworkPlotData_V3(mlnet = MLnet,lrtg_im = LRTG_im,Key=Key,Type=Type,do.check = T)
MLnet_key <- prepareMLnetworkPlotData(mlnet = MLnet,lrtg_im = LRTG_im,Key=Key,Type=Type,do.check = T)
str(MLnet_key)

colodb = pal_locuszoom(palette = "default", alpha = 0.5)(4)
names(colodb) <- nodekey
scales::show_col(colodb)

downstream <- 'Target'
gtitle <- 'TC_TAM_IL34_CSF1R'
wd <- './visualize_CCI/'
# drawMLnetworkPlot_V4(mlnet=MLnet_key,downstream=downstream,
#                      colodb=colodb,gtitle=gtitle,wd=wd,
#                      p_height = 4,p_width = 10)
drawMLnetworkPlot(mlnet=MLnet_key,downstream=downstream,
                     colodb=colodb,gtitle=gtitle,wd=wd,
                     p_height = 3.7,p_width = 6)


###############################
## AlluvialPlot: directedCCI ##
###############################

## TME-TC

wd <- './getPIM/'
files <- list.files(wd) %>% .[grep('_im_',.)] %>% .[grep('_Malignant.rds',.)]
LRTG_im_merge <- lapply(files, function(f){
  
  LRTG_im <- readRDS(paste0(wd,f))
  LRTG_im$Sender <- strsplit(f,'-|_')[[1]][4]
  LRTG_im
  # head(LRTG_im)
  
}) %>% do.call('rbind',.)

# df_MLnet_long_check <- prepareAlluviumPlotData_V2(lrtg_im = LRTG_im_merge, 
#                                                   color.by = 'Sender', # Nodekey
#                                                   do.check = TRUE)
df_MLnet_long_check <- prepareAlluviumPlotData(lrtg_im = LRTG_im_merge, 
                                                  color.by = 'Sender', # Nodekey
                                                  do.check = TRUE)
head(df_MLnet_long_check)

colodb <- c(mycolor_nt,mycolor_key,mycolor_ct)
scales::show_col(colodb)

gtitle <- 'TME_Malignant_v3'
wd = './visualize_CCI/'
drawAlluviumPlot(df_MLnet_long_check, colodb = colodb, gtitle = gtitle,
                 wd = wd,p_height=5.5, p_width=8)

## TME-TAM

wd <- './getPIM/'
files <- list.files(wd) %>% .[grep('_im_',.)] %>% .[grep('_macrophages.rds',.)]
LRTG_im_merge <- lapply(files, function(f){
  
  LRTG_im <- readRDS(paste0(wd,f))
  LRTG_im$Sender <- strsplit(f,'-|_')[[1]][4]
  LRTG_im
  # head(LRTG_im)
  
}) %>% do.call('rbind',.)

df_MLnet_long_check <- prepareAlluviumPlotData(lrtg_im = LRTG_im_merge, 
                                                  color.by = 'Sender', # Nodekey
                                                  do.check = TRUE)
head(df_MLnet_long_check)

colodb <- c(mycolor_nt,mycolor_key,mycolor_ct)
scales::show_col(colodb)

gtitle <- 'TME_Macrophages_v3'
wd = './visualize_CCI/'
drawAlluviumPlot(df_MLnet_long_check, colodb = colodb, gtitle = gtitle,
                 wd = wd,p_height=5.5, p_width=8)




