#############
## library ##
#############

library(dplyr)
library(ggsci)
library(ggplot2)
library(igraph)
library(plotrix)
library(ggraph)
library(ggalluvial)

rm(list=ls())
gc()
setwd("/home/cjy/project/giotto_merfish_dataset/giotto_merfish_dataset_layer9/")

source('../code/code.R')

###########
## color ##
###########

## load ####

load('./giotto_merfish_output.rda')

# celltype

celltype <- unique(df_anno$Cluster)

scales::show_col(pal_igv(palette = "default", alpha = 0.8)(15))
mycolors_nejm <- pal_igv(palette = "default", alpha = 0.8)(15)

mycolor_ct <- mycolors_nejm[1:length(celltype)]
names(mycolor_ct) <- celltype
scales::show_col(mycolor_ct)

# nodekey

scales::show_col(pal_locuszoom(palette = "default", alpha = 0.8)(7))
mycolors_locus <- pal_locuszoom(palette = "default", alpha = 0.8)(7)

nodekey <- c("Ligand","Receptor","TF","Target")
mycolor_key <- mycolors_locus[1:4]
names(mycolor_key) <- nodekey
scales::show_col(mycolor_key)

# nodetype

scales::show_col(pal_locuszoom(palette = "default", alpha = 0.8)(7))
mycolors_locus <- pal_locuszoom(palette = "default", alpha = 0.8)(7)

nodetype <- c("cell","Sender","Receiver")
mycolor_nt <- mycolors_locus[1:3]
names(mycolor_nt) <- nodetype
scales::show_col(mycolor_nt)

#############
## workdir ##
#############

plotdir = './visualize_CCI/'
dir.create(plotdir,recursive = T)

inputdir <- './getPIM/'
files <- list.files(inputdir)[grep('_im_',list.files(inputdir))]
CPs <- gsub('LRTG_im_clean_|LRTG_pim_clean_|\\.rds','',files)

#################
## NetworkPlot ##
#################

inputdir <- './getPIM/'

files <- list.files(inputdir)[grep('_im_',list.files(inputdir))]
files <- files[lapply(CPs, function(cp){grep(cp,files)}) %>% unlist() %>% unique()]

LRTG_detail <- lapply(files, function(f){

  LRTG_im <- readRDS(paste0(inputdir,f))
  LRTG_im <- na.omit(LRTG_im)
  c(length(unique(LRTG_im$LRpair)),length(unique(LRTG_im$Target)),
    sum(LRTG_im$IM),sum(LRTG_im$im_norm),
    sum(LRTG_im$IM)/nrow(LRTG_im),sum(LRTG_im$im_norm)/nrow(LRTG_im))
  
}) %>% do.call('rbind',.) %>% as.data.frame()
df_cellpair <- gsub('LRTG_im_clean_|.rds','',files) %>% strsplit(.,"_") %>% do.call('rbind',.) %>% as.data.frame()
LRTG_detail <- cbind(df_cellpair,LRTG_detail)
colnames(LRTG_detail) <- c('cell_from','cell_to','n_LRs','n_TGs','IM','IM_norm','mean_IM','mean_IM_norm')

colodb = mycolor_ct
scales::show_col(colodb)

for (key in colnames(LRTG_detail)[3:8]) {
  
  tmeTab <- LRTG_detail[,c('cell_from','cell_to',key)]
  colnames(tmeTab) <- c('cell_from','cell_to','n')
  
  png(paste0(plotdir,"networkPlot_",key,".png"),
      height = 7,width = 9, units = 'in', res = 300)
  DrawCellComm(tmeTab,colodb = colodb,gtitle = 'CCI')
  dev.off()
  
  pdf(paste0(plotdir,"networkPlot_",key,".pdf"),height = 7,width = 9)
  DrawCellComm(tmeTab,colodb = colodb,gtitle = 'CCI')
  dev.off()
  
}

######################
## LR activity Plot ##
######################

inputdir <- './runModel/'

for (ct in celltype) {
  
  message(paste0('running jobs: ', ct))
  files = list.files(inputdir)
  files <- files[lapply(CPs, function(cp){grep(cp,files)}) %>% unlist() %>% unique()]
  files = files[grep(paste0('_',ct,'.rds'),files)]
  
  df_LRTGscore = lapply(files, function(file){
    
    print(file)
    LRS_score = readRDS(paste0(inputdir,file))[[1]]
    LRS_score_merge = do.call('cbind',LRS_score) %>% .[,!duplicated(colnames(.))]
    
    # file <- gsub('-','_',file)
    df_LigRec <- data.frame(
      source = colnames(LRS_score_merge) %>% gsub('_.*','',.),
      target = colnames(LRS_score_merge) %>% gsub('.*_','',.),
      LRpair = colnames(LRS_score_merge),
      count = colMeans(LRS_score_merge),
      source_group = strsplit(file,'[_\\.]')[[1]][3],
      target_group = strsplit(file,'[_\\.]')[[1]][4]
    )
    
  }) %>% do.call('rbind',.)
  
  if(!is.null(df_LRTGscore)){
    
    # input
    
    df_input <- prepareEdgeBundlingPlotData(df_LRTGscore, do.check = T)
    str(df_input$df_edge)
    
    # plot
    
    colodb <- c(mycolor_nt,mycolor_ct)
    gtitle <- paste0('sender_',ct)
    drawEdgeBundlingPlot(df_input,colodb,gtitle,plotdir,p_height = 8,p_width = 7.5) 
    
  }
  
}

#################
## LR~TG score ##
#################

inputdir <- './getPIM/'

for (ct in celltype) {
  
  files = list.files(inputdir)
  files = files[grep("_im_",files)]
  files = files[grep(paste0("_",ct,".rds"),files)]
  if(length(grep("TME",files))>0) files = files[-grep("TME",files)]
  
  if(length(files)>0){
    
    LRTG_im_merge <- lapply(files, function(f){
      
      LRTG_im <- readRDS(paste0(inputdir,f))
      LRTG_im$Sender <- strsplit(f,'-|_')[[1]][4]
      LRTG_im
      # head(LRTG_im)
      
    }) %>% do.call('rbind',.)
    
    df_MLnet_long_check <- prepareAlluviumPlotData(lrtg_im = LRTG_im_merge, 
                                                      color.by = 'Sender',
                                                      do.check = TRUE)
    head(df_MLnet_long_check)
    
    colodb <- c(mycolor_nt,mycolor_key,mycolor_ct)
    gtitle <- paste0('Sender_',ct)
    drawAlluviumPlot(df_MLnet_long_check, colodb = colodb, gtitle = gtitle,
                     wd = plotdir,p_height=9, p_width=18)
    
  }
  
}

###############
## MLnetPlot ##
###############

## Astrocyte_Endothelial ####

MLnet <- readRDS('./runscMLnet/Astrocyte_Endothelial/scMLnet.rds')
MLnet$LigRec

LRTG_im <- readRDS("./getPIM/LRTG_im_clean_Astrocyte_Endothelial.rds")
head(LRTG_im[order(LRTG_im$im_norm,decreasing = T),])

Key <- c('CCK','PNOC','TRH')
Type <- 'Ligand'
MLnet_key <- prepareMLnetworkPlotData(mlnet=MLnet,lrtg_im=LRTG_im,Key=Key,Type=Type,do.check = T)
str(MLnet_key)

colodb = pal_locuszoom(palette = "default", alpha = 0.5)(4)
names(colodb) <- nodekey
scales::show_col(colodb)

downstream <- 'Target'
gtitle <- 'Astrocyte-Endothelial_Ligand'
drawMLnetworkPlot(mlnet=MLnet_key,colodb=colodb,downstream = downstream,
                     gtitle=gtitle,wd=plotdir,p_height = 4.5,p_width = 7)

Key <- c('CRHR1','CRHR2')
Type <- 'Receptor'
MLnet_key <- prepareMLnetworkPlotData(mlnet=MLnet,lrtg_im=LRTG_im,Key=Key,Type=Type,do.check = T)
str(MLnet_key)

colodb = pal_locuszoom(palette = "default", alpha = 0.5)(4)
names(colodb) <- nodekey
scales::show_col(colodb)

downstream <- 'Target'
gtitle <- 'Astrocyte-Endothelial_Receptor'
drawMLnetworkPlot(mlnet=MLnet_key,colodb=colodb,downstream = downstream,
                     gtitle=gtitle,wd=plotdir,p_height = 4.5,p_width = 7)

## Inhibitory_Endothelial ####

MLnet <- readRDS('./runscMLnet/Inhibitory_Endothelial/scMLnet.rds')
MLnet$LigRec

LRTG_im <- readRDS("./getPIM/LRTG_im_clean_Inhibitory_Endothelial.rds")
head(LRTG_im[order(LRTG_im$im_norm,decreasing = T),])

Key <- c('CRH','GAL','TAC1','OXT')
Type <- 'Ligand'
MLnet_key <- prepareMLnetworkPlotData(mlnet=MLnet,lrtg_im=LRTG_im,Key=Key,Type=Type,do.check = T)
str(MLnet_key)

colodb = pal_locuszoom(palette = "default", alpha = 0.5)(4)
names(colodb) <- nodekey
scales::show_col(colodb)

downstream <- 'Target'
gtitle <- 'Inhibitory-Endothelial_Ligand'
drawMLnetworkPlot(mlnet=MLnet_key,colodb=colodb,downstream = downstream,
                     gtitle=gtitle,wd=plotdir,p_height = 4.5,p_width = 7)

Key <- c('OPRD1','OPRK1','OPRL1')
Type <- 'Receptor'
MLnet_key <- prepareMLnetworkPlotData(mlnet=MLnet,lrtg_im=LRTG_im,Key=Key,Type=Type,do.check = T)
str(MLnet_key)

colodb = pal_locuszoom(palette = "default", alpha = 0.5)(4)
names(colodb) <- nodekey
scales::show_col(colodb)

downstream <- 'Target'
gtitle <- 'Inhibitory-Endothelial_Receptor'
drawMLnetworkPlot(mlnet=MLnet_key,colodb=colodb,downstream = downstream,
                     gtitle=gtitle,wd=plotdir,p_height = 4.5,p_width = 7)

## Excitatory_Endothelial ####

MLnet <- readRDS('./runscMLnet/Excitatory_Endothelial/scMLnet.rds')
MLnet$LigRec

LRTG_im <- readRDS("./getPIM/LRTG_im_clean_Excitatory_Endothelial.rds")
head(LRTG_im[order(LRTG_im$im_norm,decreasing = T),])

Key <- c('CRH','TRH','OXT')
Type <- 'Ligand'
MLnet_key <- prepareMLnetworkPlotData(mlnet=MLnet,lrtg_im=LRTG_im,Key=Key,Type=Type,do.check = T)
str(MLnet_key)

colodb = pal_locuszoom(palette = "default", alpha = 0.5)(4)
names(colodb) <- nodekey
scales::show_col(colodb)

downstream <- 'Target'
gtitle <- 'Excitatory-Endothelial_Ligand'
drawMLnetworkPlot(mlnet=MLnet_key,colodb=colodb,downstream = downstream,
                     gtitle=gtitle,wd=plotdir,p_height = 4.5,p_width = 7)

Key <- c('OPRL1','TACR1','TACR3')
Type <- 'Receptor'
MLnet_key <- prepareMLnetworkPlotData(mlnet=MLnet,lrtg_im=LRTG_im,Key=Key,Type=Type,do.check = T)
str(MLnet_key)

colodb = pal_locuszoom(palette = "default", alpha = 0.5)(4)
names(colodb) <- nodekey
scales::show_col(colodb)

downstream <- 'Target'
gtitle <- 'Excitatory-Endothelial_Receptor'
drawMLnetworkPlot(mlnet=MLnet_key,colodb=colodb,downstream = downstream,
                     gtitle=gtitle,wd=plotdir,p_height = 4.5,p_width = 7)
