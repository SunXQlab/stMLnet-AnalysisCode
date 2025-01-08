#############
#  library  #
#############

library(dplyr)
library(ggsci)
library(ggplot2)
library(igraph)
library(plotrix)
library(ggraph)
library(ggalluvial)
library(CellChat)
library(tidyverse)

rm(list=ls())
gc()
setwd("/home/cjy/project/giotto_seqfish_dataset/")

source('../code/code.R')

###########
## color ##
###########

load("~/cell_cell_interaction/stMLnet_cjy/apply_in_scST/giotto_seqfish_dataset/giotto_seqfish_output.rda")

# celltype

celltype <- unique(df_anno$Cluster)

scales::show_col(pal_igv(palette = "default", alpha = 0.8)(15))
mycolors_nejm <- pal_igv(palette = "default", alpha = 0.8)(15)

mycolor_ct <- mycolors_nejm[1:length(celltype)]
names(mycolor_ct) <- celltype
scales::show_col(mycolor_ct)


#############
## workdir ##
#############

plotdir = './visualize/'
dir.create(plotdir,recursive = T)

res_path <- '~/cell_cell_interaction/stMLnet_cjy/apply_in_scST/giotto_seqfish_dataset/'
inputdir <- paste0(res_path,'getPIM/')
files <- list.files(inputdir)[grep('_im_',list.files(inputdir))]
CPs <- gsub('LRTG_im_clean_|LRTG_pim_clean_|\\.rds','',files)

#################
## NetworkPlot ##
#################

res_path <- '~/cell_cell_interaction/stMLnet_cjy/apply_in_scST/giotto_seqfish_dataset/'
inputdir <- paste0(res_path,'getPIM/')

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
LRTG_detail <- na.omit(LRTG_detail)
colnames(LRTG_detail) <- c('cell_from','cell_to','n_LRs','n_TGs','IM','IM_norm','mean_IM','mean_IM_norm')

a <- c("L5-eNeuron","L5-eNeuron",0,0,0,0,0)
LRTG_detail <- rbind(LRTG_detail,a)
LRTG_detail$n_LRs <- as.numeric(LRTG_detail$n_LRs)
LRTG_detail$n_TGs <- as.numeric(LRTG_detail$n_TGs)
LRTG_detail$IM <- as.numeric(LRTG_detail$IM)
LRTG_detail$IM_norm <- as.numeric(LRTG_detail$IM_norm)
LRTG_detail$mean_IM <- as.numeric(LRTG_detail$mean_IM)
LRTG_detail$mean_IM_norm <- as.numeric(LRTG_detail$mean_IM_norm)

for (key in colnames(LRTG_detail)[3:8]) {
  
  tmeTab <- LRTG_detail[,c('cell_from','cell_to',key)] %>% spread(cell_to, key) %>% 
    column_to_rownames(.,var = "cell_from")
  tmeTab[is.na(tmeTab)] <- 0 
  tmeTab <- as.matrix(tmeTab)
  
  colordb <- mycolor_ct[rownames(tmeTab)]
  
  pdf(paste0(plotdir,"cellchat_networkPlot_",key,".pdf"),height =6,width = 6)
  netVisual_circle(tmeTab, color.use = colordb,vertex.weight = rowSums(tmeTab),alpha.edge = 0.6, 
                   weight.scale = T, label.edge= F, title.name = "Number of interactions",
                   arrow.width = 1,arrow.size = 0.3,
                   text.x = 15,text.y = 1.5)
  dev.off()
  
}

######################
## LR activity Plot ##
######################

res_path <- '~/cell_cell_interaction/stMLnet_cjy/apply_in_scST/giotto_seqfish_dataset/'
inputdir <- paste0(res_path,'runModel/')

cts <- c("OPC","Olig")

for (ct in cts) {
  
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
    res <- data.frame(ligand = df_LRTGscore$source,receptor = df_LRTGscore$target,
                      cell_from = df_LRTGscore$source_group,cell_to = df_LRTGscore$target_group, 
                      count = df_LRTGscore$count, comm_type = "growth factor")
    res <- res[order(res$count,decreasing=T),] 
    # plot
    if (dim(res)[1] > 20){
      # select top 20
      res <- res[1:20,]
    }
  
    
    colordb <- mycolor_ct[which(names(mycolor_ct) %in% c(unique(res$cell_from),unique(res$cell_to)))]
    scales::show_col(colordb)
    print(colordb)
    
    pdf(paste0(plotdir,"ChordPlot_v2_LRscore_","receiver_",ct,".pdf"),height = 4.5,width = 4.5)
    LRPlot(res,datatype='mean count',
           cell_col=colordb,
           link.arr.lwd=res$count,
           link.arr.col="#696969", # "#808080"
           link.arr.width=0.25,
           track.height_1 = uh(1,"mm"),
           track.height_2 = uh(11,"mm"),
           annotation.height_1 = 0.015,
           annotation.height_2 = 0.01,
           text.vjust = "0.5cm")
    dev.off()
    
  }
}

#############################
## LR activity Bubble Plot ##
#############################

res_path <- '~/cell_cell_interaction/stMLnet_cjy/apply_in_scST/giotto_seqfish_dataset/'
inputdir <- paste0(res_path,'runModel/')

ct <- c("OPC","Olig")

files = list.files(inputdir)
files <- files[lapply(CPs, function(cp){grep(cp,files)}) %>% unlist() %>% unique()]
files = files[grep(paste0('_',ct[1],'.rds',"|",'_',ct[2],'.rds'),files)]

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
df_LRTGscore$cellpair <- paste0(df_LRTGscore$source_group,"->",df_LRTGscore$target_group)
df_LRTGscore$LRpair <- gsub("_","-",df_LRTGscore$LRpair)

bubble_plot_LRscore(df_LRTGscore,getwd(),save_name=paste0(ct1,"_",ct2))





