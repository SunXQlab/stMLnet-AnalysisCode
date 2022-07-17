#############
## library ##
#############

library(dplyr)
library(ggsci)
library(ggplot2)
library(igraph)
library(plotrix)
library(ggraph)

rm(list=ls())
gc()
setwd("./stMLnet/apply_in_COVID19/")

source('../code/code.R')

###########
## color ##
###########

# celltype

st_covid_rctd <- readRDS("input/st_rctd.rds")
celltype <- st_covid_rctd$Cluster %>% unique()

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
df_cellpair <- gsub('LRTG_im_clean_|.rds','',files) %>% strsplit(.,"-") %>% do.call('rbind',.) %>% as.data.frame()
LRTG_detail <- cbind(df_cellpair,LRTG_detail)
colnames(LRTG_detail) <- c('cell_from','cell_to','n_LRs','n_TGs','IM','IM_norm','mean_IM','mean_IM_norm')

colodb = mycolor_ct
scales::show_col(colodb)

key <- 'n_LRs'
tmeTab <- LRTG_detail[,c('cell_from','cell_to',key)]
colnames(tmeTab) <- c('cell_from','cell_to','n')

png(paste0(plotdir,"networkPlot_",key,".png"),
    height = 7,width = 9, units = 'in', res = 300)
DrawCellComm(tmeTab,colodb = colodb,gtitle = 'CCI')
dev.off()

pdf(paste0(plotdir,"networkPlot_",key,".pdf"),height = 7,width = 9)
DrawCellComm(tmeTab,colodb = colodb,gtitle = 'CCI')
dev.off()

######################
## LR activity Plot ##
######################

inputdir <- './runModel/'

for (ct in celltype) {
  
  message(paste0('running jobs: ', ct))
  files = list.files(inputdir)
  files <- files[lapply(CPs, function(cp){grep(cp,files)}) %>% unlist() %>% unique()]
  files = files[grep(paste0('-',ct,'.rds'),files)]
  
  df_LRTGscore = lapply(files, function(file){
    
    print(file)
    LRS_score = readRDS(paste0(inputdir,file))[[1]]
    LRS_score_merge = do.call('cbind',LRS_score) %>% .[,!duplicated(colnames(.))]
    
    file <- gsub('-','_',file)
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



