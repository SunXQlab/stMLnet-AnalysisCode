#############
## library ##
#############

library(dplyr)
library(ggsci)
library(ggplot2)
library(igraph)
library(plotrix)
library(ggraph)
library(tidyverse)
library(CellChat)
library(iTALK)
library(circlize)

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
LRTG_detail <- na.omit(LRTG_detail)
colnames(LRTG_detail) <- c('cell_from','cell_to','n_LRs','n_TGs','IM','IM_norm','mean_IM','mean_IM_norm')

print(unique(LRTG_detail$cell_from))
print(unique(LRTG_detail$cell_to))

cts_less = unique(LRTG_detail$cell_from)[which(!unique(LRTG_detail$cell_from) %in% unique(LRTG_detail$cell_to))]
a = data.frame(cell_from = cts_less,
               cell_to = cts_less,
               n_LRs = rep(0,length(cts_less)),
               n_TGs = rep(0,length(cts_less)),
               IM = rep(0,length(cts_less)),
               IM_norm = rep(0,length(cts_less)),
               mean_IM = rep(0,length(cts_less)),
               mean_IM_norm = rep(0,length(cts_less)))

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
  
  pdf(paste0(plotdir,"cellchat_networkPlot_","_",key,".pdf"),height = 6,width =6)
  netVisual_circle(tmeTab, color.use = colordb,vertex.weight = rowSums(tmeTab),
                   alpha.edge = 1,edge.width.max= 8,
                   weight.scale = T, label.edge= F, title.name = "Number of interactions",
                   arrow.width = 0.5,arrow.size = 0.5,
                   text.x = 15,text.y = 1.5)
  dev.off()
  
}

######################
## LR activity Plot ##
######################

inputdir <- './runModel/'
cts <- c("Alveoli","Macrophage","Monocyte")

for (ct in cts) {
  
  message(paste0('running jobs: ', ct))
  files = list.files(inputdir)
  files <- files[lapply(CPs, function(cp){grep(cp,files)}) %>% unlist() %>% unique()]
  files = files[grep(paste0('-',ct,'.rds'),files)]
  
  df_LRTGscore = lapply(files, function(file){
    
    print(file)
    LRS_score = readRDS(paste0(inputdir,file))[[1]]
    LRS_score_merge = do.call('cbind',LRS_score) %>% .[,!duplicated(colnames(.))]
    
    cp_pair = strsplit(file,'[_\\.]')[[1]][3]
    # file <- gsub('-','_',file)
    df_LigRec <- data.frame(
      source = colnames(LRS_score_merge) %>% gsub('_.*','',.),
      target = colnames(LRS_score_merge) %>% gsub('.*_','',.),
      LRpair = colnames(LRS_score_merge),
      count = colMeans(LRS_score_merge),
      source_group = strsplit(cp_pair,'[-]')[[1]][1],
      target_group =strsplit(cp_pair,'[-]')[[1]][2]
    )
    
    df_LigRec <- df_LigRec[df_LigRec$count > 6,]
    df_LigRec <- df_LigRec[order(df_LigRec$count,decreasing=T),]
    if (dim(df_LigRec)[1] >= 8){
      # select top 20
      df_LigRec <- df_LigRec[1:5,]
    }
    # 
    # df_LigRec
  }) %>% do.call('rbind',.)
  
  if(!is.null(df_LRTGscore)){
    
    # input
    res <- data.frame(ligand = df_LRTGscore$source,receptor = df_LRTGscore$target,
                      cell_from = df_LRTGscore$source_group,cell_to = df_LRTGscore$target_group, 
                      count = df_LRTGscore$count, comm_type = "growth factor")
    res <- res[order(res$count,decreasing=T),] 
    # plot
    if (dim(res)[1] > 30){
      # select top 20
      res <- res[1:30,]
    }
    
    colordb <- mycolor_ct[which(names(mycolor_ct) %in% c(unique(res$cell_from),unique(res$cell_to)))]
    scales::show_col(colordb)
    print(colordb)
    
    pdf(paste0(plotdir,"ChordPlot_LRscore_","receiver_",ct,".pdf"),height = 4.5,width = 4.5)
    LRPlot(res,datatype='mean count',
           cell_col=colordb,
           link.arr.lwd=res$count,
           link.arr.col="#696969", # "#808080"
           link.arr.width=0.2,
           track.height_1 = uh(1,"mm"),
           track.height_2 = uh(11,"mm"),
           annotation.height_1 = 0.015,
           annotation.height_2 = 0.01,
           text.vjust = "0.5cm")
    dev.off()
    
  }
}


