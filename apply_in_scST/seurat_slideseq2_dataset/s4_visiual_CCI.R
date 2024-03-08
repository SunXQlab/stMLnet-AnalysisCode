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
library(CellChat)
library(tidyverse)
library(iTALK)
library(circlize)
library(org.Hs.eg.db)

rm(list=ls())
gc()
setwd("/home/cjy/project/seurat_slideseq2_dataset/")

source('../code/code.R')

###########
## color ##
###########

## load ####

res_path <- "/home/yll/cell_cell_interaction/stMLnet_cjy/apply_in_scST/seurat_slideseq2_dataset/"
load(paste0(res_path,'seurat_slideseq2_output.rda'))

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

inputdir <- paste0(res_path,'getPIM/')
files <- list.files(inputdir)[grep('_im_',list.files(inputdir))]
CPs <- gsub('LRTG_im_clean_|LRTG_pim_clean_|\\.rds','',files)

#################
## NetworkPlot ##
#################

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
  
  pdf(paste0(plotdir,"cellchat_networkPlot_",key,".pdf"),height =5,width = 5)
  netVisual_circle(tmeTab, color.use = colordb,vertex.weight = rowSums(tmeTab),alpha.edge = 0.6, 
                   weight.scale = T, edge.width.max= 4,label.edge= F, title.name = "Number of interactions",
                   arrow.width = 0.4,arrow.size = 0.4,
                   text.x = 15,text.y = 1.5)
  dev.off()
  
}

######################
## LR activity Plot ##
######################

inputdir <- paste0(res_path,'runModel/')

cts <- c("Interneuron","CA1-Principal-cells")

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
    
    min(res$count)
    
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
           text.vjust = "0.4cm")
    dev.off()
    
  }
}


#################
## LR~TG score ##
#################

CP <- paste0("Entorhinal-cortex","_","CA1-Principal-cells")

# load MLnet
wd <- paste0(res_path,'runscMLnet/')
MLnet <- readRDS(paste0(wd,CP,"/scMLnet.rds"))
MLnet$LigRec 

df_ligrec = data.frame(Ligand = MLnet$LigRec$source, Receptor = MLnet$LigRec$target) 
df_rectf = data.frame(Receptor = MLnet$RecTF$source, TF = MLnet$RecTF$target)
df_tftar = data.frame(TF = MLnet$TFTar$source, Target = MLnet$TFTar$target)

df_mlnet = df_ligrec %>% merge(., df_rectf, by = 'Receptor') %>% 
  merge(., df_tftar, by = 'TF') %>% 
  dplyr::select(Ligand, Receptor, TF, Target) %>% 
  arrange(Ligand, Receptor)
df_mlnet$LRpair <- paste0(df_mlnet$Ligand,"_",df_mlnet$Receptor)

# load LRTG_score
files = list.files(inputdir)
files = files[grep("_im_",files)]
file_cp = files[grep(paste0(CP,".rds"),files)]

LRTG_im_merge <- lapply(file_cp, function(f){
  
  LRTG_im <- readRDS(paste0(inputdir,f))
  LRTG_im$Sender <- strsplit(f,'-|_')[[1]][4]
  LRTG_im
  # head(LRTG_im)
  
}) %>% do.call('rbind',.)

LRTG_im_merge$LRTG <- paste0(LRTG_im_merge$LRpair,"_",LRTG_im_merge$Target)
df_mlnet$LRTG <- paste0(df_mlnet$LRpair,"_",df_mlnet$Target)
LRTG_im_merge$TF <- df_mlnet$TF[which(LRTG_im_merge$LRTG %in% df_mlnet$LRTG)]
LRTG_im_merge <- LRTG_im_merge[,-8]

df_MLnet_long1 <- prepareAlluviumPlotData(lrtg_im = LRTG_im_merge, 
                                          color.by = 'Sender',
                                          do.check = FALSE)

df_MLnet_long1$Nodekey <- factor(df_MLnet_long1$Nodekey,levels = c("Ligand", "Receptor","TF", "Target" ))

pt <-ggplot(df_MLnet_long1,
            aes(x = Nodekey, stratum = Node, alluvium = ID,
                y = Score, fill = Nodekey,label = Node)) +
  scale_x_discrete(expand = c(.1, .1)) +
  scale_fill_manual(values = mycolor_key) + 
  geom_flow() +
  geom_stratum(alpha = .5) +
  geom_text(stat = "stratum", size = 3) +
  theme_minimal()+
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 11),
    plot.title = element_text(size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5),
    panel.grid = element_blank(),
    legend.position = 'none'
  ) +
  ggtitle("Multilayer Network")
pt

## save

pdf(paste0(plotdir,'AlluviumPlot-LRTFTG_',CP,'.pdf'),height = 4,width = 5.5)
print(pt)
dev.off()

################
#  pl p-value  #
################
wd <- paset0(res_path,"/runscMLnet/")
files <- list.files(wd)[grep('TME_',list.files(wd),invert = TRUE)]

LRI_allpval <- lapply(files, function(f) {
  cat("Processing file:", f, "\n")
  
  LRI_pval <- readRDS(paste0(wd, f, "/cellpair_LRI_pval.rds"))
  print(str(LRI_pval))  # Print the structure of LRI_pval
  
  if (length(LRI_pval) > 0) {
    LRI_pval$Sender <- strsplit(f, "_")[[1]][1]
    LRI_pval$Receiver <- strsplit(f, "_")[[1]][2]
    return(LRI_pval)
  } else {
    # If LRI_pval is empty, return a placeholder or handle it as needed
    cat("Warning: Empty LRI_pval for file", f, "\n")
    return(NULL)  # or return an empty data frame, depending on your needs
  }
}) %>% do.call('rbind', .)

sender <- "Entorhinal-cortex"
LRI_pval_Inter <- LRI_allpval[LRI_allpval$Sender == sender,]
df_plot <- data.frame(cellpair = paste0(LRI_pval_Inter$Sender,"_",LRI_pval_Inter$Receiver),
                      LRpair = paste0(LRI_pval_Inter$source,"_",LRI_pval_Inter$target),
                      pval = LRI_pval_Inter$pval)

df_plot <- df_plot[order(df_plot$pval,decreasing=F),] 

# load LR signaling score
ct <- "Entorhinal-cortex"
inputdir <- paste0(res_path,'runModel/')
files = list.files(inputdir)
files = files[grep(paste0(ct),files)]  #[c(1,2,4)]

df_LRTGscore = lapply(files, function(file){
  
  print(file)
  LRS_score = readRDS(paste0(inputdir,file))[[1]]
  LRS_score_merge = do.call('cbind',LRS_score)
  if (length(unique(colnames(LRS_score_merge))) == 1){
    LRpair <- unique(colnames(LRS_score_merge))
    LRS_score_merge = LRS_score_merge[,1] 
    
    # file <- gsub('-','_',file)
    df_LigRec <- data.frame(
      source = LRpair %>% gsub('_.*','',.),
      target = LRpair %>% gsub('.*_','',.),
      LRpair = LRpair,
      count = mean(LRS_score_merge),
      source_group = strsplit(file,'[_\\.]')[[1]][3],
      target_group = strsplit(file,'[_\\.]')[[1]][4])
    
  }else{
    LRS_score_merge = do.call('cbind',LRS_score) %>% .[,!duplicated(colnames(.))]
    
    # file <- gsub('-','_',file)
    df_LigRec <- data.frame(
      source = colnames(LRS_score_merge) %>% gsub('_.*','',.),
      target = colnames(LRS_score_merge) %>% gsub('.*_','',.),
      LRpair = colnames(LRS_score_merge),
      count = colMeans(LRS_score_merge),
      source_group = strsplit(file,'[_\\.]')[[1]][3],
      target_group = strsplit(file,'[_\\.]')[[1]][4])
  }
  
}) %>% do.call('rbind',.)

df_LRTGscore$pval <- NA
for (lr in df_LRTGscore$LRpair){
  pos1 <- which(df_plot$LRpair %in% lr)
  pos2 <- which(df_LRTGscore$LRpair %in% lr)
  df_LRTGscore$pval[pos2] <- df_plot$pval[pos1]
}

df <- data.frame(cellpair = paste0(df_LRTGscore$source_group,"_",df_LRTGscore$target_group),
                 LRpair = df_LRTGscore$LRpair,
                 pval = df_LRTGscore$pval,
                 count = df_LRTGscore$count)

# colordb <- mycolor_ct[which(names(mycolor_ct) %in% c(unique(df_plot$sender),unique(df_plot$receiver)))]
# scales::show_col(colordb)
# print(colordb)

p1 <- ggplot(df, aes(x = cellpair, y = LRpair, color = pval, size = count)) +
  geom_point(pch = 16) +
  scale_color_gradient(low = "red", high = "yellow")+
  #theme_linedraw() + 
  #theme(panel.grid.major = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, size = 10,hjust= NULL, vjust = NULL),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_x_discrete(position = "bottom")
p1

pdf(paste0(plotdir,"bubble_pval_","sender_",ct,".pdf"),height = 5.5,width = 5)
p1
dev.off()
