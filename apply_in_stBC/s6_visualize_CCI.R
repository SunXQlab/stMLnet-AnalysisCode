#############
## library ##
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

setwd('E:/stMLnet/apply_in_stBC/')

source('../code/code.R')

###########
## color ##
###########

# celltype

scales::show_col(pal_nejm(palette = "default", alpha = 0.8)(8))
mycolors_nejm <- pal_nejm(palette = "default", alpha = 0.8)(8)

celltype <- c("Malignant","Macrophage","Stroma","Bcell","Endothelial","Epithelial","Tcell")
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
df_cellpair <- gsub('LRTG_im_clean_|.rds','',files) %>% strsplit(.,"-") %>% do.call('rbind',.) %>% as.data.frame()
LRTG_detail <- cbind(df_cellpair,LRTG_detail)
colnames(LRTG_detail) <- c('cell_from','cell_to','n_LRs','n_TGs','IM','IM_norm','mean_IM','mean_IM_norm')

colodb = pal_nejm(palette = "default", alpha = 0.7)(7)
names(colodb) <- celltype
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

files = list.files(inputdir)
files = files[grep("-Malignant.rds",files)]
files = files[-grep("TME",files)]

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

# input

df_input <- prepareEdgeBundlingPlotData_V2(df_LRTGscore, do.check = T)
str(df_input$df_edge)

# plot

colodb <- c(mycolor_nt,mycolor_ct)
gtitle <- 'TME_Malignant'
drawEdgeBundlingPlot(df_input,colodb,gtitle,plotdir,p_height = 9,p_width = 8.4)

#################
## LR~TG score ##
#################

inputdir <- './getPIM/'

files = list.files(inputdir)
files = files[grep("_im_",files)]
files = files[grep("-Malignant.rds",files)]
files = files[-grep("TME",files)]

LRTG_im_merge <- lapply(files, function(f){
  
  LRTG_im <- readRDS(paste0(inputdir,f))
  LRTG_im$Sender <- strsplit(f,'-|_')[[1]][4]
  LRTG_im
  # head(LRTG_im)
  
}) %>% do.call('rbind',.)

df_MLnet_long_check <- prepareAlluviumPlotData_V2(lrtg_im = LRTG_im_merge, 
                                                  color.by = 'Sender',
                                                  do.check = TRUE)
head(df_MLnet_long_check)

colodb <- c(mycolor_nt,mycolor_key,mycolor_ct)
gtitle <- 'TME_Malignant'
drawAlluviumPlot(df_MLnet_long_check, colodb = colodb, gtitle = gtitle,
                 wd = plotdir,p_height=9, p_width=18)

###############
## MLnetPlot ##
###############

MLnet <- readRDS('./runscMLnet/Tcell_Malignant/scMLnet.rds')
MLnet$LigRec

LRTG_im <- readRDS("./getPIM/LRTG_im_clean_Tcell-Malignant.rds")
head(LRTG_im)

Key <- c('SEMA4D')
Type <- 'Ligand'
MLnet_key <- prepareMLnetworkPlotData_V3(mlnet=MLnet,lrtg_im=LRTG_im,Key=Key,Type=Type,do.check = T)
str(MLnet_key)

colodb = pal_locuszoom(palette = "default", alpha = 0.5)(4)
names(colodb) <- nodekey
scales::show_col(colodb)

downstream <- 'Target'
gtitle <- 'T_TC_SEMA4D'
drawMLnetworkPlot_V4(mlnet=MLnet_key,colodb=colodb,downstream = downstream,
                     gtitle=gtitle,wd=plotdir,p_height = 4.5,p_width = 7)

Key <- c('TGFB1')
Type <- 'Ligand'
MLnet_key <- prepareMLnetworkPlotData_V3(mlnet=MLnet,lrtg_im=LRTG_im,Key=Key,Type=Type,do.check = T)
str(MLnet_key)

colodb = pal_locuszoom(palette = "default", alpha = 0.5)(4)
names(colodb) <- nodekey
scales::show_col(colodb)

downstream <- 'Target'
gtitle <- 'T_TC_TGFB1'
drawMLnetworkPlot_V4(mlnet=MLnet_key,colodb=colodb,downstream = downstream,
                     gtitle=gtitle,wd=plotdir,p_height = 4.5,p_width = 7)

#############
## Heatmap ##
#############

## load 

LRTG_pim <- readRDS("./getPIM/LRTG_pim_clean_TME-Malignant.rds")
LRTG_pim <- LRTG_pim[LRTG_pim$type == 'Receptor',]
LRTG_pim_spl <- split(LRTG_pim,LRTG_pim$regulator)


res_ORA_GO <- Perf_Enrich(LRTG_pim_spl,Type = 'ORA',DB='GO')
res_ORA_KEGG <- Perf_Enrich(LRTG_pim_spl,Type = 'ORA',DB='KEGG')
res_Enrich <- list(res_ORA_GO = res_ORA_GO,
                   res_ORA_KEGG = res_ORA_KEGG)

## ORA-KEGG

res_ORA_KEGG <- res_Enrich$res_ORA_KEGG
res_ORA_KEGG$ID %>% unique() %>% length()
df_kegg <- res_ORA_KEGG[res_ORA_KEGG$pvalue <= 0.01,]
df_kegg$geneRatio <- lapply(df_kegg$GeneRatio,function(chr){
  x = strsplit(chr,'/')[[1]][1] %>% as.numeric()
  y = strsplit(chr,'/')[[1]][2] %>% as.numeric()
  x/y
}) %>% unlist()
df_kegg$bgRatio <- lapply(df_kegg$BgRatio,function(chr){
  x = strsplit(chr,'/')[[1]][1] %>% as.numeric()
  y = strsplit(chr,'/')[[1]][2] %>% as.numeric()
  x/y
}) %>% unlist()

keep_term <- lapply(unique(df_kegg$Regulator), function(key){
  
  df_kegg[df_kegg$Regulator==key,'ID'] %>% head(.,3)
  
}) %>% unlist() %>% unique()
df_kegg$ONTOLOGY <- 'KEGG'
df_kegg <- df_kegg[df_kegg$ID %in% keep_term,c(13,1:2,5:7,10:11)]
df_kegg <- df_kegg[order(df_kegg$ONTOLOGY,df_kegg$ID),]

## ORA-GO

res_ORA_GO <- res_Enrich$res_ORA_GO
res_ORA_GO$ID %>% unique() %>% length()
df_go <- res_ORA_GO[res_ORA_GO$pvalue<=0.01,]
df_go$geneRatio <- lapply(df_go$GeneRatio,function(chr){
  x = strsplit(chr,'/')[[1]][1] %>% as.numeric()
  y = strsplit(chr,'/')[[1]][2] %>% as.numeric()
  x/y
}) %>% unlist()
df_go$bgRatio <- lapply(df_go$BgRatio,function(chr){
  x = strsplit(chr,'/')[[1]][1] %>% as.numeric()
  y = strsplit(chr,'/')[[1]][2] %>% as.numeric()
  x/y
}) %>% unlist()
df_go <- df_go[df_go$ONTOLOGY=='BP',]

keep_term <- lapply(unique(df_go$Regulator), function(key){
  
  df_go[df_go$Regulator==key,'ID'] %>% head(.,3)
  
}) %>% unlist() %>% unique()
df_go <- df_go[df_go$ID %in% keep_term,c(1:3,6:8,11:12),]
df_go <- df_go[order(df_go$ONTOLOGY,df_go$ID),]

## plot

df_plot <- rbind(df_go,df_kegg)
df_plot <- df_plot[df_plot$ONTOLOGY %in% c('KEGG','BP'),]
df_plot$Description <- factor(df_plot$Description, levels = unique(df_plot$Description))
df_plot$ONTOLOGY <- factor(df_plot$ONTOLOGY, levels = c('BP','MF','CC','KEGG'))
anno_x_loca <- lapply(unique(df_plot$ONTOLOGY), function(key){
  
  df <- df_plot[!duplicated(df_plot$ID),]
  grep(key,df$ONTOLOGY) %>% max()
  
}) %>% unlist()
anno_x_loca <- anno_x_loca+0.5
names(anno_x_loca) <- unique(df_plot$ONTOLOGY)

df_plot$Description <- gsub(
  'IgG immunoglobulin transcytosis in epithelial cells mediated by FcRn immunoglobulin receptor',
  'IgG immunoglobulin transcytosis in epithelial cells',df_plot$Description
)

pt_merge <- ggplot(data=df_plot,aes(x=Description,y=Regulator,size=geneRatio,col=p.adjust))+
  geom_hline(aes(x=Description,y=Regulator,yintercept = 1:nrow(df_plot)),size= 1.5,colour= "#E4EDF2",alpha= .5)+
  geom_vline(aes(x=Description,y=Regulator,xintercept = anno_x_loca[1]),size=0.5,linetype= "dashed")+
  geom_point()+ coord_flip() +
  scale_color_material('pink',reverse = T, alpha = 0.8) + 
  scale_size_continuous(range = c(1,4)) + 
  theme_bw() + labs(title='Function Enrichment Analysis of CCI',x='',y='') +
  theme(
    plot.title = element_text(hjust = 0.5,size = 12),
    panel.background = element_blank(),
    legend.key = element_blank(), 
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8), 
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank(), 
    legend.position = 'right',
    legend.direction = 'vertical'
  )
pt_merge

png(paste0(plotdir,"Enrichment_BP_KEGG_TEM_Malignnat.png"), 
    width = 12, height = 5, units = 'in', res = 1000)
pt_merge
dev.off()

pdf(paste0(plotdir,"Enrichment_BP_KEGG_TEM_Malignnat.pdf"), 
    width = 12, height = 5)
pt_merge
dev.off()
