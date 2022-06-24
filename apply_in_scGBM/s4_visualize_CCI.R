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

setwd('E:/stMLnet/apply_in_scGBM/')

source('../code/code.R')

###########
## color ##
###########

scales::show_col(pal_nejm(palette = "default", alpha = 0.7)(8))
mycolors_nejm <- pal_nejm(palette = "default", alpha = 0.7)(8)

celltype <- c("Malignant","macrophages","oligodendrocytes","Tcell")
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

inputdir <- 'F:/finalVersion/vaild_scRNAseq/getPIM/DEGs_from_bulk/job/'
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
      height = 7,width = 9, units = 'in', res = 300)
  DrawCellComm(tmeTab,mycolor_ct,gtitle = key)
  dev.off()
  
  pdf(paste0("./visualize_CCI/networkPlot_",key,".pdf"),height = 7,width = 9)
  DrawCellComm(tmeTab,mycolor_ct,gtitle = key)
  dev.off()
  
}

############################
## MLnetPlot: directedCCI ##
############################

## TAM_TC

wd <- "F:/finalVersion/vaild_scRNAseq/runscMLnet/DEGs_from_bulk/DEGs_from_wilcox-cpm/work_wilcox-cpm_logfc2_pval0.1/"
MLnet <- readRDS(paste0(wd,"macrophages_Malignant/scMLnet.rds"))
MLnet$LigRec

wd <- "F:/finalVersion/vaild_scRNAseq/getPIM/DEGs_from_bulk/job/"
LRTG_im <- readRDS(paste0(wd,"LRTG_im_clean_macrophages_Malignant.rds"))

Key <- c('IGF1')
Type <- 'Ligand'
MLnet_key <- prepareMLnetworkPlotData_V3(mlnet = MLnet,lrtg_im = LRTG_im,Key=Key,Type=Type,do.check = T)
str(MLnet_key)

colodb = pal_locuszoom(palette = "default", alpha = 0.5)(4)
names(colodb) <- nodekey
scales::show_col(colodb)

downstream <- 'Target'
gtitle <- 'TAM_TC_IGF1'
wd <- './visualize_CCI/'
drawMLnetworkPlot_V4(mlnet=MLnet_key,colodb=colodb,downstream = downstream,
                     gtitle=gtitle,wd=wd,p_height = 4,p_width = 7)
## TC_TAM

wd <- "F:/finalVersion/vaild_scRNAseq/runscMLnet/DEGs_from_bulk/DEGs_from_wilcox-cpm/work_wilcox-cpm_logfc2_pval0.1/"
MLnet <- readRDS(paste0(wd,"Malignant_macrophages/scMLnet.rds"))
MLnet$LigRec

wd <- "F:/finalVersion/vaild_scRNAseq/getPIM/DEGs_from_bulk/job/"
LRTG_im <- readRDS(paste0(wd,"LRTG_im_clean_Malignant_macrophages.rds"))

Key <- c('IL4R','CSF1R')
Type <- 'Receptor'
MLnet_key <- prepareMLnetworkPlotData_V3(mlnet = MLnet,lrtg_im = LRTG_im,Key=Key,Type=Type,do.check = T)
str(MLnet_key)

colodb = pal_locuszoom(palette = "default", alpha = 0.5)(4)
names(colodb) <- nodekey
scales::show_col(colodb)

downstream <- 'Target'
gtitle <- 'TC_TAM_IL4R_CSF1R'
wd <- './visualize_CCI/'
drawMLnetworkPlot_V4(mlnet=MLnet_key,downstream=downstream,
                     colodb=colodb,gtitle=gtitle,wd=wd,
                     p_height = 4,p_width = 10)

###############################
## AlluvialPlot: directedCCI ##
###############################

## TME-TC

wd <- 'F:/finalVersion/vaild_scRNAseq/getPIM/DEGs_from_bulk/job/'
files <- list.files(wd) %>% .[grep('_im_',.)] %>% .[grep('_Malignant.rds',.)]
LRTG_im_merge <- lapply(files, function(f){
  
  LRTG_im <- readRDS(paste0(wd,f))
  LRTG_im$Sender <- strsplit(f,'-|_')[[1]][4]
  LRTG_im
  # head(LRTG_im)
  
}) %>% do.call('rbind',.)

df_MLnet_long_check <- prepareAlluviumPlotData_V2(lrtg_im = LRTG_im_merge, 
                                                  color.by = 'Sender', # Nodekey
                                                  do.check = TRUE)
head(df_MLnet_long_check)

colodb <- c(mycolor_nt,mycolor_key,mycolor_ct)
scales::show_col(colodb)

gtitle <- 'TME_Malignant'
wd = './visualize_CCI/'
drawAlluviumPlot(df_MLnet_long_check, colodb = colodb, gtitle = gtitle,
                 wd = wd,p_height=9, p_width=18)

## TME-TAM

wd <- './vaild_scRNAseq/getPIM/DEGs_from_bulk/job/'
files <- list.files(wd) %>% .[grep('_im_',.)] %>% .[grep('_macrophages.rds',.)]
LRTG_im_merge <- lapply(files, function(f){
  
  LRTG_im <- readRDS(paste0(wd,f))
  LRTG_im$Sender <- strsplit(f,'-|_')[[1]][4]
  LRTG_im
  # head(LRTG_im)
  
}) %>% do.call('rbind',.)

df_MLnet_long_check <- prepareAlluviumPlotData_V2(lrtg_im = LRTG_im_merge, 
                                                  color.by = 'Sender', # Nodekey
                                                  do.check = TRUE)
head(df_MLnet_long_check)

colodb <- c(mycolor_nt,mycolor_key,mycolor_ct)
scales::show_col(colodb)

gtitle <- 'TME_Macrophages'
wd = 'F:/finalVersion/vaild_scRNAseq/visualize_CCI/AlluvialPlot/'
drawAlluviumPlot(df_MLnet_long_check, colodb = colodb, gtitle = gtitle,
                 wd = wd,p_height=9, p_width=18)

#########################
## function annotation ##
#########################
## TAM-TC ####

LRTG_im <- readRDS("F:/finalVersion/vaild_scRNAseq/getPIM/DEGs_from_bulk/job/LRTG_im_clean_macrophages_Malignant.rds")
LRTG_im_spl <- split(LRTG_im,LRTG_im$Receptor)

t1 <- Sys.time()
res_ORA_GO <- Perf_Enrich(LRTG_im_spl,Type = 'ORA',DB='GO')
res_ORA_KEGG <- Perf_Enrich(LRTG_im_spl,Type = 'ORA',DB='KEGG')
t2 <- Sys.time()
t2-t1 # 33 mins

res_Enrich <- list(res_ORA_GO = res_ORA_GO,
                   res_ORA_KEGG = res_ORA_KEGG)
saveRDS(res_Enrich,"./visualize_CCI/res_Enrich_TAM_TC_Rec.rds")

## clean

res_ORA_KEGG <- res_Enrich$res_ORA_KEGG
res_ORA_KEGG$ID %>% unique() %>% length()
df_kegg <- res_ORA_KEGG[res_ORA_KEGG$p.adjust <= 0.05,]
df_kegg$geneRatio <- lapply(df_kegg$GeneRatio,function(chr){
  # chr = df_kegg$GeneRatio[1]
  x = strsplit(chr,'/')[[1]][1] %>% as.numeric()
  y = strsplit(chr,'/')[[1]][2] %>% as.numeric()
  x/y
}) %>% unlist()
df_kegg$bgRatio <- lapply(df_kegg$BgRatio,function(chr){
  # chr = df_kegg$GeneRatio[1]
  x = strsplit(chr,'/')[[1]][1] %>% as.numeric()
  y = strsplit(chr,'/')[[1]][2] %>% as.numeric()
  x/y
}) %>% unlist()
df_kegg <- df_kegg[order(df_kegg$geneRatio,decreasing = T),]

## plot

df_kegg_ITGAV <- df_kegg[df_kegg$Regulator == 'ITGAV',]
df_kegg_ITGAV <- df_kegg_ITGAV[order(df_kegg_ITGAV$geneRatio,decreasing = T),]
df_kegg_ITGAV$Description <- factor(df_kegg_ITGAV$Description,levels = rev(df_kegg_ITGAV$Description))

pt_kegg_ITGAV <- ggplot(data=df_kegg_ITGAV[1:20,],aes(x=geneRatio,y=Description,size=Count,col=p.adjust)) +
  geom_point() + scale_color_material('pink',reverse = T, alpha = 0.5) + # scale_color_gradient(low = 'red', high = 'blue') + 
  theme_bw() + labs(title='KEGG Enrichment Analysis of ITGAV in TAM_TC',y='') +
  theme(
    plot.title = element_text(hjust = 0.5,size = 12), # 标题居中
    panel.background = element_blank(), # 去除坐标图的背景色
    legend.key = element_blank(), # 去除图例图案的背景色
    axis.text.y = element_text(size = 10),
    legend.position = 'bottom'
    # legend.direction = 'vertical' # horizontal
    # legend.box = 'vertical'
  )
pt_kegg_ITGAV

## TC-TAM ####

LRTG_im <- readRDS("./vaild_scRNAseq/getPIM/DEGs_from_bulk/job/LRTG_im_clean_Malignant_macrophages.rds")
LRTG_im_spl <- split(LRTG_im,LRTG_im$Receptor)

t1 <- Sys.time()
res_ORA_GO <- Perf_Enrich(LRTG_im_spl,Type = 'ORA',DB='GO')
res_ORA_KEGG <- Perf_Enrich(LRTG_im_spl,Type = 'ORA',DB='KEGG')
t2 <- Sys.time()
t2-t1 # 1.5 hours

res_Enrich <- list(res_ORA_GO = res_ORA_GO,
                   res_ORA_KEGG = res_ORA_KEGG)
saveRDS(res_Enrich,"./visualize_CCI/res_Enrich_TC_TAM_Rec.rds")

## clean

res_ORA_KEGG <- res_Enrich$res_ORA_KEGG
res_ORA_KEGG$ID %>% unique() %>% length()
df_kegg <- res_ORA_KEGG[res_ORA_KEGG$p.adjust<=0.05,]
df_kegg$geneRatio <- lapply(df_kegg$GeneRatio,function(chr){
  # chr = df_kegg$GeneRatio[1]
  x = strsplit(chr,'/')[[1]][1] %>% as.numeric()
  y = strsplit(chr,'/')[[1]][2] %>% as.numeric()
  x/y
}) %>% unlist()
df_kegg$bgRatio <- lapply(df_kegg$BgRatio,function(chr){
  # chr = df_kegg$GeneRatio[1]
  x = strsplit(chr,'/')[[1]][1] %>% as.numeric()
  y = strsplit(chr,'/')[[1]][2] %>% as.numeric()
  x/y
}) %>% unlist()
df_kegg <- df_kegg[order(df_kegg$geneRatio),]

## plot

# KEGG_IL4R_top

df_kegg_IL4R <- df_kegg[df_kegg$Regulator == 'IL4R',]
df_kegg_IL4R <- df_kegg_IL4R[order(df_kegg_IL4R$geneRatio,decreasing = T),]
df_kegg_IL4R$Description <- factor(df_kegg_IL4R$Description,levels = rev(df_kegg_IL4R$Description))

pt_kegg_IL4R <- ggplot(data=df_kegg_IL4R[1:20,],aes(x=geneRatio,y=Description,size=Count,col=p.adjust)) +
  geom_point() + scale_color_material('pink',reverse = T, alpha = 0.5) + # scale_color_gradient(low = 'red', high = 'blue') + 
  theme_bw() + labs(title='KEGG Enrichment Analysis of IL4R in TC-TAM',y='') +
  theme(
    plot.title = element_text(hjust = 0.5,size = 12), # 标题居中
    panel.background = element_blank(), # 去除坐标图的背景色
    legend.key = element_blank(), # 去除图例图案的背景色
    axis.text.y = element_text(size = 10),
    legend.position="bottom"
    # legend.direction = 'vertical' # horizontal
    # legend.box = 'vertical'
  )
pt_kegg_IL4R

# KEGG_CSF1R_top

df_kegg_CSF1R <- df_kegg[df_kegg$Regulator == 'CSF1R',]
df_kegg_CSF1R <- df_kegg_CSF1R[order(df_kegg_CSF1R$geneRatio,decreasing = T),]
df_kegg_CSF1R$Description <- factor(df_kegg_CSF1R$Description,levels = rev(df_kegg_IL4R$Description))
pt_kegg_CSF1R <- ggplot(data=df_kegg_CSF1R[1:20,],aes(x=geneRatio,y=Description,size=Count,col=p.adjust)) +
  geom_point() + scale_color_material('pink',reverse = T, alpha = 0.5) + # scale_color_gradient(low = 'red', high = 'blue') + 
  theme_bw() + labs(title='KEGG Enrichment Analysis of CSF1R in TC-TAM',y='') +
  theme(
    plot.title = element_text(hjust = 0.5,size = 12), # 标题居中
    panel.background = element_blank(), # 去除坐标图的背景色
    legend.key = element_blank(), # 去除图例图案的背景色
    axis.text.y = element_text(size = 10),
    legend.position = 'bottom',
    # legend.direction = 'vertical' # horizontal
    # legend.box = 'vertical'
  )
pt_kegg_CSF1R

# merge

pt_merge <- ggpubr::ggarrange(pt_kegg_IL4R,pt_kegg_CSF1R,pt_kegg_ITGAV, nrow = 1,align = 'hv')

png("./vaild_scRNAseq/visualize_CCI/enrichmentPlot/bubble_enrichment_KEGG_merge_top20.png", 
    width = 18, height = 5, units = 'in', res = 1000)
pt_merge
dev.off()

pdf("./vaild_scRNAseq/visualize_CCI/enrichmentPlot/bubble_enrichment_KEGG_merge_top20.pdf", 
    width = 18, height = 5)
pt_merge
dev.off()

###################
## GSEA analysis ##
###################
