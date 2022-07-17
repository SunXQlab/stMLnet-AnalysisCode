#############
## library ##
#############

library(dplyr)
library(ROCR)
library(ggplot2)
library(ggsci)

rm(list = ls())
gc()

setwd('./stMLnet/apply_in_stBC/')

source('../code/code.R')

###########
## color ##
###########

scales::show_col(pal_aaas(palette = "default", alpha = 0.6)(7))
mycolors_aaas <- pal_aaas(palette = "default", alpha = 0.6)(7)

mycolor_software <- mycolors_aaas[c(1:4)]
names(mycolor_software) <- c("NicheNet","stMLnet","CytoTalk","MISTy")
scales::show_col(mycolor_software)

###################
## load-software ##
###################

## NicheNet
nichenet_score <- readRDS("../other_method/NicheNet/stBC/result/tumor_LRpair_weight.rds")
nichenet_score$ligand <- as.vector(nichenet_score$ligand)
head(nichenet_score)

## CytoTalk
cytotalk_score <- readRDS("../other_method/CytoTalk/stBC/result/Malignant_LRpair_distance.rds")
cytotalk_score$from <- as.vector(cytotalk_score$from)
head(cytotalk_score)

## stMLnet
stmlnet_score <- readRDS("./getPIM/LRTG_pim_clean_TME-Malignant.rds")
stmlnet_score$pIM <- lapply(stmlnet_score$pIM,function(s){ifelse(s>=0,s,0)}) %>% unlist() %>% as.numeric()
head(stmlnet_score)

## MISTy
misty_score <- readRDS("../other_method/MISTy/stBC/results_Malignant_paraview_10/misty_score.rds")
misty_score <- as.data.frame(misty_score)
head(misty_score)

## groundtrue
degs_ls <-readRDS('../cellLine/BC_tumor_cellline/degs_ls.rds')

#######################
## evaluate-software ##
#######################

keep_datasets <- c("GSE157680_NRP1_KO","GSE15893_CXCR4_KO","GSE160990_TGFB1_Treatment",
                   "GSE54329_IL6_EFM-19_Treatment","GSE54329_IL6_HCC1428_Treatment",
                   "GSE54329_IL6_MDA-MB-134VI_Treatment","GSE54329_IL6_MDA-MB-175-VIIdsRed_Treatment")
keep_index <- match(keep_datasets,names(degs_ls))

# split

ress_ls <- lapply(keep_index, function(i){
  
  print(i)
  deg_ls <- degs_ls[[i]]
  key <- strsplit(names(degs_ls)[i],'_')[[1]][2]
  type <- ifelse(grepl('KO',names(degs_ls)[i]),'Receptor','Ligand')
  
  score_stml <- stmlnet_score %>% filter(regulator == key & type == type) 
  gene_stml <- intersect(score_stml$Target,names(deg_ls))
  label_stml <- deg_ls[gene_stml]
  pred_stml <- score_stml$pIM[match(gene_stml,score_stml$Target)]
  res_stml <- get_evaluate_metrics(pred_stml,label_stml)
  
  if(type == 'Ligand'){
    score_nich <- nichenet_score[which(nichenet_score$ligand==key),]
  }else if(type == 'Receptor' & key == 'AXL'){
    score_nich <- nichenet_score[which(nichenet_score$ligand=='GAS6'),]
  }else if(type == 'Receptor' & key == 'CXCR4'){
    score_nich <- nichenet_score[which(nichenet_score$ligand=='CXCL12'),]
  }else if(type == 'Receptor' & key == 'NRP1'){
    score_nich <- nichenet_score[which(nichenet_score$ligand=='VEGFA'),]
  }
  gene_nich <- intersect(score_nich$target,names(deg_ls))
  label_nich <- deg_ls[gene_nich]
  pred_nich <- score_nich$weight[match(gene_nich,score_nich$target)]
  res_nich <- get_evaluate_metrics(pred_nich,label_nich)
  
  score_cyto <- cytotalk_score %>% filter(from == key & from_type == type)
  gene_cyto <- intersect(score_cyto$to,names(deg_ls))
  label_cyto <- deg_ls[gene_cyto]
  pred_cyto <- score_cyto$distance[match(gene_cyto,score_cyto$to)]
  res_cyto <- get_evaluate_metrics(pred_cyto,label_cyto)
  
  score_mist <- misty_score %>% filter(regulator == regulator & type == type)
  gene_mist <- intersect(score_mist$target,names(deg_ls))
  label_mist <- deg_ls[gene_mist]
  pred_mist <- score_mist$score[match(gene_mist,score_mist$target)]
  res_mist <- get_evaluate_metrics(pred_mist,label_mist)
  
  result <- list(res_stml=res_stml, res_nich=res_nich, res_cyto=res_cyto, res_mist=res_mist)
  result
  
})
names(ress_ls) <- keep_datasets

# check

str(ress_ls,max.level = 1)

#####################
## output-software ##
#####################

## data 

evaluate_ls <- ress_ls
str(evaluate_ls,max.level = 1)
df_plot <- lapply(1:length(evaluate_ls), function(i){
  
  print(names(evaluate_ls)[i])
  res <- evaluate_ls[[i]]
  
  df <- do.call('rbind',list(res$res_stml$perf_metrics,
                             res$res_nich$perf_metrics,
                             res$res_cyto$perf_metrics,
                             res$res_mist$perf_metrics)) %>% as.data.frame()
  df$method <- c('stMLnet','NicheNet','CytoTalk','MISTy')
  df$group <- names(evaluate_ls)[i]
  df
  
}) %>% do.call('rbind',.)

## barplot

df_plot$Dataset <- strsplit(df_plot$group,'_') %>% do.call('rbind',.) %>% .[,2]
df_plot$Dataset[grep('IL6',df_plot$Dataset)] <- paste0(df_plot$Dataset[grep('IL6',df_plot$Dataset)],'_',rep(1:10,each=4))

rename_dataset <- function(x){
  switch(EXPR = x,
         'GSE54329_IL6_EFM-19_Treatment' = 'GSE54329(1)',
         'GSE54329_IL6_HCC1428_Treatment' = 'GSE54329(2)',
         'GSE54329_IL6_MDA-MB-134VI_Treatment' = 'GSE54329(3)',
         'GSE54329_IL6_MDA-MB-175-VIIdsRed_Treatment' = 'GSE54329(4)',
         'GSE160990_TGFB1_Treatment' = 'GSE160990',
         'GSE15893_CXCR4_KO' = 'GSE15893',
         'GSE157680_NRP1_KO' = 'GSE157680'
  )
}
df_plot$Dataset <- lapply(df_plot$group,rename_dataset) %>% unlist()
df_plot$ROC_AUC_label <- gsub('^0$','NA',df_plot$ROC_AUC)
df_plot$PRC_AUC_label <- gsub('^0$','NA',df_plot$PRC_AUC)

df_plot <- lapply(unique(df_plot$Dataset), function(ds){
  
  df <- df_plot[df_plot$Dataset == ds,]
  if(sum(df$ROC_AUC)==0){
    df = NULL
  }else{
    df = df
  }
  
})
delet_celllines <- which(lapply(df_plot, is.null) %>% unlist())
if(length(delet_celllines)>0) df_plot <- df_plot[-delet_celllines]
df_plot <- df_plot %>% do.call('rbind',.)
df_plot$method <-factor(df_plot$method,ordered=TRUE,levels=c('CytoTalk','NicheNet','MISTy','stMLnet'))

theme_bar <- function(..., bg='white'){
  require(grid)
  theme_classic(...) +
    theme(rect=element_rect(fill=bg),
          plot.margin=unit(rep(0.5,4), 'lines'),
          panel.background=element_rect(fill='transparent', color='black'),
          panel.border=element_rect(fill='transparent', color='transparent'),
          panel.grid=element_blank(),
          axis.title.x = element_blank(),
          axis.title.y=element_text(face = "bold",size = 14),
          axis.text = element_text(face = "bold",size = 14),
          legend.title=element_blank(),
          legend.position='top',
          legend.direction = "horizontal",
          legend.text = element_text(face = "bold",size = 12,margin = margin(r=20)),
          legend.background = element_rect( linetype="solid",colour ="black")
    )
  
}

p1 <- ggplot(df_plot, aes(x=Dataset,y=ROC_AUC,fill=method)) + 
  geom_bar(stat="identity",position=position_dodge(0.75),width=0.6) +
  scale_y_continuous(expand = c(0, 0.01)) + labs(y = 'AUCROC') +
  coord_cartesian(ylim=c(-0.02,1.1))+
  scale_fill_manual(values = mycolor_software, breaks=c('CytoTalk','NicheNet','MISTy','stMLnet')) + 
  theme_bar() + #coord_flip() +
  geom_text(aes(label=ROC_AUC_label),position=position_dodge(0.9),vjust = -0.5,size = 4) 
p1

p2 <- ggplot(df_plot, aes(x=Dataset,y=PRC_AUC,fill=method)) + 
  geom_bar(stat="identity",position=position_dodge(0.75),width=0.6) +
  scale_y_continuous(expand = c(0, 0.01)) + labs(y = 'AUCPRC') +
  coord_cartesian(ylim=c(-0.02,1.1))+
  scale_fill_manual(values = mycolor_software, breaks=c('CytoTalk','NicheNet','MISTy','stMLnet')) + 
  theme_bar() + #coord_flip() +
  geom_text(aes(label=PRC_AUC_label),position=position_dodge(0.9),vjust = -0.5,size = 4) 
p2

p11 <- p1 + theme(axis.text.x = element_blank())
p22 <- p2 + theme(legend.position = 'none')
gridExtra::grid.arrange(grobs = list(p11,p22), heights = c(1,0.9))

pdf(paste0("./output/barplot_aur_software.pdf"),
    width = 12,height = 7)
gridExtra::grid.arrange(grobs = list(p11,p22), heights = c(1,0.9))
dev.off()

png(paste0("./output/barplot_aur_software.png"),
    width = 12,height = 7, units = 'in', res = 300)
gridExtra::grid.arrange(grobs = list(p11,p22), heights = c(1,0.9))
dev.off()
