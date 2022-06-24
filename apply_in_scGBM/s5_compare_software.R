#############
## library ##
#############

library(dplyr)
library(ROCR)
library(ggplot2)
library(ggsci)

rm(list = ls())
gc()

setwd('E:/stMLnet/apply_in_scGBM/')

source('../code/code.R')

###########
## color ##
###########

scales::show_col(pal_aaas(palette = "default", alpha = 0.6)(7))
mycolors_locus <- pal_aaas(palette = "default", alpha = 0.6)(7)

mycolor_software <- mycolors_locus[c(1:3)]
names(mycolor_software) <- c("NicheNet","stMLnet","CytoTalk")
scales::show_col(mycolor_software)

##########
## load ##
##########

## groundtrue
degs_files <- "F:/finalVersion/vaild_scRNAseq/2016_science_bulk/output/DEGs_from_wilcox-cpm"
degs_file <- readRDS(paste(degs_files,'TGs_list_wilcox-cpm_logfc2_pval0.1.rds',sep = "/"))
degs <- degs_file$VehEP_TAM$degs

## stMLnet
stml_file <- './getPIM/LRTG_pim_clean_TME_TAM.rds'
stmlnet_score <- readRDS(stml_file)
stmlnet_score <- stmlnet_score %>% filter(regulator == 'CSF1R' & type == 'Receptor') %>% na.omit()

## NicheNet
nich_file <- 'E:/stMLnet/other_method/NicheNet/scGBM/result/macrophages_LRpair_weight.rds'
nichenet_score <- readRDS(nich_file)
nichenet_score <- nichenet_score[which(nichenet_score$ligand=='CSF1'),] %>% na.omit()

## CytoTalk
cyto_file <- 'E:/stMLnet/other_method/CytoTalk/scGBM/result/macrophages_LRpair_distance.rds'
cytotalk_score <- readRDS(cyto_file)
cytotalk_score$from <- as.vector(cytotalk_score$from)
cytotalk_score <- cytotalk_score[cytotalk_score$from == 'Csf1r' & cytotalk_score$from_type == 'Receptor',]

##############
## evaluate ##
##############

## stMLnet
gene_stml <- intersect(stmlnet_score$Target,rownames(degs_file$VehEP_TAM$DEGs))
pred_stml <- stmlnet_score$pIM[match(gene_stml,stmlnet_score$Target)]
label_stml <- gene_stml %in% degs_file$VehEP_TAM$degs
res_stml <- get_evaluate_metrics(pred_stml,label_stml)

## NicheNet
gene_nich <- intersect(nichenet_score$target,rownames(degs_file$VehEP_TAM$DEGs))
pred_nich <- nichenet_score$weight[match(gene_nich,nichenet_score$target)]
label_nich <- gene_nich %in% degs_file$VehEP_TAM$degs
res_nich <- get_evaluate_metrics(pred_nich,label_nich)

## CytoTalk
gene_cyto <- intersect(cytotalk_score$to,rownames(degs_file$VehEP_TAM$DEGs))
pred_cyto <- cytotalk_score$distance[match(gene_cyto,cytotalk_score$to)]
label_cyto <- gene_cyto %in% degs_file$VehEP_TAM$degs
res_cyto <- get_evaluate_metrics(pred_cyto,label_cyto)

## merge
result <- list(stMLnet = res_stml, NicheNet = res_nich, CytoTalk = res_cyto)

############
## output ##
############

## plot metrics

df_pt <- do.call('rbind',list(result$stMLnet$perf_metrics,
                              result$NicheNet$perf_metrics,
                              result$CytoTalk$perf_metrics)) %>% as.data.frame()
df_pt$method <- c('stMLnet','NicheNet','CytoTalk')

df_pt <- reshape2::melt(df_pt)
df_pt$value <- signif(df_pt$value,3)
rename_metrics <- function(x){
  switch(EXPR = x,
         'ROC_AUC' = 'AUCROC',
         'PRC_AUC' = 'AUCPR',       
         "ACC" = 'Accuracy',         
         "ERR" = 'Error',
         "PPV" = 'PPV',
         'MCC' = 'MCC'
  )
}
df_pt$variable <- df_pt$variable %>% as.character() %>% lapply(.,rename_metrics) %>% unlist()
df_pt$variable <- factor(df_pt$variable,levels = c("AUCROC","AUCPR","Accuracy","Error","PPV","MCC"))

p1 <- ggplot(df_pt, aes(x = variable, y = value, group = method, fill = method)) + 
  geom_bar(stat="identity",position=position_dodge(0.75),width=0.6) + 
  ylim(c(0,1)) + scale_fill_manual(values = mycolor_software) + theme_bw() +
  geom_text(aes(label=value),position=position_dodge(0.9),vjust = -0.5,size = 3) +
  theme(legend.background = element_rect(fill = NULL, colour = 'black'),
        legend.position = 'none',
        axis.title = element_blank(),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) 
p1

png("./output/barplot_all_metrics.png",
    width = 12, height = 4, units = 'in', res = 600)
print(p1)
dev.off()

pdf("./output/barplot_all_metrics.pdf",
    width = 12, height = 4)
print(p1)
dev.off()

## plot ROC

df_stml <- get_curve_input(res_metrics = result$stMLnet, curve_type = 'ROC')
df_nich <- get_curve_input(res_metrics = result$NicheNet, curve_type = 'ROC')
df_cyto <- get_curve_input(res_metrics = result$CytoTalk, curve_type = 'ROC')

df_ROC <- do.call('rbind',list(df_stml,df_nich,df_cyto)) %>% as.data.frame()
df_ROC$method <- c(rep('stMLnet',nrow(df_stml)),rep('NicheNet',nrow(df_nich)),rep('CytoTalk',nrow(df_cyto)))

df_pt <- df_ROC
labs <- split(df_pt,df_pt$method) %>% lapply(.,function(df){
  
  df$AUC <- paste0('AUC of ',df$method,' = ',df$AUC,'\n')
  unique(df$AUC)
  
}) %>% unlist()
df_anno <- data.frame(x = 0.55,
                      y = c(0.1,0.2,0.3),
                      lab = labs, 
                      method = names(labs))

p2 <- ggplot(df_pt, aes(x = FPR, y = TPR, group = method, color = method)) + 
  geom_line() + scale_color_manual(values = mycolor_software) + theme_bw() +
  geom_segment(aes(x=0,y=0,xend=1,yend=1),color='grey',linetype='dashed',size=0.5) +
  geom_text(data = df_anno, aes(x,y,label = lab, color = method),vjust = 1,hjust = 0, size = 4) +
  theme(legend.background = element_rect(fill = NULL, colour = 'black'),
        legend.position = 'none',
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) 
p2

png("./output/barplot_all_ROC.png",
    width = 4, height = 4, units = 'in', res = 600)
print(p2)
dev.off()

pdf("./output/barplot_all_ROC.pdf",
    width = 4, height = 4)
print(p2)
dev.off()

## plot PRC

df_stml <- get_curve_input(res_metrics = result$stMLnet, curve_type = 'PRC')
df_nich <- get_curve_input(res_metrics = result$NicheNet, curve_type = 'PRC')
df_cyto <- get_curve_input(res_metrics = result$CytoTalk, curve_type = 'PRC')

df_PRC <- do.call('rbind',list(df_stml,df_nich,df_cyto)) %>% as.data.frame()
df_PRC$method <- c(rep('stMLnet',nrow(df_stml)),rep('NicheNet',nrow(df_nich)),rep('CytoTalk',nrow(df_cyto)))

df_pt <- df_PRC
labs <- split(df_pt,df_pt$method) %>% lapply(.,function(df){
  
  df$AUC <- paste0('AUC of ',df$method,' = ',df$AUC,'\n')
  unique(df$AUC)
  
}) %>% unlist()
df_anno <- data.frame(x = 0.55,
                      y = c(0.1,0.2,0.3),
                      lab = labs, 
                      method = names(labs))

p3 <- ggplot(df_pt, aes(x = Recall, y = Precision, group = method, color = method)) + 
  geom_line() + scale_color_manual(values = mycolor_software) + theme_bw() +
  geom_segment(aes(x=0,y=0,xend=1,yend=1),color='grey',linetype='dashed',size=0.5) +
  geom_text(data = df_anno, aes(x,y,label = lab, color = method),vjust = 1,hjust = 0, size = 4) +
  theme(legend.background = element_rect(fill = NULL, colour = 'black'),
        legend.position = 'none',
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) 
p3

png("./output/barplot_all_PRC.png",
    width = 4, height = 4, units = 'in', res = 600)
print(p3)
dev.off()

pdf("./output/barplot_all_PRC.pdf",
    width = 4, height = 4)
print(p3)
dev.off()
