#############
## library ##
#############

library(dplyr)
library(ROCR)
library(ggplot2)
library(ggsci)

rm(list = ls())
gc()

setwd("~/stMLnet/apply_in_stGBM/UKF_304_T")
source('../code/code.R')

###########
## color ##
###########

scales::show_col(pal_aaas(palette = "default", alpha = 0.6)(7))
mycolors_locus <- pal_aaas(palette = "default", alpha = 0.6)(7)

mycolor_software <- mycolors_locus[c(1:4)]
names(mycolor_software) <- c("NicheNet","stMLnet","CytoTalk","MISTy")
scales::show_col(mycolor_software)

##########
## load ##
##########

## groundtrue
VehEP_TAM <- readRDS('~/cell_cell_interaction/apply_in_stGBM/2016Science/TGs_list_limma-cpm_logfc1_pval0.1.rds') # 618 degs
#VehEP_TAM <- readRDS('~/cell_cell_interaction/apply_in_stGBM/2016Science/TGs_list_wilcox-cpm_logfc1.7_pval0.1.rds') # 1437 degs
degs <- VehEP_TAM$degs

## stMLnet
#stml_file <- './getPIM/LRTG_pim_clean_TME_TAM.rds'
stml_file <- './results_TG_limma_logfc1_padj0.1/getPIM/LRTG_pim_clean_TME_TAM.rds' # 121 TGs
stmlnet_score <- readRDS(stml_file)
stmlnet_score <- stmlnet_score %>% filter(regulator == 'CSF1R' & type == 'Receptor') %>% na.omit()

## NicheNet
nich_file <- '../OtherMethods/NicheNet/result/macrophage_limma_LRpair_weight.rds'
nichenet_score <- readRDS(nich_file)
nichenet_score <- nichenet_score[which(nichenet_score$ligand=='CSF1'),] %>% na.omit()

## CytoTalk
cyto_file <- '../OtherMethods/CytoTalk/result/macrophages_LRpair_distance.rds'
cytotalk_score <- readRDS(cyto_file)
cytotalk_score$from <- as.vector(cytotalk_score$from)
cytotalk_score <- cytotalk_score[cytotalk_score$from == 'Csf1r' & cytotalk_score$from_type == 'Receptor',]

## MISTy
misty_score <- readRDS("../OtherMethods/MISTy/results_macrophages_paraview_10/misty_score.rds")
misty_score <- as.data.frame(misty_score)
head(misty_score)
misty_score <- misty_score %>% filter(regulator == 'CSF1R' & type == 'Receptor') %>% na.omit()


##############
## evaluate ##
##############

## stMLnet
gene_stml <- intersect(stmlnet_score$Target,rownames(VehEP_TAM$DEGs)) # select "CSF1R" downstream TGs
pred_stml <- stmlnet_score$pIM[match(gene_stml,stmlnet_score$Target)]
label_stml <- gene_stml %in% VehEP_TAM$degs
res_stml <- get_evaluate_metrics(pred_stml,label_stml)

## NicheNet
gene_nich <- intersect(nichenet_score$target,rownames(VehEP_TAM$DEGs))
pred_nich <- nichenet_score$weight[match(gene_nich,nichenet_score$target)]
label_nich <- gene_nich %in% VehEP_TAM$degs
res_nich <- get_evaluate_metrics(pred_nich,label_nich)

## CytoTalk
gene_cyto <- intersect(cytotalk_score$to,rownames(VehEP_TAM$DEGs))
pred_cyto <- cytotalk_score$distance[match(gene_cyto,cytotalk_score$to)]
label_cyto <- gene_cyto %in% VehEP_TAM$degs
res_cyto <- get_evaluate_metrics(pred_cyto,label_cyto)

## MISTy
gene_misty <- intersect(misty_score$target,rownames(VehEP_TAM$DEGs))
pred_misty <- misty_score$score[match(gene_misty,misty_score$target)]
label_misty <- gene_misty %in% VehEP_TAM$degs
res_misty <- get_evaluate_metrics(pred_misty,label_misty)

## merge
result <- list(stMLnet = res_stml, NicheNet = res_nich, CytoTalk = res_cyto, MISTy = res_misty)

saveRDS(result, file="./compare_results/compare_resultst_limma_logfc1_results.rds")

## plot ROC

df_stml <- get_curve_input(res_metrics = result$stMLnet, curve_type = 'ROC')
df_nich <- get_curve_input(res_metrics = result$NicheNet, curve_type = 'ROC')
df_cyto <- get_curve_input(res_metrics = result$CytoTalk, curve_type = 'ROC')
df_mist <- get_curve_input(res_metrics = result$MISTy, curve_type = 'ROC')

df_ROC <- do.call('rbind',list(df_stml,df_nich,df_cyto,df_mist)) %>% as.data.frame()
df_ROC$method <- c(rep('stMLnet',nrow(df_stml)),rep('NicheNet',nrow(df_nich)),rep('CytoTalk',nrow(df_cyto)),rep('MISTy',nrow(df_mist)))

df_pt <- df_ROC
labs <- split(df_pt,df_pt$method) %>% lapply(.,function(df){
  
  df$AUC <- paste0('AUC of ',df$method,' = ',df$AUC,'\n')
  unique(df$AUC)
  
}) %>% unlist()
df_anno <- data.frame(x = 0.55,
                      y = c(0.1,0.2,0.3,0.4),
                      lab = labs, 
                      method = names(labs))

p2 <- ggplot(df_pt, aes(x = FPR, y = TPR, group = method, color = method)) + 
  geom_line(size=1.1) + scale_color_manual(values = mycolor_software) + theme_bw() +
  geom_segment(aes(x=0,y=0,xend=1,yend=1),color='grey',linetype='dashed',size=0.8) +
  geom_text(data = df_anno, aes(x,y,label = lab, color = method),vjust = 1,hjust = 0, size = 4) +
  theme(legend.background = element_rect(fill = NULL, colour = 'black'),
        legend.position = 'none',
        axis.title = element_text(size = 14), #size=14
        axis.text = element_text(size = 12), #size=12
        panel.grid = element_blank()) 
p2

png("./compare_results/limma_lofc1_barplot_all_ROC.png",
    width = 4, height = 5, units = 'in', res = 600)
print(p2)
dev.off()

pdf("./compare_results/limma_lofc1_barplot_all_ROC_4software.pdf",
    width = 4.8, height = 4)
print(p2)
dev.off()

## plot PRC

df_stml <- get_curve_input(res_metrics = result$stMLnet, curve_type = 'PRC')
df_nich <- get_curve_input(res_metrics = result$NicheNet, curve_type = 'PRC')
df_cyto <- get_curve_input(res_metrics = result$CytoTalk, curve_type = 'PRC')
df_mist <- get_curve_input(res_metrics = result$MISTy, curve_type = 'PRC')

df_PRC <- do.call('rbind',list(df_stml,df_nich,df_cyto,df_mist)) %>% as.data.frame()
df_PRC$method <- c(rep('stMLnet',nrow(df_stml)),rep('NicheNet',nrow(df_nich)),rep('CytoTalk',nrow(df_cyto)),rep('CytoTalk',nrow(df_mist)))

df_pt <- df_PRC
labs <- split(df_pt,df_pt$method) %>% lapply(.,function(df){
  
  df$AUC <- paste0('AUC of ',df$method,' = ',df$AUC,'\n')
  unique(df$AUC)
  
}) %>% unlist()
df_anno <- data.frame(x = 0.55,
                      y = c(0.1,0.2,0.3,0.4),
                      lab = labs, 
                      method = names(labs))

p3 <- ggplot(df_pt, aes(x = Recall, y = Precision, group = method, color = method)) + 
  geom_line(size=1.1) + scale_color_manual(values = mycolor_software) + theme_bw() +
  geom_segment(aes(x=0,y=0,xend=1,yend=1),color='grey',linetype='dashed',size=0.8) +
  geom_text(data = df_anno, aes(x,y,label = lab, color = method),vjust = 1,hjust = 0, size = 4) +
  theme(legend.background = element_rect(fill = NULL, colour = 'black'),
        legend.position = 'none',
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) 
p3

png("./compare_results/limma_lofc1_barplot_all_PRC.png",
    width = 4, height = 5, units = 'in', res = 600)
print(p3)
dev.off()

pdf("./compare_results/limma_lofc1_barplot_all_PRC_4software.pdf",
    width = 4.8, height = 4)
print(p3)
dev.off()

