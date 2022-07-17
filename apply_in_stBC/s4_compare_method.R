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

mycolor_method <- mycolors_aaas[c(1:3)]
names(mycolor_method) <- c("exponent","reciprocal","constant")
scales::show_col(mycolor_method)

#################
## load-method ##
#################

## stMLnet-reci
stmlnet_score_reci <- readRDS("./getPIM/LRTG_pim_clean_TME-Malignant.rds")
stmlnet_score_reci$pIM <- lapply(stmlnet_score_reci$pIM,function(s){ifelse(s>=0,s,0)}) %>% unlist() %>% as.numeric()

## stMLnet-exp100
stmlnet_score_exp <- readRDS("./getPIM/LRTG_pim_clean_TME-Malignant-expo.rds")
stmlnet_score_exp$pIM <- lapply(stmlnet_score_exp$pIM,function(s){ifelse(s>=0,s,0)}) %>% unlist() %>% as.numeric()

## stMLnet-Constants
stmlnet_score_con <- readRDS("./getPIM/LRTG_pim_clean_TME-Malignant-constant.rds")
stmlnet_score_con$pIM <- lapply(stmlnet_score_con$pIM,function(s){ifelse(s>=0,s,0)}) %>% unlist() %>% as.numeric()

## groundtrue
degs_ls <-readRDS('../cellLine/BC_tumor_cellline/degs_ls.rds')

#####################
## evaluate-method ##
#####################

# keep
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
  
  score_reci <- stmlnet_score_reci %>% filter(regulator == key & type == type) 
  gene_reci <- intersect(score_reci$Target,names(deg_ls))
  label_reci <- deg_ls[gene_reci]
  pred_reci <- score_reci$pIM[match(gene_reci,score_reci$Target)]
  res_reci <- get_evaluate_metrics(pred_reci,label_reci)
  
  score_exp <- stmlnet_score_exp %>% filter(regulator == key & type == type) 
  gene_exp <- intersect(score_exp$Target,names(deg_ls))
  label_exp <- deg_ls[gene_exp]
  pred_exp <- score_exp$pIM[match(gene_exp,score_exp$Target)]
  res_exp <- get_evaluate_metrics(pred_exp,label_exp)
  
  score_con <- stmlnet_score_con %>% filter(regulator == key & type == type) 
  gene_con <- intersect(score_con$Target,names(deg_ls))
  label_con <- deg_ls[gene_con]
  pred_con <- score_con$pIM[match(gene_con,score_con$Target)]
  res_con <- get_evaluate_metrics(pred_con,label_con)
  
  result <- list(res_reci=res_reci, res_exp=res_exp, res_con=res_con)
  result
  
})
names(ress_ls) <- keep_datasets

# check
str(ress_ls,max.level = 1)

###################
## output-method ##
###################

## data 

evaluate_ls <- ress_ls
str(evaluate_ls,max.level = 1)
df_plot <- lapply(1:length(evaluate_ls), function(i){
  
  print(names(evaluate_ls)[i])
  res <- evaluate_ls[[i]]
  
  df <- do.call('rbind',list(res$res_reci$perf_metrics,
                             res$res_exp$perf_metrics,
                             res$res_con$perf_metrics)) %>% as.data.frame()
  df$method <- c('reciprocal','exponent','constant')
  df$group <- names(evaluate_ls)[i]
  df
  
}) %>% do.call('rbind',.)

## barplot

df_plot <- split(df_plot[,1:6],df_plot$method) %>% lapply(.,colMeans) %>% do.call('rbind',.) %>% as.data.frame()
df_plot$method <- rownames(df_plot)

df_plot_long <- reshape2::melt(df_plot,id.var=c('method'))
df_plot_long$label <- signif(df_plot_long$value,digits = 4)
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
df_plot_long$variable <- df_plot_long$variable %>% as.character() %>% lapply(.,rename_metrics) %>% unlist()
df_plot_long$variable <- factor(df_plot_long$variable,levels = c("AUCROC","AUCPR","Accuracy","PPV","Error","MCC"))

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

p1 <- ggplot(df_plot_long, aes(x=variable,y=value,fill=method)) + 
  geom_bar(stat="identity",position=position_dodge(0.75),width=0.6) +
  scale_y_continuous(expand = c(0, 0.01)) + labs(y = 'value') +
  coord_cartesian(ylim=c(-0.02,1.1))+
  scale_fill_manual(values = mycolor_method, breaks=c("constant","exponent","reciprocal")) + 
  theme_bar() + 
  geom_text(aes(label=label),position=position_dodge(0.9),vjust = -0.5,size = 4) 
p1

pdf("./output/barplot_mean_method.pdf",
    width = 12,height = 5)
p1
dev.off()

png("./output/barplot_mean_method.png",
    width = 12,height = 5, units = 'in', res = 300)
p1
dev.off()
