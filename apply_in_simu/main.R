#############
## library ##
#############

library(dplyr)
library(ranger)
library(caret)
library(doParallel)
library(parallel)
library(Metrics)
library(ROCR)
library(doSNOW)
library(ggsci)
library(ggplot2)

rm(list=ls())
gc()

setwd("./stMLnet/apply_in_simu/")

source('../code/code.R')

############################
## get LRscore and TGexpr ##
############################

## load ground true

ground <- openxlsx::read.xlsx("./data/ground_ture.xlsx", rowNames=T)
TGs = colnames(ground)
LRpairs <- lapply(1:ncol(ground), function(x){rownames(ground)})
names(LRpairs) <- TGs

## parameters

Receiver = 'CT_3'
Sender = NULL

## calculate

for(sheetID in 1:100){
  
  print(sheetID)
  
  input_ls <- prepare_input_data(sheetID = sheetID)
  exprMat <- input_ls$exprMat
  annoMat <- input_ls$annoMat
  locaMat <- input_ls$locaMat
  distMat = as.matrix(dist(locaMat))
  
  LRTG_score_reci <- calculate_LRTG_score_V2(exprMat = exprMat, distMat = distMat, annoMat = annoMat,
                                               LRpairs = LRpairs, TGs = TGs, Receiver = Receiver, Sender = Sender)
  saveRDS(LRTG_score_reci, paste0("./runModel/LRTG_score_reci_",sheetID,".rds"))
  
  LRTG_score_expo <- calculate_LRTG_score_V3(exprMat = exprMat, distMat = distMat, annoMat = annoMat,
                                             LRpairs = LRpairs, TGs = TGs, Receiver = Receiver, Sender = Sender)
  saveRDS(LRTG_score_expo, paste0("./runModel/LRTG_score_expo_",sheetID,".rds"))
  
  LRTG_score_mean <- calculate_LRTG_score_V4(exprMat = exprMat, distMat = distMat, annoMat = annoMat,
                                             LRpairs = LRpairs, TGs = TGs, Receiver = Receiver, Sender = Sender)
  saveRDS(LRTG_score_mean, paste0("./runModel/LRTG_score_mean_",sheetID,".rds")) 
  
}

####################
## train RF model ##
####################

wd <- "./runModel/"

files <-list.files(wd)
files <- files[grep('expo',files)]
for(f in files){
  
  print(f)
  label <- gsub('(LRTG_score_)|(.rds)','',f)
  LRTG_allscore <- readRDS(paste0(wd,f))
  n.TG <- length(LRTG_allscore$LRs_score)
  
  cl <- makeSOCKcluster(6)
  registerDoSNOW(cl)
  pb <- txtProgressBar(min=1, max=n.TG, style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  res_ls <- foreach(i=1:n.TG, .packages=c("dplyr","ranger",'caret'),
                    .options.snow=opts, .errorhandling = "pass"
  ) %dopar% {
    trainx = LRTG_allscore$LRs_score[[i]]
    trainy = LRTG_allscore$TGs_expr[[i]] %>% unlist()
    get_pim_auto(trainx, trainy, auto_para = TRUE, verbose = F)
  }
  names(res_ls) <- names(LRTG_allscore$LRs_score)
  close(pb)
  stopCluster(cl)
  gc() 
  
  df_im = lapply(seq(n.TG), function(i){
    
    print(paste0('gene',i))
    res = res_ls[[i]]
    im = res$df_IM
    
    if(is.null(im)|sum(im$IM,na.rm = T)==0){
      
      im = data.frame()
      
    }else{
      
      im$Ligand = stringr::str_split(im$LRpair,"_",simplify = T)[,1]
      im$Receptor = stringr::str_split(im$LRpair,"_",simplify = T)[,2]
      im$Target = names(res_ls)[i]
      im$im_norm = im$IM/sum(im$IM)
      
    }
    im
    
  })
  df_im = do.call("rbind",df_im)
  df_im = df_im[,c(2:5,1,6)]
  saveRDS(df_im, paste0('./getPIM/LRTG_im_clean_',label,'.rds'))
  
  df_pim = lapply(seq(n.TG), function(i){
    
    LRpairs = colnames(LRTG_allscore$LRs_score[[i]])
    df_LR = data.frame(LRpairs)
    df_LR$Ligand = stringr::str_split(df_LR$LRpair,"_",simplify = T)[,1]
    df_LR$Receptor = stringr::str_split(df_LR$LRpair,"_",simplify = T)[,2]
    
    res = res_ls[[i]]
    
    if(is.null(res$df_pIM)|sum(res$df_pIM$pIM,na.rm = T)==0){
      
      pim = data.frame()
      
    }else{
      
      pim = data.frame(regulator = gsub("shuffle_","",rownames(res$df_pIM)))
      pim$Target = names(res_ls)[i]
      pim$pIM = res$df_pIM$pIM
      pim$type = pim$regulator %in% df_LR$Ligand
      pim$type[pim$type==TRUE] = 'Ligand'
      pim$type[pim$type==FALSE] = 'Receptor'
      
    }
    pim
    
  })
  df_pim = do.call('rbind',df_pim)
  saveRDS(df_pim, paste0('./getPIM/LRTG_pim_clean_',label,'.rds'))
  
}

files <-list.files(wd)
files <- files[grep('mean',files)]
for(f in files){
  
  print(f)
  label <- gsub('(LRTG_score_)|(.rds)','',f)
  LRTG_allscore <- readRDS(paste0(wd,f))
  n.TG <- length(LRTG_allscore$LRs_score)
  
  cl <- makeSOCKcluster(6)
  registerDoSNOW(cl)
  pb <- txtProgressBar(min=1, max=n.TG, style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  res_ls <- foreach(i=1:n.TG, .packages=c("dplyr","ranger",'caret'),
                    .options.snow=opts, .errorhandling = "pass"
  ) %dopar% {
    trainx = LRTG_allscore$LRs_score[[i]]
    trainy = LRTG_allscore$TGs_expr[[i]] %>% unlist()
    get_pim_auto(trainx, trainy, auto_para = TRUE, verbose = F)
  }
  names(res_ls) <- names(LRTG_allscore$LRs_score)
  close(pb)
  stopCluster(cl)
  gc()
  
  df_im = lapply(seq(n.TG), function(i){
    
    print(paste0('gene',i))
    res = res_ls[[i]]
    im = res$df_IM
    
    if(is.null(im)|sum(im$IM,na.rm = T)==0){
      
      im = data.frame()
      
    }else{
      
      im$Ligand = stringr::str_split(im$LRpair,"_",simplify = T)[,1]
      im$Receptor = stringr::str_split(im$LRpair,"_",simplify = T)[,2]
      im$Target = names(res_ls)[i]
      im$im_norm = im$IM/sum(im$IM)
      
    }
    im
    
  })
  df_im = do.call("rbind",df_im)
  df_im = df_im[,c(2:5,1,6)]
  saveRDS(df_im, paste0('./getPIM/LRTG_im_clean_',label,'.rds'))
  
  df_pim = lapply(seq(n.TG), function(i){
    
    LRpairs = colnames(LRTG_allscore$LRs_score[[i]])
    df_LR = data.frame(LRpairs)
    df_LR$Ligand = stringr::str_split(df_LR$LRpair,"_",simplify = T)[,1]
    df_LR$Receptor = stringr::str_split(df_LR$LRpair,"_",simplify = T)[,2]
    
    res = res_ls[[i]]
    
    if(is.null(res$df_pIM)|sum(res$df_pIM$pIM,na.rm = T)==0){
      
      pim = data.frame()
      
    }else{
      
      pim = data.frame(regulator = gsub("shuffle_","",rownames(res$df_pIM)))
      pim$Target = names(res_ls)[i]
      pim$pIM = res$df_pIM$pIM
      pim$type = pim$regulator %in% df_LR$Ligand
      pim$type[pim$type==TRUE] = 'Ligand'
      pim$type[pim$type==FALSE] = 'Receptor'
      
    }
    pim
    
  })
  df_pim = do.call('rbind',df_pim)
  saveRDS(df_pim, paste0('./getPIM/LRTG_pim_clean_',label,'.rds'))
  
}

files <-list.files(wd)
files <- files[grep('reci',files)]
for(f in files){
  
  print(f)
  label <- gsub('(LRTG_score_)|(.rds)','',f)
  LRTG_allscore <- readRDS(paste0(wd,f))
  n.TG <- length(LRTG_allscore$LRs_score)
  
  cl <- makeSOCKcluster(6)
  registerDoSNOW(cl)
  pb <- txtProgressBar(min=1, max=n.TG, style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  res_ls <- foreach(i=1:n.TG, .packages=c("dplyr","ranger",'caret'),
                    .options.snow=opts, .errorhandling = "pass"
  ) %dopar% {
    trainx = LRTG_allscore$LRs_score[[i]]
    trainy = LRTG_allscore$TGs_expr[[i]] %>% unlist()
    get_pim_auto(trainx, trainy, ncores = 1, auto_para = TRUE, verbose = F)
  }
  names(res_ls) <- names(LRTG_allscore$LRs_score)
  close(pb)
  stopCluster(cl)
  gc() 
  
  df_im = lapply(seq(n.TG), function(i){
    
    print(paste0('gene',i))
    res = res_ls[[i]]
    im = res$df_IM
    
    if(is.null(im)|sum(im$IM,na.rm = T)==0){
      
      im = data.frame()
      
    }else{
      
      im$Ligand = stringr::str_split(im$LRpair,"_",simplify = T)[,1]
      im$Receptor = stringr::str_split(im$LRpair,"_",simplify = T)[,2]
      im$Target = names(res_ls)[i]
      im$im_norm = im$IM/sum(im$IM)
      
    }
    im
    
  })
  df_im = do.call("rbind",df_im)
  df_im = df_im[,c(2:5,1,6)]
  saveRDS(df_im, paste0('./getPIM/LRTG_im_clean_',label,'.rds'))
  
  df_pim = lapply(seq(n.TG), function(i){
    
    LRpairs = colnames(LRTG_allscore$LRs_score[[i]])
    df_LR = data.frame(LRpairs)
    df_LR$Ligand = stringr::str_split(df_LR$LRpair,"_",simplify = T)[,1]
    df_LR$Receptor = stringr::str_split(df_LR$LRpair,"_",simplify = T)[,2]
    
    res = res_ls[[i]]
    
    if(is.null(res$df_pIM)|sum(res$df_pIM$pIM,na.rm = T)==0){
      
      pim = data.frame()
      
    }else{
      
      pim = data.frame(regulator = gsub("shuffle_","",rownames(res$df_pIM)))
      pim$Target = names(res_ls)[i]
      pim$pIM = res$df_pIM$pIM
      pim$type = pim$regulator %in% df_LR$Ligand
      pim$type[pim$type==TRUE] = 'Ligand'
      pim$type[pim$type==FALSE] = 'Receptor'
      
    }
    pim
    
  })
  df_pim = do.call('rbind',df_pim)
  saveRDS(df_pim, paste0('./getPIM/LRTG_pim_clean_',label,'.rds'))
  
}

############################
## evaluate pred and true ##
############################

## ground true

ground <- openxlsx::read.xlsx("./data/ground_ture.xlsx", rowNames=T)
groundtrue <- reshape2::melt(ground)
colnames(groundtrue) <- c('TG','label')
groundtrue <- groundtrue %>% mutate(LRpair = rep(rownames(ground),ncol(ground))) %>% select(LRpair,TG,label)

## main

df_metrics <- as.data.frame(matrix(ncol = 8,
   dimnames = list(c(),c('ROC_AUC','PRC_AUC','ACC','ERR','PPV','MCC','sheet','method'))))
df_ROC <- as.data.frame(matrix(ncol = 5,
   dimnames = list(c(),c('FPR','TPR','AUC','method','sheet'))))
df_PRC <- as.data.frame(matrix(ncol = 5,
   dimnames = list(c(),c('Recall','Precision','AUC','method','sheet'))))

for(sheetID in 1:100){
  
  print(sheetID)
  
  ## reci
  reci_score <- readRDS(paste0("./getPIM/LRTG_im_clean_reci_",sheetID,".rds"))
  reci_score <- reci_score %>% select(LRpair, TG = Target, IM, im_norm)
  res_reci <- get_evaluate_metrics(pred = reci_score$im_norm, label = groundtrue$label)
  clean_reci <- clean_res_metrcs(res_metrics = res_reci, sheetID = sheetID, method = 'reci')
  
  ## expo
  expo_score <- readRDS(paste0("./getPIM/LRTG_im_clean_expo_",sheetID,".rds"))
  expo_score <- expo_score %>% select(LRpair, TG = Target, IM, im_norm)
  res_expo <- get_evaluate_metrics(pred = expo_score$im_norm, label = groundtrue$label)
  clean_expo <- clean_res_metrcs(res_metrics = res_expo, sheetID = sheetID, method = 'expo')
  
  ## mean
  mean_score <- readRDS(paste0("./getPIM/LRTG_im_clean_mean_",sheetID,".rds"))
  mean_score <- mean_score %>% select(LRpair, TG = Target, IM, im_norm)
  res_mean <- get_evaluate_metrics(pred = mean_score$im_norm, label = groundtrue$label)
  clean_mean <- clean_res_metrcs(res_metrics = res_mean, sheetID = sheetID, method = 'mean')
  
  ## merge
  df_metrics = do.call('rbind',list(df_metrics,
                                    clean_reci$df_metrics,
                                    clean_expo$df_metrics,
                                    clean_mean$df_metrics))
  df_ROC = do.call('rbind',list(df_ROC,clean_reci$df_ROC,clean_expo$df_ROC,clean_mean$df_ROC))
  df_PRC = do.call('rbind',list(df_PRC,clean_reci$df_PRC,clean_expo$df_PRC,clean_mean$df_PRC))
}

df_metrics <- df_metrics[-1,]
df_ROC <- df_ROC[-1,]
df_PRC <- df_PRC[-1,]

save(df_metrics,df_ROC,df_PRC, file = "./figure/res_sheet1-100_input.rda")
load("./figure/res_sheet1-100_input.rda")

##########
## plot ##
##########

## color

scales::show_col(pal_aaas(palette = "default", alpha = 0.6)(7))
mycolors_aaas <- pal_aaas(palette = "default", alpha = 0.6)(7)

mycolor_method <- mycolors_aaas[c(3,1,2)]
names(mycolor_method) <- c("constant","exponent","reciprocal")
scales::show_col(mycolor_method)

## plot

df_plot <- reshape2::melt(df_metrics, id.var = c('method','sheet'))
rename_methods<- function(x){
  switch(EXPR = x,
         'expo' = 'exponent',
         'reci' = 'reciprocal',
         "mean" = 'constant'
  )
}
df_plot$method <- df_plot$method %>% lapply(.,rename_methods) %>% unlist() 
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
df_plot$variable <- df_plot$variable %>% as.character() %>% lapply(.,rename_metrics) %>% unlist()
df_plot$variable <- factor(df_plot$variable,levels = c("AUCROC","AUCPR","Accuracy","PPV","Error","MCC"))

comp_list <- list(c("reciprocal","exponent"),
                  c('reciprocal','constant'),
                  c('exponent','constant'))

p5_list <- lapply(levels(df_plot$variable), function(var){
  
  pt <- ggplot(df_plot[df_plot$variable==var,], 
               aes(x = method, y = value, color = method)) + 
    scale_color_manual(values = mycolor_method) + 
    geom_violin() + geom_boxplot(width=0.1) + 
    labs(y = var,x='distance weight') + 
    scale_y_continuous(limits = c(0,1.2), breaks = seq(0,1,0.2)) + 
    scale_color_manual(values = mycolor_method, breaks=c("constant","exponent","reciprocal"))+
    scale_x_discrete(breaks=c("constant","exponent","reciprocal"),labels=c("constant","exponent","reciprocal"))+
    stat_summary(fun.y=mean, geom="point", shape=18, size=1.5, color="black") +
    ggsignif::geom_signif(comparisons = comp_list, test = 'wilcox.test', map_signif_level = F, 
                          color = 'black', step_increase = 0.15,) 
  pt <- pt + theme_classic() + theme(
    axis.title = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = 'none'
  )
  pt
  
})
str(p5_list,max.level = 1)

p5_list[[1]] <- p5_list[[1]] + theme(axis.text.x = element_blank(),axis.title.x = element_blank()) +
  scale_y_continuous(limits = c(0.4,1.2), breaks = seq(0,1,0.2))
p5_list[[2]] <- p5_list[[2]] + theme(axis.text.x = element_blank(),axis.title.x = element_blank()) +
  scale_y_continuous(limits = c(0.4,1.2), breaks = seq(0,1,0.2))
p5_list[[3]] <- p5_list[[3]] + theme(axis.text.x = element_blank(),axis.title.x = element_blank()) +
  scale_y_continuous(limits = c(0.5,1.2), breaks = seq(0,1,0.2))
p5_list[[4]] <- p5_list[[4]] + 
  scale_y_continuous(limits = c(0.5,1.2), breaks = seq(0,1,0.2))
p5_list[[5]] <- p5_list[[5]] + 
  scale_y_continuous(limits = c(0,0.5), breaks = seq(0,1,0.2))
p5_list[[6]] <- p5_list[[6]] + 
  scale_y_continuous(limits = c(0,1.2), breaks = seq(0,1,0.2))

p_merge <- ggpubr::ggarrange(p5_list[[1]], p5_list[[2]], p5_list[[3]], 
                             p5_list[[4]], p5_list[[5]], p5_list[[6]], 
                             legend = 'none',ncol = 3, nrow = 2,  align = "hv")
p_merge

pdf("./figure/violineplot_sheet1-100_merge.pdf",width = 12,height = 8)
p_merge
dev.off()

png("./figure/violineplot_sheet1-100_merge.png",width = 12,height = 8,units = 'in', res = 300)
p_merge
dev.off()
