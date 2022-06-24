#############
## library ##
#############

library(dplyr)
library(ggplot2)
library(ranger)
library(caret)
library(doParallel)
library(parallel)
library(Metrics)
library(doSNOW)

rm(list=ls())
gc()
setwd("E:/stMLnet/apply_in_scGBM/")

source('../code/code.R')

##########
## Note ##
##########

## run in centos
## use 30~50 cores
## parallel in FORK
## about 1~3h/mission

## too much TGs of scRNAseq data
## split data into two part

##########
## main ##
##########

files <- list.files("./runModel/")
for(f in  files){
  
  message(f)
  message(paste0('Start at ',as.character(Sys.time())))
  
  label <- gsub("LRTGscore_TC_TAM_euc_|.rds","",f)
  LRTG_allscore <- readRDS(paste0("./runModel/DEGs_from_bulk/",f))
  num_of_LRs <- lapply(LRTG_allscore$LRs_score, ncol) %>% unlist()
  
  index1 <- names(num_of_LRs)[num_of_LRs<=150]
  index2 <- names(num_of_LRs)[num_of_LRs>150]
  
  # part1
  
  n.TG <- length(index1)
  parallel::detectCores()
  t1 <- Sys.time()
  cl <- makeSOCKcluster(50)
  registerDoSNOW(cl)
  pb <- txtProgressBar(min=1, max=n.TG, style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  res_ls <- foreach(i=1:n.TG, .packages=c("dplyr","ranger",'Metrics','caret'),
                    .options.snow=opts, .errorhandling = "pass"
  ) %dopar% {
    tg = index1[i]
    trainx = LRTG_allscore$LRs_score[[tg]]
    trainy = LRTG_allscore$TGs_expr[[tg]]
    get_pim_auto(trainx, trainy, ncores = 1, auto_para = TRUE, verbose = F)
  }
  names(res_ls) <- index1
  close(pb)
  stopCluster(cl)
  gc()
  t2 <- Sys.time()
  t2-t1 # 15 mins
  
  # 整理
  df_im = lapply(seq(n.TG), function(i){
    
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
  df_im = na.omit(df_im)
  saveRDS(df_im, paste0('./getPIM/LRTG_im_clean_',label,'_p1.rds'))
  
  df_pim = lapply(seq(n.TG), function(i){
    
    tg = index1[i]
    LRpairs = colnames(LRTG_allscore$LRs_score[[tg]])
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
  df_pim = na.omit(df_pim)
  saveRDS(df_pim, paste0('./getPIM/LRTG_pim_clean_',label,'_p1.rds'))
  
  # part2
  
  n.TG <- length(index2)
  parallel::detectCores()
  t1 <- Sys.time()
  cl <- makeSOCKcluster(30)
  registerDoSNOW(cl)
  pb <- txtProgressBar(min=1, max=n.TG, style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  res_ls <- foreach(i=1:n.TG, .packages=c("dplyr","ranger",'Metrics','caret'),
                    .options.snow=opts, .errorhandling = "pass"
  ) %dopar% {
    tg = index2[i]
    trainx = LRTG_allscore$LRs_score[[tg]]
    trainy = LRTG_allscore$TGs_expr[[tg]]
    get_pim_auto(trainx, trainy, ncores = 1, auto_para = TRUE, verbose = F)
    
  }
  names(res_ls) <- index2
  close(pb)
  stopCluster(cl)
  gc()
  t2 <- Sys.time()
  t2-t1
  
  # 整理
  df_im = lapply(seq(n.TG), function(i){
    
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
  df_im = na.omit(df_im)
  saveRDS(df_im, paste0('./getPIM/LRTG_im_clean_',label,'_p2.rds'))
  
  df_pim = lapply(seq(n.TG), function(i){
    
    tg = index2[i]
    LRpairs = colnames(LRTG_allscore$LRs_score[[tg]])
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
  df_pim = na.omit(df_pim)
  saveRDS(df_pim, paste0('./getPIM/LRTG_pim_clean_',label,'_p2.rds'))
  
  # 合并
  df_im1 <- readRDS(paste0('./getPIM/LRTG_im_clean_',label,'_p1.rds'))
  df_im2 <- readRDS(paste0('./getPIM/LRTG_im_clean_',label,'_p2.rds'))
  df_pim1 <- readRDS(paste0('./getPIM/LRTG_pim_clean_',label,'_p1.rds'))
  df_pim2 <- readRDS(paste0('./getPIM/LRTG_pim_clean_',label,'_p2.rds'))
  
  df_im <- rbind(df_im1,df_im2) %>% as.data.frame()
  saveRDS(df_im, paste0('./getPIM/LRTG_im_clean_',label,'.rds'))
  df_pim <- rbind(df_pim1,df_pim2) %>% as.data.frame()
  saveRDS(df_pim, paste0('./getPIM/LRTG_pim_clean_',label,'.rds'))
  
  all_objs <- ls()
  rm(list = all_objs[!all_objs %in% c(keep_objs,'keep_objs')])
  gc()
  
  message(paste0('End at ',as.character(Sys.time())))
  message('####################################################')
  rstudioapi::restartSession()
  
}

