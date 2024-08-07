#############
#  library  # 
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

setwd("/home/cjy/project/giotto_merfish_dataset/giotto_merfish_dataset_layer9/")

source('../code/code.R')

#########
## run ##
#########

wd = "./runModel/"
files <- list.files(wd)

unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
unregister_dopar()

dir.create(paste0("./getPIM/"),recursive = T)
keep_objs <- c(ls(), 'keep_objs')

for (f in files) {
  
  message(paste0('running job: ',f))
  label <- gsub('(LRTG_allscore_)|(.rds)','',f)
  LRTG_allscore <- readRDS(paste0(wd,f))
  LRTG_allscore$LRs_score <- lapply(LRTG_allscore$LRs_score,function(mat){as.matrix(mat)})
  LRTG_allscore$TGs_expr <- lapply(LRTG_allscore$TGs_expr,function(vec){unlist(vec)})
  
  numofLR <- lapply(LRTG_allscore$LRs_score, ncol) %>% unlist()
  numofLR <- numofLR[numofLR>1]
  
  LRTG_allscore$LRs_score <- LRTG_allscore$LRs_score[names(numofLR)]
  LRTG_allscore$TGs_expr <- LRTG_allscore$TGs_expr[names(numofLR)]
  
  n.TG <- length(LRTG_allscore$LRs_score)
  keeps_TG <- sum(numofLR>1)
  
  if(keeps_TG==1){
    
    trainx = LRTG_allscore$LRs_score[[1]]
    trainy = LRTG_allscore$TGs_expr[[1]]
    res_ls = get_pim_auto(trainx, trainy, ncores = 10, auto_para = TRUE, verbose = F)
    
    im = res_ls$df_IM
    if(is.null(im)|sum(im$IM,na.rm = T)==0){
      
      im = data.frame()
      
    }else{
      
      im$Ligand = stringr::str_split(im$LRpair,"_",simplify = T)[,1]
      im$Receptor = stringr::str_split(im$LRpair,"_",simplify = T)[,2]
      im$Target = names(numofLR)
      im$im_norm = im$IM/sum(im$IM)
      im = im[,c(2:5,1,6)]
      
    }
    df_im = im
    saveRDS(df_im, paste0('./getPIM/LRTG_im_clean_',label,'.rds'))
    
    LRpairs = colnames(LRTG_allscore$LRs_score[[1]])
    df_LR = data.frame(LRpairs)
    df_LR$Ligand = stringr::str_split(df_LR$LRpair,"_",simplify = T)[,1]
    df_LR$Receptor = stringr::str_split(df_LR$LRpair,"_",simplify = T)[,2]
    
    if(is.null(res_ls$df_pIM)|sum(res_ls$df_pIM$pIM,na.rm = T)==0){
      
      pim = data.frame()
      
    }else{
      
      pim = data.frame(regulator = gsub("shuffle_","",rownames(res_ls$df_pIM)))
      pim$Target = names(numofLR)
      pim$pIM = res_ls$df_pIM$pIM
      pim$type = pim$regulator %in% df_LR$Ligand
      pim$type[pim$type==TRUE] = 'Ligand'
      pim$type[pim$type==FALSE] = 'Receptor'
      
    }
    df_pim = pim
    saveRDS(df_pim, paste0('./getPIM/LRTG_pim_clean_',label,'.rds'))
    
  }else if(keeps_TG>1){
    
    parallel::detectCores()
    t1 <- Sys.time()
    cl <- makeSOCKcluster(30)
    registerDoSNOW(cl)
    pb <- txtProgressBar(min=1, max=n.TG, style=3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress=progress)
    res_ls <- foreach(i=1:n.TG, 
                      .packages=c("dplyr","ranger",'Metrics','caret'),
                      .options.snow=opts,
                      .errorhandling = "pass"
    ) %dopar% {
      trainx = LRTG_allscore$LRs_score[[i]]
      trainy = LRTG_allscore$TGs_expr[[i]]
      get_pim_auto(trainx, trainy, ncores = 1, auto_para = TRUE, verbose = F)
      
    }
    names(res_ls) <- names(LRTG_allscore$LRs_score)
    close(pb)
    stopCluster(cl)
    gc()
    t2 <- Sys.time()
    
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
  
  all_objs <- ls()
  rm(list = all_objs[!all_objs %in% c(keep_objs,'keep_objs')])
  gc()
  unregister_dopar()
  
}
