#############
## library ##
#############

library(dplyr)
library(ggplot2)
library(reshape2)
library(doParallel)

rm(list = ls())
gc()

setwd("./stMLnet/prior_knowledge/")

###########
## input ##
###########

## database ####

Databases <- readRDS("./output/Databases.rds")

pkLigRec <- Databases$LigRec.DB %>% tidyr::spread(target, score, fill = 0) %>% as.data.frame()
rownames(pkLigRec) <- pkLigRec$source
pkLigRec <- pkLigRec[,-1] %>% as.matrix()
dim(pkLigRec)

pkRecTF <- Databases$RecTF.DB %>% tidyr::spread(target, score, fill = 0) %>% as.data.frame()
rownames(pkRecTF) <- pkRecTF$source
pkRecTF <- pkRecTF[,-1] %>% as.matrix()
dim(pkRecTF)

pkTFTG <- Databases$TFTG.DB %>% tidyr::spread(target, score, fill = 0) %>% as.data.frame()
rownames(pkTFTG) <- pkTFTG$source
pkTFTG <- pkTFTG[,-1] %>% as.matrix()
dim(pkTFTG)

## CellLine ####

## Ligand

clLigTG <- data.table::fread("../cellLine/clean/result/hm_CL_LigTG.csv", header = T, sep = ",", data.table = F)
rownames(clLigTG) <- clLigTG[,1]
clLigTG <- clLigTG[,-1]
dim(clLigTG)

res.ligs <- readRDS(file = "../cellLine/clean/result/res.ligs.rds")
lengths(res.ligs)

## Receptor

clRecTG <- data.table::fread("../cellLine/clean/result/hm_CL_RecTG.csv", header = T, sep = ",", data.table = F)
rownames(clRecTG) <- clRecTG[,1]
clRecTG <- clRecTG[,-1]
dim(clRecTG)

res.recs <- readRDS(file = "../cellLine/clean/result/res.recs.rds")
lengths(res.recs)

## filter

names(res.ligs)[!names(res.ligs) %in% rownames(pkLigRec)]
rownames(pkLigRec)[grep('IFN',rownames(pkLigRec))]
names(res.ligs)[names(res.ligs)=='IFNB'] <- 'IFNB1'
rownames(clLigTG)[rownames(clLigTG)=='IFNB'] <- 'IFNB1'

names(res.recs)[!names(res.recs) %in% rownames(pkRecTF)]
rownames(pkRecTF)[grep('IGF',rownames(pkRecTF))]
names(res.recs)[names(res.recs)=='IGFR1'] <- 'IGF1R'
rownames(clRecTG)[rownames(clRecTG)=='IGFR1'] <- 'IGF1R'

## save ####

save.image("./quan.ct/tunePara_input.RData")

##########
## load ##
##########

rm(list = ls())
gc()
load("./quan.ct/tunePara_input.RData")

##############
## tunePara ##
##############

## function

tunePara <- function(quan.cutoff){
  
  ############
  ## cutoff ##
  ############
  
  quan.ct <- quan.cutoff
  pkiRecTF <- pkRecTF
  rwr.ct <- quantile(unlist(as.matrix(pkiRecTF)), quan.ct)
  pkiRecTF[pkiRecTF <= rwr.ct] <- 0
  
  ################
  ## 0-1 matrix ##
  ################
  
  pkiRecTF <- pkiRecTF[colnames(pkLigRec),rownames(pkTFTG)]
  
  pkLigTG <- pkLigRec %*% pkiRecTF %*% pkTFTG
  pkLigTG[pkLigTG != 0] <- 1
  
  pkRecTG <- pkiRecTF %*% pkTFTG
  pkRecTG[pkRecTG != 0] <- 1
  
  ##############
  ## keep dim ##
  ##############
  
  LigsInBoth <- intersect(rownames(pkLigTG), rownames(clLigTG))
  TGsInBoth <- intersect(colnames(pkLigTG), colnames(clLigTG))
  
  pkfLigTG <- pkLigTG[LigsInBoth,TGsInBoth]
  clLigTG <- clLigTG[LigsInBoth,TGsInBoth]
  
  RecsInBoth <- intersect(rownames(pkRecTG), rownames(clRecTG))
  TGsInBoth <- intersect(colnames(pkRecTG), colnames(clRecTG))
  
  pkfRecTG <- pkRecTG[RecsInBoth,TGsInBoth]
  clRecTG <- clRecTG[RecsInBoth,TGsInBoth]
  
  ######################
  ## objective matrix ##
  ######################
  
  LigTG <- as.matrix(pkfLigTG) * as.matrix(clLigTG)
  # table(unlist(LigTG))
  
  RecTG <- as.matrix(pkfRecTG) * as.matrix(clRecTG)
  # table(unlist(RecTG))
  
  ########################
  ## objective function ##
  ########################
  
  a=sum(pkiRecTF!=0)/(nrow(pkiRecTF)*ncol(pkiRecTF)) 
  
  b=sum(pkLigTG)/(nrow(pkLigTG)*ncol(pkLigTG)) 
  c=sum(pkRecTG)/(nrow(pkRecTG)*ncol(pkRecTG)) 
  
  d=sum(clLigTG)/(nrow(clLigTG)*ncol(clLigTG))
  e=sum(clRecTG)/(nrow(clRecTG)*ncol(clRecTG))
  
  f=sum(LigTG)/(nrow(LigTG)*ncol(LigTG))
  g=sum(RecTG)/(nrow(RecTG)*ncol(RecTG))
  
  h <- f/d
  i <- g/e
  
  J <- abs((h+i)-(b+c))
 
  ############
  ## return ##
  ############
  
  return(c(a,b,c,d,e,f,g,h,i,J))
  
}

## input

quan.cutoff1 <- c(seq(0.1,1,by = 0.05))
quan.cutoff2 <- c(seq(0.9,0.995,by = 0.005))

## test

t1 <- Sys.time()
message(paste0('Start at ',as.character(t1)))
cl <- makeCluster(6)
registerDoParallel(cl)
res1 <- foreach(i = quan.cutoff1,
                .export = c('pkLigRec','pkRecTF','pkTFTG','clLigTG','clRecTG'),
                .packages = c('dplyr'),
                .combine = 'rbind') %dopar% tunePara(i) %>% as.data.frame(.)
colnames(res1) <- c('a','b','c','d','e','f','g','h','i','J')
res1$quan.cutoff <- quan.cutoff1
stopCluster(cl)
t2 <- Sys.time()
message(paste0('End at ',as.character(t2)))
t2-t1 

t1 <- Sys.time()
message(paste0('Start at ',as.character(t1)))
cl <- makeCluster(6)
registerDoParallel(cl)
res2 <- foreach(i = quan.cutoff2,
               .export = c('pkLigRec','pkRecTF','pkTFTG','clLigTG','clRecTG'),
               .packages = c('dplyr'),
               .combine = 'rbind') %dopar% tunePara(i) %>% as.data.frame(.)
colnames(res2) <- c('a','b','c','d','e','f','g','h','i','J')
res2$quan.cutoff <- quan.cutoff2
stopCluster(cl)
t2 <- Sys.time()
message(paste0('End at ',as.character(t2)))
t2-t1 

saveRDS(res1, file = "./quan.ct/result/result_res1.rds")
saveRDS(res2, file = "./quan.ct/result/result_res1.rds")

## check

res1$quan.cutoff[which.max(res1$J1.2)]
res2$quan.cutoff[which.max(res2$J1.2)]
res <- rbind(res1,res2)

## plot

df.res <- res[res$quan.cutoff %in% seq(0,1,0.05), c('J1.2','quan.cutoff')]
p1 <-ggplot(df.res, aes(x=quan.cutoff, y=J1.2)) +
  geom_line() + geom_point() + labs(y='Score') + 
  geom_vline(aes(xintercept=0.95), colour="grey", linetype="dashed", size=1) + 
  scale_x_continuous(limits = c(0.1,1), breaks = seq(0,1,0.1)) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12)
  )
p1

df.res <- res[res$quan.cutoff %in% seq(0.9,0.999,by = 0.005), c('J1.2','quan.cutoff')]
p2 <-ggplot(df.res, aes(x=quan.cutoff, y=J1.2)) +
  geom_line() + geom_point() + labs(y='Score') + 
  geom_vline(aes(xintercept=0.98), colour="grey", linetype="dashed", size=1) + 
  scale_x_continuous(breaks = seq(0.9,1,0.01)) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12)
  )
p2

pdf('./quan.ct/result/plot_result.pdf', width = 10,height = 4.5)
gridExtra::grid.arrange(p1,p2,nrow=1)
dev.off()

png('./quan.ct/result/plot_result.png', width = 10,height = 4.5, units = 'in', res = 300)
gridExtra::grid.arrange(p1,p2,nrow=1)
dev.off()

