#############
## library ##
#############

library(dplyr)
library(Giotto)
library(Seurat)
library(SeuratWrappers)

rm(list=ls())
gc()

setwd("./stMLnet/apply_in_stBC/")

source('../code/code.R')

######################
## compare software ##
######################
## color ####

library(RColorBrewer)
cols <- brewer.pal(3, "Set1")
scales::show_col(cols)
names(cols) <- c("stMLnet","Random","CytoTalk")
cols[2:3] <- c('grey','DimGrey')

## load ####

gio_bc <- readRDS("./input/giotto_bc.rds")
annoMat <- data.frame(Barcode=gio_bc@cell_metadata$cell_ID %>% as.character(),
                      Cluster=gio_bc@cell_metadata$celltype %>% as.character())

locaMat <- data.frame(gio_bc@spatial_locs[,1:2])
rownames(locaMat) <- gio_bc@cell_metadata$cell_ID
distMat <- as.matrix(dist(locaMat))

seed <- 4321
exprMat = gio_bc@norm_expr
exprMat.Impute <- run_Imputation(exprMat, use.seed = T,seed = seed)

## main ####

t1 <- Sys.time()
celltype <- c('Bcell','Tcell','Stroma','Epithelial','Macrophage','Malignant')
for (Receiver in celltype) {
  
  print(Receiver)
  Sender = NULL
  
  ###########
  ## input ##
  ###########
  
  ## predicted LR~TG ####
  
  wd <- "./runscMLnet/"
  files = list.files(wd)
  files_tumor = files[grep(paste0("_",Receiver,"$"),files)]
  mulNetList = lapply(files_tumor, function(file_tumor){
    
    readRDS(paste0(wd,file_tumor,"/scMLnet.rds"))
    
  })
  names(mulNetList) = stringr::str_split(files_tumor,pattern = "_|\\.",simplify = T)[,1]
  mulNetList = mulNetList[!unlist(lapply(mulNetList, function(mulnet){nrow(mulnet$LigRec)==0}))]
  
  mulNet_tab = lapply(mulNetList, function(mlnet){
    
    ligrec = data.frame(Ligand = mlnet$LigRec$source, Receptor = mlnet$LigRec$target) 
    rectf = data.frame(Receptor = mlnet$RecTF$source, TF = mlnet$RecTF$target)
    tftg = data.frame(TF = mlnet$TFTar$source, Target = mlnet$TFTar$target)
    
    res = ligrec %>% merge(., rectf, by = 'Receptor') %>% 
      merge(., tftg, by = 'TF') %>% 
      dplyr::select(Ligand, Receptor, TF, Target) %>% 
      arrange(Ligand, Receptor)
    
  })
  mulNet_tab = do.call("rbind", mulNet_tab)
  
  LRpairs = by(mulNet_tab, as.character(mulNet_tab$Target), function(x){paste(x$Ligand, x$Receptor, sep = "_")})
  LRpairs = lapply(LRpairs, function(lrtg){lrtg[!duplicated(lrtg)]})
  TGs = names(LRpairs)
  res_pred <- list(LRpairs=LRpairs,TGs=TGs)
  
  ## cytotalk LR~TG ####
  
  cytotalk_wd <- '../other_method/CytoTalk/stBC/result/'
  cytotalk_score <- readRDS(paste0(cytotalk_wd,Receiver,"_LRpair_distance.rds"))
  
  folders <- list.files(cytotalk_wd)
  folders <- folders[grep(paste0("_",Receiver),folders)]
  cytotalk_crosstalk <- lapply(1:length(folders),function(i){
    
    wd <- paste0(cytotalk_wd,folders[i],"/IllustratePCSF/")
    files <- list.files(wd)
    
    corsstalkedge <- read.table(paste0(wd,files[1]), header = T, sep = "\t")
    corsstalkedge <- corsstalkedge[grep(Receiver,corsstalkedge$Receptor),]
    corsstalkedge$Ligand <- stringr::str_split(corsstalkedge$Ligand," ",simplify = T)[,1]
    corsstalkedge$Receptor <- stringr::str_split(corsstalkedge$Receptor," ",simplify = T)[,1]
    corsstalkedge[,1:2]
    
  }) %>% do.call('rbind',.) %>% .[!duplicated(.),]
  cytotalk_crosstalk$LRpair <- paste(cytotalk_crosstalk$Ligand,cytotalk_crosstalk$Receptor,sep="_")
  
  TGs <- cytotalk_score$to %>% unique() %>% as.character()
  LRpairs <- lapply(TGs, function(tg){
    
    Recs <-  cytotalk_score$from[cytotalk_score$to == tg & cytotalk_score$from_type == 'Receptor']
    LRpairs <- cytotalk_crosstalk$LRpair[cytotalk_crosstalk$Receptor %in% Recs]
    unique(LRpairs)
    
  })
  names(LRpairs) <- TGs
  res_cyto <- list(LRpairs=LRpairs,TGs=TGs) 
  
  ## simulated LR~TG ####
  
  Databases <- readRDS('../prior_knowledge/output/Databases.rds')
  LigRec.DB <- Databases$LigRec.DB
  LigRec.DB <- LigRec.DB[LigRec.DB$source %in% rownames(exprMat) &
                           LigRec.DB$target %in% rownames(exprMat),]
  LigRec <- paste(LigRec.DB$source,LigRec.DB$target,sep = "_")
  
  seed <- 1
  res_rand <- res_pred
  set.seed(seed)
  res_rand$TGs <- sample(rownames(exprMat), size=length(res_rand$TGs),replace = T)
  res_rand$LRpairs <- lapply(1:length(res_rand$LRpairs),function(i){
    
    set.seed(seed)
    sample(LigRec, size=length(res_rand$LRpairs[[i]]),replace = T)
    
  })
  names(res_rand$LRpairs) <- res_rand$TGs
  
  ##########
  ## main ##
  ##########
  
  ## get pred LRscore ####
  
  res_LRscore_pred <- calculate_LRTG_score_V2(exprMat = exprMat.Impute, distMat = distMat, annoMat = annoMat, 
                                              LRpairs = res_pred$LRpairs, TGs = res_pred$TGs, 
                                              Receiver = Receiver, Sender = Sender)
  saveRDS(res_LRscore_pred,paste0("./runCor/software/LRscore_pred_",Receiver,".rds"))
  
  mi_LRscore_pred <- getMI(res_LRscore_pred)
  pcc_LRscore_pred <- getPCC(res_LRscore_pred)
  write.csv(mi_LRscore_pred,paste0("./runCor/software/mi_LRscore_pred_",Receiver,".csv"))
  write.csv(pcc_LRscore_pred,paste0("./runCor/software/pcc_LRscore_pred_",Receiver,".csv"))
  
  ## get simu LRscore ####
  
  res_LRscore_simu <- calculate_LRTG_score_V2(exprMat = exprMat.Impute, distMat = distMat, annoMat = annoMat, 
                                              LRpairs = res_rand$LRpairs, TGs = res_rand$TGs, 
                                              Receiver = Receiver, Sender = Sender)
  saveRDS(res_LRscore_simu,paste0("./runCor/software/LRscore_simu_",Receiver,".rds"))
  
  mi_LRscore_simu <- getMI(res_LRscore_simu)
  pcc_LRscore_simu <- getPCC(res_LRscore_simu)
  write.csv(mi_LRscore_simu,paste0("./runCor/software/mi_LRscore_simu_",Receiver,".csv"))
  write.csv(pcc_LRscore_simu,paste0("./runCor/software/pcc_LRscore_simu_",Receiver,".csv"))
  
  ## get cyto LRscore ####
  
  res_LRscore_cyto <- calculate_LRTG_score_V2(exprMat = exprMat.Impute, distMat = distMat, annoMat = annoMat, 
                                              LRpairs = res_cyto$LRpairs, TGs = res_cyto$TGs, 
                                              Receiver = Receiver, Sender = Sender)
  saveRDS(res_LRscore_cyto,paste0("./runCor/software/LRscore_cyto_",Receiver,".rds"))
  
  mi_LRscore_cyto <- getMI(res_LRscore_cyto)
  pcc_LRscore_cyto <- getPCC(res_LRscore_cyto)
  write.csv(mi_LRscore_cyto,paste0("./runCor/software/mi_LRscore_cyto_",Receiver,".csv"))
  write.csv(pcc_LRscore_cyto,paste0("./runCor/software/pcc_LRscore_cyto_",Receiver,".csv"))
  
}
t2 <- Sys.time()
t2-t1

## plot ####

mi_list <- list()
pcc_list <- list()
pval_list <- list()
for (Receiver in celltype) {
  
  print(Receiver)
  Sender = NULL
  
  ##########
  ## main ##
  ##########
  
  ## get pred LRscore ####
  
  mi_LRscore_pred <- read.csv(paste0("./runCor/software/mi_LRscore_pred_",Receiver,".csv"),row.names = 1)
  pcc_LRscore_pred <- read.csv(paste0("./runCor/software/pcc_LRscore_pred_",Receiver,".csv"),row.names = 1)
  
  ## get simu LRscore ####
  
  mi_LRscore_simu <- read.csv(paste0("./runCor/software/mi_LRscore_simu_",Receiver,".csv"),row.names = 1)
  pcc_LRscore_simu <- read.csv(paste0("./runCor/software/pcc_LRscore_simu_",Receiver,".csv"),row.names = 1)
  
  ## get cyto LRscore ####
  
  mi_LRscore_cyto <- read.csv(paste0("./runCor/software/mi_LRscore_cyto_",Receiver,".csv"),row.names = 1)
  pcc_LRscore_cyto <- read.csv(paste0("./runCor/software/pcc_LRscore_cyto_",Receiver,".csv"),row.names = 1)
  
  ##########
  ## plot ##
  ##########
  
  ## plot-MI ####
  
  df_plot = rbind(mi_LRscore_pred,
                  mi_LRscore_cyto,
                  mi_LRscore_simu) %>% as.data.frame()
  df_plot$group = c(rep('stMLnet',nrow(mi_LRscore_pred)),
                    rep('CytoTalk',nrow(mi_LRscore_cyto)),
                    rep('Random',nrow(mi_LRscore_simu)))
  
  p <- ggplot(df_plot, aes(x = MI, group = group, color = group)) + 
    geom_histogram(aes(y = ..density..), fill = 'white', binwidth = 0.01) + 
    theme_bw() + labs(title = Receiver) + 
    scale_color_manual(values = cols) + 
    theme(legend.position = c(0.8,0.7),
          legend.background = element_rect(fill = NULL, colour = 'black'),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12)) +
    geom_density(alpha= 0.2, size = 1)
  
  p1 <- p + rotate()+ theme_classic() + theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )
  
  comp_list <- list(c('stMLnet','CytoTalk'),c('stMLnet','Random'),c('Random','CytoTalk'))
  p2 <- ggplot(df_plot, aes(x = group, y = MI, color = group)) + 
    scale_color_manual(values = cols) + geom_boxplot() + 
    stat_summary(fun.y=mean, geom="point", shape=18, size=3, color="black") +
    geom_signif(comparisons = comp_list, test = t.test, map_signif_level = T, color = 'black', step_increase = 0.1) 
  p2 <- p2 + theme_classic() + theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45,hjust = 1),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )
  
  p_merge <- ggarrange(p2,p1, widths = c(1.5,2.5), legend = 'bottom',
                       ncol = 2, nrow = 1,  align = "hv", 
                       common.legend = TRUE)
  mi_list[[Receiver]] <- p_merge
  
  ## plot-PCC ####
  
  df_plot = rbind(pcc_LRscore_pred,
                  pcc_LRscore_cyto,
                  pcc_LRscore_simu) %>% as.data.frame()
  df_plot$group = c(rep('stMLnet',nrow(pcc_LRscore_pred)),
                    rep('CytoTalk',nrow(pcc_LRscore_cyto)),
                    rep('Random',nrow(pcc_LRscore_simu)))
  df_plot$R = abs(df_plot$R)
  
  p <- ggplot(df_plot, aes(x = R, group = group, color = group)) + 
    geom_histogram(aes(y = ..density..), fill = 'white', binwidth = 0.01) + 
    # xlim(c(0,1)) + 
    theme_bw() + labs(title = Receiver) + 
    scale_color_manual(values = cols) + 
    theme(legend.position = c(0.8,0.7),
          legend.background = element_rect(fill = NULL, colour = 'black'),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12)) +
    geom_density(alpha= 0.2, size = 1)
  
  p1 <- p + rotate()+ theme_classic() + theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )
  
  comp_list <- list(c('stMLnet','CytoTalk'),c('stMLnet','Random'),c('Random','CytoTalk'))
  p2 <- ggplot(df_plot, aes(x = group, y = R, color = group)) + 
    scale_color_manual(values = cols) + geom_boxplot() + 
    # scale_y_continuous(limits = c(0,1.1), breaks = seq(0,1,0.25)) +
    stat_summary(fun.y=mean, geom="point", shape=18, size=3, color="black") +
    geom_signif(comparisons = comp_list, test = t.test, map_signif_level = T, color = 'black', step_increase = 0.1) 
  p2 <- p2 + theme_classic() + theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45,hjust = 1),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )
  
  p_merge <- ggarrange(p2,p1, widths = c(1.5,2.5), legend = 'bottom',
                       ncol = 2, nrow = 1,  align = "hv", 
                       common.legend = TRUE)
  pcc_list[[Receiver]] <- p_merge
  
  ## plot-Pval ####
  
  df_plot = rbind(pcc_LRscore_pred,
                  pcc_LRscore_cyto,
                  pcc_LRscore_simu) %>% as.data.frame()
  df_plot$group = c(rep('stMLnet',nrow(pcc_LRscore_pred)),
                    rep('CytoTalk',nrow(pcc_LRscore_cyto)),
                    rep('Random',nrow(pcc_LRscore_simu)))
  df_plot$pval = -log10(df_plot$pval)
  df_plot$pval[is.infinite(df_plot$pval)] = max(df_plot$pval[!is.infinite(df_plot$pval)])+1
  
  p <- ggplot(df_plot, aes(x = pval, group = group, color = group)) + 
    geom_histogram(aes(y = ..density..), fill = 'white', binwidth = 0.01) + 
    # xlim(c(0,1)) + 
    theme_bw() + labs(title = Receiver) +  
    scale_color_manual(values = cols) + 
    theme(legend.position = c(0.8,0.7),
          legend.background = element_rect(fill = NULL, colour = 'black'),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12)) +
    geom_density(alpha= 0.2, size = 1)
  
  p1 <- p + rotate()+ theme_classic() + theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )
  
  comp_list <- list(c('stMLnet','CytoTalk'),c('stMLnet','Random'),c('Random','CytoTalk'))
  p2 <- ggplot(df_plot, aes(x = group, y = pval, color = group)) + 
    scale_color_manual(values = cols) + geom_boxplot() + 
    stat_summary(fun.y=mean, geom="point", shape=18, size=3, color="black") +
    ggpubr::geom_signif(comparisons = comp_list, test = t.test, map_signif_level = T, color = 'black', step_increase = 0.1) 
  p2 <- p2 + theme_classic() + theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45,hjust = 1),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )
  
  p_merge <- ggarrange(p2,p1, widths = c(1.5,2.5), legend = 'bottom',
                       ncol = 2, nrow = 1,  align = "hv", 
                       common.legend = TRUE)
  pval_list[[Receiver]] <- p_merge
  
}

str(mi_list,max.level = 1)
p_merge <- ggpubr::ggarrange(
  mi_list[[1]],mi_list[[2]],mi_list[[3]],
  mi_list[[4]],mi_list[[5]],mi_list[[6]],
  ncol = 3, nrow = 2, align = 'hv')
pdf(paste0("./runCor/software/merge_MI.pdf"),
    width = 13,height = 8)
print(p_merge)
dev.off()
png(paste0("./runCor/software/merge_MI.png"),
    width = 13,height = 8, units = 'in', res = 300)
print(p_merge)
dev.off()

str(pcc_list,max.level = 1)
p_merge <- ggpubr::ggarrange(
  pcc_list[[1]],pcc_list[[2]],pcc_list[[3]],
  pcc_list[[4]],pcc_list[[5]],pcc_list[[6]],
  ncol = 3, nrow = 2, align = 'hv')
pdf(paste0("./runCor/software/-merge_PCC.pdf"),
    width = 13,height = 8)
print(p_merge)
dev.off()

png(paste0("./runCor/software/merge_PCC.png"),
    width = 13,height = 8, units = 'in', res = 300)
print(p_merge)
dev.off()

str(pval_list,max.level = 1)
p_merge <- ggpubr::ggarrange(
  pval_list[[1]],pval_list[[2]],pval_list[[3]],
  pval_list[[4]],pval_list[[5]],pval_list[[6]],
  ncol = 3, nrow = 2, align = 'hv')
pdf(paste0("./runCor/software/merge_Pval.pdf"),
    width = 13,height = 8)
print(p_merge)
dev.off()
png(paste0("./runCor/software/merge_Pval.png"),
    width = 13,height = 8, units = 'in', res = 300)
print(p_merge)
dev.off()

######################
## compare distance ##
######################
## color ####

library(RColorBrewer)
cols <- brewer.pal(3, "Set1")[-3]
scales::show_col(cols)
names(cols) <- c("close","far")

## load ####

gio_bc <- readRDS("./input/giotto_bc.rds")
annoMat <- data.frame(Barcode=gio_bc@cell_metadata$cell_ID %>% as.character(),
                      Cluster=gio_bc@cell_metadata$celltype %>% as.character())

locaMat <- data.frame(gio_bc@spatial_locs[,1:2])
rownames(locaMat) <- gio_bc@cell_metadata$cell_ID
distMat <- as.matrix(dist(locaMat))

seed <- 4321
exprMat = gio_bc@norm_expr
exprMat.Impute <- run_Imputation(exprMat, use.seed = T,seed = seed)

wd <- "./runscMLnet/"
files = list.files(wd)
files_tumor = files[grep(paste0("_",Receiver,"$"),files)]
mulNetList = lapply(files_tumor, function(file_tumor){
  
  readRDS(paste0(wd,file_tumor,"/scMLnet.rds"))
  
})
names(mulNetList) = stringr::str_split(files_tumor,pattern = "_|\\.",simplify = T)[,1]
mulNetList = mulNetList[!unlist(lapply(mulNetList, function(mulnet){nrow(mulnet$LigRec)==0}))]

mulNet_tab = lapply(mulNetList, function(mlnet){
  
  ligrec = data.frame(Ligand = mlnet$LigRec$source, Receptor = mlnet$LigRec$target) 
  rectf = data.frame(Receptor = mlnet$RecTF$source, TF = mlnet$RecTF$target)
  tftg = data.frame(TF = mlnet$TFTar$source, Target = mlnet$TFTar$target)
  
  res = ligrec %>% merge(., rectf, by = 'Receptor') %>% 
    merge(., tftg, by = 'TF') %>% 
    dplyr::select(Ligand, Receptor, TF, Target) %>% 
    arrange(Ligand, Receptor)
  
})
mulNet_tab = do.call("rbind", mulNet_tab)

LRpairs = by(mulNet_tab, as.character(mulNet_tab$Target), function(x){paste(x$Ligand, x$Receptor, sep = "_")})
LRpairs = lapply(LRpairs, function(lrtg){lrtg[!duplicated(lrtg)]})
TGs = names(LRpairs)

## parameter ####

Receiver = 'Malignant'
Sender = NULL

close.ct = 0.25
far.ct = 0.75

## main ####

## close group

res_LRscore_close <- calculate_LRTG_score_V2(exprMat = exprMat.Impute, distMat = distMat, annoMat = annoMat, 
                                             group = 'close', LRpairs = LRpairs, TGs = TGs, 
                                             Receiver = Receiver, Sender = Sender, close.ct = close.ct)
mi_LRscore_close <- getMI(res_LRscore_close)
pcc_LRscore_close <- getPCC(res_LRscore_close)
write.csv(mi_LRscore_close,paste0("./runCor/distance/mi_LRscore_close",close.ct,".csv"))
write.csv(pcc_LRscore_close,paste0("./runCor/distance/pcc_LRscore_close",close.ct,".csv"))

## far group

res_LRscore_far <- calculate_LRTG_score_V2(exprMat = exprMat.Impute, distMat = distMat, annoMat = annoMat, 
                                           group = 'far', LRpairs = LRpairs, TGs = TGs, 
                                           Receiver = Receiver, Sender = Sender, far.ct = far.ct)
mi_LRscore_far <- getMI(res_LRscore_far)
pcc_LRscore_far <- getPCC(res_LRscore_far)
write.csv(mi_LRscore_far,paste0("./runCor/distance/mi_LRscore_far",far.ct,".csv"))
write.csv(pcc_LRscore_far,paste0("./runCor/distance/pcc_LRscore_far",far.ct,".csv"))

## plot ####

## plot-MI

df_plot = rbind(mi_LRscore_close,
                mi_LRscore_far) %>% as.data.frame()
df_plot$group = c(rep('close',nrow(mi_LRscore_close)),
                  rep('far',nrow(mi_LRscore_far)))

p <- ggplot(df_plot, aes(x = MI, group = group, color = group)) + 
  geom_histogram(aes(y = ..density..), fill = 'white', bins = 50) + 
  # xlim(c(0,1)) + 
  theme_bw() + labs(title = Receiver) + 
  scale_color_manual(values = cols) + 
  theme(legend.position = c(0.8,0.7),
        legend.background = element_rect(fill = NULL, colour = 'black'),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  geom_density(alpha= 0.2, size = 1)

p1 <- p + rotate()+ theme_classic() + theme(
  axis.title.y = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.line.y = element_blank(),
  axis.title = element_text(size = 14),
  axis.text = element_text(size = 12)
)

comp_list <- list(c('close','far'))
p2 <- ggplot(df_plot, aes(x = group, y = MI, color = group)) + 
  scale_color_manual(values = cols) + geom_boxplot() + 
  stat_summary(fun.y=mean, geom="point", shape=18, size=3, color="black") +
  geom_signif(comparisons = comp_list, test = t.test, map_signif_level = T, color = 'black', step_increase = 0.1) 
p2 <- p2 + theme_classic() + theme(
  axis.title.x = element_blank(),
  axis.text.x = element_text(angle = 45,hjust = 1),
  axis.title = element_text(size = 14),
  axis.text = element_text(size = 12)
)

pt_mi <- ggarrange(p2,p1, widths = c(2,3), legend = 'none',
                   ncol = 2, nrow = 1,  align = "hv")
pt_mi

## plot-PCC

df_plot = rbind(pcc_LRscore_close,
                pcc_LRscore_far) %>% as.data.frame()
df_plot$group = c(rep('close',nrow(pcc_LRscore_close)),
                  rep('far',nrow(pcc_LRscore_far)))
df_plot$R = abs(df_plot$R)

p <- ggplot(df_plot, aes(x = R, group = group, color = group)) + 
  geom_histogram(aes(y = ..density..), fill = 'white', bins = 50) + 
  # xlim(c(0,1)) + 
  theme_bw() + labs(title = Receiver) + 
  scale_color_manual(values = cols) + 
  theme(legend.position = c(0.8,0.7),
        legend.background = element_rect(fill = NULL, colour = 'black'),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  geom_density(alpha= 0.2, size = 1)

p1 <- p + rotate()+ theme_classic() + theme(
  axis.title.y = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.line.y = element_blank(),
  axis.title = element_text(size = 14),
  axis.text = element_text(size = 12)
)

comp_list <- list(c('close','far'))
p2 <- ggplot(df_plot, aes(x = group, y = R, color = group)) + 
  scale_color_manual(values = cols) + geom_boxplot() + 
  stat_summary(fun.y=mean, geom="point", shape=18, size=3, color="black") +
  geom_signif(comparisons = comp_list, test = t.test, map_signif_level = T, color = 'black', step_increase = 0.1) 
p2 <- p2 + theme_classic() + theme(
  axis.title.x = element_blank(),
  axis.text.x = element_text(angle = 45,hjust = 1),
  axis.title = element_text(size = 14),
  axis.text = element_text(size = 12)
)

pt_pcc <- ggarrange(p2,p1, widths = c(2,3), legend = 'none',
                    ncol = 2, nrow = 1,  align = "hv")
pt_pcc

## plot-Pval

df_plot = rbind(pcc_LRscore_close,
                pcc_LRscore_far) %>% as.data.frame()
df_plot$group = c(rep('close',nrow(pcc_LRscore_close)),
                  rep('far',nrow(pcc_LRscore_far)))
df_plot$pval = -log10(df_plot$pval)
df_plot$pval[is.infinite(df_plot$pval)] = max(df_plot$pval[!is.infinite(df_plot$pval)])+1

p <- ggplot(df_plot, aes(x = pval, group = group, color = group)) + 
  geom_histogram(aes(y = ..density..), fill = 'white', bins = 50) + 
  # xlim(c(0,1)) + 
  theme_bw() + labs(title = Receiver) + 
  scale_color_manual(values = cols) + 
  theme(legend.position = c(0.8,0.7),
        legend.background = element_rect(fill = NULL, colour = 'black'),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  geom_density(alpha= 0.2, size = 1)

p1 <- p + rotate()+ theme_classic() + theme(
  axis.title.y = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.line.y = element_blank(),
  axis.title = element_text(size = 14),
  axis.text = element_text(size = 12)
)

comp_list <- list(c('close','far'))
p2 <- ggplot(df_plot, aes(x = group, y = pval, color = group)) + 
  scale_color_manual(values = cols) + geom_boxplot() + labs(y = '-log10P') +
  # scale_y_continuous(limits = c(0,1.1), breaks = seq(0,1,0.25)) +
  stat_summary(fun.y=mean, geom="point", shape=18, size=3, color="black") +
  geom_signif(comparisons = comp_list, test = t.test, map_signif_level = T, color = 'black', step_increase = 0.1) 
p2 <- p2 + theme_classic() + theme(
  axis.title.x = element_blank(),
  axis.text.x = element_text(angle = 45,hjust = 1),
  axis.title = element_text(size = 14),
  axis.text = element_text(size = 12)
)

pt_pval <- ggarrange(p2,p1, widths = c(2,3), legend = 'none',
                     ncol = 2, nrow = 1,  align = "hv")
pt_pval

## merge

pt_merge <- ggarrange(pt_pcc,pt_pval,pt_mi,legend = 'none',ncol = 3, nrow = 1,  align = "hv")
pt_merge

pdf(paste0("./runCor/distance/merge_all_close",close.ct,"_far",far.ct,".pdf"),
    width = 13,height = 4)
print(pt_merge)
dev.off()

png(paste0("./runCor/distance/merge_all_close",close.ct,"_far",far.ct,".png"),
    width = 13,height = 4, units = 'in', res = 300)
print(pt_merge)
dev.off()
