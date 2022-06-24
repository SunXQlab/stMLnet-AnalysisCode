#############
## library ##
#############

library(dplyr)
library(ggsci)
library(ggplot2)
library(igraph)
library(plotrix)
library(ggraph)
library(org.Hs.eg.db)
library(clusterProfiler)     
library(ggalluvial)

rm(list=ls())
gc()
setwd("E:/stMLnet/apply_in_COVID19/")

source('../code/code.R')

###########
## color ##
###########

# celltype

st_covid_rctd <- readRDS("input/st_rctd.rds")
celltype <- st_covid_rctd$Cluster %>% unique()

scales::show_col(pal_igv(palette = "default", alpha = 0.8)(15))
mycolors_nejm <- pal_igv(palette = "default", alpha = 0.8)(15)

mycolor_ct <- mycolors_nejm[1:length(celltype)]
names(mycolor_ct) <- celltype
scales::show_col(mycolor_ct)

# nodekey

scales::show_col(pal_locuszoom(palette = "default", alpha = 0.8)(7))
mycolors_locus <- pal_locuszoom(palette = "default", alpha = 0.8)(7)

nodekey <- c("Ligand","Receptor","TF","Target")
mycolor_key <- mycolors_locus[1:4]
names(mycolor_key) <- nodekey
scales::show_col(mycolor_key)

# nodetype

scales::show_col(pal_locuszoom(palette = "default", alpha = 0.8)(7))
mycolors_locus <- pal_locuszoom(palette = "default", alpha = 0.8)(7)

nodetype <- c("cell","Sender","Receiver")
mycolor_nt <- mycolors_locus[1:3]
names(mycolor_nt) <- nodetype
scales::show_col(mycolor_nt)

################
## get detail ##
################
## get cell pairs #### 

inputdir <- './getPIM/'
files <- list.files(inputdir)[grep('_im_',list.files(inputdir))]
df_cellpair <- gsub('LRTG_im_clean_|.rds','',files) %>% strsplit(.,"-") %>% do.call('rbind',.) %>% as.data.frame()
colnames(df_cellpair) <- c('Sender','Receiver') 
df_cellpair$CP1 <- paste(df_cellpair$Sender,df_cellpair$Receiver,sep ='_')
df_cellpair$CP2 <- paste(df_cellpair$Receiver,df_cellpair$Sender,sep ='_')

celltypes <- unique(c(df_cellpair$Sender,df_cellpair$Receiver))
cp_of_inter <- data.frame(matrix(ncol = 4,dimnames = list(c(),c('ct1','ct2','keys_of_ct1','keys_of_ct2'))))
for (i in 1:length(celltypes)) {
  
  ct1 <- celltypes[i]
  for (j in (i+1):length(celltypes)) {
    
    ct2 <- celltypes[j]
    cp1 <- paste0(ct1,'_',ct2)
    cp2 <- paste0(ct2,'_',ct1)
    if(cp1 %in% list.files('./runscMLnet/') & cp2 %in% list.files('./runscMLnet/')){
      
      cat('check in ',cp1,' and ',cp2,'\n')
      mlnet1 <- readRDS(paste0('./runscMLnet/',cp1,'/scMLnet.rds'))
      mlnet2 <- readRDS(paste0('./runscMLnet/',cp2,'/scMLnet.rds'))
      
      ct2_tgs <- unique(mlnet1$TFTar$target)
      ct1_tgs <- unique(mlnet2$TFTar$target)
      
      ct1_ligs <- unique(mlnet1$LigRec$source)
      ct2_ligs <- unique(mlnet2$LigRec$source)
      
      ct1_keys <- intersect(ct1_ligs,ct1_tgs)
      ct2_keys <- intersect(ct2_ligs,ct2_tgs)
      
      if(length(ct2_keys)>0|length(ct1_keys)>0){
        cp_of_inter <- rbind(cp_of_inter,c(ct1,ct2,length(ct1_keys),length(ct2_keys)))
      }
      
    }
    
  }
  
}

cp_of_inter <- na.omit(cp_of_inter)
cp_of_inter <- cp_of_inter[cp_of_inter$keys_of_ct1 != 0 & cp_of_inter$keys_of_ct2 != 0,]
rownames(cp_of_inter) <- 1:nrow(cp_of_inter)

## get genes list ####

key_of_inter <- list()
for (k in 1:nrow(cp_of_inter)) {
  
  ct1 <- cp_of_inter$ct1[k]
  ct2 <- cp_of_inter$ct2[k]
  
  mlnet1 <- readRDS(paste0('./runscMLnet/',paste(ct1,ct2,sep = '_'),'/scMLnet.rds'))
  mlnet2 <- readRDS(paste0('./runscMLnet/',paste(ct2,ct1,sep = '_'),'/scMLnet.rds'))
  
  ct2_tgs <- unique(mlnet1$TFTar$target)
  ct1_tgs <- unique(mlnet2$TFTar$target)
  
  ct1_ligs <- unique(mlnet1$LigRec$source)
  ct2_ligs <- unique(mlnet2$LigRec$source)
  
  ct1_keys <- intersect(ct1_ligs,ct1_tgs)
  ct2_keys <- intersect(ct2_ligs,ct2_tgs)
  
  key_of_inter[[paste(ct1,ct2,sep = '_')]] <- list(
    ct1_keys = ct1_keys,
    ct2_keys = ct2_keys
  )
  
}

key_of_inter$Alveoli_Macrophage
key_of_inter$Alveoli_Monocyte
key_of_inter$Macrophage_Monocyte

## transform to table ####

fbloop <- lapply(1:length(key_of_inter), function(i){
  
  ls <- key_of_inter[[i]]
  cp <- names(key_of_inter)[i]
  ct1 <- strsplit(cp,'_')[[1]][1]
  ct2 <- strsplit(cp,'_')[[1]][2]
  ct1_keys <- paste(ls$ct1_keys,collapse = ' ')
  ct2_keys <- paste(ls$ct2_keys,collapse = ' ')
  rbind(c(ct1,ct2,ct2_keys),c(ct2,ct1,ct1_keys))
  
}) %>% do.call('rbind',.) %>% as.data.frame()
colnames(fbloop) <- c('Sender','Receiver','Siganl')

cts_of_interest <- c('Alveoli','Macrophage','Monocyte')
fbloop <- fbloop[fbloop$Sender %in% cts_of_interest & fbloop$Receiver %in% cts_of_interest,]

write.csv(fbloop,'./feedback/feddback_loop_detail.csv')

##############
## plot cor ##
##############
## calculate cor ####

pt_list <- list()
for (i in 1:nrow(fbloop)) {
  
  ## prepare
  
  sender <- fbloop$Sender[i]
  receiver <- fbloop$Receiver[i]
  tgs <- unlist(strsplit(fbloop$Siganl[i],' '))
  
  res_LRscore <- readRDS(paste0('./runModel/LRTG_allscore_',sender,'-',receiver,'.rds'))
  res_LRscore$LRs_score <- res_LRscore$LRs_score[tgs]
  res_LRscore$TGs_expr <- res_LRscore$TGs_expr[tgs]
  
  ## calculate
  
  mi_LRscore_tgs <- as.data.frame(matrix(ncol = 3,dimnames = list(c(),c("LRpair","TG","MI"))))
  pcc_LRscore_tgs <- as.data.frame(matrix(ncol = 4,dimnames = list(c(),c("LRpair","TG","R","pval"))))
  for (j in 1:length(tgs)) {
    
    print(j)
    res_LRscore_tg <- list(
      LRs_score = list(
        res_LRscore$LRs_score[[tgs[j]]]
      ),
      TGs_expr = list(
        res_LRscore$TGs_expr[[tgs[j]]]
      )
    )
    names(res_LRscore_tg$LRs_score) <- tgs[j]
    names(res_LRscore_tg$TGs_expr) <- tgs[j]
    
    mi_LRscore_tg <- getMI(res_LRscore_tg)
    pcc_LRscore_tg <- getPCC(res_LRscore_tg)
    write.csv(mi_LRscore_tg,paste0("./feedback/result/mi_",paste(sender,receiver,tgs[j],sep = '_'),".csv"))
    write.csv(pcc_LRscore_tg,paste0("./feedback/result/pcc_",paste(sender,receiver,tgs[j],sep = '_'),".csv"))
    
    mi_LRscore_tgs <- rbind(mi_LRscore_tgs,mi_LRscore_tg)
    pcc_LRscore_tgs <- rbind(pcc_LRscore_tgs,pcc_LRscore_tg)
    
  }
  mi_LRscore_tgs <- na.omit(mi_LRscore_tgs)
  pcc_LRscore_tgs <- na.omit(pcc_LRscore_tgs)
  
  ## plot-PCC
  
  df_plot <- pcc_LRscore_tgs
  
  pt <- ggplot(df_plot, aes(x = TG, y = R, color = TG)) + 
    geom_violin() + geom_boxplot(width=0.1) + 
    labs(title = paste(sender,receiver,sep = '-'))
  pt <- pt + theme_classic() + theme(
    legend.position = 'none',
    plot.title=element_text(hjust=0.5,size = 14),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = -45,hjust = 0),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )
  pt_list[[paste(sender,receiver,sep = '-')]] <- pt
  
}

## plot cor ####

p1 <- ggarrange(pt_list$`Macrophage-Alveoli`,
                pt_list$`Alveoli-Macrophage`,
                pt_list$`Monocyte-Macrophage`,
                nrow = 1,align = 'hv',widths = c(2,1,1))
p2 <- ggarrange(pt_list$`Monocyte-Alveoli`,
                pt_list$`Alveoli-Monocyte`,
                pt_list$`Macrophage-Monocyte`,
                nrow = 1,align = 'hv',widths = c(2,1,1))
p3 <- ggarrange(p1,p2,nrow = 2, align = 'hv')

pdf(paste0("./feedback/figure/merge_PCC.pdf"),
    width = 10,height = 6)
print(p3)
dev.off()

png(paste0("./feedback/figure/merge_PCC.png"),
    width = 10,height = 6, units = 'in', res = 300)
print(p3)
dev.off()

##################
## plot heatmap ##
##################

R.utils::setOption( "clusterProfiler.download.method",'auto' )

celltypes <- unique(fbloop$Sender)
for (receiver in celltypes) {
  
  ## prepare ####

  LRTG_pim <- lapply(celltypes[celltypes!=receiver], function(sender){

    LRTG_pim <- readRDS(paste0("./getPIM/LRTG_pim_clean_",sender,'-',receiver,".rds"))
    LRTG_pim <- LRTG_pim[LRTG_pim$type == 'Receptor',]
    LRTG_pim <- LRTG_pim[,1:2]

  }) %>% do.call('rbind',.) %>% .[!duplicated(.),]
  LRTG_pim_spl <- split(LRTG_pim,LRTG_pim$regulator)

  ## calculate ####

  res_ORA_GO <- Perf_Enrich(LRTG_pim_spl,Type = 'ORA',DB='GO')
  res_ORA_KEGG <- Perf_Enrich(LRTG_pim_spl,Type = 'ORA',DB='KEGG')
  res_Enrich <- list(res_ORA_GO = res_ORA_GO,
                     res_ORA_KEGG = res_ORA_KEGG)
  saveRDS(res_Enrich, paste0("./feedback/result/enrich_Sender-",receiver,".rds"))
  
  ## clean ####
  
  res_ORA_KEGG <- res_Enrich$res_ORA_KEGG
  res_ORA_KEGG$ID %>% unique() %>% length()
  df_kegg <- res_ORA_KEGG[res_ORA_KEGG$pvalue <= 0.01,]
  df_kegg$geneRatio <- lapply(df_kegg$GeneRatio,function(chr){
    x = strsplit(chr,'/')[[1]][1] %>% as.numeric()
    y = strsplit(chr,'/')[[1]][2] %>% as.numeric()
    x/y
  }) %>% unlist()
  df_kegg$bgRatio <- lapply(df_kegg$BgRatio,function(chr){
    x = strsplit(chr,'/')[[1]][1] %>% as.numeric()
    y = strsplit(chr,'/')[[1]][2] %>% as.numeric()
    x/y
  }) %>% unlist()
  keep_term <- lapply(unique(df_kegg$Regulator), function(key){
    
    df_kegg[df_kegg$Regulator==key,'ID'] %>% head(.,10)
    
  }) %>% unlist() %>% unique()
  df_kegg$ONTOLOGY <- 'KEGG'
  df_kegg <- df_kegg[df_kegg$ID %in% keep_term,c(13,1:2,5:7,10:11)]
  df_kegg <- df_kegg[order(df_kegg$ONTOLOGY,df_kegg$ID),]
  
  res_ORA_GO <- res_Enrich$res_ORA_GO
  res_ORA_GO$ID %>% unique() %>% length()
  df_go <- res_ORA_GO[res_ORA_GO$pvalue<=0.01,]
  df_go$geneRatio <- lapply(df_go$GeneRatio,function(chr){
    x = strsplit(chr,'/')[[1]][1] %>% as.numeric()
    y = strsplit(chr,'/')[[1]][2] %>% as.numeric()
    x/y
  }) %>% unlist()
  df_go$bgRatio <- lapply(df_go$BgRatio,function(chr){
    x = strsplit(chr,'/')[[1]][1] %>% as.numeric()
    y = strsplit(chr,'/')[[1]][2] %>% as.numeric()
    x/y
  }) %>% unlist()
  df_go <- df_go[df_go$ONTOLOGY=='BP',]
  keep_term <- lapply(unique(df_go$Regulator), function(key){
    
    df_go[df_go$Regulator==key,'ID'] %>% head(.,10)
    
  }) %>% unlist() %>% unique()
  df_go <- df_go[df_go$ID %in% keep_term,c(1:3,6:8,11:12),]
  df_go <- df_go[order(df_go$ONTOLOGY,df_go$ID),]
  
  df_plot <- rbind(df_go,df_kegg)
  df_plot <- df_plot[df_plot$ONTOLOGY %in% c('KEGG','BP'),]
  df_plot$Description <- factor(df_plot$Description, levels = unique(df_plot$Description))
  df_plot$ONTOLOGY <- factor(df_plot$ONTOLOGY, levels = c('BP','KEGG'))
  anno_x_loca <- lapply(unique(df_plot$ONTOLOGY), function(key){
    
    df <- df_plot[!duplicated(df_plot$ID),]
    grep(key,df$ONTOLOGY) %>% max()
    
  }) %>% unlist()
  anno_x_loca <- anno_x_loca+0.5
  names(anno_x_loca) <- unique(df_plot$ONTOLOGY)
  
  ## plot ####
  
  pt_merge <- ggplot(data=df_plot,aes(x=Description,y=Regulator,size=geneRatio,col=p.adjust))+
    geom_hline(aes(x=Description,y=Regulator,yintercept = 1:nrow(df_plot)),size= 1.5,colour= "#E4EDF2",alpha= .5)+
    geom_vline(aes(x=Description,y=Regulator,xintercept = anno_x_loca[1]),size=0.5,linetype= "dashed")+
    geom_point()+ coord_flip() +
    scale_color_material('pink',reverse = T, alpha = 0.8) +
    scale_size_continuous(range = c(1,4)) + 
    theme_bw() + labs(title=paste0('Function Enrichment Analysis (Sender to ',receiver,')'),x='',y='') +
    theme(
      plot.title = element_text(hjust = 0.5,size = 12),
      panel.background = element_blank(), 
      legend.key = element_blank(), 
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank(),
      legend.position = 'bottom',
      legend.direction = 'vertical'
    )
  pt_merge
  
  pdf(paste0("./feedback/figure/Enrich_Sender-",receiver,".pdf"),
      width = 8,height = 5)
  print(pt_merge)
  dev.off()
  
  png(paste0("./feedback/figure/Enrich_Sender-",receiver,".png"),
      width = 8,height = 5, units = 'in', res = 300)
  print(pt_merge)
  dev.off()
  
}

###################
## plot Alluvium ##
###################

inputdir <- './getPIM/'
files <- list.files(inputdir)[grep('_im_',list.files(inputdir))]
CPs <- gsub('LRTG_im_clean_|LRTG_pim_clean_|\\.rds','',files)
plotdir = './feedback/figure/'

celltypes <- unique(fbloop$Sender)
pt_list <- list()
for (receiver in celltypes) {
  
  for (sender in celltypes[celltypes!=receiver]) {
    
    files <- list.files(inputdir)
    files <- files[lapply(CPs, function(cp){grep(cp,files)}) %>% unlist() %>% unique()]
    files <- files[grep("_im_",files)]
    files <- files[grep(paste0(sender,"-",receiver),files)]
    
    LRTG_im_merge <- lapply(files, function(f){
      
      LRTG_im <- readRDS(paste0(inputdir,f))
      LRTG_im <- na.omit(LRTG_im)
      LRTG_im$Sender <- sender
      LRTG_im
      
    }) %>% do.call('rbind',.)
    
    if(!is.null(LRTG_im_merge)){
      
      df_MLnet_long_check <- prepareAlluviumPlotData_V2(lrtg_im = LRTG_im_merge, 
                                                        color.by = 'Sender',
                                                        do.check = TRUE)
      head(df_MLnet_long_check)
      
      colodb <- c(mycolor_nt,mycolor_key,mycolor_ct)
      gtitle <- paste0(sender,"-",receiver)
      pt <- drawAlluviumPlot(df_MLnet_long_check, colodb = colodb, gtitle = gtitle,
                             wd = plotdir,p_height=6, p_width=5)
      
      pt_list[[gtitle]] <- pt
      
    }
    
  }
  
}

names(pt_list)
pt_merge <- ggarrange(pt_list[[1]],pt_list[[2]],pt_list[[3]],
                      pt_list[[6]],pt_list[[4]],pt_list[[5]],
                      nrow = 2,ncol = 3,align = 'hv')
pt_merge

pdf(paste0(plotdir,"AlluviumPlot-merge-feedback.pdf"),
    width = 12,height = 12)
print(pt_merge)
dev.off()

png(paste0(plotdir,"AlluviumPlot-merge-feedback.png"),
    width = 12,height = 12, units = 'in', res = 300)
print(pt_merge)
dev.off()

################
## plot MLnet ##
################

df_key <- list(
  Alveoli = c("C3","GPC3"),
  Macrophage = c("FN1","TGM2"),
  Monocyte = c('NID1','VCAN')
)

celltypes <- unique(fbloop$Sender)
pt_list <- list()
for (receiver in celltypes) {
  
  ## prepare
  
  senders <- celltypes[celltypes!=receiver]
  
  MLnet_merge <- list(
    LigRec = as.data.frame(matrix(ncol = 2,dimnames = list(c(),c('source','target')))),
    RecTF =  as.data.frame(matrix(ncol = 2,dimnames = list(c(),c('source','target')))),
    TFTar =  as.data.frame(matrix(ncol = 2,dimnames = list(c(),c('source','target'))))
  )
  for (sender in senders) {
    
    MLnet <- readRDS(paste0('./runscMLnet/',sender,'_',receiver,'/scMLnet.rds'))
    MLnet_merge$LigRec <- rbind(MLnet_merge$LigRec,MLnet$LigRec) %>% as.data.frame() %>% na.omit()
    MLnet_merge$RecTF <- rbind(MLnet_merge$RecTF,MLnet$RecTF) %>% as.data.frame() %>% na.omit()
    MLnet_merge$TFTar <- rbind(MLnet_merge$TFTar,MLnet$TFTar) %>% as.data.frame() %>% na.omit()
    
  }
  
  LRTG_im_merge <- as.data.frame(matrix(ncol = 6,dimnames = list(c(),c("LRpair","Ligand","Receptor","Target","IM","im_norm"))))
  for (sender in senders) {
    
    LRTG_im <- readRDS(paste0('./getPIM//LRTG_im_clean_',sender,'-',receiver,'.rds'))
    LRTG_im_merge <- rbind(LRTG_im_merge,LRTG_im) %>% as.data.frame() %>% na.omit()
    
  }
  
  ## data
  
  Key <- df_key[senders] %>% unlist()
  Type <- 'Ligand'
  MLnet_key <- prepareMLnetworkPlotData_V3(mlnet=MLnet_merge,lrtg_im=LRTG_im_merge,Key=Key,Type=Type,do.check = T)
  str(MLnet_key)
  
  colodb = pal_locuszoom(palette = "default", alpha = 0.5)(4)
  names(colodb) <- nodekey
  scales::show_col(colodb)
  
  ## plot
  
  downstream <- 'Target'
  gtitle <- paste0('Sender-',receiver)
  plotdir <- "./feedback/figure/"
  drawMLnetworkPlot_V4(mlnet=MLnet_key,colodb=colodb,downstream = downstream,
                       gtitle=gtitle,wd=plotdir,p_height = 4.5,p_width = 7)
  
  
}
