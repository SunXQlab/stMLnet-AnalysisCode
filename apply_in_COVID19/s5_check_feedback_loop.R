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
setwd('/home/yll/cell_cell_interaction/apply_in_COVID19')

source('/home/yll/cell_cell_interaction/apply_in_COVID19/feedback_loop_function.R')

###########
## color ##
###########

# celltype
res_path <- '~/cell_cell_interaction/stMLnet_cjy/apply_in_COVID19/'
st_covid_rctd <- readRDS(paste0(res_path,"input/st_rctd.rds"))

celltype <- st_covid_rctd$Cluster %>% unique()

scales::show_col(pal_igv(palette = "default", alpha = 0.6)(15))
mycolors_nejm <- pal_igv(palette = "default", alpha = 0.6)(15)

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

inputdir <- paste0(res_path,'getPIM/')
files <- list.files(inputdir)[grep('_im_',list.files(inputdir))]
df_cellpair <- gsub('LRTG_im_clean_|.rds','',files) %>% strsplit(.,"-") %>% do.call('rbind',.) %>% as.data.frame()
colnames(df_cellpair) <- c('Sender','Receiver') 
df_cellpair$CP1 <- paste(df_cellpair$Sender,df_cellpair$Receiver,sep ='_')
df_cellpair$CP2 <- paste(df_cellpair$Receiver,df_cellpair$Sender,sep ='_')

celltypes <- unique(c(df_cellpair$Sender,df_cellpair$Receiver))
cp_of_inter <- data.frame(matrix(ncol = 4,dimnames = list(c(),c('ct1','ct2','keys_of_ct1','keys_of_ct2'))))

wd_path <- res_path
for (i in 1:length(celltypes)) {
  
  ct1 <- celltypes[i]
  for (j in (i+1):length(celltypes)) {
    
    ct2 <- celltypes[j]
    cp1 <- paste0(ct1,'_',ct2)
    cp2 <- paste0(ct2,'_',ct1)
    if(cp1 %in% list.files(paste0(wd_path,'/runscMLnet/')) & cp2 %in% list.files(paste0(wd_path,'/runscMLnet/'))){
      
      cat('check in ',cp1,' and ',cp2,'\n')
      mlnet1 <- readRDS(paste0(wd_path,'/runscMLnet/',cp1,'/scMLnet.rds'))
      mlnet2 <- readRDS(paste0(wd_path,'/runscMLnet/',cp2,'/scMLnet.rds'))
      
      ct2_tgs <- unique(mlnet1$TFTar$target)
      ct1_tgs <- unique(mlnet2$TFTar$target)
      
      ct1_ligs <- unique(mlnet1$LigRec$source)
      ct2_ligs <- unique(mlnet2$LigRec$source)
      
      ct1_keys <- intersect(ct1_ligs,ct1_tgs)
      ct2_keys <- intersect(ct2_ligs,ct2_tgs)
      
      # add by yll 2024-12-04
      ct1_keys_filter <- c()
      ct2_keys_filter <- c()
      for (key in ct1_keys) {
        
        res_ct1_key1 <- process_ml_net(mlnet1, lig_key = key)
        res_ct1_loop <- res_ct1_key1[res_ct1_key1$Target %in% ct2_keys, ]
        
        res_ct2_key2 <- process_ml_net(mlnet2, lig_key = ct2_keys)
        res_ct2_loop <- res_ct2_key2[res_ct2_key2$Target %in% key, ]
        
        if (nrow(res_ct2_loop) > 0) {
          ct1_keys_filter <- c(ct1_keys_filter, key)
        }
        
      }
      
      for (key in ct2_keys) {
        
        res_ct1_key1 <- process_ml_net(mlnet1, lig_key = ct1_keys)
        res_ct1_loop <- res_ct1_key1[res_ct1_key1$Target %in% key, ]
        
        res_ct2_key2 <- process_ml_net(mlnet2, lig_key = key)
        res_ct2_loop <- res_ct2_key2[res_ct2_key2$Target %in% ct1_keys, ]
        
        if (nrow(res_ct1_loop) > 0) {
          ct2_keys_filter <- c(ct2_keys_filter, key)
        }
      }
      
      ct1_keys <- ct1_keys_filter
      ct2_keys <- ct2_keys_filter
      
      
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
  
  mlnet1 <- readRDS(paste0(wd_path, '/runscMLnet/', paste(ct1, ct2, sep = '_'), '/scMLnet.rds'))
  mlnet2 <- readRDS(paste0(wd_path, '/runscMLnet/', paste(ct2, ct1, sep = '_'), '/scMLnet.rds'))
  
  ct2_tgs <- unique(mlnet1$TFTar$target)
  ct1_tgs <- unique(mlnet2$TFTar$target)
  ct1_ligs <- unique(mlnet1$LigRec$source)
  ct2_ligs <- unique(mlnet2$LigRec$source)
  
  ct1_keys <- intersect(ct1_ligs, ct1_tgs)
  ct2_keys <- intersect(ct2_ligs, ct2_tgs)
  
  # add by yll 2024-12-04
  ct1_keys_filter <- c()
  ct2_keys_filter <- c()
  for (key in ct1_keys) {
    
    res_ct1_key1 <- process_ml_net(mlnet1, lig_key = key)
    res_ct1_loop <- res_ct1_key1[res_ct1_key1$Target %in% ct2_keys, ]
    
    res_ct2_key2 <- process_ml_net(mlnet2, lig_key = ct2_keys)
    res_ct2_loop <- res_ct2_key2[res_ct2_key2$Target %in% key, ]
    
    if (nrow(res_ct2_loop) > 0) {
      ct1_keys_filter <- c(ct1_keys_filter, key)
    }
    
  }
  
  for (key in ct2_keys) {
    
    res_ct1_key1 <- process_ml_net(mlnet1, lig_key = ct1_keys)
    res_ct1_loop <- res_ct1_key1[res_ct1_key1$Target %in% key, ]
    
    res_ct2_key2 <- process_ml_net(mlnet2, lig_key = key)
    res_ct2_loop <- res_ct2_key2[res_ct2_key2$Target %in% ct1_keys, ]
    
    if (nrow(res_ct1_loop) > 0) {
      ct2_keys_filter <- c(ct2_keys_filter, key)
    }
  }
  
  ct1_keys <- ct1_keys_filter
  ct2_keys <- ct2_keys_filter
  
  key_of_inter[[paste(ct1, ct2, sep = '_')]] <- list(
    ct1_keys = ct1_keys,
    ct2_keys = ct2_keys
  )
}

key_of_inter$Alveoli_Macrophage
key_of_inter$Alveoli_Monocyte
key_of_inter$Macrophage_Monocyte

saveRDS(cp_of_inter,file = './feedback/cp_of_inter_yll_methods.rds')
saveRDS(key_of_inter,file = './feedback/key_of_inter_yll_methods.rds')

## transform to table ####

fbloop <- lapply(1:length(key_of_inter), function(i){
  
  ls <- key_of_inter[[i]]
  cp <- names(key_of_inter)[i]
  ct1 <- strsplit(cp,'_')[[1]][1]
  ct2 <- strsplit(cp,'_')[[1]][2]
  ct1_keys <- paste(ls$ct1_keys,collapse = ' ')
  ct2_keys <- paste(ls$ct2_keys,collapse = ' ')
  rbind(c(ct1,ct2,ct1_keys),c(ct2,ct1,ct2_keys))
  
}) %>% do.call('rbind',.) %>% as.data.frame()
colnames(fbloop) <- c('Sender','Receiver','Siganl')

cts_of_interest <- c('Alveoli','Macrophage','Monocyte')
fbloop <- fbloop[fbloop$Sender %in% cts_of_interest & fbloop$Receiver %in% cts_of_interest,]

# write.csv(fbloop,'./feedback/feddback_loop_detail.csv')

# plot feedback loop network 
key_of_plot <- key_of_inter[names(key_of_inter) %in% c('Alveoli_Macrophage','Alveoli_Monocyte','Macrophage_Monocyte')]
#key_genes <- c("C3/GPC3", "FN1/TGM2", "C3/GPC3", "NID1/VCAN", "FN1/TGM2", "NID1/VCAN")

for (cp in names(key_of_plot)){
  ct1 <- stringr::str_split(cp,"_",simplify = T)[,1]
  ct2 <- stringr::str_split(cp,"_",simplify = T)[,2]
  
  ct1_key <- key_of_plot[[cp]]$ct1_keys
  ct2_key <- key_of_plot[[cp]]$ct2_keys
  
  if (length(ct1_key)>1|length(ct1_key)>1){
    for (j in ct1_key){
      for (k in ct2_key){
        plot_feedback_network(wd_path = wd_path,ct1,ct2,j,k,vertex.size=38 
                              ,rescale=FALSE,do.check=T)
      }
    }
  }else{
    plot_feedback_network(wd_path = wd_path,ct1,ct2,ct1_key,ct2_key,vertex.size=38
                          ,rescale=FALSE,do.check=T)
    
  }
}

# plot feedback loop between cell type

plot_feedback_loop(floop,wd_path=getwd())




