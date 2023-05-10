
# ENV  -------------------------------------------------------------------

setwd('/home/cjy/project/giotto_merfish_dataset')
rm(list = ls())
gc()

# Compare --------------------------------------------------------------------

layers <- c('layer3','layer6','layer9','layer12')
floders <- list.files(path = './',pattern = paste(layers,collapse = '|'))

## summary

Ligs_summ <- list()
Recs_summ <- list()
TGs_summ <- list()
ct_summ <- list()
for (layer in layers) {
  
  # layer = 'layer3'
  load(paste0('giotto_merfish_dataset_',layer,'/giotto_merfish_output.rda'))
  
  cts <- unique(df_anno$Cluster) %>% sort()
  Ligs_summ[[layer]] <- lengths(Ligs_expr_list)[cts]
  Recs_summ[[layer]] <- lengths(Recs_expr_list)[cts]
  TGs_summ[[layer]] <- lengths(ICGs_list)[cts]
  ct_summ[[layer]] <- split(df_anno$Cluster,df_anno$Cluster) %>% lengths() %>% .[cts]

}

str(Ligs_summ)
str(Recs_summ)
str(TGs_summ)
str(ct_summ)

df_Ligs_summ <- do.call('cbind',Ligs_summ)
df_Recs_summ <- do.call('cbind',Recs_summ)
df_TGs_summ <- do.call('cbind',TGs_summ)
df_ct_summ <- do.call('cbind',ct_summ)

## visualize

cts <- rownames(df_TGs_summ)
TGs_list <- list()
for (layer in layers) {
  
  # layer = 'layer3'
  load(paste0('giotto_merfish_dataset_',layer,'/giotto_merfish_output.rda'))
  for (ct in cts) {
    TGs_list[[ct]][[layer]] <- unlist(ICGs_list[ct])
  }

}
str(TGs_list)

library(ggvenn)
ct = 'Astrocyte'
ggvenn(TGs_list[[ct]],layers)
ct = 'Excitatory'
ggvenn(TGs_list[[ct]],layers)
ct = 'Inhibitory'
ggvenn(TGs_list[[ct]],layers)

