
# ENV -------------------------------------------------------------------

library(Giotto)
library(plotly)

setwd('/home/cjy/project/')
rm(list = ls())
gc()

# Giotto-MERFISH -----------------------------------------------------------------


# choose your directory
my_working_dir = paste0(getwd(),'/giotto_merfish_dataset/')

# avoid certification issues with wget
getSpatialDataset(dataset = 'merfish_preoptic', directory = my_working_dir, method = 'wget', extra = '--no-check-certificate')

# check data
files <- list.files(my_working_dir)
files

# location
df_loca <- read.table(paste0(my_working_dir, "merFISH_3D_data_cell_locations.txt"),header = T)
head(df_loca)
plot_ly(x = df_loca$x, y = df_loca$y, z = df_loca$z, showscale = FALSE) #%>% add_surface()
unique(df_loca$z)
split(df_loca,df_loca$z) %>% lapply(.,nrow) %>% unlist()

# expression
df_count <- read.table(paste0(my_working_dir, "merFISH_3D_data_expression.txt.gz"),header = T)
df_count[1:4,1:4]
dim(df_count)
range(df_count)

# meta
df_meta <- read.table(paste0(my_working_dir, "merFISH_3D_metadata.txt"),header = T,sep = '\t')
head(df_meta)
table(df_meta$orig_cell_types)
unique(df_meta$orig_cell_types)

# Giotto --------------------------------------------------------------

# 1. Giotto object
# to automatically save figures in save_dir set save_plot to TRUE
temp_dir = paste0(getwd(),'/giotto_merfish_dataset/')
instrs = createGiottoInstructions(save_dir = temp_dir,
                                  save_plot = TRUE,
                                  show_plot = TRUE)

# create giotto object from provided paths
expr_path = paste0(my_working_dir, "merFISH_3D_data_expression.txt.gz")
loc_path = paste0(my_working_dir, "merFISH_3D_data_cell_locations.txt")
meta_path = paste0(my_working_dir, "merFISH_3D_metadata.txt")

# create Giotto object
merFISH_test <- createGiottoObject(raw_exprs = expr_path,
                                   spatial_locs = loc_path,
                                   instructions = instrs)


# add additional metadata if wanted
metadata = data.table::fread(meta_path)
merFISH_test = addCellMetadata(merFISH_test, new_metadata = metadata$layer_ID, vector_name = 'layer_ID')
merFISH_test = addCellMetadata(merFISH_test, new_metadata = metadata$orig_cell_types, vector_name = 'orig_cell_types')

# 2. Processing
# 2.1. pre-test filter parameters
filterDistributions(merFISH_test, detection = 'genes',
                    save_param = list(save_name = '2_a_distribution_genes'))
filterDistributions(merFISH_test, detection = 'cells',
                    save_param = list(save_name = '2_b_distribution_cells'))
filterCombinations(merFISH_test,
                   expression_thresholds = c(0,1e-6,1e-5),
                   gene_det_in_min_cells = c(500, 1000, 1500),
                   min_det_genes_per_cell = c(1, 5, 10),
                   save_param = list(save_name = '2_c_filter_combos'))

# 2.2. filter data
merFISH_test <- filterGiotto(gobject = merFISH_test,
                             gene_det_in_min_cells = 0,
                             min_det_genes_per_cell = 0)
## 2.3. normalize
merFISH_test <- normalizeGiotto(gobject = merFISH_test, scalefactor = 10000, verbose = T)
merFISH_test <- addStatistics(gobject = merFISH_test)
merFISH_test <- adjustGiottoMatrix(gobject = merFISH_test, expression_values = c('normalized'),
                                   batch_columns = NULL, covariate_columns = c('nr_genes', 'total_expr'),
                                   return_gobject = TRUE,
                                   update_slot = c('custom'))

# 2.4. save according to giotto instructions
# 2D
spatPlot(gobject = merFISH_test, point_size = 1.5,
         save_param = list(save_name = '2_d_spatial_locations2D'))
# 3D
spatPlot3D(gobject = merFISH_test, point_size = 2.0, axis_scale = 'real',
           save_param = list(save_name = '2_e_spatial_locations3D'))

# 3. Dimension Reduction
# only 155 genes, use them all (default)
merFISH_test <- runPCA(gobject = merFISH_test, genes_to_use = NULL, scale_unit = FALSE, center = TRUE)
screePlot(merFISH_test, save_param = list(save_name = '3_a_screeplot'))

merFISH_test <- runUMAP(merFISH_test, dimensions_to_use = 1:8, n_components = 3, n_threads = 4)

plotUMAP_3D(gobject = merFISH_test, point_size = 1.5,
            save_param = list(save_name = '3_b_UMAP_reduction'))

# 4. Clustering
## sNN network (default)
merFISH_test <- createNearestNetwork(gobject = merFISH_test, dimensions_to_use = 1:8, k = 15)
## Leiden clustering
merFISH_test <- doLeidenCluster(gobject = merFISH_test, resolution = 0.2, n_iterations = 200,
                                name = 'leiden_0.2')
plotUMAP_3D(gobject = merFISH_test, cell_color = 'leiden_0.2', point_size = 1.5, show_center_label = F,
            save_param = list(save_name = '4_a_UMAP_leiden'))

# 5. Co-Visualization
spatDimPlot3D(gobject = merFISH_test, show_center_label = F,
              cell_color = 'leiden_0.2', dim3_to_use = 3,
              axis_scale = 'real', spatial_point_size = 2.0,
              save_param = list(save_name = '5_a_covis_leiden'))

spatPlot2D(gobject = merFISH_test, point_size = 1.5,
           cell_color = 'leiden_0.2',
           group_by = 'layer_ID', cow_n_col = 2, group_by_subset = c(260, 160, 60, -40, -140, -240),
           save_param = list(save_name = '5_b_leiden_2D'))

spatDimPlot3D(gobject = merFISH_test, show_center_label = F,
              cell_color = 'orig_cell_types', dim3_to_use = 3,
              axis_scale = 'real', spatial_point_size = 2.0,
              save_param = list(save_name = '5_a_covis_celltype'))

spatPlot2D(gobject = merFISH_test, point_size = 1.5,
           cell_color = 'orig_cell_types',
           group_by = 'layer_ID', cow_n_col = 2, group_by_subset = c(260, 160, 60, -40, -140, -240),
           save_param = list(save_name = '5_b_celltype_2D'))

table(merFISH_test@cell_metadata$orig_cell_types,merFISH_test@cell_metadata$leiden_0.2)

# 6. Cell-Type Marker Gene Detection
markers = findMarkers_one_vs_all(gobject = merFISH_test,
                                 method = 'gini',
                                 expression_values = 'normalized',
                                 cluster_column = 'leiden_0.2',
                                 min_genes = 1, rank_score = 2)
markers[, head(.SD, 2), by = 'cluster']

ICGs_list = split(markers$genes,markers$cluster)
str(ICGs_list)

# violinplot
topgini_genes = unique(markers[, head(.SD, 2), by = 'cluster']$genes)
violinPlot(merFISH_test, genes = topgini_genes, cluster_column = 'leiden_0.2', strip_position = 'right',
           save_param = c(save_name = '6_a_violinplot'))

topgini_genes = unique(markers[, head(.SD, 6), by = 'cluster']$genes)
plotMetaDataHeatmap(merFISH_test, expression_values = 'scaled',
                    metadata_cols = c('leiden_0.2'),
                    selected_genes = topgini_genes,
                    save_param = c(save_name = '6_b_clusterheatmap_markers'))

# 7. Cell-Type Annotation and Visualization
# known markers and DEGs
selected_genes = c('Myh11', 'Klf4', 'Fn1', 'Cd24a', 'Cyr61', 'Nnat', 'Trh', 'Selplg', 'Pou3f2', 'Aqp4', 'Traf4',
                   'Pdgfra', 'Opalin', 'Mbp', 'Ttyh2', 'Fezf1', 'Cbln1', 'Slc17a6', 'Scg2', 'Isl1', 'Gad1')
cluster_order = c(6, 11, 10, 12, 4, 9, 8, 5, 13, 3, 7, 1, 2)

plotMetaDataHeatmap(merFISH_test, expression_values = 'scaled',
                    metadata_cols = c('leiden_0.2'),
                    selected_genes = selected_genes,
                    custom_gene_order = rev(selected_genes),
                    custom_cluster_order = cluster_order,
                    save_param = c(save_name = '7_a_clusterheatmap_markers'))

## name clusters
clusters_cell_types_hypo = c('Inhibitory', 'Ambiguous', 'Excitatory', 'Astrocyte','OD Mature', 'Endothelial',
                             'Inhibitory', 'OD Mature', 'OD Immature', 'Ependymal', 'Endothelial', 'Microglia', 'OD Mature')
names(clusters_cell_types_hypo) = as.character(sort(cluster_order))
merFISH_test = annotateGiotto(gobject = merFISH_test, annotation_vector = clusters_cell_types_hypo,
                              cluster_column = 'leiden_0.2', name = 'cell_types')

## show heatmap
plotMetaDataHeatmap(merFISH_test, expression_values = 'scaled',
                    metadata_cols = c('cell_types'),
                    selected_genes = selected_genes,
                    custom_gene_order = rev(selected_genes),
                    custom_cluster_order = clusters_cell_types_hypo,
                    save_param = c(save_name = '7_b_clusterheatmap_markers_celltypes'))

# 7.1. Visualization

## color
celltype <- unique(merFISH_test@cell_metadata$cell_types)
scales::show_col(pal_igv(palette = "default", alpha = 0.8)(15))
mycolors_nejm <- pal_igv(palette = "default", alpha = 0.8)(15)

mycolor_ct <- mycolors_nejm[1:length(celltype)]
names(mycolor_ct) <- c('OD Mature','Astrocyte','Endothelial','Inhibitory','OD Immature',
                       'Excitatory','Microglia','Ependymal','Ambiguous')
scales::show_col(mycolor_ct)

## plot
pt <- spatPlot2D(gobject = merFISH_test, point_size = 1.0,
           cell_color = 'cell_types', cell_color_code = mycolor_ct,
           group_by = 'layer_ID', cow_n_col = 2, group_by_subset = c(seq(260, -290, -100)),
           save_param = c(save_name = '7_e_spatPlot2D_cell_types_all_setColor'))

pdf("./giotto_merfish_dataset/spatPlot2D_celltype_setColor.pdf",width = 14,height = 16)
print(pt)
dev.off()

png("./giotto_merfish_dataset/spatPlot2D_celltype_setColor.png",res = 80,width = 700, height = 800)
print(pt)
dev.off()
  
## save
saveRDS(merFISH_test, file = paste0(my_working_dir, "giotto_merfish_object.rds"))

# Inputs --------------------------------------------------------------------

layerIDs <- unique(merFISH_test@cell_metadata$layer_ID)
for (i in 1:length(layerIDs)) {
  
  layerID <- layerIDs[i]
  keep_cells = merFISH_test@cell_metadata$cell_ID[merFISH_test@cell_metadata$layer_ID == layerID]
  merFISH_filter <- subsetGiotto(merFISH_test, cell_ids = keep_cells)

  keep_cells = merFISH_filter@cell_metadata$cell_ID[merFISH_filter@cell_metadata$cell_types!='Ambiguous']
  keep_genes = slot(merFISH_filter, 'gene_ID')
  merFISH_filter = subsetGiotto(merFISH_filter, cell_ids = keep_cells, gene_ids = keep_genes)
  
  tmp_cell_types <- merFISH_filter@cell_metadata$cell_types
  tmp_cell_types[tmp_cell_types == 'OD Immature'] <- 'OD-Immature'
  tmp_cell_types[tmp_cell_types == 'OD Mature'] <- 'OD-Mature'
  table(tmp_cell_types)
  merFISH_filter@cell_metadata$cell_types <- tmp_cell_types
  
  table(merFISH_filter@cell_metadata$cell_types)
  df_anno = data.frame(Barcode=merFISH_filter@cell_metadata$cell_ID,
                       Cluster=merFISH_filter@cell_metadata$cell_types)
  head(df_anno)
  
  df_count = merFISH_filter@raw_exprs
  df_count[1:4,1:4]
  
  df_norm = merFISH_filter@norm_expr
  df_norm[1:4,1:4]
  
  head(merFISH_filter@spatial_locs)
  df_loca = data.frame(merFISH_filter@spatial_locs[,1:2])
  rownames(df_loca) = merFISH_filter@spatial_locs$cell_ID
  colnames(df_loca) = c('dim.x','dim.y')
  head(df_loca)
  
  ICGs_list = findMarkers_one_vs_all(gobject = merFISH_filter,
                                     method = 'gini',
                                     expression_values = 'normalized',
                                     cluster_column = 'cell_types',
                                     min_genes = 1, rank_score = 2)
  ICGs_list = split(ICGs_list$gene, ICGs_list$cluster)
  str(ICGs_list)
  
  Databases <- readRDS('../prior_knowledge/output/Databases.rds')
  ligs_in_db <- Databases$LigRec.DB$source %>% unique() %>% stringr::str_to_title()
  ligs_in_db <- intersect(ligs_in_db, slot(merFISH_single, 'gene_ID'))
  recs_in_db <- Databases$LigRec.DB$target %>% unique() %>% stringr::str_to_title()
  recs_in_db <- intersect(recs_in_db, slot(merFISH_single, 'gene_ID'))
  
  expr.ct <- 0.5
  pct.ct <- 0.1
  data <- as.matrix(merFISH_filter@norm_expr)
  clusters <- df_anno$Cluster %>% as.character() %>% unique()
  
  abundant.cutoff = expr.ct
  all_mean <- rowMeans(data)
  hist(log10(all_mean), breaks=100, main="", col="grey80",
       xlab=expression(Log[10]~"average count"))
  abline(v=log10(abundant.cutoff), col="red", lwd=2, lty=2)
  
  meanExpr_of_LR <- lapply(clusters, function(cluster){
    
    cluster.ids <- df_anno$Barcode[df_anno$Cluster == cluster]
    source_mean <- rowMeans(data[,cluster.ids])
    names(source_mean) <- rownames(data)
    source_mean
    
  }) %>% do.call('cbind',.) %>% as.data.frame()
  colnames(meanExpr_of_LR) <- clusters
  
  pct_of_LR <- lapply(clusters, function(cluster){
    
    cluster.ids <- df_anno$Barcode[df_anno$Cluster == cluster]
    dat <- data[,cluster.ids]
    pct <- rowSums(dat>0)/ncol(dat)
    names(pct) <- rownames(data)
    pct
    
  }) %>% do.call('cbind',.) %>% as.data.frame()
  colnames(pct_of_LR) <- clusters
  
  Recs_expr_list <- lapply(clusters, function(cluster){
    
    recs <- rownames(data)[meanExpr_of_LR[,cluster] >= expr.ct & pct_of_LR[,cluster] >= pct.ct]
    intersect(recs, recs_in_db)
    
  })
  names(Recs_expr_list) <- clusters
  str(Recs_expr_list)
  
  Ligs_expr_list <- lapply(clusters, function(cluster){
    
    ligs <- rownames(data)[meanExpr_of_LR[,cluster] >= expr.ct & pct_of_LR[,cluster] >= pct.ct]
    intersect(ligs, ligs_in_db)
    
  })
  names(Ligs_expr_list) <- clusters
  str(Ligs_expr_list)
  
  rownames(df_count) <- toupper(rownames(df_count))
  rownames(df_norm) <- toupper(rownames(df_norm))
  ICGs_list <- lapply(ICGs_list, toupper)
  Ligs_expr_list <- lapply(Ligs_expr_list, toupper)
  Recs_expr_list <- lapply(Recs_expr_list, toupper)
  
  my_working_dir <- paste0(getwd(),'/giotto_merfish_dataset_layer',i,'/')
  dir.create(my_working_dir,recursive = T)
  
  save(df_anno,df_count,df_norm,df_loca,
       ICGs_list,Ligs_expr_list,Recs_expr_list,
       file = paste0(my_working_dir,"/giotto_merfish_output.rda"))
  
  setwd('/home/cjy/project/giotto_merfish_dataset')
  
}