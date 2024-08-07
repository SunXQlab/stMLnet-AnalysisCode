#############
## library ##
#############

library(Giotto)
library(plotly)

setwd('/home/cjy/project/')
rm(list = ls())
gc()

source('../code/code.R')

# Giotto-seqFISH+ --------------------------------------------------------------

# choose your directory
my_working_dir = paste0(getwd(),'/giotto_seqfish_dataset/')

# download
getSpatialDataset(dataset = 'seqfish_SS_cortex', directory = my_working_dir, method = 'wget')

# check data
files <- list.files(my_working_dir)
files

df_loca <- read.table(paste0(my_working_dir,'cortex_svz_centroids_coord.txt'),header = T)
head(df_loca)
plot(df_loca$X,df_loca$Y)

df_count <- read.table(paste0(my_working_dir,'cortex_svz_expression.txt'),header = T)
df_count[1:4,1:4]
dim(df_count)
range(df_count)

df_meta <- read.table(paste0(my_working_dir,'cortex_svz_centroids_annot.txt'),header = T,sep = '\t')
head(df_meta)
table(df_meta$FOV)
table(df_meta$cell_types)

# Giotto --------------------------------------------------------------

# 1. Giotto object

# 1.1. (optional) set Giotto instructions
temp_dir = paste0(getwd(),'/giotto_seqfish_dataset/')
instrs = createGiottoInstructions(save_plot = TRUE,
                                  show_plot = TRUE,
                                  save_dir = temp_dir)

# 1.2. create giotto object from provided paths
expr_path = paste0(my_working_dir, "cortex_svz_expression.txt")
loc_path = paste0(my_working_dir, "cortex_svz_centroids_coord.txt")
meta_path = paste0(my_working_dir, "cortex_svz_centroids_annot.txt")

# 1.3. This dataset contains multiple field of views which need to be stitched together
## first merge location and additional metadata
SS_locations = data.table::fread(loc_path)
cortex_fields = data.table::fread(meta_path)
SS_loc_annot = data.table::merge.data.table(SS_locations, cortex_fields, by = 'ID')
SS_loc_annot[, ID := factor(ID, levels = paste0('cell_',1:913))]
data.table::setorder(SS_loc_annot, ID)

## 1.4. create file with offset information
my_offset_file = data.table::data.table(field = c(0, 1, 2, 3, 4, 5, 6),
                                        x_offset = c(0, 1654.97, 1750.75, 1674.35, 675.5, 2048, 675),
                                        y_offset = c(0, 0, 0, 0, -1438.02, -1438.02, 0))

## 1.5. create a stitch file
stitch_file = stitchFieldCoordinates(location_file = SS_loc_annot,
                                     offset_file = my_offset_file,
                                     cumulate_offset_x = T,
                                     cumulate_offset_y = F,
                                     field_col = 'FOV',
                                     reverse_final_x = F,
                                     reverse_final_y = T)
stitch_file = stitch_file[,.(ID, X_final, Y_final)]
my_offset_file = my_offset_file[,.(field, x_offset_final, y_offset_final)]

# 2. Preprocessing

## 2.1. create Giotto object
SS_seqfish <- createGiottoObject(raw_exprs = expr_path,
                                 spatial_locs = stitch_file,
                                 offset_file = my_offset_file,
                                 instructions = instrs)

## 2.2. add additional annotation if wanted
SS_seqfish = addCellMetadata(SS_seqfish,
                             new_metadata = cortex_fields,
                             by_column = T,
                             column_cell_ID = 'ID')

## 2.3. subset data to the cortex field of views
cell_metadata = pDataDT(SS_seqfish)
cortex_cell_ids = cell_metadata[FOV %in% 0:4]$cell_ID
SS_seqfish = subsetGiotto(SS_seqfish, cell_ids = cortex_cell_ids)
# 523 cells x 10000 genes

## 2.4. filter
SS_seqfish <- filterGiotto(gobject = SS_seqfish,
                           expression_threshold = 1,
                           gene_det_in_min_cells = 10,
                           min_det_genes_per_cell = 10,
                           expression_values = c('raw'),
                           verbose = T)

## 2.5. normalize
SS_seqfish <- normalizeGiotto(gobject = SS_seqfish, scalefactor = 6000, verbose = T)

## 2.6. add gene & cell statistics
SS_seqfish <- addStatistics(gobject = SS_seqfish)

## 2.7. adjust expression matrix for technical or known variables
SS_seqfish <- adjustGiottoMatrix(gobject = SS_seqfish, expression_values = c('normalized'),
                                 batch_columns = NULL, covariate_columns = c('nr_genes', 'total_expr'),
                                 return_gobject = TRUE,
                                 update_slot = c('custom'))

## 2.8. visualize
spatPlot(gobject = SS_seqfish,
         save_param = list(save_name = '2_a_spatplot'))

# 3. Dimension Reduction

## highly variable genes (HVG)
SS_seqfish <- calculateHVG(gobject = SS_seqfish, method = 'cov_loess', difference_in_cov = 0.1,
                           save_param = list(save_name = '3_a_HVGplot', base_height = 5, base_width = 5))

## select genes based on HVG and gene statistics, both found in gene metadata
gene_metadata = fDataDT(SS_seqfish)
featgenes = gene_metadata[hvg == 'yes' & perc_cells > 4 & mean_expr_det > 0.5]$gene_ID

## run PCA on expression values (default)
SS_seqfish <- runPCA(gobject = SS_seqfish, genes_to_use = featgenes, scale_unit = F, center = F)
screePlot(SS_seqfish, save_param = list(save_name = '3_b_screeplot'))

plotPCA(gobject = SS_seqfish,
        save_param = list(save_name = '3_c_PCA_reduction'))

## run UMAP and tSNE on PCA space (default)
SS_seqfish <- runUMAP(SS_seqfish, dimensions_to_use = 1:15, n_threads = 10)
plotUMAP(gobject = SS_seqfish,
         save_param = list(save_name = '3_d_UMAP_reduction'))

SS_seqfish <- runtSNE(SS_seqfish, dimensions_to_use = 1:15)
plotTSNE(gobject = SS_seqfish,
         save_param = list(save_name = '3_e_tSNE_reduction'))

# 4. Clustering

## sNN network (default)
SS_seqfish <- createNearestNetwork(gobject = SS_seqfish, dimensions_to_use = 1:15, k = 15)

## Leiden clustering
SS_seqfish <- doLeidenCluster(gobject = SS_seqfish, resolution = 0.4, n_iterations = 1000)
plotUMAP(gobject = SS_seqfish,
         cell_color = 'leiden_clus', show_NN_network = T, point_size = 2.5,
         save_param = list(save_name = '4_a_UMAP_leiden'))

## Leiden subclustering for specified clusters
SS_seqfish = doLeidenSubCluster(gobject = SS_seqfish, cluster_column = 'leiden_clus',
                                resolution = 0.2, k_neighbors = 10,
                                hvg_param = list(method = 'cov_loess', difference_in_cov = 0.1),
                                pca_param = list(expression_values = 'normalized', scale_unit = F),
                                nn_param = list(dimensions_to_use = 1:5),
                                selected_clusters = c(5, 6, 7),
                                name = 'sub_leiden_clus_select')

## set colors for clusters
subleiden_order = c( 1.1, 5.1, 5.2,  2.1, 3.1,
                     4.1, 6.2, 6.1,
                     7.1,  7.2, 9.1, 8.1)
subleiden_colors = Giotto:::getDistinctColors(length(subleiden_order))
names(subleiden_colors) = subleiden_order

plotUMAP(gobject = SS_seqfish,
         cell_color = 'sub_leiden_clus_select', cell_color_code = subleiden_colors,
         show_NN_network = T, point_size = 2.5, show_center_label = F,
         legend_text = 12, legend_symbol_size = 3,
         save_param = list(save_name = '4_b_UMAP_leiden_subcluster'))

## show cluster relationships
showClusterHeatmap(gobject = SS_seqfish, cluster_column = 'sub_leiden_clus_select',
                   save_param = list(save_name = '4_c_heatmap', units = 'cm'),
                   row_names_gp = grid::gpar(fontsize = 9), column_names_gp = grid::gpar(fontsize = 9))

showClusterDendrogram(SS_seqfish, h = 0.5, rotate = T, cluster_column = 'sub_leiden_clus_select',
                      save_param = list(save_name = '4_d_dendro', units = 'cm'))

# 5. Visualize the Spatial and Expression Space

# expression and spatial
spatDimPlot(gobject = SS_seqfish, cell_color = 'sub_leiden_clus_select',
            cell_color_code = subleiden_colors,
            dim_point_size = 2, spat_point_size = 2,
            save_param = list(save_name = '5_a_covis_leiden'))

# selected groups and provide new colors
groups_of_interest = c(6.1, 6.2, 7.1, 7.2)
group_colors = c('red', 'green', 'blue', 'purple'); names(group_colors) = groups_of_interest

spatDimPlot(gobject = SS_seqfish, cell_color = 'sub_leiden_clus_select',
            dim_point_size = 2, spat_point_size = 2,
            select_cell_groups = groups_of_interest, cell_color_code = group_colors,
            save_param = list(save_name = '5_b_covis_leiden_selected'))

# 6. Cell-Type Marker Gene Expression

## gini ##
gini_markers_subclusters = findMarkers_one_vs_all(gobject = SS_seqfish,
                                                  method = 'gini',
                                                  expression_values = 'normalized',
                                                  cluster_column = 'sub_leiden_clus_select',
                                                  min_genes = 20,
                                                  min_expr_gini_score = 0.5,
                                                  min_det_gini_score = 0.5)
topgenes_gini = gini_markers_subclusters[, head(.SD, 2), by = 'cluster']

# violinplot
violinPlot(SS_seqfish, genes = unique(topgenes_gini$genes), cluster_column = 'sub_leiden_clus_select',
           strip_text = 8, strip_position = 'right', cluster_custom_order = unique(topgenes_gini$cluster),
           save_param = c(save_name = '6_a_violinplot_gini', base_width = 5, base_height = 10))

# cluster heatmap
topgenes_gini2 = gini_markers_subclusters[, head(.SD, 6), by = 'cluster']
plotMetaDataHeatmap(SS_seqfish, selected_genes = unique(topgenes_gini2$genes),
                    custom_gene_order = unique(topgenes_gini2$genes),
                    custom_cluster_order = unique(topgenes_gini2$cluster),
                    metadata_cols = c('sub_leiden_clus_select'), x_text_size = 10, y_text_size = 10,
                    save_param = c(save_name = '6_b_metaheatmap_gini'))

# 7. Cell-Type Annotation

## general cell types
# create vector with names
clusters_cell_types_cortex = c('L6 eNeuron', 'L4 eNeuron', 'L2/3 eNeuron', 'L5 eNeuron',
                               'Lhx6 iNeuron', 'Adarb2 iNeuron',
                               'endothelial', 'mural',
                               'OPC','Olig',
                               'astrocytes', 'microglia')

names(clusters_cell_types_cortex) = c(1.1, 2.1, 3.1, 4.1,
                                      5.1, 5.2,
                                      6.1, 6.2,
                                      7.1, 7.2,
                                      8.1, 9.1)

SS_seqfish = annotateGiotto(gobject = SS_seqfish, annotation_vector = clusters_cell_types_cortex,
                            cluster_column = 'sub_leiden_clus_select', name = 'cell_types')

## color
celltype <- unique(SS_seqfish@cell_metadata$cell_types)

scales::show_col(pal_igv(palette = "default", alpha = 0.8)(15))
mycolors_nejm <- pal_igv(palette = "default", alpha = 0.8)(15)

mycolor_ct <- mycolors_nejm[1:length(celltype)]
names(mycolor_ct) <- celltype
scales::show_col(mycolor_ct)

## plot
pt <- spatDimPlot(gobject = SS_seqfish, cell_color = 'cell_types', cell_color_code = mycolor_ct,
            dim_point_size = 2, spat_point_size = 2, dim_show_cluster_center = F, dim_show_center_label = T,
            save_param = c(save_name = '7_b_covisualization_setColor', base_width = 6, base_height = 8))

pdf("./giotto_seqfish_dataset/spatDimPlot_celltype_setColor.pdf",width = 6,height = 8)
print(pt)
dev.off()

png("./giotto_seqfish_dataset/spatDimPlot_celltype_setColor.png",res = 80,width = 600, height = 800)
print(pt)
dev.off()

## save
saveRDS(SS_seqfish, file = paste0(my_working_dir, "giotto_seqfish_object.rds"))

# Inputs --------------------------------------------------------------------

# filtering
table(SS_seqfish@cell_metadata$cell_types)
keep_cells = slot(SS_seqfish, 'cell_ID')
keep_genes = slot(SS_seqfish, 'gene_ID')
seqfish_filter = subsetGiotto(SS_seqfish, cell_ids = keep_cells, gene_ids = keep_genes)

tmp_cell_types <- seqfish_filter@cell_metadata$cell_types
tmp_cell_types[tmp_cell_types == 'Adarb2 iNeuron'] <- 'Adarb2-iNeuron'
tmp_cell_types[tmp_cell_types == 'L2/3 eNeuron'] <- 'L2-L3-eNeuron'
tmp_cell_types[tmp_cell_types == 'L4 eNeuron'] <- 'L4-eNeuron'
tmp_cell_types[tmp_cell_types == 'L5 eNeuron'] <- 'L5-eNeuron'
tmp_cell_types[tmp_cell_types == 'L6 eNeuron'] <- 'L6-eNeuron'
tmp_cell_types[tmp_cell_types == 'Lhx6 iNeuron'] <- 'Lhx6-iNeuron'
table(tmp_cell_types)
seqfish_filter@cell_metadata$cell_types <- tmp_cell_types

# annotation
table(seqfish_filter@cell_metadata$cell_types)
df_anno = data.frame(Barcode=seqfish_filter@cell_metadata$cell_ID,
                     Cluster=seqfish_filter@cell_metadata$cell_types)
head(df_anno)

# expression
df_count = seqfish_filter@raw_exprs
df_count[1:4,1:4]

df_norm = seqfish_filter@norm_expr
df_norm[1:4,1:4]

# location
head(seqfish_filter@spatial_locs)
df_loca = data.frame(seqfish_filter@spatial_locs[,1:2])
rownames(df_loca) = seqfish_filter@spatial_locs$cell_ID
colnames(df_loca) = c('dim.x','dim.y')
head(df_loca)

# Signals
assayobj = CreateAssayObject(data = df_norm)
BarCluTable = df_anno

ICGs_list <- lapply(unique(df_anno$Cluster), function(Clu){

  BarListResult <- getBarList(Clu, df_norm, df_anno)
  Clus.1 <- BarListResult[[1]]
  Clus.2 <- BarListResult[[2]]
  
  DEGs <- FindMarkers(object = assayobj, cells.1 = Clus.1, cells.2 = Clus.2,
                      logfc.threshold = 0.25, min.pct = 0.1, verbose = T)
  DEGs <- DEGs %>% filter(p_val_adj <= 0.05)
  rownames(DEGs)
  
})
names(ICGs_list) <- unique(df_anno$Cluster)
ICGs_list

Databases <- readRDS('../prior_knowledge/output/Databases.rds')
ligs_in_db <- Databases$LigRec.DB$source %>% unique() %>% stringr::str_to_title()
ligs_in_db <- intersect(ligs_in_db, slot(SS_seqfish, 'gene_ID'))
recs_in_db <- Databases$LigRec.DB$target %>% unique() %>% stringr::str_to_title()
recs_in_db <- intersect(recs_in_db, slot(SS_seqfish, 'gene_ID'))

expr.ct <- 1
pct.ct <- 0.1
data <- as.matrix(SS_seqfish@norm_expr)
clusters <- tmp_cell_types %>% as.character() %>% unique()

abundant.cutoff = expr.ct
all_mean <- rowMeans(data)
hist(log10(all_mean), breaks=100, main="", col="grey80",
     xlab=expression(Log[10]~"average count"))
abline(v=log10(abundant.cutoff), col="red", lwd=2, lty=2)

meanExpr_of_LR <- lapply(clusters, function(cluster){
  
  cluster.ids <- df_anno$Barcode[tmp_cell_types == cluster]
  source_mean <- rowMeans(data[,cluster.ids])
  names(source_mean) <- rownames(data)
  source_mean
  
}) %>% do.call('cbind',.) %>% as.data.frame()
colnames(meanExpr_of_LR) <- clusters

pct_of_LR <- lapply(clusters, function(cluster){
  
  cluster.ids <- df_anno$Barcode[tmp_cell_types == cluster]
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

# save
save(df_anno,df_count,df_norm,df_loca,
     ICGs_list,Ligs_expr_list,Recs_expr_list,
     file = paste0(my_working_dir, "giotto_seqfish_output.rda"))
