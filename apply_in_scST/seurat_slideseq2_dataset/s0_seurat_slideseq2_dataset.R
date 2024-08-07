#############
#  library  #
#############

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

setwd('/home/cjy/project/seurat_slideseq2_dataset/')
rm(list = ls())
gc()

# Seurat-Slide-seq2 --------------------------------------------------------------

# choose your directory
my_working_dir = paste0(getwd(),'/seurat_slideseq2_dataset/')

InstallData("ssHippo")
slide.seq <- LoadData("ssHippo")

# 1. Data preprocessing

plot1 <- VlnPlot(slide.seq, features = "nCount_Spatial", pt.size = 0, log = TRUE) + NoLegend()
slide.seq$log_nCount_Spatial <- log(slide.seq$nCount_Spatial)
plot2 <- SpatialFeaturePlot(slide.seq, features = "log_nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

slide.seq <- SCTransform(slide.seq, assay = "Spatial", ncells = 3000, verbose = FALSE)
slide.seq <- RunPCA(slide.seq)
slide.seq <- RunUMAP(slide.seq, dims = 1:30)
slide.seq <- FindNeighbors(slide.seq, dims = 1:30)
slide.seq <- FindClusters(slide.seq, resolution = 0.3, verbose = FALSE)

plot1 <- DimPlot(slide.seq, reduction = "umap", label = TRUE)
plot2 <- SpatialDimPlot(slide.seq, stroke = 0)
plot1 + plot2

SpatialDimPlot(slide.seq, 
               cells.highlight = CellsByIdentities(
                 object = slide.seq, idents = c(2,5,12)
                 ), 
               facet.highlight = TRUE)

# 2. Integration with a scRNA-seq reference

ref <- readRDS(paste0(my_working_dir,'mouse_hippocampus_reference.rds'))
anchors <- FindTransferAnchors(reference = ref, query = slide.seq, normalization.method = "SCT",
                               npcs = 50)
predictions.assay <- TransferData(anchorset = anchors, refdata = ref$celltype, prediction.assay = TRUE,
                                  weight.reduction = slide.seq[["pca"]], dims = 1:50)
slide.seq[["predictions"]] <- predictions.assay

DefaultAssay(slide.seq) <- "predictions"
SpatialFeaturePlot(slide.seq, features = c("Dentate Principal cells", "CA3 Principal cells", "Entorhinal cortex",
                                           "Endothelial tip", "Ependymal", "Oligodendrocyte"), alpha = c(0.1, 1))

slide.seq$predicted.id <- GetTransferPredictions(slide.seq)
Idents(slide.seq) <- "predicted.id"
SpatialDimPlot(slide.seq, 
               cells.highlight = CellsByIdentities(
                 object = slide.seq, 
                 idents = c("CA3 Principal cells", 
                            "Dentate Principal cells", "Endothelial tip")
                 ), 
               facet.highlight = TRUE)

# 3. Identification of Spatially Variable Features

DefaultAssay(slide.seq) <- "SCT"
slide.seq <- FindSpatiallyVariableFeatures(slide.seq, 
                                           assay = "SCT", 
                                           slot = "scale.data", 
                                           features = VariableFeatures(slide.seq)[1:1000],
                                           selection.method = "moransi", 
                                           x.cuts = 100, y.cuts = 100)

SpatialFeaturePlot(slide.seq, 
                   features = head(
                     SpatiallyVariableFeatures(slide.seq, selection.method = "moransi"),6), ncol = 3, alpha = c(0.1, 1), max.cutoff = "q95")

# 4. annotation

## celltypes
tmp_cell_types <- slide.seq$predicted.id
tmp_cell_types[tmp_cell_types == 'CA1 Principal cells'] <- 'CA1-Principal-cells'
tmp_cell_types[tmp_cell_types == 'CA3 Principal cells'] <- 'CA3-Principal-cells'
tmp_cell_types[tmp_cell_types == 'Dentate hilum'] <- 'Dentate-hilum'
tmp_cell_types[tmp_cell_types == 'Dentate Principal cells'] <- 'Dentate-Principal-cells'
tmp_cell_types[tmp_cell_types == 'Endothelial stalk'] <- 'Endothelial-stalk'
tmp_cell_types[tmp_cell_types == 'Endothelial tip'] <- 'Endothelial-tip'
tmp_cell_types[tmp_cell_types == 'Entorhinal cortex'] <- 'Entorhinal-cortex'
tmp_cell_types[tmp_cell_types == 'Neurogenesis (SGZ)'] <- 'Neurogenesis'
tmp_cell_types[tmp_cell_types == 'Resident macrophage'] <- 'Resident-macrophage'
slide.seq$predicted.id.new <- tmp_cell_types

## color
celltype <- unique(slide.seq$predicted.id.new)

scales::show_col(pal_igv(palette = "default", alpha = 0.8)(20))
mycolors_nejm <- pal_igv(palette = "default", alpha = 0.8)(20)

mycolor_ct <- mycolors_nejm[1:length(celltype)]
names(mycolor_ct) <- celltype
scales::show_col(mycolor_ct)
mycolor_ct

## plot
pt1 <- SpatialDimPlot(slide.seq, group.by = 'predicted.id.new', pt.size.factor = 0.8,stroke = 0,cols = mycolor_ct)
pt2 <- SpatialDimPlot(slide.seq, group.by = 'predicted.id.new', pt.size.factor = 0.8,stroke = 0,cols = mycolor_ct) + NoLegend()
pt2

pdf("./seurat_slideseq2_dataset/SpatialDimPlot_celltype_setColor.pdf",width = 4,height = 4)
print(pt2)
dev.off()

png("./seurat_slideseq2_dataset/SpatialDimPlot_celltype_setColor.png",width = 400, height = 400,res = 10)
print(pt2)
dev.off()

# save
saveRDS(slide.seq, file = paste0(my_working_dir, "seurat_slideseq2_object.rds"))

# Inputs --------------------------------------------------------------------

## filtering
table(slide.seq$predicted.id.new)
slide.seq_filter = subset(slide.seq, subset = predicted.id.new != "Unassigned")
slide.seq_filter[['SCT']]
slide.seq_filter[['Spatial']]

range(slide.seq_filter[['Spatial']]@counts)

# annotation
table(slide.seq_filter$predicted.id.new)
df_anno = data.frame(Barcode=colnames(slide.seq_filter),
                     Cluster=slide.seq_filter$predicted.id.new)
head(df_anno)
unique(tmp_cell_types)

# expression
df_count = slide.seq_filter@assays$SCT@counts
df_count[1:4,1:4]
dim(df_count)

df_norm = slide.seq_filter@assays$SCT@data
df_norm[1:4,1:4]
dim(df_norm)

# location
head(slide.seq_filter@images$image@coordinates)
df_loca = data.frame(slide.seq_filter@images$image@coordinates[,1:2])
rownames(df_loca) = colnames(slide.seq_filter)
colnames(df_loca) = c('dim.x','dim.y')
head(df_loca)

# signals
st_markers <- lapply(unique(tmp_cell_types),function(clu){
  
  markers <- FindMarkers(slide.seq_filter, ident.1 = clu,logfc.threshold = 0.5, min.pct = 0.05) 
  markers$ident.1 <- clu
  markers$gene <- rownames(markers)
  markers
  
}) %>% do.call('rbind',.) %>% as.data.frame()
table(st_markers$p_val_adj <= 0.05,st_markers$ident.1)
df_markers <- st_markers[st_markers$p_val_adj<=0.05,]
ICGs_list <- split(df_markers$gene,df_markers$ident.1)
str(ICGs_list)

Databases <- readRDS('../prior_knowledge/output/Databases.rds')
ligs_in_db <- Databases$LigRec.DB$source %>% unique() #%>% stringr::str_to_title()
ligs_in_db <- intersect(ligs_in_db, rownames(slide.seq_filter))
recs_in_db <- Databases$LigRec.DB$target %>% unique() #%>% stringr::str_to_title()
recs_in_db <- intersect(recs_in_db, rownames(slide.seq_filter))

expr.ct <- 0.01
pct.ct <- 0.05
data <- as.matrix(df_norm)
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
my_working_dir <- '/home/cjy/project/seurat_slideseq2_dataset/'
save(ICGs_list,Ligs_expr_list,Recs_expr_list,
     file = paste0(my_working_dir, "seurat_slideseq2_output.rda"))
