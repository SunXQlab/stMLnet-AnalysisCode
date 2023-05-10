#############
## library ##
#############

library(Seurat)
library(clusterProfiler)
library(enrichplot)
library(GSEABase) 
library(ggsci)
library(Giotto)

rm(list = ls())
gc()

setwd("D:/AA-luluyan-phd/code/data_analysis/GBM_ST_Analysis/stMLnet_apply_in_stGBM")

###########
## color ##
###########

mycols <- pal_lancet("lanonc", alpha = 0.7)(9)
scales::show_col(mycols)

##############
## function ##
##############

runGSEA <- function(tgList){
  
  ## genelist
  
  geneList <- tgList
  geneList <- geneList[order(geneList, decreasing = T)]
  geneList <- geneList[geneList!=0] 
  geneList <- geneList[!is.na(geneList)] 
  
  ## msigDB
  
  gmtfile ='./msigDB/c5.go.bp.v7.5.1.symbols.gmt'
  geneset <- read.gmt(gmtfile)
  geneset$term <- geneset$term %>% tolower()
  length(unique(geneset$term))
  # 7481 genesets
  
  ## RUN ############
  
  gsea_gobp <- GSEA(geneList, TERM2GENE=geneset, 
                    minGSSize = 1, pvalueCutoff = 1, verbose=FALSE, seed = 10)
  
  return(gsea_gobp)
}

############
## TC_TAM ##
############

## get logfc ####

gio_stGBM <- readRDS("./giotto_stGBM.rds")

GCMat = gio_stGBM@norm_expr
BarCluTable = data.frame(Barcode=gio_stGBM@cell_metadata$cell_ID,
                         Cluster=gio_stGBM@cell_metadata$celltype)
clusters <- BarCluTable$Cluster %>% as.character() %>% unique()
rownames(BarCluTable) = gio_stGBM@cell_metadata$cell_ID

seur <- CreateSeuratObject(counts = GCMat, meta.data = BarCluTable, 
                           assay = "Spatial")
Idents(seur) <- seur$Cluster

clusters <- seur@active.ident %>% as.character() %>% unique()
df_markers <- lapply(clusters, function(cluster){
  
  df <- FindMarkers(seur, ident.1 = cluster, min.pct = 0, logfc.threshold = 0,assay = "Spatial")
  df$gene <- rownames(df)
  df$ident.1 <- cluster
  df
  
}) %>% do.call('rbind',.)
str(df_markers)

saveRDS(df_markers,file="./df_markers_logfc0.05.rds")

LRTG_im <- readRDS("./getPIM/LRTG_im_clean_Malignant_macrophages.rds")

df_markers_macro <- df_markers[df_markers$ident.1 == 'macrophages',]
gl_macro <- df_markers_macro$avg_log2FC
names(gl_macro) = df_markers_macro$gene
gl_macro = gl_macro[unique(LRTG_im$Target)] %>% na.omit()

## run ####

res_go_logfc <- runGSEA(gl_macro)
result_go_logfc <- res_go_logfc@result
result_go_logfc <- result_go_logfc[result_go_logfc$pvalue <= 0.01 ,]
result_go_logfc$ID %>% gsub('_',' ',.)

saveRDS(res_go_logfc,file="./visualize_CCI/res_GSEA_macro_Rec.rds")

## plot ####

hit_terms <- result_go_logfc$ID %>% gsub('_',' ',.)
hit_terms[grep("myeloid|macro",hit_terms,ignore.case = T)]
hit_terms[grep("neuron|vascu",hit_terms,ignore.case = T)]

hits <- result_go_logfc$ID[grep("myeloid|macro",result_go_logfc$ID,ignore.case = T)]
pdf("./visualize_CCI/gsea_go_logfc_tam_v6.pdf", width = 18, height = 12)
gseaplot2(res_go_logfc, geneSetID = hits, pvalue_table = T, base_size = 24, color = mycols[1:length(hits)])
dev.off()

png("./visualize_CCI/gsea_go_logfc_tam_pval0.01.png", width = 18, height = 12,units = 'in',res = 300)
gseaplot2(res_go_logfc, geneSetID = hits, pvalue_table = T, base_size = 18, color = mycols[1:length(hits)])
dev.off()

############
## TAM_TC ##
############

## get logfc ####
LRTG_im <- readRDS("./getPIM/LRTG_im_clean_macrophages_Malignant.rds")

df_markers_mali <- df_markers[df_markers$ident.1 == 'Malignant',]
gl_mali <- df_markers_mali$avg_log2FC
names(gl_mali) <- df_markers_mali$gene
gl_mali <- gl_mali[unique(LRTG_im$Target)]
gl_mali <- na.omit(gl_mali)

## run ####

res_go_logfc <- runGSEA(gl_mali)
result_go_logfc <- res_go_logfc@result
result_go_logfc <- result_go_logfc[result_go_logfc$pvalue <= 0.5 ,]
result_go_logfc$ID %>% gsub('_',' ',.)

saveRDS(res_go_logfc,file="./visualize_CCI/res_GSEA_maligant_Rec.rds")
## plot ####

hit_terms <- result_go_logfc$ID %>% gsub('_',' ',.)
hit_terms[grep("neuron|vascu",hit_terms,ignore.case = T)] 

hits <- result_go_logfc$ID[grep("neuron|vascu",hit_terms,ignore.case = T)][c(1,2,3,4)]
hits <- result_go_logfc$ID[c(1,2,3,4)]
pdf("./visualize_CCI/gsea_go_logfc_tc_v7.pdf", width = 18, height = 12)
gseaplot2(res_go_logfc, geneSetID = hits, pvalue_table = T, base_size = 24, color = mycols[1:length(hits)])
dev.off()

png("./visualize_CCI/gsea_go_logfc_tc_pval0.01_v5.png", width = 18, height = 12, units = 'in', res=300)
gseaplot2(res_go_logfc, geneSetID = hits, pvalue_table = T, base_size = 18, color = mycols[1:length(hits)])
dev.off()
