#############
## library ##
#############

library(Seurat)
library(clusterProfiler)
library(enrichplot)
library(GSEABase) 
library(ggsci)

rm(list = ls())
gc()

setwd('E:/stMLnet/apply_in_scGBM/')

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
  
  gmtfile ='./msigDB/c5.go.bp.v7.4.symbols.gmt'
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

seur <- readRDS("./vaild_scRNAseq/2019_cell_scRNAseq/result/seur.rds")
DefaultAssay(seur) <- 'alra'

clusters <- seur@active.ident %>% as.character() %>% unique()
df_markers <- lapply(clusters, function(cluster){
  
  df <- FindMarkers(seur, ident.1 = cluster, min.pct = 0.05, logfc.threshold = 0.1)
  df$gene <- rownames(df)
  df$ident.1 <- cluster
  df
  
}) %>% do.call('rbind',.)
str(df_markers)

df_markers_macro <- df_markers[df_markers$ident.1 == 'macrophages',]
gl_macro <- df_markers_macro$avg_log2FC
names(gl_macro) = df_markers_macro$gene
gl_macro = gl_macro[unique(LRTG_im$Target)]

## run ####

res_go_logfc <- runGSEA(gl_macro)
result_go_logfc <- res_go_logfc@result
result_go_logfc <- result_go_logfc[result_go_logfc$pvalue <= 0.01 ,]
result_go_logfc$ID %>% gsub('_',' ',.)

## plot ####

hit_terms <- result_go_logfc$ID %>% gsub('_',' ',.)
hit_terms[grep("myeloid|macro",hit_terms,ignore.case = T)]

hits <- result_go_logfc$ID[grep("myeloid|macro",result_go_logfc$ID,ignore.case = T)]
pdf("./visualize_CCI/gsea_go_logfc_tam_pval0.01.pdf", width = 18, height = 12)
gseaplot2(res_go_logfc, geneSetID = hits, pvalue_table = T, base_size = 18, color = mycols[1:length(hits)])
dev.off()

png("./visualize_CCI/gsea_go_logfc_tam_pval0.01.png", width = 18, height = 12,units = 'in',res = 300)
gseaplot2(res_go_logfc, geneSetID = hits, pvalue_table = T, base_size = 18, color = mycols[1:length(hits)])
dev.off()

############
## TAM_TC ##
############

## get logfc ####

seur <- readRDS("./vaild_scRNAseq/2019_cell_scRNAseq/result/seur.rds")
DefaultAssay(seur) <- 'alra'

clusters <- seur@active.ident %>% as.character() %>% unique()
df_markers <- lapply(clusters, function(cluster){
  
  df <- FindMarkers(seur, ident.1 = cluster, min.pct = 0.05, logfc.threshold = 0.1)
  df$gene <- rownames(df)
  df$ident.1 <- cluster
  df
  
}) %>% do.call('rbind',.)
str(df_markers)

df_markers_mali <- df_markers[df_markers$ident.1 == 'Malignant',]
gl_mali <- df_markers_mali$avg_log2FC
names(gl_mali) <- df_markers_mali$gene
gl_mali <- gl_mali[unique(LRTG_im$Target)]
gl_mali <- na.omit(gl_mali)

## run ####

res_go_logfc <- runGSEA(gl_mali)
result_go_logfc <- res_go_logfc@result
result_go_logfc <- result_go_logfc[result_go_logfc$pvalue <= 0.01 ,]
result_go_logfc$ID %>% gsub('_',' ',.)

## plot ####

hit_terms <- result_go_logfc$ID %>% gsub('_',' ',.)
hit_terms[grep("neuron|vascu",hit_terms,ignore.case = T)] 

hits <- result_go_logfc$ID[grep("neuron|vascu",result_go_logfc$ID,ignore.case = T)][c(1,2,5,7)]
pdf("./visualize_CCI/gsea_go_logfc_tc_pval0.01.pdf", width = 18, height = 12)
gseaplot2(res_go_logfc, geneSetID = hits, pvalue_table = T, base_size = 18, color = mycols[1:length(hits)])
dev.off()

png("./visualize_CCI/gsea_go_logfc_tc_pval0.01.png", width = 18, height = 12, units = 'in', res=300)
gseaplot2(res_go_logfc, geneSetID = hits, pvalue_table = T, base_size = 18, color = mycols[1:length(hits)])
dev.off()
