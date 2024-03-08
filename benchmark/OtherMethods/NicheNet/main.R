library(nichenetr)
library(tidyverse)
library(Seurat)
options(stringsAsFactors = F)

rm(list = ls())
gc()

setwd("E:/stMLnet/other_method/NicheNet/scGBM/")

##############################################################################################

## prepare

# 载入表达矩阵和细胞注释
seur <- readRDS("E:/stMLnet/apply_in_scGBM/input/seur.rds")

# 表达矩阵
expression = seur@assays$alra@data
expression <- t(as.matrix(expression))

# 细胞注释
sample_info <- data.frame(
  cellType = as.vector(seur$orig.annotation),
  cellID = as.vector(rownames(seur@meta.data))
)

# 先验配体~靶基因模型
# download.file(url = "https://zenodo.org/record/3260758/files/ligand_target_matrix.rds",
#               destfile = './data/ligand_target_matrix.rds')
ligand_target_matrix = readRDS('../data/ligand_target_matrix.rds')
ligand_target_matrix[1:5,1:5]

## 定义表达基因(配体，受体)
receiver_ids = sample_info %>% filter(cellType == 'macrophages') %>% pull(cellID)
sender_ids = sample_info %>% filter(cellType != 'macrophages') %>% pull(cellID)

## 定义基因集(靶基因)
# 细胞相互作用潜在影响的靶基因
TGs_list <- readRDS('F:/finalVersion/vaild_scRNAseq/2016_science_bulk/output/DEGs_from_wilcox-cpm/TGs_list_wilcox-cpm_logfc2_pval0.1.rds')
geneset =  TGs_list$RebEP_TC$degs

# 单细胞数据定义高表达基因：平均值
expressed_genes_sender = (colSums(expression[sender_ids,]>0)/length(sender_ids)) > 0.1
expressed_genes_sender = names(expressed_genes_sender)[expressed_genes_sender]
expressed_genes_receiver = (colSums(expression[receiver_ids,]>0)/length(receiver_ids)) > 0.1
expressed_genes_receiver = names(expressed_genes_receiver)[expressed_genes_receiver]

## 背景基因
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
head(background_expressed_genes)

## 定义潜在先验配体(根据高表达配体，受体和先验数据库)
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
head(lr_network)

ligands = lr_network %>% pull(from) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_sender)

receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_receiver)

lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) 
head(lr_network_expressed)

potential_ligands = lr_network_expressed %>% pull(from) %>% unique()
head(potential_ligands)

## 配体活性分析(ligand activity analysis)
ligand_activities = predict_ligand_activities(geneset = geneset, 
                                              background_expressed_genes = background_expressed_genes, 
                                              ligand_target_matrix = ligand_target_matrix, 
                                              potential_ligands = potential_ligands)

ligand_activities = ligand_activities %>% arrange(-pearson) 

# 可视化配体活性(pearson)
p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) +
  geom_histogram(color="black", fill="darkorange")  +
  geom_density(alpha=.1, fill="orange") +
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(20, pearson) %>% pull(pearson))), 
             color="red", linetype="dashed", size=1) +
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()
p_hist_lig_activity

# 选取所有配体作为下游分析目标方便验证
best_upstream_ligands = ligand_activities %>% arrange(-pearson) %>% pull(test_ligand)
head(best_upstream_ligands)


# 获取激活配体-靶基因连接 Get the active ligand-target links
active_ligand_target_links_df = best_upstream_ligands %>% 
  lapply(get_weighted_ligand_target_links,
         geneset = geneset, 
         ligand_target_matrix = ligand_target_matrix, 
         n = nrow(ligand_target_matrix)) %>% 
  bind_rows()
active_ligand_target_links_df
active_ligand_target_links_df %>% group_by(ligand) %>% summarize(count=n())

## save
saveRDS(active_ligand_target_links_df, paste0("./result/macrophages_LRpair_weight.rds"))