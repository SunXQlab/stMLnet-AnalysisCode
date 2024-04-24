Select_Lig_Rec_TGs <- function(ExprMat, AnnoMat, LocaMat, Databases, python_path, min.pct, expr.ct, pct.ct){
  
  # workdir
  de_count = ExprMat
  de_coords = LocaMat
  de_cell_type = AnnoMat
  my_python_path = python_path
  myinstr = createGiottoInstructions(python_path = my_python_path)
  # Create a Giotto object 
  gio_bc <- createGiottoObject(raw_exprs = de_count,
                               spatial_locs = de_coords,
                               instructions = myinstr)
  
  gio_bc <- addCellMetadata(gobject = gio_bc, 
                            new_metadata = de_cell_type$Cluster, 
                            vector_name = "celltype")
  
  # normalize
  gio_bc <- normalizeGiotto(gio_bc)
  gio_bc <- addStatistics(gobject = gio_bc)
  
  # create a spatial Delaunay network (default)
  gio_bc = createSpatialNetwork(gobject = gio_bc, method = "Delaunay")
  
  # select top 25th highest expressing genes
  gene_metadata <- fDataDT(gio_bc)
  lapply(seq(0.1,1,0.05), quantile, x = gene_metadata$mean_expr_det, na.rm = T) %>% unlist()
  high_expressed_genes = gene_metadata[mean_expr_det > 0.75]$gene_ID  # mean_expr_det > 0.75
  
  # identify ICGs
  CPGscoresHighGenes =  findICG(gobject = gio_bc, 
                                selected_genes = high_expressed_genes,
                                spatial_network_name = 'Delaunay_network',
                                cluster_column = 'celltype',
                                diff_test = 'permutation', 
                                adjust_method = 'fdr',
                                nr_permutations = 500,
                                do_parallel = T, cores = 6)
  
  # filter ICGs
  CPGscoresFilt = filterICG(CPGscoresHighGenes, direction = "both")
  table(CPGscoresFilt$CPGscores$spec_int[CPGscoresFilt$CPGscores$cell_type=="Malignant"])
  ICGs_list = lapply(unique(CPGscoresFilt$CPGscores$cell_type), function(x){
    y=CPGscoresFilt$CPGscores[CPGscoresFilt$CPGscores$cell_type==x,]
    z=lapply(unique(y$int_cell_type), function(t){
      y$genes[y$int_cell_type==t]
    })
    names(z)=unique(y$int_cell_type)
    z
  })
  names(ICGs_list) = unique(CPGscoresFilt$CPGscores$cell_type)
  str(ICGs_list)
  
  # calculate ligand list
  rownames(de_cell_type) <- de_cell_type$Barcode
  st_bc_A_RCTD <- CreateSeuratObject(counts = de_count, meta.data = de_cell_type, assay = 'Spatial')
  st_bc_A_RCTD <- SCTransform(st_bc_A_RCTD, assay = 'Spatial')
  Idents(st_bc_A_RCTD) <- st_bc_A_RCTD@meta.data$Cluster
  
  ligs_in_db <- Databases$LigRec.DB$source %>% unique()
  ligs_in_db <- intersect(ligs_in_db, rownames(st_bc_A_RCTD))
  
  clusters <- st_bc_A_RCTD@active.ident %>% as.character() %>% unique()
  df_markers_ligs <- lapply(clusters, function(cluster){
    
    df <- FindMarkers(st_bc_A_RCTD, ident.1 = cluster, features = ligs_in_db, only.pos = T, 
                      min.pct = min.pct) 
    df$gene <- rownames(df)
    df$ident.1 <- cluster
    df
    
  }) %>% do.call('rbind',.)
  
  Ligs_up_list <- split(df_markers_ligs$gene,df_markers_ligs$ident.1)
  str(Ligs_up_list)
  
  # calculate receptor list
  data <- gio_bc@norm_expr
  BarCluTable <- gio_bc@cell_metadata[,1:2]
  colnames(BarCluTable) <- c('Barcode','Cluster')
  
  recs_in_db <- Databases$LigRec.DB$target %>% unique()
  
  clusters <- BarCluTable$Cluster %>% as.character() %>% unique()
  
  meanExpr_of_LR <- lapply(clusters, function(cluster){
    
    cluster.ids <- BarCluTable$Barcode[BarCluTable$Cluster == cluster]
    source_mean <- rowMeans(data[,cluster.ids])
    names(source_mean) <- rownames(data)
    source_mean
    
  }) %>% do.call('cbind',.) %>% as.data.frame()
  colnames(meanExpr_of_LR) <- clusters
  
  pct_of_LR <- lapply(clusters, function(cluster){
    
    cluster.ids <- BarCluTable$Barcode[BarCluTable$Cluster == cluster]
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
  
  ex_inputs <- list(exprMat = data, annoMat = de_cell_type, locaMat = de_coords,
                    ligs_of_inter = Ligs_up_list, recs_of_inter = Recs_expr_list,
                    tgs_of_inter = ICGs_list)
  
  return(ex_inputs)
}
