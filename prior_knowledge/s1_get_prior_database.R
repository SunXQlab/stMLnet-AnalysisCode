#############
## library ##
#############

library(dplyr)
library(tidyr)
library(Hmisc)

options(connectionObserver = NULL)
library(org.Hs.eg.db)

library(UpSetR)
library(VennDiagram)
library(ggplot2)
library(RColorBrewer)

library(graphite) 
library(doParallel)
library(foreach)

rm(list = ls())
gc()

setwd("./stMLnet/prior_knowledge/")

##############
## LigRecDB ##
##############

## CellChat ####

if(F){
  
  # download.file(url = 'https://github.com/sqjin/CellChat/raw/master/data/CellChatDB.human.rda',
  #               destfile = './LigRec/human_CellChatDB.rda', method = 'internal')
  load("./LigRec/human_CellChatDB.rda")
  
  hm.cc.interaction <- CellChatDB.human$interaction
  hm.cc.complex <- CellChatDB.human$complex
  hm.cc.cofactor <- CellChatDB.human$cofactor
  hm.cc.geneInfo <- CellChatDB.human$geneInfo
  
  
  hm.cc.complex$combine <- apply(hm.cc.complex,1,function(x){
    x = x[nchar(x)!=0]
    paste(x,collapse = "_")
  })
  
  hm.cc <- hm.cc.interaction[,c("ligand","receptor","evidence")]
  hm.cc <- hm.cc[!duplicated(hm.cc[,1:2]),]
  hm.cc$ligand <- lapply(hm.cc$ligand,function(x){
    x = ifelse(x %in% rownames(hm.cc.complex),hm.cc.complex[x,'combine'],x)
  })
  hm.cc$receptor <- lapply(hm.cc$receptor,function(x){
    x = ifelse(x %in% rownames(hm.cc.complex),hm.cc.complex[x,'combine'],x)
  })
  
  hm.cc <- hm.cc %>% mutate(ligand = strsplit(as.character(ligand), "_")) %>% unnest(ligand)
  hm.cc <- hm.cc %>% mutate(receptor = strsplit(as.character(receptor), "_")) %>% unnest(receptor)
  hm.cc <- hm.cc[!duplicated(hm.cc[,1:2]),]
  
  hm.cc <- hm.cc %>% mutate(primary.source = rep("CellChat",nrow(hm.cc))) %>% 
    mutate(secondary.source = evidence) %>%
    dplyr::select(source = ligand, target = receptor, 
                  primary.source, secondary.source, note = evidence)
  
  ### KEGG > PMID/PMC
  sort(unique(hm.cc$secondary.source))
  ### KEGG -> Pathwaydb
  unique(hm.cc$secondary.source[grep('KEGG',hm.cc$secondary.source,ignore.case = T)])
  hm.cc$secondary.source[grep('KEGG',hm.cc$secondary.source,ignore.case = T)] <- 'KEGG'
  ### PMID/PMC -> PMID
  unique(hm.cc$secondary.source[grep('PMID|PMC',hm.cc$secondary.source,ignore.case = T)])
  hm.cc$secondary.source[grep('PMID|PMC',hm.cc$secondary.source,ignore.case = T)] <- 'PMID'
  
  g2s <- toTable(org.Hs.egSYMBOL)
  filter_source <- hm.cc$source[!hm.cc$source%in%g2s$symbol]
  filter_target <- hm.cc$target[!hm.cc$target%in%g2s$symbol]
  filter_gene <- unique(c(filter_source, filter_target))
  filter_interaction_id <- unique(c(which(hm.cc$source %in% filter_gene),
                                    which(hm.cc$target %in% filter_gene)))
  filter_interaction <- hm.cc[filter_interaction_id,]
  hm.cc <- hm.cc[-na.omit(filter_interaction_id),]
  
  unique(hm.cc$source)
  unique(hm.cc$target)
  table(!duplicated(hm.cc[,1:2]))

  saveRDS(hm.cc, "./output/hm.LigRec.cc.rds")
  
}

hm.cc <- readRDS("./output/hm.LigRec.cc.rds")

## connectomeDB2020 ####

if(F){
  
  
  # download.file(url = 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-18873-z/MediaObjects/41467_2020_18873_MOESM4_ESM.xlsx',
  #               destfile = './LigRec/human_mouse_connectomeDB2020.xlsx', method = 'wininet')
  hm.connectomeDB <- xlsx::read.xlsx("./LigRec/human_mouse_connectomeDB2020.xlsx", 
                                     sheetName = 'literature_support')
  hm.connectomeDB <- hm.connectomeDB[,c(2,3,4,7)]
  hm.connectomeDB <- hm.connectomeDB[!duplicated(hm.connectomeDB[,3:4]),]
  
  hm.cn <- hm.connectomeDB %>% 
    mutate(primary.source = rep("connectomeDB2020",nrow(hm.connectomeDB))) %>%
    dplyr::select(source = Ligand.gene.symbol, target = Receptor.gene.symbol,
                  primary.source, secondary.source = Source, note = PMID.support) %>%
    as_tibble()
  
  hm.cn$secondary.source[grep('Ramilowski_2015',hm.cn$secondary.source)] <- "Ramilowski_2015"
  hm.cn$secondary.source[grep('Hou et al. 2020',hm.cn$secondary.source)] <- "Hou_2020"
  hm.cn$secondary.source[grep('RNA-Magnet',hm.cn$secondary.source)] <- "RNA-Magnet"
  hm.cn$secondary.source[grep('CellphoneDB',hm.cn$secondary.source)] <- "CellphoneDB"
  hm.cn$secondary.source[grep('SingleCellSignalR',hm.cn$secondary.source)] <- "SingleCellSignalR"
  hm.cn$secondary.source[grep('ICELLNET',hm.cn$secondary.source)] <- "ICELLNET"
  
  hm.cn$note <- paste0('PMID:',hm.cn$note)
  table(hm.cn$secondary.source)
  
  unique(hm.cn$source)
  hm.cn$source <- as.character(hm.cn$source)
  unique(hm.cn$target)
  hm.cn$target <- as.character(hm.cn$target)
  
  saveRDS(hm.cn, "./output/hm.LigRec.cn.rds")
}

hm.cn <- readRDS("./output/hm.LigRec.cn.rds")

## iTALK ####

if(F){
  
  # download.file(url = 'https://github.com/Coolgenome/iTALK/blob/master/data/LR_database.rda',
  #               destfile = './LigRec/human_mouse_iTALK.rda', method = 'internal')
  load("./LigRec/human_mouse_iTALK.rda")
  
  database <- na.omit(database)
  database <- database[!duplicated(database[,c(2,4)]),]
  
  hm.it <- tibble(database)
  hm.it <- hm.it[,c(1,2,4,6)]
  hm.it <- hm.it %>% mutate(primary.source = rep("ITALK",nrow(hm.it))) %>%
    mutate(secondary.source = rep("PMID",nrow(hm.it))) %>%
    dplyr::select(key = Pair.Name, source = Ligand.ApprovedSymbol, 
                  target = Receptor.ApprovedSymbol, 
                  primary.source, secondary.source, note = Classification)
  
  hm.cn.exc <- xlsx::read.xlsx("./LigRec/human_mouse_connectomeDB2020.xlsx", sheetName = 'excluded pairs')
  hm.cn.exc <- hm.cn.exc[hm.cn.exc$Primary.source == 'Ramilowski_2015',]
  hm.cn.exc.key <- paste(hm.cn.exc$Ligand.Approved.Symbol,hm.cn.exc$Receptor.Approved.Symbol,sep = "_")
  hm.it <- hm.it[!hm.it$key %in% hm.cn.exc.key,]
  table(!duplicated(hm.it$key))
  
  head(hm.it)
  unique(hm.it$source)
  hm.it$source <- as.character(hm.it$source)
  unique(hm.it$target)
  hm.it$target <- as.character(hm.it$target)
  hm.it <- hm.it[,-1]
  
  saveRDS(hm.it, "./output/hm.LigRec.it.rds")
  
}

hm.it <- readRDS("./output/hm.LigRec.it.rds")

## NicheNet ####

if(F){
  
  # download.file(url = 'https://github.com/saeyslab/nichenetr/blob/master/data/lr_network.rda',
  #               destfile = './LigRec/human_mouse_NicheNet_lr_network.rda', method = 'wininet')
  load("./LigRec/human_mouse_NicheNet_lr_network.rda")
  
  hm.nn <- lr_network
  table(hm.nn$database)
  table(hm.nn$source)
  
  hm.nn <- hm.nn %>% mutate(primary.source = rep("NicheNet",nrow(hm.nn))) %>% 
    mutate(key = paste(from, to, sep = "_")) %>%
    dplyr::select(key, source = from, target = to, 
                  primary.source, secondary.source = database, note = source)
  
  table(hm.nn$secondary.source)
  hm.nn <- hm.nn %>% filter(secondary.source != 'ppi_prediction') %>% 
    filter(secondary.source != 'ppi_prediction_go')
  
  sort(unique(hm.nn$secondary.source))
  ### kegg -> KEGG
  unique(hm.nn$secondary.source[grep('kegg',hm.nn$secondary.source,ignore.case = T)])
  hm.nn$secondary.source[grep('kegg',hm.nn$secondary.source,ignore.case = T)] <- 'KEGG'
  ### ramilowski -> Ramilowski_2015
  unique(hm.nn$secondary.source[grep('ramilowski',hm.nn$secondary.source,ignore.case = T)])
  hm.nn$secondary.source[grep('ramilowski',hm.nn$secondary.source,ignore.case = T)] <- 'Ramilowski_2015'
  ### guide2pharmacology -> Guide2PHARMACOLOGY
  unique(hm.nn$secondary.source[grep('guide2pharmacology',hm.nn$secondary.source,ignore.case = T)])
  hm.nn$secondary.source[grep('guide2pharmacology',hm.nn$secondary.source,ignore.case = T)] <- 'Guide2PHARMACOLOGY'
  
  hm.cn.exc <- xlsx::read.xlsx("./LigRec/human_mouse_connectomeDB2020.xlsx", sheetName = 'excluded pairs')
  hm.cn.exc <- hm.cn.exc[hm.cn.exc$Primary.source == 'Ramilowski_2015',]
  hm.cn.exc.key <- paste(hm.cn.exc$Ligand.Approved.Symbol,hm.cn.exc$Receptor.Approved.Symbol,sep = "_")
  hm.nn <- hm.nn[!(hm.nn$key %in% hm.cn.exc.key & hm.nn$secondary.source == 'Ramilowski_2015'),]
  table(!duplicated(hm.nn$key))
  
  head(hm.nn)
  unique(hm.nn$source)
  unique(hm.nn$target)
  unique(hm.nn$secondary.source)
  hm.nn <- hm.nn[,-1]

  saveRDS(hm.nn, "./output/hm.LigRec.nn.rds")
  
}

hm.nn <- readRDS("./output/hm.LigRec.nn.rds")

## merge ####

hm.LigRec <- do.call('rbind', list(hm.LigRec.cn, hm.LigRec.it, hm.LigRec.nn, hm.LigRec.cc))
save(hm.LigRec, file = "./output/scMLnet.human.LigRec.rda")


if(T){
  
  load("./output/scMLnet.human.LigRec.rda")
  
  # veen for source
  
  table(hm.LigRec$source)
  listinput <- list(iTALK = unique(hm.LigRec$source[hm.LigRec$primary.source == 'ITALK']),
                    NicheNet = unique(hm.LigRec$source[hm.LigRec$primary.source == 'NicheNet']),
                    connectomeDB2020 = unique(hm.LigRec$source[hm.LigRec$primary.source == 'connectomeDB2020']),
                    CellChat = unique(hm.LigRec$source[hm.LigRec$primary.source == 'CellChat']))
  
  p <- venn.diagram(listinput, resolution = 1500, imagetype = "tiff",
                    main.cex = 2, main.pos = c(0.5, 1.05), cat.cex = 1.4, cex = 1.4,
                    fill=c('Tomato','SkyBlue','PaleGreen2','MediumPurple'), alpha=c(0.5,0.5,0.5,0.5),
                    main="Intersection of Ligands \ncollected from different databsources",
                    filename = NULL)
  pdf("./figure/veen_Ligand_from_datasources.pdf",width = 7,height = 7)
  grid.draw(p)
  dev.off()
  
  # veen for target
  
  table(hm.LigRec$target)
  listinput <- list(ITALK = unique(hm.LigRec$target[hm.LigRec$primary.source == 'ITALK']),
                    NicheNet = unique(hm.LigRec$target[hm.LigRec$primary.source == 'NicheNet']),
                    connectomeDB2020 = unique(hm.LigRec$target[hm.LigRec$primary.source == 'connectomeDB2020']),
                    CellChat = unique(hm.LigRec$target[hm.LigRec$primary.source == 'CellChat']))
  
  p <- venn.diagram(listinput, resolution = 1500, imagetype = "tiff",
                    main.cex = 2, main.pos = c(0.5, 1.05), cat.cex = 1.4, cex = 1.4,
                    fill=c('Tomato','SkyBlue','PaleGreen2','MediumPurple'), alpha=c(0.5,0.5,0.5,0.5),
                    main="Intersection of Receptors \ncollected from different datasources",
                    filename = NULL)
  pdf("./figure/veen_Receptor_from_datasources.pdf",width = 7,height = 7)
  grid.draw(p)
  dev.off()
  
  # doupie for intersection of primary.source and secondary.source
  
  hm.count.ps <- as.data.frame(table(hm.LigRec$primary.source, hm.LigRec$secondary.source))
  hm.count.ps$primary.source <- as.vector(hm.count.ps$Var1)
  hm.count.ps$secondary.source <- as.vector(hm.count.ps$Var2)
  hm.count.ps <- hm.count.ps[order(hm.count.ps$primary.source, hm.count.ps$Freq, decreasing = T),]
  hm.count.ps$ymax <- unlist(lapply(seq(nrow(hm.count.ps)), function(i){sum(hm.count.ps$Freq[seq(i)])}))
  hm.count.ps$ymin <- c(0,hm.count.ps$ymax[-nrow(hm.count.ps)])
  hm.count.ps <- hm.count.ps[,c(4,5,3,6,7)]
  mycols <- c('Tomato','SkyBlue','PaleGreen2','MediumPurple',
              brewer.pal(length(unique(c(hm.LigRec$secondary.source))), 'Set3'))
  names(mycols) <- unique(c(hm.LigRec$primary.source, hm.LigRec$secondary.source))
  
  pdf(file='./figure/doupie_LigRec_source.pdf', width = 9, height = 6)
  p <- ggplot(hm.count.ps) + 
    geom_rect(aes(fill = secondary.source, ymax = ymax, ymin = ymin, xmax = 4, xmin = 3)) +
    geom_rect(aes(fill = primary.source, ymax = ymax, ymin = ymin, xmax = 3, xmin = 0)) +
    scale_fill_manual(name="data.source", values = mycols) +
    xlim(c(0, 4)) + theme(aspect.ratio = 1) + 
    coord_polar(theta = 'y') + theme_bw() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.title = element_text(size = 12)
    )
  p
  dev.off()
  
  # node and edge
  
  hm.count <- hm.count.ps[,1:3]
  hm.count <- hm.count[!hm.count$Freq==0,]
  hm.count <- apply(hm.count, 1, function(hm.count.x){
    
    hm.LigRec.x <- hm.LigRec[hm.LigRec$primary.source == hm.count.x[1] & 
                               hm.LigRec$secondary.source == hm.count.x[2],]
    
    nr_nodes <- length(unique(c(hm.LigRec.x$source, hm.LigRec.x$target)))
    nr_source_nodes <- length(unique(hm.LigRec.x$source))
    nr_target_nodes <- length(unique(hm.LigRec.x$target))
    
    c(nr_nodes, nr_source_nodes, nr_target_nodes)
    
  }) %>% t() %>% cbind(hm.count,.)
  colnames(hm.count)[-c(1:2)] <- c('nr_edges', 'nr_nodes', 'nr_source_nodes', 'nr_target_nodes')
  write.csv(hm.count,"./output/hm.LigRec.detail.csv", row.names = F)
  
}

############
## TFTGDB ##
############
## TRRUST ####

if(F){
  
  # download.file(url = 'https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv', 
  #               destfile = './TFTG/human_TRRUST.tsv', method = 'internal')
  hm.tt <- read.table("./TFTG/human_TRRUST.tsv",header = F,sep = "\t")
  colnames(hm.tt) <- c('TF', 'Target', 'Mode of Regulation', 'PMID')
  
  hm.tt$key <- paste(hm.tt$TF, hm.tt$Target, sep = "_")
  sort(table(hm.tt$key), decreasing = T)
  hm.tt <- hm.tt[order(hm.tt$key,hm.tt$`Mode of Regulation`),]
  hm.tt <- hm.tt[!duplicated(hm.tt[,1:2]),]
  
  hm.tt <- as_tibble(hm.tt)
  hm.tt$PMID <- paste0("PMID:",hm.tt$PMID)
  hm.tt <- hm.tt %>% mutate(primary.source = rep("TRRUST",nrow(hm.tt))) %>% 
    mutate(secondary.source = rep("PMID",nrow(hm.tt))) %>% 
    dplyr::select(source = TF, target = Target, primary.source, secondary.source, 
                  note = PMID, note2 = `Mode of Regulation`)
  
  saveRDS(hm.tt, "./output/hm.TFTG.tt.rds")
  
}

hm.tt <- readRDS("./output/hm.TFTG.tt.rds")

## HTRIdb ####

if(F){

  ### download from http://www.lbbc.ibb.unesp.br/htri
  hm.htridb <- read.table("./TFTG/human_HTRIdb.txt", 
                          header = F, skip = 1, sep = "\t", stringsAsFactors = F)
  hm.htridb <- hm.htridb[,seq(1,14,by = 2)]
  colnames(hm.htridb) <- read.table("./TFTG/human_HTRIdb.txt", 
                                    header = F, nrows = 1,sep = "\t")[-8] %>% unlist()
  head(hm.htridb)
  
  hm.ht <- hm.htridb
  hm.ht <- hm.ht %>% mutate(primary.source = rep("HTRIdb",nrow(hm.ht))) %>% 
    mutate(secondary.source = rep("PMID", nrow(hm.ht))) %>%
    mutate(note = paste0("PMID:", PUBMED_ID)) %>%
    dplyr::select(source = SYMBOL_TF, target = SYMBOL_TG, primary.source, secondary.source, note) %>%
    as_tibble()
  
  hm.ht$key <- paste(hm.ht$source, hm.ht$target, sep = "_")
  table(duplicated(hm.ht$key))
  hm.ht <- hm.ht[order(hm.ht$key),]
  hm.ht <- hm.ht[!duplicated(hm.ht[,1:2]),]
  hm.ht <- hm.ht[,-6]
  
  saveRDS(hm.ht, "./output/hm.TFTG.ht.rds")
  
}

hm.ht <- readRDS("./output/hm.TFTG.ht.rds")

## RegNetwork ####

if(F){
  
  ### download from RegNet: http://www.regnetworkweb.org/home.jsp
  hm.regnet.pred <- read.csv("./TFTG/RegNetwork/human_RegNetwork_predicted.csv", stringsAsFactors = F)
  hm.regnet.expe <- read.csv("./TFTG/RegNetwork/human_RegNetwork_experimental.csv", stringsAsFactors = F)
  hm.regnet <- rbind(hm.regnet.pred, hm.regnet.expe)
  head(hm.regnet)
  
  hm.rn <- hm.regnet[,c(1,3,5,6)]
  hm.rn <- hm.rn %>% mutate(primary.source = rep("RegNetwork",nrow(hm.rn))) %>%
    dplyr::select(source = regulator_symbol, target = target_symbol, 
                  primary.source, secondary.source = database, note = evidence) %>% 
    as_tibble()
  
  hm.rn.dbs <- unique(hm.rn$secondary.source) %>% strsplit(.,',') %>% 
    unlist() %>% table()
  ### Pathwaydb > Regdb > PPI > Annodb
  ### KEGG -> Pathwaydb
  unique(hm.rn$secondary.source[grep('kegg',hm.rn$secondary.source,ignore.case = T)])
  hm.rn$secondary.source[grep('kegg',hm.rn$secondary.source,ignore.case = T)] <- 'KEGG'
  ### MicroT/miRanda/miRBase/miRecords/miRTarBase/PicTar/Tarbase/TargetScan/TransmiR/tred -> Regdb
  unique(hm.rn$secondary.source[grep('mir|micro|tar|tred',hm.rn$secondary.source,ignore.case = T)])
  hm.rn$secondary.source[grep('mir|micro|tar|tred',hm.rn$secondary.source,ignore.case = T)] <- 'Regdb'
  ### HPRD(human) -> PPI
  unique(hm.rn$secondary.source[grep('hprd',hm.rn$secondary.source,ignore.case = T)])
  hm.rn$secondary.source[grep('hprd',hm.rn$secondary.source,ignore.case = T)] <- 'PPI'
  ### UCSC/Ensembl -> Annodb
  unique(hm.rn$secondary.source[grep('ucsc|ensembl',hm.rn$secondary.source,ignore.case = T)])
  hm.rn$secondary.source[grep('ucsc|ensembl',hm.rn$secondary.source,ignore.case = T)] <- 'Annodb'
  table(hm.rn$secondary.source)
  
  hm.rn$key <- paste(hm.rn$source, hm.rn$target, sep = "_")
  hm.rn <- hm.rn[!duplicated(hm.rn[,c('key','secondary.source')]),]
  hm.rn <- hm.rn[,1:5]
  
  filter_TFTG <- grepl('hsa-mir',hm.rn$source,ignore.case = T) | grepl('hsa-mir',hm.rn$target,ignore.case = T)
  table(filter_TFTG)
  hm.rn <- hm.rn[!filter_TFTG,]
  
  saveRDS(hm.rn, "./output/hm.TFTG.rn.rds")
  
}

hm.rn <- readRDS("./output/hm.TFTG.rn.rds")

## GTRD ####

if(F){
  
  # download.file(url = 'http://gtrd20-06.biouml.org/downloads/20.06/intervals/target_genes/Homo%20sapiens.tar.gz',
  #               destfile = './prior_knowledge/TFTG/GTRD/Homo sapiens.tar.gz', method = "internal")
  ### from website: http://gtrd.biouml.org/#!table/gtrd_current.tf_links/Brief
  
  hm.gd.tf <- read.csv("./TFTG/GTRD/human_GTRD_TFs.csv")
  hm.gd.tf <- hm.gd.tf %>% dplyr::select(symbol = Gene, GTRD = GTRD) %>%  as_tibble()
  
  hm.gd = lapply(seq(nrow(hm.gd.tf)), function(i){
    
    cat(paste0(i,"\n"))
    gtrd = hm.gd.tf$GTRD[i] %>% as.character()
    symbol = hm.gd.tf$symbol[i] %>% as.character()
    
    dl.url <- paste0('./TFTG/GTRD/Homo sapiens/transcripts promoter[-1000,+100]/',gtrd,'.txt')
    error_tag = file.exists(dl.url)
    if(error_tag){
      
      gtrd.tg <- read.table(dl.url, header = T, sep = "\t")
      gtrd.tg <- gtrd.tg %>% mutate(source = symbol) %>% dplyr::select(source, target = ID)
      
    }else{
      
      gtrd.tg <- data.frame(source = NA, target = NA)
      
    }
    gtrd.tg
    
  }) %>% do.call('rbind',.) %>% 
    na.omit() %>%  as_tibble()
  
  et2ent <- toTable(org.Hs.egENSEMBLTRANS)
  et2sy <- toTable(org.Hs.egSYMBOL)
  ent2sy <- merge(et2sy, et2ent, by = 'gene_id') %>% na.omit()
  saveRDS(ent2sy,'./TFTG/org.Hs.en2sy.3.10.rds')
  
  hm.gd$symbol <- ent2sy$symbol[match(hm.gd$target,ent2sy$trans_id)]
  hm.gd <- hm.gd %>% dplyr::select(source, target = symbol) %>% na.omit() %>% .[!duplicated(.),]
  table(hm.gd$source) %>% .[order(.)]
  
  hm.gd <- hm.gd %>% mutate(primary.source = rep("GTRD",nrow(hm.gd))) %>%
    mutate(secondary.source = rep("meta analysis",nrow(hm.gd)))
  
  table(hm.gd$source) %>% .[order(.)]
  hm.gd.tf[grep('^T$',hm.gd.tf$symbol),]
  hm.gd$source[grep('^T$',hm.gd$source)] <- 'TBXT'
  
  hm.gd$key <- paste(hm.gd$source, hm.gd$target, sep = "_")
  table(duplicated(hm.gd$key))
  hm.gd <- hm.gd[order(hm.gd$key),]
  hm.gd <- hm.gd[!duplicated(hm.gd[,1:2]),]
  hm.gd <- hm.gd[,-5]
  
  saveRDS(hm.gd, "./output/hm.TFTG.gd.rds")
  
}

hm.gd <- readRDS("./output/hm.TFTG.gd.rds")

## merge ####

filter_TFTG <- grepl('hsa-mir',hm.TFTG$source,ignore.case = T) | grepl('hsa-mir',hm.TFTG$target,ignore.case = T)
table(filter_TFTG)
hm.TFTG <- hm.TFTG[!filter_TFTG,]

save(hm.TFTG, file = "./output/scMLnet.human.TFTG.rda")

if(T){
  
  load("./output/scMLnet.human.TFTG.rda")
  
  # veen for source
  
  table(hm.TFTG$source)
  listinput <- list(GTRD = unique(hm.TFTG$source[hm.TFTG$primary.source == 'GTRD']),
                    HTRIdb = unique(hm.TFTG$source[hm.TFTG$primary.source == 'HTRIdb']),
                    RegNetwork = unique(hm.TFTG$source[hm.TFTG$primary.source == 'RegNetwork']),
                    TRRUST = unique(hm.TFTG$source[hm.TFTG$primary.source == 'TRRUST']))
  
  p <- venn.diagram(listinput, resolution = 1500, imagetype = "tiff",
                    main.cex = 2, main.pos = c(0.5, 1.05), cat.cex = 1.4, cex = 1.4,
                    fill=c('Tomato','SkyBlue','PaleGreen2','MediumPurple'), alpha=c(0.5,0.5,0.5,0.5),
                    main="Intersection of TFs \ncollected from different datasources",
                    filename = NULL)
  pdf("./figure/veen_TF_from_datasources.pdf",width = 7,height = 7)
  grid.draw(p)
  dev.off()
  
  # veen for target
  
  table(hm.TFTG$target)
  listinput <- list(GTRD = unique(hm.TFTG$target[hm.TFTG$primary.source == 'GTRD']),
                    HTRIdb = unique(hm.TFTG$target[hm.TFTG$primary.source == 'HTRIdb']),
                    RegNetwork = unique(hm.TFTG$target[hm.TFTG$primary.source == 'RegNetwork']),
                    TRRUST = unique(hm.TFTG$target[hm.TFTG$primary.source == 'TRRUST']))
  
  p <- venn.diagram(listinput, resolution = 1500, imagetype = "tiff",
                    main.cex = 2, main.pos = c(0.5, 1.05), cat.cex = 1.4, cex = 1.4,
                    fill=c('Tomato','SkyBlue','PaleGreen2','MediumPurple'), alpha=c(0.5,0.5,0.5,0.5),
                    main="Intersection of TGs \ncollected from different datasources",
                    filename = NULL)
  pdf("./figure/veen_TG_from_datasources.pdf",width = 7,height = 7)
  grid.draw(p)
  dev.off()
  
  # upset for intersection of primary.source and secondary.source
  
  hm.count.ps <- as.data.frame(table(hm.TFTG$primary.source, hm.TFTG$secondary.source))
  hm.count.ps$primary.source <- as.vector(hm.count.ps$Var1)
  hm.count.ps$secondary.source <- as.vector(hm.count.ps$Var2)
  hm.count.ps <- hm.count.ps[order(hm.count.ps$primary.source, hm.count.ps$Freq, decreasing = T),]
  hm.count.ps$ymax <- unlist(lapply(seq(nrow(hm.count.ps)), function(i){sum(hm.count.ps$Freq[seq(i)])}))
  hm.count.ps$ymin <- c(0,hm.count.ps$ymax[-nrow(hm.count.ps)])
  hm.count.ps <- hm.count.ps[,c(4,5,3,6,7)]
  mycols <- c('Tomato','SkyBlue','PaleGreen2','MediumPurple',
              brewer.pal(length(unique(c(hm.TFTG$secondary.source))), 'Set3'))
  names(mycols) <- unique(c(hm.TFTG$primary.source, hm.TFTG$secondary.source))
  
  pdf(file='./figure/doupie_TFTG_source.pdf', width = 9, height = 6)
  p <- ggplot(hm.count.ps) + 
    geom_rect(aes(fill = secondary.source, ymax = ymax, ymin = ymin, xmax = 4, xmin = 3)) +
    geom_rect(aes(fill = primary.source, ymax = ymax, ymin = ymin, xmax = 3, xmin = 0)) +
    scale_fill_manual(name="data.source", values = mycols) +
    xlim(c(0, 4)) + theme(aspect.ratio = 1) + 
    coord_polar(theta = 'y') + theme_bw() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.title = element_text(size = 12)
    )
  p
  dev.off()
  
  # node and edge
  
  hm.count <- hm.count.ps[,1:3]
  hm.count <- hm.count[!hm.count$Freq==0,]
  hm.count <- apply(hm.count, 1, function(hm.count.x){
    
    hm.TFTG.x <- hm.TFTG[hm.TFTG$primary.source == hm.count.x[1] & 
                               hm.TFTG$secondary.source == hm.count.x[2],]
    
    nr_nodes <- length(unique(c(hm.TFTG.x$source, hm.TFTG.x$target)))
    nr_source_nodes <- length(unique(hm.TFTG.x$source))
    nr_target_nodes <- length(unique(hm.TFTG.x$target))
    
    c(nr_nodes, nr_source_nodes, nr_target_nodes)
    
  }) %>% t() %>% cbind(hm.count,.)
  colnames(hm.count)[-c(1:2)] <- c('nr_edges', 'nr_nodes', 'nr_source_nodes', 'nr_target_nodes')
  write.csv(hm.count,"./output/hm.TFTG.detail.csv", row.names = F)
  
}

#############
## RecTFDB ##
#############
## graphite ####

if(F){
  
  pathwayDatabases()
  human.dbs <- pathwayDatabases() %>% 
    dplyr::filter(species == "hsapiens") %>% 
    dplyr::select(database) %>% unlist()

  options(Ncpus = 4)
  hm.dbs <- lapply(human.dbs, function(h.db){
    
    # h.db <- human.dbs[1]
    message(h.db)
    pws <- pathways(species = "hsapiens", database = h.db)
    pws.sym <- convertIdentifiers(pws, "SYMBOL")
    
    pws.dat <- lapply(pws.sym, function(pw.sym){
      # pw.sym <- pws.sym[[2]]
      pw.dat <- edges(pw.sym)
      if(nrow(pw.dat)){
        pw.dat <- pw.dat %>% mutate(numnodes = pw.dat %>% dplyr::select(src,dest) %>% unlist() %>% unique() %>% length(),
                                    numedges = nrow(pw.dat),
                                    id = pw.sym@id, 
                                    pathway = pw.sym@title, 
                                    database = pw.sym@database,
                                    species = pw.sym@species)
      }
      pw.dat
    })
    pws.dat <- do.call(rbind,pws.dat)
    pws.dat
    
  })
  hm.pw <- do.call(rbind,hm.dbs) %>% as_tibble()
  
  ## check
  
  table(hm.pw$src_type)
  table(hm.pw$dest_type)
  
  table(hm.pw$direction)
  table(hm.pw$type)
  table(hm.pw$type,hm.pw$direction)
  
  table(hm.pw$database)
  table(hm.pw$species)
  
  ## binding / undirected interaction according to detail
  
  add_pairs <- hm.pw[hm.pw$direction == "undirected" | 
                       hm.pw$type == "Process(binding/association)",] 
  add_pairs <- add_pairs[,c(3:4,1:2,5:12)]
  hm.pw <- rbind(hm.pw, add_pairs)
  
  ## database
  
  unique(hm.pw$database)
  ### BioCarta	KEGG	NCI	PANTHER	Pathbank	PharmGKB	Reactome	SMPDB
  
  ### biocarta -> BioCarta
  unique(hm.pw$database[grep('biocarta',hm.pw$database,ignore.case = T)])
  hm.pw$database[grep('biocarta',hm.pw$database,ignore.case = T)] <- 'BioCarta'
  
  ### nci -> NCI
  unique(hm.pw$database[grep('nci',hm.pw$database,ignore.case = T)])
  hm.pw$database[grep('nci',hm.pw$database,ignore.case = T)] <- 'NCI'
  
  ### panther -> PANTHER
  unique(hm.pw$database[grep('panther',hm.pw$database,ignore.case = T)])
  hm.pw$database[grep('panther',hm.pw$database,ignore.case = T)] <- 'PANTHER'
  
  ### Pathbank -> PathBank
  unique(hm.pw$database[grep('Pathbank',hm.pw$database,ignore.case = T)])
  hm.pw$database[grep('Pathbank',hm.pw$database,ignore.case = T)] <- 'PathBank'
  
  ## filter
  
  hm.pw <- hm.pw %>% mutate(primary.source = rep("graphite",nrow(hm.pw))) %>%
    filter(src_type != "CHEBI" & dest_type != "CHEBI") %>% 
    dplyr::select(source = src, target = dest, primary.source, secondary.source = database,
                  Ed.type = type, PW.names = pathway, PW.id = id)
  
  ## save
  
  save(hm.pw, file = "./output/scMLnet.human.pw.rda")
  
}


if(F){
  
  ## load
  
  load("./output/scMLnet.human.pw.rda")
  
  ## color
  
  mycols <- brewer.pal(length(unique(hm.pw$secondary.source)), 'Set2')
  names(mycols) <- unique(hm.pw$secondary.source)
  
  ## barplot for pathway in each database
  
  hm.pw.ct<- hm.pw %>% group_by(secondary.source, PW.names) %>% 
    summarise(., count = n()) %>% ungroup() %>% group_by(secondary.source) %>% 
    count(secondary.source) %>% ungroup() %>% arrange(desc(n))
  hm.pw.ct$secondary.source = factor(hm.pw.ct$secondary.source, levels = rev(hm.pw.ct$secondary.source))
  hm.pw.ct$log2n <- log2(hm.pw.ct$n)
  
  pdf(file='./figure/barplot_pathway_source.pdf', width = 7, height = 6)
  png(file='./figure/barplot_pathway_source.png', 
      width = 7, height = 6, units = 'in', res = 300)
  p1 <- ggplot(data = hm.pw.ct, mapping=aes(x=secondary.source, y=n, fill=secondary.source))+
    geom_bar(stat="identity",position=position_dodge()) +
    scale_fill_manual(name="data.source", values = mycols) +
    labs(x=NULL,y='number of pathways',fill=NULL) +   
    coord_cartesian(ylim = c(0,2000)) + 
    theme_classic() + 
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 14),
      axis.text.y = element_text(size = 12),
      axis.text.x = element_text(size = 10),
      legend.position = 'none'
    )
  p2 <- ggplot(data = hm.pw.ct, mapping=aes(x=secondary.source, y=n, fill=secondary.source))+
    geom_bar(stat="identity",position=position_dodge()) +
    scale_fill_manual(name="data.source", values = mycols) +
    labs(x=NULL,y=NULL,fill=NULL) +   
    coord_cartesian(ylim = c(48500,48600)) +  
    scale_y_continuous(breaks = c(48500,48600, 10)) +
    theme_classic() + 
    theme(
      legend.position = 'none',
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 12),
      axis.line.x.bottom = element_line(size=0, colour = 'white')
    )
  ggpubr::ggarrange(p2,p1,heights=c(1/4, 3/4),ncol = 1, nrow = 2, align = "v") 
  dev.off()
  
}

## detail ####

if(F){
  
  ### load interaction
  
  load("./output/scMLnet.human.pw.rda")
  load("./output/scMLnet.human.LigRec.rda")
  load("./output/scMLnet.human.TFTG.rda")
  
  Ligands <- hm.LigRec$source %>% unique() %>% as.character()
  Receptors <- hm.LigRec$target %>% unique() %>% as.character()
  TFs <- hm.TFTG$source %>% unique() %>% as.character()
  TGs <- hm.TFTG$target %>% unique() %>% as.character()
  
  ## get hm.db.count
  
  hm.db.count <- as.data.frame(table(hm.pw$primary.source,hm.pw$secondary.source))
  hm.db.count <- apply(hm.db.count, 1, function(hm.db.count.x){
    
    hm.pw.x <- hm.pw[hm.pw$primary.source == hm.db.count.x[1] & 
                       hm.pw$secondary.source == hm.db.count.x[2],]
    
    nr_nodes <- length(unique(c(hm.pw.x$source, hm.pw.x$target)))
    nr_source_nodes <- length(unique(hm.pw.x$source))
    nr_target_nodes <- length(unique(hm.pw.x$target))
    nr_pws <- length(unique(hm.pw.x$PW.names))
    
    c(nr_nodes, nr_source_nodes, nr_target_nodes, nr_pws)
    
  }) %>% t() %>% cbind(hm.db.count,.)
  colnames(hm.db.count) <- c('primary.source', 'secondary.source', 'nr_edges', 
                             'nr_nodes', 'nr_source_nodes', 'nr_target_nodes', 'nr_pws')
  
  hm.db.count <- lapply(seq(nrow(hm.db.count)), function(x){
    
    hm.pw.db.x <- hm.pw %>% filter(secondary.source == hm.db.count[x,2])
    nodes <- unique(c(hm.pw.db.x$source, hm.pw.db.x$target))
    nr_Ligands <- length(intersect(Ligands,nodes))
    nr_Receptors <- length(intersect(Receptors,nodes))
    nr_TFs <- length(intersect(TFs,nodes))
    
    c(nr_Ligands,nr_Receptors,nr_TFs)
    
  }) %>% 
    do.call('rbind',.) %>% cbind(hm.db.count,.)
  colnames(hm.db.count)[8:10] <- c('nr_Ligands','nr_Receptors','nr_TFs')
  write.csv(hm.db.count, "./output/hm.RecTF.db.detail.csv", row.names = F)
  
}