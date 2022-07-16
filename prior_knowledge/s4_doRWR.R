#############
## library ##
#############

library(igraph)
library(dplyr)
library(parallel)
library(foreach)
library(doParallel)
library(doSNOW)

rm(list = ls())
gc()

setwd("E:/stMLnet/prior_knowledge/")

##############
## function ##
##############

## RWR算法
doRWR <- function(mygraph, restart=0.75, 
                  mode = c('directed','undirected'), 
                  normalise.affM = c("none","quantile"),
                  setSeeds=NULL, SeedInRow = TRUE, verbose=T)
{
  
  ## function
  ####################################################
  
  ## A function to make sure the sum of elements in each steady probability vector is one
  sum2one <- function(PTmatrix){
    col_sum <- apply(PTmatrix, 2, sum)
    #col_sum <- Matrix::colSums(PTmatrix, sparseResult=F)
    col_sum_matrix <- matrix(rep(col_sum, nrow(PTmatrix)), ncol=ncol(PTmatrix), nrow=nrow(PTmatrix), byrow =T)
    res <- as.matrix(PTmatrix)/col_sum_matrix
    res[is.na(res)] <- 0
    return(res)
  }
  
  ## A function to normalize columns of a matrix to have the same Quantiles
  normalizeQuantiles <- function (A, ties=TRUE) {
    n <- dim(A)
    if(is.null(n)) return(A)
    if(n[2] == 1) return(A)
    O <- S <- array(, n)
    nobs <- rep(n[1], n[2])
    i <- (0:(n[1] - 1))/(n[1] - 1)
    for(j in 1:n[2]){
      Si <- sort(A[, j], method = "quick", index.return = TRUE)
      nobsj <- length(Si$x)
      if (nobsj < n[1]){
        nobs[j] <- nobsj
        isna <- is.na(A[, j])
        S[, j] <- stats::approx((0:(nobsj - 1))/(nobsj - 1), Si$x, 
                                i, ties = "ordered")$y
        O[!isna, j] <- ((1:n[1])[!isna])[Si$ix]
      }else{
        S[, j] <- Si$x
        O[, j] <- Si$ix
      }
    }
    m <- rowMeans(S)
    for (j in 1:n[2]){
      if(ties) r<-rank(A[, j])
      if(nobs[j] < n[1]){
        isna <- is.na(A[, j])
        if(ties){ 
          A[!isna, j] <- stats::approx(i, m, (r[!isna] - 1)/(nobs[j] - 1), ties = "ordered")$y
        }else{ 
          A[O[!isna, j], j] <- stats::approx(i, m, (0:(nobs[j] - 1))/(nobs[j] - 1), ties = "ordered")$y
        }
      }else{
        if(ties){
          A[, j] <- stats::approx(i, m, (r - 1)/(n[1] - 1), ties = "ordered")$y
        }else{
          A[O[, j], j] <- m
        }
      }
    }
    
    return(A)
  }
  
  ## A function to perform RWR algorithm on each seed
  runRWR = function(Seed, G, restart) {
    
    # prepare preference vector P(the initial probability distribution), only the seeds have values different from zero in P.
    P = rep(0, times = length(igraph::V(G)))
    P[match(Seed,igraph::V(G)$name)] = 1
    
    PTvector = igraph::page_rank(G, algo = c("prpack"), vids = igraph::V(G)$name,
                                 directed = ifelse(mode=="directed",TRUE,FALSE), weights = NULL,
                                 damping = restart, personalized = P) %>% .$vector
    
    return(PTvector)
    
  }
  
  ####################################################
  
  ## start
  startT <- Sys.time()
  if(verbose){
    message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
    message("", appendLF=T)
  }
  
  ## check
  mode <- match.arg(mode)
  
  normalise.affM <- match.arg(normalise.affM)
  
  if(is.null(restart) || is.na(restart) || restart<0 || restart>100){
    r <- 0.75
  }else if(restart>1 && restart<100){
    r <- restart/100
  }else{
    r <- restart
  }
  
  if (class(mygraph) != "igraph"){
    stop("The function must apply to either 'igraph' object.\n")
  }
  
  if(is.null(setSeeds)){
    Seeds <- V(mygraph)$name ## 每组seeds包括至少一个节点, 默认所有节点
  }else{
    
    ## check mapping between input and graph
    ind <- match(setSeeds, V(mygraph)$name)
    if(length(ind)==0){
      stop("All setSeeds do not included in the input graph.\n")
    }else if(length(ind)!=length(setSeeds)){
      Seeds <- V(mygraph)$name[ind]
      warning("At least one input of setSeeds do not included in the input graph.\n")
    }else{
      Seeds <- setSeeds
    }
    
  }
  
  ## main
  
  if(verbose){
    now <- Sys.time()
    message(sprintf("First, RWR using %1.2f restart probability (%s) ...", r, as.character(now)), appendLF=T)
  }
  
  PTmatrix = lapply(Seeds, runRWR, G=mygraph, restart=r) %>% do.call("cbind",.)
  
  if(verbose){
    now <- Sys.time()
    message(sprintf("Second, rescale steady probability vector (%s) ...", as.character(now)), appendLF=T)
  }
  
  PTmatrix <- sum2one(PTmatrix) # input/output: full matrix
  # make sure the sum of elements in each steady probability vector is one
  
  if(ncol(PTmatrix) == 1){
    normalise.affM <- "none"
  }
  if(normalise.affM=="quantile"){
    PTmatrix <- normalizeQuantiles(PTmatrix)
    if(verbose){
      now <- Sys.time()
      message(sprintf("Finally, normalise the affinity matrix using %s normalisation (%s) ...", normalise.affM, as.character(now)), appendLF=T)
      
    }
  }
  
  ## should keep seed in row（row->col） 
  if(SeedInRow){
    
    PTmatrix <- t(PTmatrix)
    colnames(PTmatrix) <- V(mygraph)$name 
    rownames(PTmatrix) <- Seeds
    PTmatrix <- Matrix::Matrix(PTmatrix, sparse=T)
    
  }else{
    
    colnames(PTmatrix) <- Seeds
    rownames(PTmatrix) <- V(mygraph)$name 
    PTmatrix <- Matrix::Matrix(PTmatrix, sparse=T)
    
  }
  
  ## end
  
  endT <- Sys.time()
  if(verbose){
    message(paste(c("\nFinish at ",as.character(endT)), collapse=""), appendLF=T)
  }
  
  runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
  if(verbose){
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=T)
  }
  
  invisible(PTmatrix)
}

## 划分数据集
splitData <- function(data, k=3){
  
  getCVgroup <- function(data,k=3){
    cvlist <- list()
    datasize <- nrow(data)
    n <- rep(1:k,ceiling(datasize/k))[1:datasize]    #将数据分成K份，并生成的完成数据集n
    temp <- sample(n,datasize)   #把n打乱
    dataseq <- 1:datasize
    cvlist <- lapply(1:k,function(x) dataseq[temp==x])  #dataseq中随机生成k个随机有序数据列
    return(cvlist)
  }
  cvlist <- getCVgroup(data)
  
  datas <- lapply(1:k, function(i){
    
    index <- cvlist[[i]]
    train <- data[-index, ]
    rownames(train) <- 1:nrow(train)
    test <- data[index, ]
    rownames(test) <- 1:nrow(test)
    
    list(train=train,
         test=test)
    
  })
  return(datas)
  
}

## 构建图
buildNet <- function(datas, mode = c('directed','undirected'), 
                     weight = c('weighted','unweighted'), 
                     train.model = FALSE, verbose=TRUE){	
  
  mode <- match.arg(mode)
  weight <- match.arg(weight)
  
  if(verbose){
    message(sprintf("\tConstructing a %s %s graph!", mode, weight), appendLF=T)
  }
  
  if(train.model){
    
    V <- c(datas$train$source,datas$train$target,datas$test$source,datas$test$target) %>% unique() %>% sort()
    edge <- datas$train
    
  }else{
    
    V <- c(datas$source,datas$target) %>% unique() %>% sort()
    edge <- datas
    
  }
  
  edge$from <- match(edge$source,V)
  edge$to <- match(edge$target,V)
  
  adjM <- matrix(0,nrow = length(V),ncol = length(V), dimnames = list(V,V))
  if(weight == "unweighted"){
    
    for(i in 1:nrow(edge)){
      
      if(mode=="directed"){
        
        adjM[edge$from[i],edge$to[i]] <- 1
        
      }else{
        
        adjM[edge$from[i],edge$to[i]] <- 1
        adjM[edge$to[i],edge$from[i]] <- 1
        
      }
    }
    adjM <- as(adjM, "dgCMatrix")
    
  }else{
    
    # check weight
    if(!'weight' %in% colnames(edge)){
      
      edge$weight = 1
      warning("data do not have a weight column, all the edges have the same weight")
    }
    
    for(i in 1:nrow(edge)){
      
      if(mode=="directed"){
        
        adjM[edge$from[i],edge$to[i]] <- edge$weight[i]
        
      }else{
        
        adjM[edge$from[i],edge$to[i]] <- edge$weight[i]
        adjM[edge$to[i],edge$from[i]] <- edge$weight[i]
        
      }
    }
    adjM <- as(adjM, "dgCMatrix")
    
  }
  
  if(weight == "unweighted"){
    
    # directed: The graph will be directed and a matrix element gives the number of edges between two vertices.
    # undirected: This is exactly the same as max, for convenience. Note that it is not checked whether the matrix is symmetric.
    net <- igraph::graph_from_adjacency_matrix(adjM, mode=mode, weighted=NULL)
    
  }else{
    
    # directed: The graph will be directed and a matrix element gives the edge weights.
    # undirected: First we check that the matrix is symmetric. It is an error if not. Then only the upper triangle is used to create a weighted undirected graph.
    net <- igraph::graph_from_adjacency_matrix(adjM, mode=mode, weighted=TRUE)
    
  }
  
  # delete loop and multiple edges
  net <- igraph::simplify(net,remove.multiple = TRUE, remove.loops = TRUE)
  
  # delete isolated nodes
  isolates <- which(igraph::degree(net, mode = c("all")) == 0) - 1
  net <- delete.vertices(net, names(isolates))
  
  return (net) 	
  
}


##########
## load ##
##########

## load LigRec.DB

load("./output/scMLnet.human.LigRec.rda")
LigRec.DB <- hm.LigRec %>%
  .[!duplicated(.[,1:4]),] %>%
  group_by(source, target) %>% 
  summarise(score = n()) %>% 
  ungroup()
head(LigRec.DB)

## load TFTG.DB

load("./output/scMLnet.human.TFTG.rda")
TFTG.DB <- hm.TFTG %>%
  .[!duplicated(.[,c(1,2,4)]),] %>%
  group_by(source, target) %>% 
  summarise(score = n()) %>% 
  ungroup()
head(TFTG.DB)

## load pw.DB

load("./output/scMLnet.human.pw.rda")
pw.DB <- hm.pw %>% 
  .[!duplicated(.[,1:4]),] %>%
  group_by(source, target) %>% 
  summarise(weight = n()) %>% 
  ungroup() %>% 
  as.data.frame()
head(pw.DB)

## filter Rec,TF

LigRec.DB <- LigRec.DB[LigRec.DB$target %in% c(pw.DB$source, pw.DB$target),]
TFTG.DB <- TFTG.DB[TFTG.DB$source %in% c(pw.DB$source, pw.DB$target),]

##################
## get RecTF.DB ##
##################

## input

r = 0.61 # s3_tunepara_database.R

Receptors <- LigRec.DB %>% dplyr::select(target) %>% unlist() %>% unique() 
TFs <- TFTG.DB %>% dplyr::select(source) %>% unlist() %>% unique()

## doRWR

net <- buildNet(pw.DB, mode = 'directed', weight = 'weighted')
RecSiganltab <- doRWR(net, restart=r, setSeeds=Receptors, verbose=T, SeedInRow = T)
saveRDS(RecSiganltab,"./output/RecSiganltab.rds")

## filter

# keep TF in columns
RecTFtab <- RecSiganltab[,colnames(RecSiganltab) %in% TFs]

# delect Receptor without downstream TFs
filter_Recs <- Matrix::rowSums(RecTFtab) %>% .[.==0] %>% names()
RecTFtab <- RecTFtab[!rownames(RecTFtab) %in% filter_Recs,]
fivenum(RecTFtab)

# delect Receptor without downstream TFs
LigRec.DB <- LigRec.DB[!LigRec.DB$target %in% filter_Recs,]

###########
## clean ##
###########

## save

RecTF.DB <- RecTFtab %>% as.matrix() %>% reshape2::melt() %>% 
  dplyr::select(source=Var1, target=Var2, score=value) %>%
  as_tibble()
RecTF.DB$source <- as.character(RecTF.DB$source)
RecTF.DB$target <- as.character(RecTF.DB$target)

Databases <- list(LigRec.DB = LigRec.DB,
                  RecTF.DB = RecTF.DB,
                  TFTG.DB = TFTG.DB)

saveRDS(Databases, "./output/Databases.rds")



