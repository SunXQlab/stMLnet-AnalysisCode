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

setwd("./stMLnet/prior_knowledge/")

##############
## function ##
##############

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
    Seeds <- V(mygraph)$name 
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

splitData <- function(data, k=3){
  
  getCVgroup <- function(data,k=3){
    cvlist <- list()
    datasize <- nrow(data)
    n <- rep(1:k,ceiling(datasize/k))[1:datasize]   
    temp <- sample(n,datasize)
    dataseq <- 1:datasize
    cvlist <- lapply(1:k,function(x) dataseq[temp==x])  
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
    
     net <- igraph::graph_from_adjacency_matrix(adjM, mode=mode, weighted=NULL)
    
  }else{
    
    net <- igraph::graph_from_adjacency_matrix(adjM, mode=mode, weighted=TRUE)
    
  }
  
  # delete loop and multiple edges
  net <- igraph::simplify(net,remove.multiple = TRUE, remove.loops = TRUE)
  
  # delete isolated nodes
  isolates <- which(degree(net, mode = c("all")) == 0) - 1
  net <- delete.vertices(net, names(isolates))
  
  return (net) 	
  
}

###############
## parameter ##
###############

gs1 <- data.frame(r = sample(seq(0.01,0.99,0.01),n.iter,replace = F))
saveRDS(gs1, "./RWR/gs1.rds")

gs2 <- data.frame(ds = rep(1:3,each=30), r = rep(gs1$r,time3=3))
gs2 <- gs2[order(gs2$r,decreasing = T),]
rownames(gs2) <- 1:nrow(gs2)
saveRDS(gs2, "./RWR/gs2.rds")

#################
## run stMLnet ##
#################

## load stMLnet database

load(file = "./output/scMLnet.human.pw.rda")
RecTF_from_stMLnet = hm.pw %>% 
  .[!duplicated(.[,1:4]),] %>%
  group_by(source, target) %>% 
  summarise(weight = n()) %>% 
  ungroup()
data_from_stMLnet <- RecTF_from_stMLnet %>% as.data.frame()

## parameter

n.iter <- 30
gs2 <- readRDS("./RWR/gs2.rds")

## support function

mod_new <- function(i) {
  
  r=gs2$r[i]
  ds=gs2$ds[i]
  
  g <- netlist[[ds]]
  AffM <- doRWR(g, restart=r, verbose=F)
  
  results <- list(dataset = ds, 
                  restart_probability = r,
                  AffM = AffM)
  
  return(results)
  
}
eva_new1 <- function(i){
  
  s = (i-1)*3+1
  e = i*3
  AffMlist <- res_stMLnet_RWR[s:e]
  r <- AffMlist[[1]]$restart_probability
  
  AUClist <- lapply(1:length(AffMlist), function(i){
    evaluateAUC(AffMlist[[i]]$AffM, datalist[[i]])
  }) %>% unlist()
  AUClist <- list(mean_AUC = mean(AUClist), AUC = AUClist)
  
  results <- list(Parameter = list(restart_probability = r),
                  AUC = AUClist
  )
  
  return(results)
  
} 
eva_new2 <- function(i){
  
  s = (i-1)*3+1
  e = i*3
  AffMlist <- res_stMLnet_RWR[s:e]
  r <- AffMlist[[1]]$restart_probability
  
  Precisionlist <- lapply(1:length(AffMlist), function(i){
    evaluatePrecision(AffMlist[[i]]$AffM, datalist[[i]])
  }) %>% unlist()
  Precisionlist <- list(mean_Precision = mean(Precisionlist), Precision = Precisionlist)
  
  results <- list(Parameter = list(restart_probability = r),
                  Precision = Precisionlist
  )
  
  return(results)
  
} 

## split data 

datalist <- splitData(data_from_stMLnet)

## built net

netlist <- lapply(datalist, buildNet, mode = "directed", weight = "weighted", train.model = TRUE)

## do rwr

parallel::detectCores()
t1 <- Sys.time()
message("Start at ",as.character(t1))
cl <- makeSOCKcluster(30)
registerDoSNOW(cl)
pb <- txtProgressBar(min=1, max=90, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)
res_stMLnet_RWR <- foreach(i=1:90, .packages=c("dplyr","igraph"), 
                           .options.snow=opts, .errorhandling = 'pass'
) %dopar% {
  res <- mod_new(i)
  gc()
  res 
} 
close(pb)
stopCluster(cl)
t2 <- Sys.time()
message("End at ",as.character(t2))
t2-t1 # 1.8 hour
gc()

## get AUC

num <- length(res_stMLnet_RWR)/3
parallel::detectCores()
cl <- makeSOCKcluster(15)
t1 <- Sys.time()
message("Start at ",as.character(t1))
pb <- txtProgressBar(min=1, max=num, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)
res_stMLnet_RWR_1 <- foreach(i=1:num, 
                             .packages=c("dplyr","Matrix"), 
                             .options.snow=opts, 
                             .errorhandling = 'pass'
) %do% {
  res <- eva_new1(i)
  gc()
  res
} # 
close(pb)
stopCluster(cl)
t2 <- Sys.time()
message("End at ",as.character(t2))
t2-t1 # 23 mins
saveRDS(res_stMLnet_RWR_1, "./RWR/res_stMLnet_RWR_1.rds")

## get precision

num <- length(res_stMLnet_RWR)/3
parallel::detectCores()
t1 <- Sys.time()
cl <- makeSOCKcluster(10)
registerDoSNOW(cl)
pb <- txtProgressBar(min=1, max=num, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)
res_stMLnet_RWR_2 <- foreach(i=1:num, 
                             .packages=c("dplyr"), 
                             .options.snow=opts,
                             .errorhandling = 'pass'
) %do% {
  res <- eva_new2(i)
  gc()
  res
} 
close(pb)
stopCluster(cl)
t2 <- Sys.time()
t2-t1 # 2.9 hours
saveRDS(res_stMLnet_RWR_2, "./RWR/res_stMLnet_RWR_2.rds")
gc()

## clean

res_stMLnet_RWR_1 <- readRDS("./RWR/res_stMLnet_RWR_1.rds")
res_stMLnet_RWR_2 <- readRDS("./RWR/res_stMLnet_RWR_2.rds")

df_metrics_st <- lapply(1:length(res_stMLnet_RWR_1), function(i){
  
  r <- res_stMLnet_RWR_1[[i]]$Parameter$restart_probability
  AUC <- res_stMLnet_RWR_1[[i]]$AUC$mean_AUC
  Precision <- res_stMLnet_RWR_2[[i]]$Precision$mean_Precision
  
  c(r,AUC,Precision)
  
}) %>% do.call("rbind",.) %>% as.data.frame()
colnames(df_metrics_st) <- c('r','mean_AUC','mean_Precision')

saveRDS(df_metrics_st, "./RWR/df_metrics_st.rds")

##################
## run Omnipath ##
##################

## load Omnipath database

hm.omni <- readRDS(file = "../other_method/Omnipath/cleaned_data/RecTF_from_interactions.rds")
RecTF_from_Omnipath <- hm.omni %>% filter(weight > 0)
data_from_Omnipath <- RecTF_from_Omnipath %>% as.data.frame()

## parameter

n.iter <- 30
gs1 <- readRDS("./RWR/gs1.rds")

## support function

mod <- function(i) {
  
  r=gs1$r[i]
  
  AffMlist <- lapply(netlist, function(g){
    doRWR(g, restart=r, setSeeds=NULL, verbose=F)
  })
  
  AUClist <- lapply(1:length(AffMlist), function(i){
    evaluateAUC(AffMlist[[i]], datalist[[i]])
  }) %>% unlist()
  AUClist <- list(mean_AUC = mean(AUClist), AUC = AUClist)
  
  Precisionlist <- lapply(1:length(AffMlist), function(i){
    evaluatePrecision(AffMlist[[i]], datalist[[i]])
  }) %>% unlist()
  Precisionlist <- list(mean_Precision = mean(Precisionlist), Precision = Precisionlist)
  
  results <- list(Parameter = list(restart_probability = r),
                  AUC = AUClist,
                  Precision = Precisionlist)
  
  return(results)
  
}

## split data 

datalist <- splitData(data_from_Omnipath)

## built net

netlist <- lapply(datalist, buildNet, mode = "directed", weight = "weighted", train.model = TRUE)

## do rwr + get metrics

parallel::detectCores()
t1 <- Sys.time()
message("Start at ",as.character(t1))
cl <- makeSOCKcluster(30)
registerDoSNOW(cl)
pb <- txtProgressBar(min=1, max=n.iter, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)
res_Omnipath_RWR <- foreach(i=1:n.iter, .packages=c("dplyr","igraph"), 
                            .options.snow=opts, .errorhandling = "pass"
) %dopar% mod(i) 
close(pb)
stopCluster(cl)
t2 <- Sys.time()
message("End at ",as.character(t2))
t2-t1

## save 

saveRDS(res_Omnipath_RWR, "./RWR/res_Omnipath_RWR.rds")

## clean

df_metrics_om <- lapply(1:length(res_Omnipath_RWR), function(i){
  
  r <- res_Omnipath_RWR[[i]]$Parameter$restart_probability
  AUC <- res_Omnipath_RWR[[i]]$AUC$mean_AUC
  Precision <- res_Omnipath_RWR[[i]]$Precision$mean_Precision
  
  c(r,AUC,Precision)
  
}) %>% do.call("rbind",.) %>% as.data.frame()
colnames(df_metrics_om) <- c('r','mean_AUC','mean_Precision')

saveRDS(df_metrics_om, "./RWR/df_metrics_om.rds")

##################
## run NicheNet ##
##################

## load NicheNet database

load("../other_method/NicheNet/data/nichenet.human.ppi.rda")
hm.ppi.Unique <- hm.ppi %>% dplyr::distinct(source, target)

RecTF_from_NicheNet <- hm.ppi
RecTF_from_NicheNet = RecTF_from_NicheNet %>% 
  group_by(source, target) %>% 
  summarise(weight = n()) %>% 
  ungroup()
data_from_NicheNet <- RecTF_from_NicheNet %>% as.data.frame()#[1:500,]

## parameter

n.iter <- 30
gs2 <- readRDS("./RWR/gs2.rds")

## support function

mod_new <- function(i) {
  
  r=gs2$r[i]
  ds=gs2$ds[i]
  
  g <- netlist[[ds]]
  AffM <- doRWR(g, restart=r, setSeeds=NULL, verbose=F, SeedInRow = FALSE)
  
  results <- list(dataset = ds, 
                  restart_probability = r,
                  AffM = AffM)
  
  return(results)
  
}
eva_new1 <- function(i){
  
  s = (i-1)*3+1
  e = i*3
  AffMlist <- res_NicheNet_RWR[s:e]
  r <- AffMlist[[1]]$restart_probability
  
  AUClist <- lapply(1:length(AffMlist), function(i){
    evaluateAUC(AffMlist[[i]]$AffM, datalist[[i]], FromInRow = FALSE)
  }) %>% unlist()
  AUClist <- list(mean_AUC = mean(AUClist), AUC = AUClist)
  
  results <- list(Parameter = list(restart_probability = r),
                  AUC = AUClist
  )
  
  return(results)
  
} 
eva_new2 <- function(i){
  
  s = (i-1)*3+1
  e = i*3
  AffMlist <- res_NicheNet_RWR[s:e]
  r <- AffMlist[[1]]$restart_probability
  
  Precisionlist <- lapply(1:length(AffMlist), function(i){
    evaluatePrecision(AffMlist[[i]]$AffM, datalist[[i]], FromInRow = FALSE)
  }) %>% unlist()
  Precisionlist <- list(mean_Precision = mean(Precisionlist), Precision = Precisionlist)
  
  results <- list(Parameter = list(restart_probability = r),
                  Precision = Precisionlist
  )
  
  return(results)
  
} 

## split data 

datalist <- splitData(data_from_NicheNet)

## bulit net

netlist <- lapply(datalist, buildNet, mode = "directed", weight = "weighted", train.model = TRUE)

## part1

if(T){
  
  ## do rwr
  
  parallel::detectCores()
  t1 <- Sys.time()
  message('Start at ',as.character(t1))
  cl <- makeSOCKcluster(63)
  registerDoSNOW(cl)
  pb <- txtProgressBar(min=1, max=63, style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  res_NicheNet_RWR <- foreach(i=10:72, .packages=c("dplyr","igraph"), 
                              .options.snow=opts, .errorhandling = 'pass'
  ) %dopar% {
    res <- mod_new(i)
    gc()
    res
  } # 10.36 hour
  close(pb)
  stopCluster(cl)
  t2 <- Sys.time()
  message('End at ',as.character(t2))
  t2-t1
  
  ## get AUC
  
  num <- length(res_NicheNet_RWR)/3
  parallel::detectCores()
  cl <- makeSOCKcluster(15)
  t1 <- Sys.time()
  message("Start at ",as.character(t1))
  pb <- txtProgressBar(min=1, max=num, style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  res_NicheNet_RWR_1 <- foreach(i=1:num, 
                                .packages=c("dplyr","Matrix"), 
                                .options.snow=opts, 
                                .errorhandling = 'pass'
  ) %do% {
    res <- eva_new1(i)
    gc()
    res
  } 
  close(pb)
  stopCluster(cl)
  t2 <- Sys.time()
  message("End at ",as.character(t2))
  t2-t1 
  
  saveRDS(res_NicheNet_RWR_1, "./RWR/res_NicheNet_RWR_1.rds")
  
  ## get precision
  
  num <- length(res_NicheNet_RWR)/3
  parallel::detectCores()
  t1 <- Sys.time()
  message("Start at ",as.character(t1))
  cl <- makeSOCKcluster(2)
  registerDoSNOW(cl)
  pb <- txtProgressBar(min=1, max=num, style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  res_NicheNet_RWR_2 <- foreach(i=1:num, 
                                .packages=c("dplyr"), 
                                .options.snow=opts,
                                .errorhandling = 'pass'
  ) %do% {
    res <- eva_new2(i)
    gc()
    res
  } # 
  close(pb)
  stopCluster(cl)
  t2 <- Sys.time()
  message("End at ",as.character(t2))
  t2-t1 
  
  saveRDS(res_NicheNet_RWR_2, "./RWR/res_NicheNet_RWR_2.rds") 
  gc()
  
}

## part2

if(T){
  
  ## support function
  eva_new1 <- function(i){
    
    s = (i-1)*3+1
    e = i*3
    AffMlist <- res_NicheNet_RWR2[s:e]
    r <- AffMlist[[1]]$restart_probability
    
    AUClist <- lapply(1:length(AffMlist), function(i){
      evaluateAUC(AffMlist[[i]]$AffM, datalist[[i]], FromInRow = FALSE)
    }) %>% unlist()
    AUClist <- list(mean_AUC = mean(AUClist), AUC = AUClist)
    
    results <- list(Parameter = list(restart_probability = r),
                    AUC = AUClist
    )
    
    return(results)
    
  } 
  eva_new2 <- function(i){
    
    s = (i-1)*3+1
    e = i*3
    AffMlist <- res_NicheNet_RWR2[s:e]
    r <- AffMlist[[1]]$restart_probability
    
    Precisionlist <- lapply(1:length(AffMlist), function(i){
      evaluatePrecision(AffMlist[[i]]$AffM, datalist[[i]], FromInRow = FALSE)
    }) %>% unlist()
    Precisionlist <- list(mean_Precision = mean(Precisionlist), Precision = Precisionlist)
    
    results <- list(Parameter = list(restart_probability = r),
                    Precision = Precisionlist
    )
    
    return(results)
    
  } 
  
  ## do rwr
  
  parallel::detectCores()
  t1 <- Sys.time()
  message("Start at ",as.character(t1))
  cl <- makeSOCKcluster(27)
  registerDoSNOW(cl)
  pb <- txtProgressBar(min=1, max=27, style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  res_NicheNet_RWR2 <- foreach(i=c(1:9,73:90), .packages=c("dplyr","igraph"), 
                               .options.snow=opts, .errorhandling = 'pass'
  ) %dopar% {
    res <- mod_new(i)
    gc()
    res
  }
  close(pb)
  stopCluster(cl)
  gc()
  t2 <- Sys.time()
  message("End at ",as.character(t2))
  t2-t1
  
  ## get AUC
  
  num <- length(res_NicheNet_RWR2)/3
  parallel::detectCores()
  cl <- makeSOCKcluster(9)
  t1 <- Sys.time()
  message("Start at ",as.character(t1))
  pb <- txtProgressBar(min=1, max=num, style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  res_NicheNet_RWR2_1 <- foreach(i=1:num, 
                                 .packages=c("dplyr","Matrix"), 
                                 .options.snow=opts, 
                                 .errorhandling = 'pass'
  ) %do% {
    res <- eva_new1(i)
    gc()
    res
  } 
  close(pb)
  stopCluster(cl)
  gc()
  t2 <- Sys.time()
  message("End at ",as.character(t2))
  t2-t1
  
  saveRDS(res_NicheNet_RWR2_1, "./RWR/res_NicheNet_RWR2_1.rds")
  
  ## get precision
  
  num <- length(res_NicheNet_RWR2)/3
  parallel::detectCores()
  t1 <- Sys.time()
  message("Start at ",as.character(t1))
  cl <- makeSOCKcluster(2)
  registerDoSNOW(cl)
  pb <- txtProgressBar(min=1, max=num, style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  res_NicheNet_RWR2_2 <- foreach(i=1:num, 
                                 .packages=c("dplyr"), 
                                 .options.snow=opts,
                                 .errorhandling = 'pass'
  ) %do% {
    res <- eva_new2(i)
    gc()
    res
  } # 
  close(pb)
  stopCluster(cl)
  gc()
  t2 <- Sys.time()
  message("End at ",as.character(t2))
  t2-t1 
  
  saveRDS(res_NicheNet_RWR2_2, "./RWR/res_NicheNet_RWR2_2.rds") 
  gc()
  
}

## clean

res_NicheNet_RWR_1 <- readRDS("./RWR/res_NicheNet_RWR_1.rds")
res_NicheNet_RWR_2 <- readRDS("./RWR/res_NicheNet_RWR_2.rds")
res_NicheNet_RWR2_1 <- readRDS("./RWR/res_NicheNet_RWR2_1.rds")
res_NicheNet_RWR2_2 <- readRDS("./RWR/res_NicheNet_RWR2_2.rds")

df_metrics_nn1 <- lapply(1:length(res_NicheNet_RWR_1), function(i){
  
  r <- res_NicheNet_RWR_1[[i]]$Parameter$restart_probability
  AUC <- res_NicheNet_RWR_1[[i]]$AUC$mean_AUC
  Precision <- res_NicheNet_RWR_2[[i]]$Precision$mean_Precision
  
  c(r,AUC,Precision)
  
}) %>% do.call("rbind",.) %>% as.data.frame()
df_metrics_nn2 <- lapply(1:length(res_NicheNet_RWR2_1), function(i){
  
  r <- res_NicheNet_RWR2_1[[i]]$Parameter$restart_probability
  AUC <- res_NicheNet_RWR2_1[[i]]$AUC$mean_AUC
  Precision <- res_NicheNet_RWR2_2[[i]]$Precision$mean_Precision
  
  c(r,AUC,Precision)
  
}) %>% do.call("rbind",.) %>% as.data.frame()
df_metrics_nn <- rbind(df_metrics_nn1,df_metrics_nn2) %>% as.data.frame()
colnames(df_metrics_nn) <- c('r','mean_AUC','mean_Precision')
df_metrics_nn <- df_metrics_nn[order(df_metrics_nn$mean_AUC,decreasing = T),]
head(df_metrics_nn)

saveRDS(df_metrics_nn, "./RWR/df_metrics_nn.rds")

##########
## plot ##
##########

## load 

df_metrics_st <- readRDS("./RWR/df_metrics_st.rds")[,1:3]
df_metrics_nn <- readRDS("./RWR/df_metrics_nn.rds")
df_metrics_om <- readRDS("./RWR/df_metrics_om.rds")[,1:3]

df_metrics <- do.call('rbind', list(df_metrics_st,df_metrics_om,df_metrics_nn))
df_metrics$software <- rep(c('stMLnet','Omnipath','NicheNet'),each=30)

df_plot <- df_metrics %>% reshape2::melt(id = c('r','software'))

## check

df <- df_plot[df_plot$software == 'stMLnet' & df_plot$variable == 'mean_Precision',]
df[order(df$value,decreasing = T),] %>% head(.,5)

df <- df_plot[df_plot$software == 'stMLnet' & df_plot$variable == 'mean_AUC',]
df[order(df$value,decreasing = T),] %>% head(.,5)

## color

library(RColorBrewer)
cols <- brewer.pal(3, "Set1")
names(cols) <- c("stMLnet","NicheNet","Omnipath")

# plot 

p1 <- ggplot(df_plot[df_plot$variable == 'mean_AUC',], aes(x=r,y=value,group=software,color=software)) + 
  labs(x='restart probability',y='AUC') + xlim(c(0,1)) + 
  scale_y_continuous(breaks=seq(0.65,1,0.05)) +
  geom_point() + geom_line() + theme_bw() + 
  scale_color_manual(values = cols) + 
  geom_vline(aes(xintercept=0.61), colour="grey", linetype="dashed") +
  geom_hline(aes(yintercept=df_plot$value[df_plot$software=='stMLnet' &
                                            df_plot$r==0.61 &
                                            df_plot$variable=='mean_AUC']), 
             colour="grey", linetype="dashed") +
  theme(
    legend.position = 'bottom',
    legend.title = element_blank(),
    legend.direction = 'horizontal',
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) 
p1

p2 <- ggplot(df_plot[df_plot$variable == 'mean_Precision',], aes(x=r,y=value,group=software,color=software)) + 
  labs(x='restart probability',y='Precision') + xlim(c(0,1)) + 
  scale_y_continuous(breaks=seq(0,0.6,0.1)) +
  geom_point() + geom_line() + theme_bw() + 
  scale_color_manual(values = cols) + 
  geom_vline(aes(xintercept=0.61), colour="grey", linetype="dashed") +
  geom_hline(aes(yintercept=df_plot$value[df_plot$software=='stMLnet' &
                                            df_plot$r==0.61 &
                                            df_plot$variable=='mean_Precision']), 
             colour="grey", linetype="dashed") +
  theme(
    legend.position = 'bottom',
    legend.title = element_blank(),
    legend.direction = 'horizontal',
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) 
p2

pdf("./RWR/plot_metrics.pdf",height = 4,width = 7.5)
png("./RWR/plot_metrics.png",height = 4,width = 7.5,units = 'in',res = 300)  
gridExtra::grid.arrange(p1,p2, nrow=1)
dev.off()
