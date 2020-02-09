crisp <- function(u, method, tv){
  if(missing(u))
    stop("Missing membership matrix for crisp labeling")
  if(missing(method) & missing(tv))
    method <- "max"
  if(missing(method) & !missing(tv))
    method <- "threshold"
  if(!is.element(method, c("max", "threshold")))
    stop("Invalid method. The valid options are 'max' or 'threshold'")
  if(!missing(tv)){
    if(!is.numeric(tv))
      stop("Threshold value for crisp labelling should be a number")
    if((tv < 0) || (tv > 1))
      stop("Threshold value should be between 0 and 1")
  }else{
    tv <- 0.50
  }
  clabels <- numeric(nrow(u))
  if(method == "max") 
    clabels <- apply(u, 1, function(x) which.max(x))
  else
    clabels <- apply(u, 1, function(x) ifelse(all(x < tv), 0, which.max(x)))
  return(clabels)
}

.crisp2 <- function(u, method, tv){
  if(missing(u))
    stop("Missing membership matrix for crisp labeling")
  if(missing(method) & missing(tv))
    method <- "max"
  if(missing(method) & !missing(tv))
    method <- "threshold"
  method <- pmatch(method, c("max", "threshold"))
  if(is.na(method))
    stop("Invalid method. The valid options are 'max' or 'threshold'")
  if(method == -1) 
    stop("Unapplied method in this package. The valid options are 'max' or 'threshold'")
  if(!missing(tv)){
    if(!is.numeric(tv))
      stop("Threshold value for crisp labelling should be a number")
    if((tv < 0) || (tv > 1))
      stop("Threshold value should be between 0 and 1")
  }else{
    tv <- 0.50
  }
  clabels <- numeric(nrow(u))
  if(method == "max"){
    idx <- 1
    for(i in 1:nrow(u)){
      maxc <- 0
      for(j in 1:ncol(u)){
        if(u[i,j] > maxc){
          maxc <- u[i,j]
          idx <- j
        }
      }
      clabels[i] <- idx
    }
  }else {
    for(i in 1:nrow(u)){
      clusts <- which(u[i,] >= tv)
      if(length(clusts) == 0){
        idx <- 0
      }else{
        maxc <- 0
        for(j in clusts){
          if(u[i,j] > maxc){
            maxc <- u[i,j]
            idx <- j
          }
        }
      }
      clabels[i] <- idx
    }
  }
  return(clabels)
}

.compdist <- function(a, b, dmetric="euclidean", pw=2){
  if(missing(a))
    stop("Missing data arguments to compute the distances")
  if(missing(b))
    b <- a
  if(!is.numeric(a) || !is.numeric(b))
    stop("Input data arguments must be numeric to compute the distances")
  dmetrics <- c("euclidean", "sqeuclidean", "manhattan", "minkowski", "chebyshev", 
    "pearsonchi", "neymanchi", "sqchi", "divergence", "addsymchi", "prosymchi", "clark", 
    "canberra", "sorensen", "lorentzian", "sqchord", "cosine", "correlation")
  dmetric <- match.arg(dmetric, dmetrics)
  if(dmetric=="euclidean")
    distance <- sqrt(sum(t(a-b) * (a-b)))
  else if (dmetric=="sqeuclidean")
    distance <- sum(t(a-b) * (a-b))
  else if (dmetric=="manhattan")
    distance <- sum(abs(a-b))
  else if (dmetric=="minkowski")
    distance <- sum(abs(a-b)^pw)^(1/pw)
  else if (dmetric=="chebyshev")
    distance <- max(abs(a-b))
  else if (dmetric=="pearsonchi")
    distance <- sum((t(a-b) * (a-b) / b))
  else if (dmetric=="neymanchi")
    distance <- sum(t(a-b) * (a-b) / a)
  else if (dmetric=="sqchi")
    distance <- sum(t(a-b) * (a-b) / (a+b))
  else if (dmetric=="addsymchi")
    distance <- sum((t(a-b) * (a-b)) * ((a+b) / (a*b)))
  else if (dmetric=="prosymchi")
    distance <- 2 * sum((t(a-b) * (a-b)) / (a+b))
  else if (dmetric=="divergence")
    distance <- 2 * sum((t(a-b) * (a-b)) / (t(a+b) * (a+b)))
  else if (dmetric=="clark")
    distance <- sqrt(sum(abs((a-b) / (a+b))^2))
  else if (dmetric=="sqchord")
    distance <- sum(sqrt(a)-sqrt(b))^2
  else if (dmetric=="canberra")
    distance <- sum(abs(a-b)/(a+b))
  else if (dmetric=="sorensen")
    distance <- sum(abs(a-b))/sum(a+b)
  else if (dmetric=="lorentzian")
    distance <- sum(log(1+abs(a-b), exp(1)))
  else if (dmetric=="cosine")
    distance <- sum(a%*%b)/(sqrt(sum(a*a)) * sqrt(sum(b*b)))
  else if (dmetric=="correlation")
    distance <- (1-cor(a,b))/2
  return(distance)
}

.get.dtype <- function(dmetric="euclidean", pw=2){
  dmetrics <- c("euclidean", "sqeuclidean", "manhattan", "minkowski", "sqchord", "chebyshev", 
    "pearsonchi", "neymanchi", "sqchi", "divergence", "addsymchi", "prosymchi", "clark", 
    "canberra", "sorensen", "lorentzian", "cosine", "correlation")
  dmetric <- match.arg(dmetric, dmetrics)
  if(dmetric=="euclidean")
    dtype <- 2
  else if (dmetric=="sqeuclidean")
    dtype <- 2
  else if (dmetric=="manhattan")
    dtype <- 2
  else if (dmetric=="minkowski")
    dtype <- pw
  else if (dmetric=="chebyshev")
    dtype <- 2
  else if (dmetric=="pearsonchi")
    dtype <- 2
  else if (dmetric=="neymanchi")
    dtype <- 2
  else if (dmetric=="sqchi")
    dtype <- 2
   else if (dmetric=="addsymchi")
    dtype <- 2
  else if (dmetric=="prosymchi")
    dtype <- 2
  else if (dmetric=="divergence")
    dtype <- 2
  else if (dmetric=="clark")
    dtype <- 2
  else if (dmetric=="sqchord")
    dtype <- 2
  else if (dmetric=="canberra")
    dtype <- 1
   else if (dmetric=="sorensen")
    dtype <- 1
  else if (dmetric=="lorentzian")
    dtype <- 1
  else if (dmetric=="cosine")
    dtype <- 2
  else if (dmetric=="correlation")
    dtype <- 2
  return(dtype)
}

get.dmetrics <- function(dmt="all"){
    dmt <- match.arg(dmt, c("all", "l1","l2", "li", "lp", "sl2", "sc"))
    dmetrics <- c("manhattan", "canberra", "sorensen", "lorentzian", "euclidean", 
      "minkowski", "chebyshev", "sqeuclidean", "pearsonchi", "neymanchi", "sqchi", 
      "divergence", "addsymchi", "prosymchi", "clark", "sqchord", "cosine", "correlation")
    names(dmetrics) <- c("Manhattan", "Canberra", "Sorensen", "Lorentzian", "Euclidean", 
      "Minkowski",  "Chebyshev", "Squared Euclidean", "Pearson Chi", "Neyman Chi", "Squared Chi",
      "Divergence", "Additive Symmetric", "Probabilistic Symmetric", "Clark", "Squared Chord", 
      "Cosine", "Correlation") 
    dmtypes <- c("l1", "l1", "l1", "l1", "l2", "lp", "li", "sl2", "sl2", "sl2", "sl2", "sl2", "sl2", "sl2", "sl2", "sc","cos","cor") 
    if(dmt=="all"){
      smetrics <- dmetrics
      smetricst <- dmtypes
    }
    else {
      smetrics <- dmetrics[dmtypes==dmt]
      smetricst <- dmtypes[dmtypes==dmt]
    }
    dmlist <- matrix(nrow=length(smetrics), ncol=3)
    for(i in 1:nrow(dmlist)){
      dmlist[i,1] <- smetrics[i]
      dmlist[i,2] <- names(smetrics[i])
      dmlist[i,3] <- smetricst[i]
    }
    colnames(dmlist) <- c("dmetric", "Description", "Norm")
    dmlist <- dmlist[order(dmlist[,1]),]
    return(dmlist)
}

.sumsqr <- function(x, v, clusters){
  sumsqr <- function(x) sum(scale(x, scale = FALSE)^2)
  bwss <- sumsqr(v[clusters,])
  wss <- sapply(split(as.data.frame(x), clusters), sumsqr)
  twss <- sum(wss)
  tss <- bwss + twss
  ss <- list(bwss, wss, twss, tss)
  names(ss) <- c("between.ss", "within.ss", "tot.within.ss", "tot.ss")
  return(ss)
}

comp.omega <- function(d, u, m=2,  pco=NULL, K=1){
  if(!missing(pco)){
    if(pco$algorithm=="PCM")
      u <- pco$t
    if(pco$algorithm=="FCM")
      u <- pco$u
    d <- pco$d
    m <- pco$m
  }
  if(missing(u))
    stop("Missing typicality matrix")
  if(is.null(u))
    stop("Memberships matrix contains no data")
  if(is(u, "data.frame"))
    u <- as.matrix(u)
  if(!is(u, "matrix"))
    stop("Memberships data must be numeric data frame or matrix")
  if(any(is.na(u)))
    stop("Memberships data should not contain NA values")
  if(!is.numeric(u))
   stop("Memberships data must be a numeric data frame or matrix")
  if(missing(d))
    stop("Missing distance matrix")
  if(is.null(d))
    stop("Distance matrix contains no data")
  if(is(d, "data.frame"))
    d <- as.matrix(d)
  if(!is(d,"matrix"))
    stop("Distance data must be numeric data frame or matrix")
  if(any(is.na(d)))
    stop("Distance data should not contain NA values")
  if(!is.numeric(m))
    stop("Argument 'm' must be a numeric value")
  if(!is.numeric(K))
    stop("Argument 'K' must be a numeric value")
  omega <- numeric(ncol(u))
  for(i in 1:ncol(u))
    omega[i] <- K * ((u[,i]^m) %*% d[,i])/sum(u[,i]^m)
  return(omega)
}

ppclust2 <- function(objx, otype, ...){
  if(!inherits(objx, "ppclust")) 
    stop("Input object should be an instance of 'ppclust' class")
  otypes <- c("fanny", "fclust", "kmeans", "vegclust")
  otype <- match.arg(otype, otypes)
  items <- list()
  if(otype=="vegclust"){
    vmethod <- objx$algorithm
    if(objx$algorithm=="HCM"){
      vmethod <- "KM"
    }
    if(!is.null(objx$t)){
      memb <- objx$t
      expc <- objx$eta
    }
    else{
      memb <- objx$u 
      expc <- objx$m
    }  
    items$mode <- "raw"
    items$method <- vmethod
    items$m <- expc
    items$dnoise <- NULL
    items$omega <- objx$omega
    items$alpha <- objx$inpargs$con.val
    items$memb <- memb
    colnames(items$memb) <- paste0("M",1:ncol(memb))
    items$mobileCenters <- objx$centers
    items$fixedCenters <- NULL
    items$dist2clusters <- objx$d
    colnames(items$dist2clusters) <- paste0("M",1:ncol(memb))
    items$withinss <- objx$sumsqrs$withinss
    items$size <- objx$csize
    names(items$size) <- paste0("M",1:length(objx$csize))
    items$functional <- objx$func.val[objx$best.start]
    items$iter <- objx$iter[objx$best.start]
    class(items) <- c("vegclust")
  }
  if(otype=="fclust"){
    if(!is.null(objx$t)){
      memb <- objx$t
      expc <- objx$eta
    }
    else{
      memb <- objx$u 
      expc <- objx$m
    } 
    items$U <- memb
    items$H <- objx$d
    items$F <- objx$f
    items$clus <- objx$cluster
    items$medoid <- NULL
    items$value <- objx$func.val
    items$cput <- objx$comp.time
    items$iter <- objx$iter
    items$k <- objx$k
    items$m <- expc
    items$ent <- NULL
    items$b <- NULL
    items$vp <- NULL
    items$delta <- NULL
    items$stand <- objx$inpargs$stand
    if(objx$inpargs$stand==TRUE)
      items$Xca <- scale(objx$x, center = TRUE, scale = TRUE)[, ]
    else
      items$Xca <- objx$x
    items$X <- objx$x
    items$call <- objx$call
    class(items) <- c("fclust")
  }
  if(otype=="kmeans"){
    items$cluster <- unname(objx$cluster)
    items$centers <- objx$v
    items$totss <- objx$sumsqrs$tot.ss
    items$withinss <- unname(objx$sumsqrs$withinss)
    items$tot.withinss <- objx$sumsqrs$tot.withinss
    items$betweenss <- objx$sumsqrs$betweenss
    items$size <- objx$csize
    items$iter <- objx$iter[objx$best.start]
    items$ifault <- 0L
    class(items) <- c("kmeans")
  }
  if(otype=="fanny"){
    if(!is.null(objx$u))
      memb <- objx$u
    else
      memb <- objx$t   
    items$membership <- memb
    dunnfk <- numeric(2) 
    dunnfk[1] <- sum(objx$u^2) / nrow(objx$x)
    dunnfk[2] <- (dunnfk[1]-1/objx$k) / (1-1/objx$k)
    names(dunnfk) <- c("dunn_coeff", "normalized")
    items$coeff <- as.vector(dunnfk)
    items$memb.exp <- objx$m
    items$clustering <- unname(as.integer(objx$cluster))
    items$k.crisp <- objx$k
    objective <- numeric(2)
    names(objective) <- c("objective","tolerance") 
    objective[1] <- objx$func.val[objx$best.start]
    objective[2] <- objx$inpargs$con.val
    items$objective <- objective
    convergence <- integer(3)
    names(convergence) <- c("iterations", "converged", "maxit")
    convergence[1] <- objx$iter[objx$best.start]
    convergence[2] <- 1
    convergence[3] <- objx$inpargs$iter.max
    items$convergence <- convergence
    items$diss <- NULL
    items$call <- NULL
    items$silinfo <- NULL
    items$data <- as.matrix(objx$x)
    class(items) <- c("fanny","partition")
  }
  return(items)
}

as.ppclust <- function(objx, ...){
  if(is.null(objx))
    stop("Input an object to be converted to 'ppclust'")
  otype <- class(objx)[1]
  otypes <- c("fanny", "fclust", "kmeans", "vegclust")
  otype <- match.arg(otype, otypes)
  items = list()
  if(otype=="fclust"){
    items$u <- objx$U
    items$t <- objx$T
    items$v <- objx$H
    items$f <- objx$F
    items$cluster <- objx$clus
    items$func.val <- objx$func
    items$comp.time <- objx$cput
    items$iter <- objx$iter
    items$k <- objx$k
    items$m <- objx$m
    items$ent <- objx$ent
    items$b <- objx$b
    items$delta <- objx$delta
    items$inpargs$stand <- objx$stand
    items$x <- objx$X
    items$call <- objx$call
  }
  if(otype=="kmeans"){
    items$cluster <- objx$cluster
    items$v <- objx$centers
    items$sumsqrs$tot.ss <- objx$totss
    items$sumsqrs$withinss <- objx$withinss
    items$sumsqrs$tot.withinss <- objx$tot.withinss
    items$sumsqrs$betweenss <- objx$betweenss
    items$csize <- objx$size
    items$iter <- objx$iter
  }
  if(otype=="vegclust"){
    items$u <-objx$memb
    items$t <-objx$memb
    items$v <-objx$mobileCenters
    items$d <-objx$dist2clusters
    items$func.val <- objx$functional
    items$iter <- objx$iter
    items$k <- objx$k
    items$m <- objx$m
    items$omega <- objx$eta
    items$csize <- objx$size
    items$algorithm <- objx$method
  }
  if(otype=="fanny"){
    items$u <- objx$membership
    items$m <- objx$memb.exp
    items$cluster <- objx$clustering
    items$k <- objx$k.crisp
    items$func.val <- objx$objective[1]
    items$inpargs$con.val <- objx$objective[2]
    items$inpargs$iter.max <- objx$convergence[3]
    items$iter <- objx$convergence[1]
    items$x <- objx$data
    items$algorithm <- "FCM"
  }
  class(items) <- c("ppclust")
  return(items)
}

is.ppclust <- function(objx){
  class(objx)=="ppclust"
}

summary.ppclust <- function (object, ...){
  pco <- object
  if(missing(pco))
    stop("Missing input object")
  if(!inherits(pco, "ppclust")) 
    stop("The input object must be an instance of the class 'ppclust'")
  talgs <- c("FPCM", "PCA", "PCM", "PCMR", "PFCM", "UPFC")
  if(is.element(pco$algorithm, talgs))
    u <- pco$t
  else 
    u <- pco$u
  if(is.null(pco$k))
    k <- ncol(u)
  else 
    k <- pco$k
  if(is.null(pco$cluster))
    clusters <- crisp(u)
  else
    clusters <- pco$cluster
  cat(paste0("Summary for '", (as.list(match.call())[-1])[1],"'\n"))
  cat("\nNumber of data objects: " , nrow(pco$x), "\n")
  cat("\nNumber of clusters: ", k, "\n")
  cat("\nCrisp clustering vector:\n")
  print(unname(pco$cluster))
  if(!is.null(pco$v0)){
    cat("\nInitial cluster prototypes:\n")
    print(round(pco$v0,9))
  }
  cat("\nFinal cluster prototypes:\n")
  print(round(pco$v,9))
  cat("\nDistance between the final cluster prototypes\n")
  print(round(dv <- dist(pco$v)^2,9))
  
  if(!is.null(pco$v0)){
    cat("\nDifference between the initial and final cluster prototypes\n")
    print(round(pco$v-pco$v0, 9))
    cat("\nRoot Mean Squared Deviations (RMSD):", sqrt(sum((pco$v-pco$v0)^2)/k),"\n")
    cat("Mean Absolute Deviation (MAD):", sum(abs(pco$v-pco$v0))/nrow(pco$v)*ncol(pco$v),"\n")
  }
  cat("\nMembership degrees matrix (top and bottom 5 rows): \n")
  print(round(head(u, 5),9))
  cat("...\n")
  print(round(tail(u, 5),9))
  if(pco$algorithm != "KM"){
    cat("\nDescriptive statistics for the membership degrees by clusters\n")
    nc <- minmemdeg <- q1memdeg <- meanmemdeg <- medmemdeg <- q3memdeg <- maxmemdeg <- rep(0, k)
    for(i in 1:k){
      nc[i] <- length(u[clusters==i,i])
      quantiles <- quantile(u[clusters==i,i], names=FALSE)
      minmemdeg[i] <- quantiles[1]
      q1memdeg[i] <- quantiles[2]
      medmemdeg[i] <- quantiles[3]
      q3memdeg[i] <- quantiles[4]
      maxmemdeg[i] <- quantiles[5]
      meanmemdeg[i] <- mean(u[clusters==i,i])
    }
    df <- cbind(nc, minmemdeg, q1memdeg, meanmemdeg, medmemdeg, q3memdeg, maxmemdeg)
    colnames(df) <- c("Size", "Min", "Q1", "Mean", "Median", "Q3", "Max")
    rownames(df) <- paste0("Cluster ", 1:k)
    print(df)
  }
  dunnfk <- c() 
  dunnfk[1] <- sum(u^2) / nrow(u)
  dunnfk[2] <- (dunnfk[1]-1/k) / (1-1/k)
  names(dunnfk) <- c("dunn_coeff", "normalized")
  cat("\nDunn's Fuzziness Coefficients:\n")
  print(dunnfk)
  cat("\nWithin cluster sum of squares by cluster:\n")
  print(pco$sumsqrs$within.ss)
  cat("(between_SS / total_SS = ", paste0(round(pco$sumsqrs$between.ss/pco$sumsqrs$tot.ss*100,2),"%)"),"\n")
  cat("\nAvailable components: \n")
  print(names(pco))
  invisible(pco)
}

.colpal <- function(pno){
  paletmat <- matrix(nrow=5, ncol=20)
  paletmat[1,] <- c("dodgerblue", "firebrick3", "green4", "mediumorchid1", 
                    "navyblue","orange", "seagreen4","maroon1",
                    "yellow4", "lightsteelblue4","orangered1","mediumpurple3",
                    "khaki4","gray47", "darkseagreen3","darksalmon",
                    "deepskyblue1","blue1","chocolate4","brown1")
  paletmat[2,] <- c("dodgerblue", "burlywood4", "aquamarine4", "red", 
                    "chocolate","azure4", "chocolate4", "coral4", "blue",
                    "blueviolet", "brown4", "darkcyan", "darkgoldenrod",
                    "darkgreen", "bisque4", "darkorchid", "chartreuse4",
                    "indianred1", "mediumslateblue", "peru")
  paletmat[3,] <- c("dodgerblue", "lightblue4", "lightcoral", "lightcyan4", 
                   "lightgoldenrod4", "lightpink4", "lightsalmon", "lightsalmon4", 
                   "lightskyblue4", "lightyellow4", "limegreen", "olivedrab4", 
                   "magenta2", "magenta4", "maroon2", "maroon4", "hotpink1", 
                   "hotpink3", "khaki3", "khaki4") 
  paletmat[4,] <- c("seashell4", "sienna", "skyblue", "skyblue4", "slateblue",
                   "slateblue4", "slategray", "slategray4", "tan3", "tan4", 
                   "thistle", "thistle4", "tomato", "tomato4", "turquoise1", 
                   "turquoise4", "violet", "violetred", "violetred4", "yellowgreen") 
  paletmat[5,] <- c("dodgerblue", "deepskyblue", "deepskyblue1", "deepskyblue2", "deepskyblue3", 
                   "deepskyblue4", "lightskyblue", "lightskyblue1", "lightskyblue2",
                   "lightskyblue3", "lightskyblue4", "skyblue" ,"skyblue1",
                   "skyblue2", "skyblue3", "skyblue4","aquamarine1", "aquamarine2", "aquamarine3", "aquamarine4")
  return(paletmat[pno,])
}

plotcluster <- function(objx, mt, cm, tv, cp=1, pt=19, trans=FALSE){
  if(!is.ppclust(objx))
    stop("The input should be an object of ppclust")
  if(missing(mt)){
    if(is.null(objx$t))
      mt <- "u"
    else 
      mt <- "t"
  }
  if(is.na(pmatch(mt, c("u","t"))))
    stop("Invalid membership type with 'mt' argument. Use 'u' for fuzzy memberships or 't' for typicalities")
  if(mt=="u"){
    if(is.null(objx$u))
      stop("Null membership matrix with 'mt' argument. Use 't' for typicalities matrix")
    else
      memb <- objx$u
  }else{
     if(is.null(objx$t))
      stop("Null typicality matrix with 'mt' argument. Use 'u' for memberships matrix")
    else
      memb <- objx$t
  }
  if(!is.logical(trans))
    stop("Invalid transparency option. Use 'FALSE' for solid colors or 'TRUE' for transparent colors")
  if(missing(cm) && missing(tv))
    cm <- "max"
  if(missing(cm) && !missing(tv))
    cm <- "threshold"
  if(is.na(pmatch(cm, c("max","threshold"))))
    stop("Invalid crisping method. Use 'max' for the maximum memberships or 'threshold' for the memberships exceeding tv, selected threshold value")
  if(missing(tv))
    tv <- 0.5
  if(!is.numeric(tv))
    stop("The argument 'tv', threshold value must be 0 and 1")
  if(tv < 0 || tv > 1)
    stop("The argument 'tv', threshold value must be 0 and 1")
  if(trans)
    cols <- adjustcolor(.colpal(cp), alpha.f = 0.7)
  else
    cols <- .colpal(cp)
  cluster <- .crisp2(u=memb, method=cm, tv=tv)
  if(pt=="#")
    pt <- cluster + 48
  cluster[cluster == 0] <- NA
  pairs(objx$x, col = cols[cluster], pch = pt, cex=1.1)
}
