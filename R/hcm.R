hcm <- function(x, centers, dmetric="euclidean", pw=2, alginitv="kmpp", 
       nstart=1, iter.max=1000, con.val=1e-09, stand=FALSE, numseed){
  if(missing(x))
    stop("Missing data set")
  if(is.null(x))
    stop("Data set is null")
  if((is.data.frame(x)) || (is.vector(x)))
    x <- as.matrix(x)
  if(!is.matrix(x))
    stop("Data set must be a vector, data frame or matrix")
  if(any(is.na(x)))
    stop("Data set should not contain NA values. Remove NAs and try again")
  if(!is.numeric(x))
   stop("Data set must be a numeric vector, data frame or matrix")
  n <- nrow(x) ; p <- ncol(x)
  if(!missing(numseed))
    set.seed(numseed)
  if(missing(centers))
    stop("Missing argument 'centers'")
  if(is.data.frame(centers))
    centers <- as.matrix(centers)
  if(is.matrix(centers)){
    if(is.null(centers))
      stop("Centers contains no data")
    if(any(is.na(centers)))
      stop("Centers should not contain NA values")
    if(!is.numeric(centers))
     stop("Centers should be numeric")
    else
      k <- nrow(centers)
    if(anyDuplicated(centers) > 0)
      stop("Identical centers are not allowed. Check the centers data.")
  }
  else{
    if(!is.numeric(centers))
      stop("Centers should be an integer for the number of clusters or a numeric prototypes matrix")
    k <- ceiling(centers)
  }
  if((k < 1) || (k > n))
    stop(paste0("k, number of clusters should be between 1 and ", n, ". Re-enter the 'centers' argument"))

  if(!missing(numseed)){
    if(!is.numeric(numseed))
      stop("Argument numseed should be a number")
   else
      set.seed(numseed)
  }
  alginitv <- match.arg(alginitv, inaparc::get.algorithms("prototype"))
  compd <- parse(text = paste0(".compdist(x[i,], v[j,], dmetric='", dmetric, "', p=", pw, ")")) 
  compv <- parse(text = paste0("inaparc::", alginitv, "(x, k)$v"))
  if(!is.matrix(centers))
    centers <- matrix(nrow=k, ncol=p, eval(compv))
  if(!is.logical(stand))
    stop("Value of argument 'stand' should be TRUE or FALSE")
  if(stand){
    x <- scale(x, center = TRUE, scale = TRUE)[, ]
    centers <- scale(centers, center = TRUE, scale = TRUE)[, ]
  }
  if(!is.numeric(nstart))
    stop("Number of starts must be integer")
  if(nstart < 1)
    stop("Number of starts cannot be less than 1")
  if(nstart%%ceiling(nstart) > 0)
    nstart <- ceiling(nstart)
  if(!is.numeric(iter.max))
    stop("Maximum number of iterations must be a positive integer")
  else
    iter.max <- ceiling(iter.max)
  if(iter.max <= 1)
    stop("Maximum number of iterations must be equal to or greater than 1")
  func.val <- numeric(nstart)
  comp.time <- numeric(nstart)
  iter.num <- numeric(nstart)
  best.func <- Inf
  for(start.idx in 1:nstart){
    if(start.idx > 1){
      set.seed(as.integer(Sys.time()) + start.idx)
      centers <- eval(compv)
      if(stand)
        centers <- scale(centers, center = TRUE, scale = TRUE)[, ]
      if(!missing(numseed))
        set.seed(numseed)
    }
    d <- u <- matrix(nrow = n, ncol = k, 0)
    v0 <- v <- centers
    prevv <- v + 2*con.val
    iter <- 0
    cputime <- system.time(
      while((con.val < abs(prevv - v)) && (iter < iter.max)){
        prevv <- v
        iter <- iter + 1
        for(j in 1:k)
          for(i in 1:n)
            d[i,j] <- eval(compd)
        for(j in 1:k){
          for(i in 1:n){
            u[i,] <- 0
            u[i, which.min(d[i,])] <- 1
          }
        }
        clusters <- crisp(u)
        ssqs <- .sumsqr(x, v, clusters)
        wss <- ssqs$within.ss
        bwss <- ssqs$between.ss
        twss <- ssqs$tot.within.ss
        tss <- ssqs$tot.ss
        for(j in 1:k){
          xc <- as.matrix(x[u[,j]==1,])
          if(ncol(xc)==1)
            v[j,]  <- mean(xc)
          else
            v[j,] <- apply(xc, 2, mean)         
        }
      }
    )
    comp.time[start.idx] <- cputime[1]
    iter.num[start.idx] <- iter
    obj.func <- twss
    func.val[start.idx] <- obj.func
    if(obj.func < best.func){
      best.func <- obj.func
      best.v0 <- v0
      best.u <- u
      best.d <- d
      best.v <- v
      best.clusters <- clusters
      best.betweenss <- bwss
      best.withinss <- wss
      best.tot.withinss <- twss
      best.tot.ss <- tss    
      best.start <- start.idx
    }
  }
  clabels <- best.clusters
  csize <- numeric(k)
  for(i in 1:k)
    csize[i] <- sum(clabels==i)
  csumsqrs <- list(best.betweenss, best.withinss, best.tot.withinss, best.tot.ss)
  names(csumsqrs) <- c("between.ss", "within.ss", "tot.within.ss", "tot.ss")
  if(is.null(rownames(x)))
    rnames <- paste(1:n)
  else
    rnames <- rownames(x)
  if(is.null(colnames(x)))
    cnames <- paste("p", 1:p, sep = "")
  else
    cnames <- colnames(x)
  rownames(x) <- rnames; colnames(x) <- cnames
  rownames(best.v) <- paste("Cluster", 1:k, sep = " ")
  colnames(best.v) <- cnames
  colnames(u) <- rownames(best.v)
  rownames(u) <- rownames(x)
  colnames(v0) <- colnames(best.v)
  rownames(v0) <- rownames(best.v)
  rownames(best.d) <- rnames
  colnames(best.d) <- rownames(best.v)
  names(clabels) <- paste(1:n, sep = " ")
  names(csize) <- paste(1:k, sep = " ")
  result = list()
    result$u <- best.u
    result$t <- NULL
    result$v <- best.v
    result$v0 <- best.v0
    result$d <- best.d
    result$f <- NULL
    result$x <- x
    result$cluster <- clabels
    result$csize <- csize
    result$sumsqrs <- csumsqrs
    result$k <- k
    result$m <- NULL
    result$eta <- NULL
    result$b <- NULL
    result$beta <- NULL
    result$delta <- NULL
    result$gamma <- NULL
    result$omega <- NULL
    result$ent <- NULL
    result$iter <- iter.num
    result$best.start <- best.start
    result$func.val <- func.val
    result$comp.time <-  comp.time
    result$inpargs <- list()
    result$inpargs[1] <- as.integer(iter.max)
    result$inpargs[2] <- con.val
    result$inpargs[3] <- dmetric	
    result$inpargs[4] <- alginitv	
    result$inpargs[5] <- NA	
    result$inpargs[6] <- NA
    result$inpargs[7] <- NA
    result$inpargs[8] <- stand
    names(result$inpargs) <- c("iter.max", "con.val", "dmetric", "alginitv", "alginitu", "fixcent", "fixmemb", "stand")
    result$algorithm <- "HCM"
    result$call <- match.call()
  class(result) <- c("ppclust")
  return(result)
}
