ekm <- function (x, centers, dmetric="euclidean", alginitv="hartiganwong",  
       nstart=1, iter.max=1000, stand=FALSE, numseed){
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
  }
  else{
    if(!is.numeric(centers))
      stop("Centers should be an integer for the number of clusters or a numeric prototypes matrix")
    k <- ceiling(centers)
  }
  if((k < 1) || (k > n))
    stop(paste0("k, number of clusters should be between 1 and ", n, ". Check the value of 'centers' argument"))

  if(!missing(numseed)){
    if(!is.numeric(numseed))
      stop("Argument numseed should be a number")
   else
      set.seed(numseed)
  }
  alginitv <- match.arg(alginitv, inaparc::get.algorithms("prototype"))
  compv <- parse(text = paste0("inaparc::", alginitv, "(x, k)$v"))
  if(!is.matrix(centers))
    centers <- matrix(nrow=k, ncol=p, eval(compv))
  dmetrics <- c("euclidean")
  dmetric <- match.arg(dmetric, dmetrics)
  compd <- parse(text = paste0(".compdist(x[i,], v[j,], dmetric='", dmetric, "')"))
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
  d <- matrix(nrow = n, ncol = k, 0)
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
    cputime <- system.time(
      res.km <- kmeans(x = x, centers = centers, iter.max = iter.max)
    )
    comp.time[start.idx] <- cputime[1]
    iter.num[start.idx] <- res.km$iter
    obj.func <- res.km$tot.withinss
    func.val[start.idx] <- obj.func
    if(obj.func < best.func){
      best.func <- obj.func
      v0 <- centers
      best.v <- v <- res.km$centers
      for(i in 1:n)
        for(j in 1:k)
          d[i,j] <- eval(compd)
       best.d <- d
       betweenss <- res.km$betweenss
       withinss <- res.km$withinss
       tot.withinss <- res.km$tot.withinss
       tot.ss <- res.km$totss     
       cluster <- res.km$cluster
       size <- res.km$size
       best.start <- start.idx
    }
  }
  u <- matrix(nrow=n, ncol=k, 0)
  for(i in 1:n)
    for(j in 1:k)
      u[i, cluster[i]] <- 1
  csumsqrs <- list(betweenss, withinss, tot.withinss, tot.ss)
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
  names(cluster) <- paste(1:n, sep = " ")
  names(size) <- paste(1:k, sep = " ")
  result = list()
    result$u <- u
    result$t <- NULL
    result$v <- best.v
    result$v0 <- v0
    result$d <- best.d
    result$f <- NULL
    result$x <- x
    result$cluster <- cluster
    result$csize <- size
    result$sumsqrs <- csumsqrs
    result$k <- k
    result$m <- NULL
    result$eta <- NULL
    result$a <- NULL
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
    result$inpargs[2] <- NA
    result$inpargs[3] <- dmetric	
    result$inpargs[4] <- alginitv	
    result$inpargs[5] <- NA	
    result$inpargs[6] <- NA
    result$inpargs[7] <- NA
    result$inpargs[8] <- stand
    names(result$inpargs) <- c("iter.max", "con.val", "dmetric", "alginitv", "alginitu", "fixcent", "fixmemb", "stand")
    result$algorithm <- "KM"
    result$call <- match.call()
  class(result) <- c("ppclust")
  return(result)
}
