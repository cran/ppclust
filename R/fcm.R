fcm <- function(x, centers, memberships, m=2, 
       dmetric="sqeuclidean", pw=2, alginitv="kmpp", 
       alginitu="imembrand", nstart=1, iter.max=1000, con.val=1e-09, 
       fixcent=FALSE, fixmemb=FALSE, stand=FALSE, numseed){
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
    stop(paste("k, number of clusters should be between 1 and", n, ". Check the value of 'centers' argument"))

  if(!missing(numseed)){
    if(!is.numeric(numseed))
      stop("Argument numseed should be a number")
   else
      set.seed(numseed)
  }
  alginitv <- match.arg(alginitv, inaparc::get.algorithms("prototype"))
  compv <- parse(text = paste0("inaparc::", alginitv, "(x,", k, ")$v"))
  algsinitu <- match.arg(alginitu, inaparc::get.algorithms("membership"))
  alginitu <- match.arg(alginitu, algsinitu)
  if(alginitu=="imembrand")
    compu <- parse(text = paste0("inaparc::", alginitu, "(n,", k, ")$u"))
  else
    compu <- parse(text = paste0("inaparc::", alginitu, "(x,", k, ")$u"))
  compd <- parse(text = paste0(".compdist(x[i,], v[j,], dmetric='", dmetric, "', p=", pw, ")"))
  dt <- .get.dtype(dmetric, pw=pw)
  if(!is.matrix(centers))
    centers <- matrix(nrow=k, ncol=p, eval(compv))
  if(!missing(memberships)){
    if(class(memberships) == "data.frame")
      memberships <- as.matrix(memberships)
    if(class(memberships) != "matrix")
      stop("The initial membership degrees matrix is not a numeric data.frame or matrix")
  }else{
    memberships <- matrix(nrow = n, ncol = k, eval(compu))
  }
  if(is.null(memberships))
    stop("The initial membership matrix cannot be empty")
  if(any(!is.numeric(memberships)))
    stop("The initial membership matrix is not a numeric data.frame or matrix")
  if(any(is.na(memberships)))
    stop("The initial membership matrix should not contain NAs")
  if(n != nrow(memberships))
    stop("The number of rows of initial membership matrix is different from that of data set")
  if(k != ncol(memberships))
    stop("The number of columns of initial membership matrix is not equal to k, number of clusters")
  if(sum(memberships) != n)
    memberships <- memberships/apply(memberships, 1, sum)
  if(!is.numeric(m))
    stop("The fuzziness exponent (m) should be a number")
  if(m < 1)
    stop("The fuzziness exponent (m) should be a number equals to or greater than 1")
  if(!is.numeric(nstart))
    stop("Number of starts must be integer")
  if(nstart < 1)
    stop("Number of starts cannot be less than 1")
  if(nstart%%ceiling(nstart) > 0)
    nstart <- ceiling(nstart)
  if(!is.numeric(con.val))
    stop("Convergence value must be a number")
  if(con.val <= 0)
    stop("Convergence value can not be 0 or a negatif value")
  if(!is.numeric(iter.max))
    stop("Maximum number of iteration must be a positive integer")
  else
   iter.max <- ceiling(iter.max)
  if(iter.max <= 1)
    stop("Maximum number of iterations must be equal to or greater than 1")
  if(!is.logical(fixcent))
    stop("Argument 'fixcent' should be a TRUE or FALSE")
  if(!is.logical(fixmemb))
    stop("Argument 'fixmemb' should be a TRUE or FALSE")
  if(fixcent && fixmemb)
    stop("Arguments 'fixcent' and 'fixmemb' should not be a TRUE at the same time")
  if(!is.logical(stand))
    stop("Value of argument 'stand' should be TRUE or FALSE")
  if(stand){
    x <- scale(x, center = TRUE, scale = TRUE)[, ]
    centers <- scale(centers, center = TRUE, scale = TRUE)[, ]
  }
  func.val <- numeric(nstart)
  comp.time <- numeric(nstart)
  iter.num <- numeric(nstart)
  best.func <- Inf
  for(start.idx in 1:nstart){
    if(start.idx > 1){
      set.seed(as.integer(Sys.time()) + start.idx)
      if(!fixcent)
        centers <- eval(compv)
      if(!fixmemb) 
        memberships <- eval(compu)
      if(stand)
        centers <- scale(centers, center = TRUE, scale = TRUE)[, ]
      if(!missing(numseed))
        set.seed(numseed)
    }
    v0 <- v <- centers
    u <- memberships
    d <- matrix(nrow = n, ncol = k, 0)
    prevv<- v + 2*con.val
    iter <- 0
    cputime <- system.time(
      while((iter < iter.max) && (sum(abs(prevv - v)) > con.val)){
        iter <- iter + 1
        prevv <- v
        for(j in 1:k)
          for(i in 1:n)
            d[i,j] <- eval(compd)
        for(j in 1:k)
          for(i in 1:n)
            if(any(d[i,] == 0))
              u[i,] <- rep(1/k, k)
            else    
              u[i, j] <- 1/(sum((d[i, j]/d[i, ])^((2/dt)/(m - 1))))
             v <- t(u^m) %*% x / colSums(u^m)
      }
    )
    comp.time[start.idx] <- cputime[1]
    iter.num[start.idx] <- iter
    obj.func <- sum(d * (u^m))
    func.val[start.idx] <- obj.func
    if(obj.func < best.func){
      best.func <- obj.func
      best.u <- u
      best.v <- v
      best.d <- d
      best.start <- start.idx
    }
  }
  clabels <- crisp(best.u)
  csize <- numeric(k)
  for(i in 1:k)
    csize[i] <- sum(clabels==i)
  csumsqrs <- .sumsqr(x, best.v, clabels)
  if(is.null(rownames(x)))
    rnames <- paste(1:n)
  else
    rnames <- rownames(x)
  if(is.null(colnames(x)))
    cnames <- paste("p", 1:p, sep = "-")
  else
    cnames <- colnames(x)
  rownames(x) <- rnames; colnames(x) <- cnames
  rownames(best.v) <- paste("Cluster", 1:k, sep = " ")
  colnames(best.v) <- cnames
  colnames(v0) <- colnames(best.v)
  rownames(v0) <- rownames(best.v)
  rownames(best.u) <- rnames
  colnames(best.u) <- rownames(best.v)
  rownames(best.d) <- rnames
  colnames(best.d) <- rownames(best.v)
  names(clabels) <- paste(1:n, sep = " ")
  names(csize) <- paste(c(1:k), sep = " ")
  result = list()
    result$u <- best.u
    result$t <- NULL
    result$v <- best.v
    result$v0 <- v0
    result$d <- best.d
    result$f <- NULL
    result$x <- x
    result$cluster <- clabels
    result$csize <- csize
    result$sumsqrs <- csumsqrs
    result$k <- k
    result$m <- m
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
    result$inpargs[2] <- con.val
    result$inpargs[3] <- dmetric	
    result$inpargs[4] <- alginitv	
    result$inpargs[5] <- alginitu	
    result$inpargs[6] <- fixcent
    result$inpargs[7] <- fixmemb
    result$inpargs[8] <- stand
    names(result$inpargs) <- c("iter.max", "con.val", "dmetric", "alginitv", "alginitu", "fixcent", "fixmemb", "stand")
    result$algorithm <- "FCM"
    result$call <- match.call()
  class(result) <- c("ppclust")
  return(result)
}
