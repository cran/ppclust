fpppcm <- function(x, centers, memberships, m=2, eta=2, K=1, omega, 
          dmetric="sqeuclidean", pw=2, alginitv="kmpp", alginitu="imembrand", 
          nstart=1, iter.max=1000, con.val=1e-09, 
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
  if(!is.matrix(centers))
    centers <- matrix(nrow=k, ncol=p, eval(compv))
  algsinitu <- match.arg(alginitu, inaparc::get.algorithms("membership"))
  alginitu <- match.arg(alginitu, algsinitu)
  if(alginitu=="imembrand")
    compu <- parse(text = paste0("inaparc::", alginitu, "(n,", k, ")$u"))
  else
    compu <- parse(text = paste0("inaparc::", alginitu, "(x,", k, ")$u"))
  compd <- parse(text = paste0(".compdist(x[i,], v[j,], dmetric='", dmetric, "', p=", pw, ")"))
  if(!is.matrix(centers))
    centers <- matrix(nrow=k, ncol=p, eval(compv))
  if(!missing(memberships)){
    if(is(memberships,"data.frame"))
      memberships <- as.matrix(memberships)
    if(!is(memberships, "matrix"))
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
    stop("The fuzziness exponent (m) should be a number greater than 1")
  if(!is.numeric(eta))
    stop("The typicality exponent (eta) should be a number")
  if(eta < 1)
    stop("The typicality exponent (eta) should be a number greater than 1")
  if(K < 1)
    stop("The coefficient K should be a greater than or equal to 1")
  if(!is.numeric(K))
    stop("The coefficient K should be a number")
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
  if(stand){
    x <- scale(x, center = TRUE, scale = TRUE)[, ]
    centers <- scale(centers, center = TRUE, scale = TRUE)[, ]
  }
  d <- matrix(nrow = n, ncol = k, 0)
  if(!missing(omega)){
    if(!is.vector(omega)||!is.numeric(omega))
      stop("Argument 'omega' must be a numeric vector")
    if(length(omega)!= p)
      stop(paste0("Argument 'omega' must have ", p, " elements"))
  }
  else{
    res.fcm <- fcm(x, centers)
    memberships <- res.fcm$u
    centers <- res.fcm$v
    omega <- comp.omega(d=res.fcm$d, u=memberships, m=eta, K=K)
  }
  omega.recomp <- FALSE
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
      omega.recomp <- TRUE
    }
    t <- matrix(nrow = n, ncol = k, 0)
    v0 <- v <- centers
    u <- memberships   
    prevu <- u + 2*con.val
    iter <- 0
    cputime <- system.time(
      while((iter < iter.max) && (con.val < sum(abs(prevu-u)))){
        iter <- iter + 1
        prevu <- u
        for(i in 1:n)
          for(j in 1:k)
            d[i,j] <- eval(compd)
        if(omega.recomp){
          res.fcm <- fcm(x, centers)
          memberships <- res.fcm$u
          centers <- res.fcm$v
          omega <- comp.omega(d=res.fcm$d, u=memberships, m=eta, K=K)
        }
        for(j in 1:k){
          for(i in 1:n){
            if(min(d[i,]) == 0){
              t[i,] <- rep(0, k)
              t[i,j] <- 1
            }
            else{
              t[i,j] <- 1/(1 + (d[i,j]/omega[j])^(1/(eta-1)))
            }
          }
        }
        for(j in 1:k){
          for(i in 1:n){
            if(min(d[i,]) == 0){
              u[i,] <- rep(0, k)
              u[i,j] <- 1
            }
            else{
              nom <- t[i,j]^eta * d[i,j] + omega[j] * (1-t[i,j])^eta
              denom <- 0
              for(l in 1:k)
                denom <- denom + t[i,l]^eta * d[i,l] + omega[l] * (1-t[i,l])^eta
              u[i,j] <- (1/nom^(1/(m-1))) / (1/denom^(1/(m-1)))
            }
          }
        }
        v <- t(u^m * t^eta) %*% x / colSums(u^m * t^eta)
      }
    )
    comp.time[start.idx] <- round(cputime[1],9)
    iter.num[start.idx] <- iter
    obj.func <- 0
    for(j in 1:k)
      obj.func <- obj.func + sum(u[,j]^m * (t[,j]^eta * d[,j] + (1-t[,j])^eta * omega[j]))
    func.val[start.idx] <- obj.func
    if(obj.func < best.func){
      best.func <- obj.func
      best.t <- t
      best.v <- v
      best.d <- d
      best.start <- start.idx
    }
  }
  clabels <- crisp(best.t)
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
  rownames(best.d) <- rnames
  colnames(best.d) <- rownames(best.v)
  names(clabels) <- paste(1:n, sep = " ")
  names(csize) <- paste(c(1:k), sep = " ")
  result = list()
    result$u <- NULL
    result$t <- best.t
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
    result$eta <- eta
    result$a <- NULL
    result$b <- NULL
    result$beta <- NULL
    result$delta <- NULL
    result$gamma <- NULL
    result$omega <- omega
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
    result$algorithm <- "PFCM"
    result$call <- match.call()
  class(result) <- c("ppclust")
  return(result)
}
