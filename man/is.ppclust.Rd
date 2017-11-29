\name{is.ppclust}
\alias{is.ppclust}
\title{
Check the class of object for \sQuote{ppclust}
}
\description{
Checks the class of given object whether it is an instance of the \code{ppclust} class or not. 
}
\usage{
is.ppclust(objx)
}

\arguments{
  \item{objx}{an object to be checked for its class.}
}

\value{\code{TRUE} if \code{objx} is a valid \code{ppclust} object and \code{FALSE} for the other types of object classes.
}

\author{
Zeynel Cebeci
}

\seealso{
 \code{\link{as.ppclust}},
 \code{\link{ppclust2}},
 \code{\link{summary.ppclust}}
}

\examples{
data(iris)

# Run FCM for 3 clusters
res.fcm <- fcm(x=iris[,1:4], centers=2)

# Test for a ppclust object returned by the fcm function
is.ppclust(res.fcm)

# Test for a matrix object
x <- matrix(nrow=3, ncol=2, c(1,4,2,5,7,8))
is.ppclust(x)
}
