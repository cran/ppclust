\name{summary.ppclust}
\alias{summary.ppclust}
\title{
Summarize the clustering results
}
\description{
Summarizes the clustering results for an object which is an instance of \sQuote{ppclust} class.}

\usage{
\method{summary}{ppclust}(object, \dots)
}

\arguments{
  \item{object}{an object of \code{ppclust} class to be summarized.}
  \item{\dots}{additional arguments for S3 method summary.}
}

\value{summary of the clustering results from the object of \code{ppclust} class.}

\author{
Zeynel Cebeci
}

\seealso{
 \code{\link{as.ppclust}},
 \code{\link{is.ppclust}},
 \code{\link{ppclust2}}
}

\examples{
data(iris)
# Run FCM for three clusters
res.fcm <- fcm(x=iris[,1:4], centers=3)

# Summarize the result
summary(res.fcm)
}

\concept{possibilistic c-means clustering}
\concept{prototype-based clustering}
\concept{partitioning clustering}
\concept{cluster analysis}

\keyword{cluster}