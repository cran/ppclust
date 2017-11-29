\name{crisp}
\alias{crisp}
\title{
Crisp the fuzzy membership degrees
}
\description{
Crisps the fuzzy and possibilistic membership degrees from the fuzzy or possibilistic clustering algorithms.
}
\usage{
crisp(u, method, tv)
}

\arguments{
  \item{u}{a numeric matrix containing the fuzzy or possibilistic membership degrees of the data objects.}
  \item{method}{a string for selection of the crisping method. The default is \option{max} that assigns the data object to the cluster in which the object has maximum membership. The alternative is \option{threshold} that assigns the objects to a cluster if its maximum membership degree is greater than \code{tv}, a threshold value.}
  \item{tv}{a number for the threshold membership degree. The default is 0.5 with the \code{method} is \option{threshold} if it is not speficied by the user.}
}

\value{
  \item{cluster}{a numeric vector containing the indexes (labels) of clusters for the maximum membership of the objects.}
}

\details{
The function \code{crisp} produces the crisp or hard membership degrees of the objects in order to place them into only one cluster. 
}

\author{
Zeynel Cebeci
}

\examples{
data(iris)
x <- iris[,1:4]

# Run FCM 
res.fcm <- fcm(x, centers=3)

# Crisp the fuzzy memberships degrees and plot the crisp memberships
cllabels <- crisp(res.fcm$u)
plot(x, col=cllabels)
}

\concept{fuzzy c-means clustering}
\concept{prototype-based clustering}
\concept{partitioning clustering}
\concept{cluster analysis}

\keyword{cluster}
