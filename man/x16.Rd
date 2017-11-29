\name{x16}
\alias{x16}
\docType{data}

\title{
Synthetic data set of two variables forming two clusters 
}

\description{
A synthetic data set described by Yang & Wu (2006). It consists of two continous variables forming two well-separated clusters in addition to two noise points. 
}

\usage{data(x16)}
\format{A data frame with 16 rows and 2 numeric variables:
 \describe{
  \item{p1}{a numeric variable ranging from 50 to 150}
  \item{p2}{a numeric variable ranging from 145 to 200}
  \item{cl}{a numeric variable ranging from 1 to 4}
 }
}
\note{
The data set \code{x16} is recommended to test the performances of the possibilistic and noise clustering algorithms.
}

\references{
Yang, M. S. & Wu, K. L. (2006). Unsupervised possibilistic clustering. \emph{Pattern Recognition}, 39(1): 5-21. <doi:10.1016/j.patcog.2005.07.005>
}

\examples{
data(x16)
x <- x16[,-3]
# descriptive statistics of the variables
summary(x)
# scatter plots for the variable pairs
pairs(x, col=x16$cl, pch=20, cex=2)
}

\concept{synthetic datasets}
\concept{clustering}
\concept{cluster analysis}

\keyword{datasets}
