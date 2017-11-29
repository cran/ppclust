\name{x12}
\alias{x12}
\docType{data}

\title{
Synthetic data set of two variables
}

\description{
A synthetic data set described by Pal et al (2005). It consists of two continous variables forming two well-separated clusters in addition to two noise. 
}

\usage{data(x12)}
\format{A data frame with 12 rows and 2 numeric variables:
 \describe{
  \item{p1}{a numeric variable ranging from -5.0 to 5.0}
  \item{p2}{a numeric variable ranging from -1.67 to 10.0}
 }
}
\note{
The data set \code{x12} is recommended to test the performances of the possibilistic and probabilistic clustering algorithms.
}

\references{
Pal, N. R., Pal, K. & Bezdek, J. C. (2005). A possibilistic fuzzy c-means clustering algorithm. \emph{IEEE
Trans. Fuzzy Systems}, 13 (4): 517-530. <doi:10.1109/TFUZZ.2004.840099>
}

\examples{
data(x12)
# Descriptive statistics of the data set
summary(x12)
# Plot the data set
pairs(x12, pch=19, cex=2)
}

\concept{synthetic datasets}
\concept{clustering}
\concept{cluster analysis}

\keyword{datasets}
