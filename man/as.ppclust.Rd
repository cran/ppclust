\name{as.ppclust}
\alias{as.ppclust}
\title{
Convert object to \sQuote{ppclust} class
}
\description{
Converts an object of the classes \code{\link[cluster]{fanny.object}}, \code{\link[fclust]{summary.fclust}}, \code{\link[stats]{kmeans}} or \code{\link[vegclust]{vegclust}} to \sQuote{ppclust} class.
}
\usage{
as.ppclust(objx, ...)
}

\arguments{
  \item{objx}{an object to be converted to an instance of \code{ppclust} class.}
  \item{...}{additional arguments.}
}

\value{an object of \code{ppclust} class.}

\author{
Zeynel Cebeci
}

\seealso{
 \code{\link{is.ppclust}},
 \code{\link{ppclust2}},
 \code{\link{summary.ppclust}}
}

\examples{
data(iris)
# Create an fclust object
ofc <- fclust::FKM(X=iris[,1:4], k=3)

# Test the class of object 'ofc' 
class(ofc)

# Convert 'ofc' to ppclust object
opc <- as.ppclust(ofc)

# Test the class of 'opc' object
class(opc)
}
