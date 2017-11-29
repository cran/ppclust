\name{ppclust2}
\alias{ppclust2}
\title{
Convert \sQuote{ppclust} objects to the other types of cluster objects
}
\description{
Converts an object of \sQuote{ppclust} class to the other types of cluster objects.
}
\usage{
ppclust2(objx, otype, \dots)
}

\arguments{
  \item{objx}{an object of \code{ppclust} class.}
  \item{otype}{target object class type for conversion.}
  \item{\dots}{additional arguments.}
}

\value{an object of \code{\link[cluster]{fanny.object}}, \code{\link[fclust]{summary.fclust}}, \code{\link[stats]{kmeans}} or \code{\link[vegclust]{vegclust}} class.
}

\author{
Zeynel Cebeci
}

\seealso{
 \code{\link{as.ppclust}},
 \code{\link{is.ppclust}},
 \code{\link{summary.ppclust}}
}

\examples{
data(iris)
# Create a object of ppclust
opc <- fcm(x=iris[,1:4], centers=3)

# Test the class of opc object
is.ppclust(opc)

# Convert ppclust object 'opc' to the fanny object
ofc <- ppclust2(opc, otype="fanny")

# Test the class of 'ofc' for ppclust
is.ppclust(ofc)

# Test the class of 'ofc'
class(ofc)

# Convert ppclust object 'opc' to fclust object
ofc <- ppclust2(opc, otype="fclust")

# Test the class of 'ofc'
class(ofc)
}
