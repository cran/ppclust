\name{get.dmetrics}
\alias{get.dmetrics}
\title{
List the names of distance metrics
}
\description{
Displays the distance metric options for calculation of the distances between the data objects and the cluster centers.
}
\usage{
get.dmetrics(dmt="all")
}
\arguments{
  \item{dmt}{a string for the type of distance metrics. The default is \option{all} for the list of all of the distance metrics which are available in the package. The other options are \option{l1}, \option{l2}, \option{lp}, \option{sl2} and \option{sc} for L1, L2, Lp, squared L2 and squared Chord, respectively.}
}

\value{
  \item{dmlist}{a two-column matrix containing the distance metrics and their descriptions.}
}

\references{
Podani, J. (2000). \emph{Introduction to the Exploration of Multivariate Biological Data}. Backhuys Publishers, Leiden, The Netherlands. 407 pages. <\href{https://isbnsearch.org/isbn/9057820676}{ISBN 90-5782-067-6}>

McCune, B., & Grace, J. B. (2002). \emph{Analysis of ecological communities}. Gleneden Beach, Oregon: MjM Software Design. <ISBN:0972129006>
}

\seealso{
 \code{\link{ppclust-package}}
}

\author{
Zeynel Cebeci
}

\examples{
# Get all distance metrics
dmlist <- get.dmetrics(dmt="all")
dmlist

# Get only L2 type distance metrics
dmlist <- get.dmetrics(dmt="l2")
dmlist
}

\concept{distance}
\concept{partitioning clustering}
\concept{cluster analysis}

\keyword{cluster}
