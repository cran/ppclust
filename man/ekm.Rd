\name{ekm}
\alias{ekm}
\title{
K-Means Clustering Using Different Seeding Techniques
}
\description{
The function \code{ekm} partitions a numeric data set by using the K-means clustering algorithm. It is a wrapper function of the standard \code{\link[stats]{kmeans}} function with more initialization (seeding) techniques and output values obtained in the multiple starts of the algorithm. 
}
\usage{
ekm(x, centers, dmetric="euclidean", alginitv="hartiganwong",  
    nstart=1, iter.max=1000, stand=FALSE, numseed)
}

\arguments{
  \item{x}{a numeric vector, data frame or matrix.}
  \item{centers}{an integer specifying the number of clusters or a numeric matrix containing the initial cluster centers.}
  \item{dmetric}{a string for the distance calculation method. The default is \option{euclidean} for the Euclidean distances.}
  \item{alginitv}{a string for the initialization of cluster prototypes matrix. The default is \option{hartiganwong} for Hartigan-Wong seeding method. See \code{\link[inaparc]{get.algorithms}} for the alternative options.}
  \item{nstart}{an integer for the number of starts for clustering. The default is 1.}
  \item{iter.max}{an integer for the maximum number of iterations allowed. The default is 1000.}
  \item{stand}{a logical flag to standardize data. Its default value is \code{FALSE}. If its value is \code{TRUE}, the data matrix \code{x} is standardized.}
  \item{numseed}{a seeding number to set the seed of R's random number generator.}
}

\details{
K-Means (KM) clustering algorithm partitions a data set of \eqn{n}{n} objects into \eqn{k}{k}, a pre-defined number of clusters. It is an iterative relocation algorithm, and in each iteration step the objects are assigned to the nearest cluster by using the Euclidean distances. The objective of KM is to minimize total intra-cluster variance (or the sum of squared errors):

\eqn{J_{KM}(\mathbf{X}; \mathbf{V})=\sum\limits_{i=1}^n d^2(\vec{x}_i, \vec{v}_j)}{J_{KM}(\mathbf{X}; \mathbf{V})=\sum\limits_{i=1}^n d^2(\vec{x}_i, \vec{v}_j)}

In the above equation for \eqn{J_{KM}}{J_{KM}}:

\eqn{d^2(\vec{x}_i, \vec{v}_j)}{d^2(\vec{x}_i, \vec{v}_j)} is the distance measure between the object \eqn{\vec{x}_j}{\vec{x}_j} and cluster prototype \eqn{\vec{v}_i}{\vec{v}_i}.The Euclidean distance metric is usually employed with the implementations of K-means.

See \code{\link[stats]{kmeans}} and \code{\link{ppclust-package}} for more details about the terms of the objective function \eqn{J_{KM}}{J_{KM}}.

The update equation of the cluster prototypes:
  
\eqn{\vec{v}_{j} =\frac{1}{n_j} \sum\limits_{i=1}^{n_j} x_{ij} \;;\; 1 \leq j \leq k}{\vec{v}_{j} =\frac{1}{n_j} \sum\limits_{i=1}^{n_j} x_{ij} \;;\; 1 \leq j \leq k}
}

\value{an object of class \sQuote{ppclust}, which is a list consists of the following items:
   \item{x}{a numeric matrix containing the processed data set.}
   \item{v}{a numeric matrix containing the final cluster prototypes (centers of clusters).}
   \item{u}{a numeric matrix containing the crisp membership degrees of the data objects.}
   \item{d}{a numeric matrix containing the distances of objects to the final cluster prototypes.}
   \item{k}{an integer for the number of clusters.}
   \item{cluster}{a numeric vector containing the cluster labels found by defuzzying the fuzzy membership degrees of the objects.}
   \item{csize}{a numeric vector containing the number of objects in the clusters.}
   \item{iter}{an integer vector for the number of iterations in each start of the algorithm.}
   \item{best.start}{an integer for the index of start with the minimum objective functional.}
   \item{func.val}{a numeric vector for the objective function values in each start of the algorithm.}
   \item{comp.time}{a numeric vector for the execution time in each start of the algorithm.}
   \item{stand}{a logical value, \code{TRUE} shows that \code{x} data set contains the standardized values of raw data.}
   \item{wss}{a number for the within-cluster sum of squares for each cluster.}
   \item{bwss}{a number for the between-cluster sum of squares.}
   \item{tss}{a number for the total within-cluster sum of squares.}
   \item{twss}{a number for the total sum of squares.}
   \item{algorithm}{a string for the name of partitioning algorithm. It is \sQuote{KM} with this function.}
   \item{call}{a string for the matched function call generating this \sQuote{ppclust} object.}
}

\author{
Zeynel Cebeci, Figen Yildiz & Hasan Onder
}

\references{
MacQueen, J.B. (1967). Some methods for classification and analysis of multivariate observations. 
In \emph{Proc. of 5th Berkeley Symp. on Mathematical Statistics and Probability}, Berkeley, University of California Press, 1: 281-297. <\url{http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.308.8619&rep=rep1&type=pdf}>
}

\seealso{
 \code{\link[stats]{kmeans}},
 \code{\link{fcm}},
 \code{\link{fcm2}},
 \code{\link{fpcm}},
 \code{\link{fpppcm}},
 \code{\link{gg}},
 \code{\link{gk}},
 \code{\link{gkpfcm}},
 \code{\link{hcm}},
 \code{\link{pca}},
 \code{\link{pcm}},
 \code{\link{pcmr}},
 \code{\link{pfcm}},
 \code{\link{upfc}}
}

\examples{
# Load dataset iris 
data(iris)
x <- iris[,-5]

# Run EKM for 3 clusters
res.ekm <- ekm(x, centers=3)

# Print and plot the clustering result
print(res.ekm$cluster)
plot(x, col=res.ekm$cluster, pch=16)
}

\concept{hard c-means clustering}
\concept{k-means clustering}
\concept{prototype-based clustering}
\concept{partitioning clustering}
\concept{cluster analysis}

\keyword{cluster}
