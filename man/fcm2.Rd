\name{fcm2}
\alias{fcm2}
\title{
Type-2 Fuzzy C-Means Clustering
}
\description{
Partitions a numeric data set by using the Type-2 Fuzzy C-Means (FCM2) clustering algorithm (Rhee & Hwang, 2001). It has been reported that it is effective for spherical clusters in the data sets, but it fails when the data sets contain non-spherical and complex structures (Gosain & Dahiya, 2016).
}
\usage{
fcm2(x, centers, memberships, m=2, dmetric="sqeuclidean", pw = 2, 
    alginitv="kmpp", alginitu="imembrand", 
    nstart=1, iter.max=1000, con.val=1e-09, 
    fixcent=FALSE, fixmemb=FALSE, stand=FALSE, numseed)
}

\arguments{
  \item{x}{a numeric vector, data frame or matrix.}
  \item{centers}{an integer specifying the number of clusters or a numeric matrix containing the initial cluster centers.}
  \item{memberships}{a numeric matrix containing the initial membership degrees. If missing, it is internally generated.}
  \item{m}{a number greater than 1 to be used as the fuzziness exponent or fuzzifier. The default is 2.}
  \item{dmetric}{a string for the distance metric. The default is \option{sqeuclidean} for the squared Euclidean distances. See \code{\link{get.dmetrics}} for the alternative options.}
  \item{pw}{a number for the power of Minkowski distance calculation. The default is 2 if the \code{dmetric} is \option{minkowski}.}
  \item{alginitv}{a string for the initialization of cluster prototypes matrix. The default is \option{kmpp} for K-means++ initialization method (Arthur & Vassilvitskii, 2007). For the list of alternative options, see \code{\link[inaparc]{get.algorithms}}.}
  \item{alginitu}{a string for the initialization of memberships degrees matrix. The default is \option{imembrand} for random sampling of initial membership degrees.}
  \item{nstart}{an integer for the number of starts for clustering. The default is 1.}
  \item{iter.max}{an integer for the maximum number of iterations allowed. The default is 1000.}
  \item{con.val}{a number for the convergence value between the iterations. The default is 1e-09.}
  \item{fixcent}{a logical flag to make the initial cluster centers not changed along the different starts of the algorithm. The default is \code{FALSE}. If it is \code{TRUE}, the initial centers are not changed in the successive starts of the algorithm when the \code{nstart} is greater than 1.}
  \item{fixmemb}{a logical flag to make the initial membership degrees not changed along the different starts of the algorithm. The default is \code{FALSE}. If it is \code{TRUE}, the initial memberships are not changed in the successive starts of the algorithm when the \code{nstart} is greater than 1.}
  \item{stand}{a logical flag to standardize data. Its default value is \code{FALSE}. If its value is \code{TRUE}, the data matrix \code{x} is standardized.}
  \item{numseed}{an optional seeding number to set the seed of R's random number generator.}
}

\details{
In the Type-2 Fuzzy C-Means (T2FCM) clustering algorithm proposed by (Rhee & Hwang, 2001), the idea is that all data points should not have the same contribution in computing the cluster prototypes, instead the data points with higher membership value should contribute much more in the cluster prototypes (Gosain & Dahiya, 2016). Based on this idea, a type-2 membership value is computed by using the following equation:

\eqn{a_{ij} =  u_{ij}^m - \frac{1 - u_{ij}^m}{2}}{a_{ij} =  u_{ij}^m - \frac{1 - u_{ij}^m}{2}}

In the above equation, \eqn{u_{ij}}{u_{ij}} and \eqn{a_{ij}}{a_{ij}} respectively stand for the Type-1 (standard)  and Type-2 membership values. \eqn{a_{ij}}{a_{ij}} is only used to update the cluster prototypes.

The objective function of FCM2 and other terms are the same with those of FCM except the following update equation for the cluster prototypes:

\eqn{\vec{v}_{j} =\frac{\sum\limits_{i=1}^n a_{ij}^m \vec{x}_i}{\sum\limits_{i=1}^n a_{ij}^m} \;\;; 1 \leq j\leq k}{\vec{v}_{j} =\frac{\sum\limits_{i=1}^n a_{ij}^m \vec{x}_i}{\sum\limits_{i=1}^n a_{ij}^m} \;\;; 1 \leq j\leq k}

}

\value{an object of class \sQuote{ppclust}, which is a list consists of the following items:
   \item{x}{a numeric matrix containing the processed data set.}
   \item{v}{a numeric matrix containing the final cluster prototypes (centers of clusters).}
   \item{u}{a numeric matrix containing the fuzzy memberships degrees of the data objects.}
   \item{d}{a numeric matrix containing the distances of objects to the final cluster prototypes.}
   \item{k}{an integer for the number of clusters.}
   \item{m}{a number for the fuzzifier.}
   \item{cluster}{a numeric vector containing the cluster labels found by defuzzying the fuzzy membership degrees of the objects.}
   \item{csize}{a numeric vector containing the number of objects in the clusters.}
   \item{iter}{an integer vector for the number of iterations in each start of the algorithm.}
   \item{best.start}{an integer for the index of start that produced the minimum objective functional.}
   \item{func.val}{a numeric vector for the objective function values in each start of the algorithm.}
   \item{comp.time}{a numeric vector for the execution time in each start of the algorithm.}
   \item{stand}{a logical value, \code{TRUE} shows that data set \code{x} contains the standardized values of raw data.}
   \item{wss}{a number for the within-cluster sum of squares for each cluster.}
   \item{bwss}{a number for the between-cluster sum of squares.}
   \item{tss}{a number for the total within-cluster sum of squares.}
   \item{twss}{a number for the total sum of squares.}
   \item{algorithm}{a string for the name of partitioning algorithm. It is \sQuote{FCM} with this function.}
   \item{call}{a string for the matched function call generating this \sQuote{ppclust} object.}
}

\author{
Zeynel Cebeci & Alper Tuna Kavlak
}

\references{
Arthur, D. & Vassilvitskii, S. (2007). K-means++: The advantages of careful seeding, in \emph{Proc. of the 18th Annual ACM-SIAM Symposium on Discrete Algorithms}, p. 1027-1035. <\url{http://ilpubs.stanford.edu:8090/778/1/2006-13.pdf}>

Rhee, F.C.H. & Hwang, C. (2001). A type-2 fuzzy c-means clustering algorithm. In <\emph{IEEE 9th IFSA World Congress and 20th NAFIPS Conf. 2001.}, 4:1926-1929. <doi:10.1109/NAFIPS.2001.944361>

Gosain, A. & Dahiya, S. (2016). Performance analysis of various fuzzy clustering algorithms: A review. \emph{Procedia Comp. Sci.}, 79:100-111. <doi:10.1016/j.procs.2016.03.014>
}

\seealso{
 \code{\link{ekm}},
 \code{\link{fcm}},
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
# Load dataset X12
data(x12)

# Initialize the prototype matrix using K-means++ algorithm
v <- inaparc::kmpp(x12, k=2)$v
# Initialize the membership degrees matrix 
u <- inaparc::imembrand(nrow(x12), k=2)$u

# Run FCM2 with the initial prototypes and memberships
fcm2.res <- fcm2(x12, centers=v, memberships=u, m=2)

# Show the fuzzy membership degrees for the top 5 objects
head(fcm2.res$u, 5)
}

\concept{type-2 fuzzy c-means clustering}
\concept{prototype-based clustering}
\concept{partitioning clustering}
\concept{cluster analysis}

\keyword{cluster}