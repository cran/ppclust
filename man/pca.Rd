\name{pca}
\alias{pca}
\title{
Possibilistic Clustering Algorithm
}
\description{
Partitions a numeric data set by using the Possibilistic Clustering Algorithm (PCA) which has been proposed by Yang & Wu (2006).
}
\usage{
pca(x, centers, memberships, m=2, eta=2, 
   dmetric="sqeuclidean", pw=2, alginitv="kmpp", alginitu="imembrand", 
   nstart=1, iter.max=1000, con.val=1e-09, 
   fixcent=FALSE, fixmemb=FALSE, stand=FALSE, numseed)
}

\arguments{
  \item{x}{a numeric vector, data frame or matrix.}
  \item{centers}{an integer specifying the number of clusters or a numeric matrix containing the initial cluster centers.}
  \item{memberships}{a numeric matrix containing the initial membership degrees. If missing, it is internally generated.}
  \item{m}{a number greater than 1 to be used as the fuzziness exponent. The default is 2.}
  \item{eta}{a number greater than 1 to be used as the typicality exponent. The default is 2.}
  \item{dmetric}{a string for the distance metric. The default is \option{sqeuclidean} for the squared Euclidean distances. See \code{\link{get.dmetrics}} for the alternative options.}
  \item{pw}{a number for the power of Minkowski distance calculation. The default is 2 if the \code{dmetric} is \option{minkowski}.}
  \item{alginitv}{a string for the initialization of cluster prototypes matrix. The default is \option{kmpp} for K-means++ initialization method (Arthur & Vassilvitskii, 2007). For the list of alternative options see \code{\link[inaparc]{get.algorithms}}.}
  \item{alginitu}{a string for the initialization of memberships degrees matrix. The default is \option{imembrand} for random sampling of initial membership degrees.}
  \item{nstart}{an integer for the number of starts for clustering. The default is 1.}
  \item{iter.max}{an integer for the maximum number of iterations allowed. The default is 1000.}
  \item{con.val}{a number for the convergence value between the iterations. The default is 1e-09.}
  \item{fixcent}{a logical flag to fix the initial cluster centers. The default is \code{FALSE}. If it is \code{TRUE}, the initial centers are not changed in the successive starts of the algorithm when the \code{nstart} is greater than 1.}
  \item{fixmemb}{a logical flag to fix the initial membership degrees. The default is \code{FALSE}. If it is \code{TRUE}, the initial memberships are not changed in the successive starts of the algorithm when the \code{nstart} is greater than 1.}
  \item{stand}{a logical flag to standardize data. Its default value is \code{FALSE}. If its value is \code{TRUE}, the data matrix \code{x} is standardized.}
  \item{numseed}{a seeding number to set the seed of R's random number generator.}
}

\details{
Unlike the Possibilistic C-Means (PCM) algorithm requiring the results of a previous run of Fuzzy C-Means (FCM) clustering in order to calculate the parameter \eqn{\Omega}{\Omega}, Possibilistic Clustering Algorithm (PCA) is based on the FCM objective function, the partition coefficient (PC) and partition entropy (PE) validity indexes. So that PCA directly computes the typicality values and needs not run FCM beforehand to compute this parameter. The resulting membership becomes the exponential function, hence, it is reported that it is robust to noise and outliers (Yang & Wu, 2006). However, Wu et al (2010) reported that PCA is very sensitive to initializations and sometimes generates coincident clusters.

The objective function of PCA is:

\eqn{J_{PCA}(\mathbf{X}; \mathbf{V}, \mathbf{T})=\sum\limits_{j=1}^k \sum\limits_{i=1}^n t_{ij}^m \; d^2(\vec{x}_i, \vec{v}_j) + \frac{\beta}{m^2\sqrt{k}} \sum\limits_{j=1}^k \sum\limits_{i=1}^n (t_{ij}^m \; log \; t_{ij}^m - t_{ij}^m)}{J_{PCA}(\mathbf{X}; \mathbf{V}, \mathbf{T})=\sum\limits_{j=1}^k \sum\limits_{i=1}^n t_{ij}^m \; d^2(\vec{x}_i, \vec{v}_j) + \frac{\beta}{m^2\sqrt{k}} \sum\limits_{j=1}^k \sum\limits_{i=1}^n (t_{ij}^m \; log \; t_{ij}^m - t_{ij}^m)}

Where:

\eqn{t_{ij} = exp\Big(- \frac{m \sqrt{k} \; d^2(\vec{x}_i, \vec{v}_j)}{\beta}\Big) \;\;; {1\leq i\leq n},\; {1\leq j\leq k}}{t_{ij} = exp\Big(- \frac{m \sqrt{k} \; d^2(\vec{x}_i, \vec{v}_j)}{\beta}\Big) \;\;; {1\leq i\leq n},\; {1\leq j\leq k}}

The update equation for cluster prototypes:

\eqn{\vec{v}_{j} =\frac{\sum\limits_{i=1}^n t_{ij}^m \; \vec{x}_i}{\sum\limits_{i=1}^n t_{ij}^m} \;\;; {1\leq j\leq k}}{\vec{v}_{j} =\frac{\sum\limits_{i=1}^n t_{ij}^m \; \vec{x}_i}{\sum\limits_{i=1}^n t_{ij}^m} \;\;; {1\leq j\leq k}}

Where:

\eqn{\beta = \frac{\sum\limits_{i=1}^n \; d^2(\vec{x}_i, \overline{x})}{n}} with \eqn{\overline{x}=\frac{\sum\limits_{i=1}^n \vec{x}_i}{n}}
}

\value{an object of class \sQuote{ppclust}, which is a list consists of the following items:
   \item{v}{a numeric matrix containing the final cluster prototypes.}
   \item{u}{a numeric matrix containing the fuzzy membership degrees of the data objects.}
   \item{t}{a numeric matrix containing the typicality degrees of the data objects.}
   \item{d}{a numeric matrix containing the distances of objects to the final cluster prototypes.}
   \item{x}{a numeric matrix containing the processed data set.}
   \item{cluster}{a numeric vector containing the cluster labels found by defuzzifying the typicality degrees of the objects.}
   \item{csize}{a numeric vector containing the number of objects in the clusters.}
   \item{k}{an integer for the number of clusters.}
   \item{m}{a number for the fuziness exponent.}
   \item{eta}{a number for the typicality exponent.}
   \item{omega}{a numeric vector of reference distances.}
   \item{iter}{an integer vector for the number of iterations in each start of the algorithm.}
   \item{best.start}{an integer for the index of start that produced the minimum objective functional.}
   \item{func.val}{a numeric vector for the objective function values in each start of the algorithm.}
   \item{comp.time}{a numeric vector for the execution time in each start of the algorithm.}
   \item{stand}{a logical value, \code{TRUE} shows that data set \code{x} contains the standardized values of raw data.}
   \item{wss}{a number for the within-cluster sum of squares for each cluster.}
   \item{bwss}{a number for the between-cluster sum of squares.}
   \item{tss}{a number for the total within-cluster sum of squares.}
   \item{twss}{a number for the total sum of squares.}
   \item{algorithm}{a string for the name of partitioning algorithm. It is \sQuote{PCM} with this function.}
   \item{call}{a string for the matched function call generating this \sQuote{ppclust} object.}
}

\author{
Zeynel Cebeci, Alper Tuna Kavlak
}
\references{
Arthur, D. & Vassilvitskii, S. (2007). K-means++: The advantages of careful seeding, in \emph{Proc. of the 18th Annual ACM-SIAM Symposium on Discrete Algorithms}, p. 1027-1035. <\url{http://ilpubs.stanford.edu:8090/778/1/2006-13.pdf}>

Yang, M. S. & Wu, K. L. (2006). Unsupervised possibilistic clustering. \emph{Pattern Recognition}, 39(1): 5-21. <doi:10.1016/j.patcog.2005.07.005>
}

\seealso{
 \code{\link{ekm}},
 \code{\link{fcm}},
 \code{\link{fcm2}},
 \code{\link{fpcm}},
 \code{\link{fpppcm}},
 \code{\link{gg}},
 \code{\link{gk}},
 \code{\link{gkpfcm}},
 \code{\link{hcm}},
 \code{\link{pcm}},
 \code{\link{pcmr}},
 \code{\link{pfcm}},
 \code{\link{upfc}}
}

\examples{
# Load dataset X16
data(x16)
x <- x16[,-3]

# Initialize the prototype matrix using K-means++
v <- inaparc::kmpp(x, k=2)$v
# Initialize the membership degrees matrix 
u <- inaparc::imembrand(nrow(x), k=2)$u

# Run PCA
pca.res <- pca(x, centers=v, memberships=u, m=2, eta=2)

# Display the fuzzy membership degrees 
print(round(pca.res$u,2))

# Display the typicality degrees
print(round(pca.res$t,2))
}

\concept{possibilistic c-means clustering}
\concept{prototype-based clustering}
\concept{partitioning clustering}
\concept{cluster analysis}

\keyword{cluster}