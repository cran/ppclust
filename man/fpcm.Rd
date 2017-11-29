\name{fpcm}
\alias{fpcm}
\title{
Fuzzy Possibilistic C-Means Clustering
}
\description{
Partitions a numeric data set by using the Fuzzy and Possibilistic C-Means (FPCM) clustering algorithm (Pal et al, 1997).
}
\usage{
fpcm(x, centers, memberships, m=2, eta=2,  
     dmetric="sqeuclidean", pw=2, alginitv="kmpp", alginitu="imembrand", 
     nstart=1, iter.max=1000, con.val=1e-09, 
     fixcent=FALSE, fixmemb=FALSE, stand=FALSE, numseed)
}

\arguments{
  \item{x}{a numeric vector, data frame or matrix.}
  \item{centers}{an integer specifying the number of clusters or a numeric matrix containing the initial cluster centers.}
  \item{memberships}{a numeric matrix containing the initial membership degrees. If missing, it is internally generated.}
  \item{m}{a number greater than 1 to be used as the fuzziness exponent or fuzzifier. The default is 2.}
  \item{eta}{a number greater than 1 to be used as the typicality exponent. The default is 3.}
  \item{dmetric}{a string for the distance metric. The default is \option{sqeuclidean} for the squared Euclidean distances. See \code{\link{get.dmetrics}} for the alternative options.}
  \item{pw}{a number for the power of Minkowski distance calculation. The default is 2 if the \code{dmetric} is \option{minkowski}.}
  \item{alginitv}{a string for the initialization of cluster prototypes matrix. The default is \option{kmpp} for K-means++ initialization method (Arthur & Vassilvitskii, 2007). For the list of alternative options see \code{\link[inaparc]{get.algorithms}}.}
  \item{alginitu}{a string for the initialization of memberships degrees matrix. The default is \option{imembrand} for random sampling of initial membership degrees.}
  \item{nstart}{an integer for the number of starts for clustering. The default is 1.}
  \item{iter.max}{an integer for the maximum number of iterations allowed. The default is 1000.}
  \item{con.val}{a number for the convergence value between the iterations. The default is 1e-09.}
  \item{fixcent}{a logical flag to make the initial cluster centers not changed along the different starts of the algorithm. The default is \code{FALSE}. If it is \code{TRUE}, the initial centers are not changed in the successive starts of the algorithm when the \code{nstart} is greater than 1.}
  \item{fixmemb}{a logical flag to make the initial membership degrees not changed along the different starts of the algorithm. The default is \code{FALSE}. If it is \code{TRUE}, the initial memberships are not changed in the successive starts of the algorithm when the \code{nstart} is greater than 1.}
  \item{stand}{a logical flag to standardize data. Its default value is \code{FALSE}. If its value is \code{TRUE}, the data matrix \code{x} is standardized.}
  \item{numseed}{a seeding number to set the seed of R's random number generator.}
}

\details{
Fuzzy and Possibilistic C Means (FPCM) algorithm which has been proposed by Pal et al (1997) indended to combine the characteristics of FCM and PCM, and hence, was also so-called \dfn{Mixed C-Means} (MCM) algorithm.

The objective function of FPCM is:

\eqn{J_{FPCM}(\mathbf{X}; \mathbf{V}, \mathbf{U}, \mathbf{T})=\sum\limits_{i=1}^n (u_{ij}^m + t_{ij}^\eta) \; d^2(\vec{x}_i, \vec{v}_j)}{J_{FPCM}(\mathbf{X}; \mathbf{V}, \mathbf{U}, \mathbf{T})=\sum\limits_{i=1}^n (u_{ij}^m + t_{ij}^\eta) \; d^2(\vec{x}_i, \vec{v}_j)}

In the above equation:

\eqn{\mathbf{X} = \{\vec{x}_1, \vec{x}_2,\dots, \vec{x}_n\} \subseteq\Re^p}{\mathbf{X} = \{\vec{x}_1, \vec{x}_2,\dots, \vec{x}_n\} \subseteq\Re^p} is the data set for \eqn{n}{n} objects in the \emph{p}-dimensional data space \eqn{\Re}{\Re}, 

\eqn{\mathbf{V} = \{\vec{v}_1, \vec{v}_2, \dots, \vec{v}_k\} \subseteq\Re^n}{\mathbf{V} = \{\vec{v}_1, \vec{v}_2, \dots, \vec{v}_k\} \subseteq\Re^n} is the protoype matrix of the clusters,

\eqn{\mathbf{U} = \{u_{ij}\}}{\mathbf{U} = \{u_{ij}\}} is the matrix for a fuzzy partition of \eqn{\mathbf{X}}{\mathbf{X}}, 

\eqn{\mathbf{T} = \{t_{ij}\}}{\mathbf{T} = \{t_{ij}\}} is the matrix for a possibilistic partition of \eqn{\mathbf{X}}{\mathbf{X}}, 

\eqn{d^2(\vec{x}_i, \vec{v}_j)}{d^2(\vec{x}_i, \vec{v}_j)} is the squared Euclidean distance between the object \eqn{\vec{x}_j}{\vec{x}_j} and cluster prototype \eqn{\vec{v}_i}{\vec{v}_i}. 

\eqn{d^2(\vec{x}_i , \vec{v}_j) = ||\vec{x}_i - \vec{v}_j||^2 = (\vec{x}_i - \vec{v}_j)^T (\vec{x}_i - \vec{v}_j)}{d^2(\vec{x}_i , \vec{v}_j) = ||\vec{x}_i - \vec{v}_j||^2 = (\vec{x}_i - \vec{v}_j)^T (\vec{x}_i - \vec{v}_j)}

\eqn{m}{m} is the fuzzifier to specify the amount of fuzziness for the clustering; \eqn{1\leq m\leq \infty}. It is usually chosen as 2. 

\eqn{\eta}{\eta} is the typicality exponent to specify the amount of typicality for the clustering; \eqn{1\leq \eta\leq \infty}. It is usually chosen as 2. 

FPCM must satisfy the following constraints:

\eqn{\sum\limits_{j=1}^k u_{ij} = 1 \;\;;\; 1 \leq i\leq n}{\sum\limits_{j=1}^k u_{ij} = 1 \;\;;\; 1 \leq i\leq n}

\eqn{\sum\limits_{i=1}^n t_{ij} = 1 \;\;;\; 1 \leq j\leq k}{\sum\limits_{i=1}^n t_{ij} = 1 \;\;;\; 1 \leq j\leq k}

The objective function of FPCM is minimized by using the following update equations:

\eqn{u_{ij} =\Bigg[\sum\limits_{j=1}^k \Big(\frac{d^2(\vec{x}_i, \vec{v}_j)}{d^2(\vec{x}_i, \vec{v}_l)}\Big)^{1/(m-1)} \Bigg]^{-1} \;\;; 1 \leq i \leq n,\; 1 \leq l \leq k}{u_{ij} =\Bigg[\sum\limits_{j=1}^k \Big(\frac{d^2(\vec{x}_i, \vec{v}_j)}{d^2(\vec{x}_i, \vec{v}_l)}\Big)^{1/(m-1)} \Bigg]^{-1} \;\;; 1\leq i \leq n,\; 1 \leq l \leq k}

\eqn{t_{ij} =\Bigg[\sum\limits_{l=1}^n \Big(\frac{d^2(\vec{x}_i, \vec{v}_j)}{d^2(\vec{x}_i, \vec{v}_l)}\Big)^{1/(\eta-1)} \Bigg]^{-1} \;\;; 1 \leq i \leq n, \; 1 \leq j \leq k}{t_{ij} =\Bigg[\sum\limits_{l=1}^n \Big(\frac{d^2(\vec{x}_i, \vec{v}_j)}{d^2(\vec{x}_i, \vec{v}_l)}\Big)^{1/(\eta-1)} \Bigg]^{-1} \;\;; 1\leq i \leq n, \; 1 \leq j \leq k}

\eqn{\vec{v}_{j} =\frac{\sum\limits_{i=1}^n (u_{ij}^m + t_{ij}^\eta) \vec{x}_i}{\sum\limits_{i=1}^n (u_{ij}^m + t_{ij}^\eta)} \;\;; {1\leq j\leq k}}{\vec{v}_{j} =\frac{\sum\limits_{i=1}^n (u_{ij}^m + t_{ij}^\eta) \vec{x}_i}{\sum\limits_{i=1}^n (u_{ij}^m + t_{ij}^\eta)} \;\;; {1\leq j\leq k}}
}

\value{an object of class \sQuote{ppclust}, which is a list consists of the following items:
   \item{x}{a numeric matrix containing the processed data set.}
   \item{v}{a numeric matrix containing the final cluster prototypes (centers of clusters).}
   \item{u}{a numeric matrix containing the fuzzy memberships degrees of the data objects.}
   \item{d}{a numeric matrix containing the distances of objects to the final cluster prototypes.}
   \item{k}{an integer for the number of clusters.}
   \item{m}{a number for the fuzzifier.}
   \item{eta}{a number for the typicality exponent.}
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
Zeynel Cebeci, Alper Tuna Kavlak & Figen Yildiz
}
\references{
Arthur, D. & Vassilvitskii, S. (2007). K-means++: The advantages of careful seeding, in \emph{Proc. of the 18th Annual ACM-SIAM Symposium on Discrete Algorithms}, p. 1027-1035. <\url{http://ilpubs.stanford.edu:8090/778/1/2006-13.pdf}>

Pal, N.R., Pal, K., & Bezdek, J.C. (1997). A mixed c-means clustering model. In \emph{Proc. of the 6th IEEE Int. Conf. on Fuzzy Systems}, 1, pp. 11-21. <doi:10.1109/FUZZY.1997.616338>
}

\seealso{
 \code{\link{ekm}},
 \code{\link{fcm}},
 \code{\link{fcm2}},
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

# Initialize the prototype matrix using K-means++
v <- inaparc::kmpp(x, k=3)$v

# Initialize the memberships degrees matrix 
u <- inaparc::imembrand(nrow(x), k=3)$u

# Run FPCM with the initial prototypes and memberships
fpcm.res <- fpcm(x, centers=v, memberships=u, m=2, eta=2)

# Show the fuzzy membership degrees for the top 5 objects
head(fpcm.res$u, 5)

# Show the possibilistic membership degrees for the top 5 objects
head(fpcm.res$t, 5)
}

\concept{fuzzy possibilistic c-means clustering}
\concept{possibilistic c-means clustering}
\concept{mixed c-means clustering}
\concept{prototype-based clustering}
\concept{partitioning clustering}
\concept{cluster analysis}

\keyword{cluster}