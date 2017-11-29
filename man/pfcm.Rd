\name{pfcm}
\alias{pfcm}
\title{
Possibilistic Fuzzy C-Means Clustering Algorithm
}
\description{
Partitions a numeric data set by using the Possibilistic Fuzzy C-Means (PFCM) clustering algorithm proposed by Pal et al (2005).
}
\usage{
pfcm(x, centers, memberships, m=2, eta=2, K=1, omega, a, b,
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
  \item{a}{a number for the relative importance of the fuzzy part of the objective function. The default is 1.}
  \item{b}{a number for the relative importance of the possibilistic part of the objective function. The default is 1.}
  \item{K}{a number greater than 0 to be used as the weight of penalty term. The default is 1.}
  \item{omega}{a numeric vector of reference distances. If missing, it is internally generated.}
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
In FPCM, the constraint corresponding to the sum of all the typicality values of all data objects to a cluster must be equal to one causes problems; particularly for a big data set. In order to avoid this problem Pal et al (2005) proposed Possibilistic Fuzzy C-Means (PFCM) clustering algorithm with following objective function:

\eqn{J_{PFCM}(\mathbf{X}; \mathbf{V}, \mathbf{U}, \mathbf{T})=\sum\limits_{j=1}^k \sum\limits_{i=1}^n (a \; u_{ij}^m + b \; t_{ij}^\eta) \; d^2(\vec{x}_i, \vec{v}_j) + \sum\limits_{j=1}^k \Omega_j \sum\limits_{i=1}^n (1-t_{ij})^\eta}{J_{PFCM}(\mathbf{X}; \mathbf{V}, \mathbf{U}, \mathbf{T})=\sum\limits_{j=1}^k \sum\limits_{i=1}^n (a \; u_{ij}^m + b \; t_{ij}^\eta) \; d^2(\vec{x}_i, \vec{v}_j) + \sum\limits_{j=1}^k \Omega_j \sum\limits_{i=1}^n (1-t_{ij})^\eta}

The fuzzy membership degrees in the probabilistic part of the objective function \eqn{J_{PFCM}}{J_{PFCM}} is calculated in the same way as in FCM, as follows:

\eqn{u_{ij} =\Bigg[\sum\limits_{j=1}^k \Big(\frac{d^2(\vec{x}_i, \vec{v}_j)}{d^2(\vec{x}_i, \vec{v}_l)}\Big)^{1/(m-1)} \Bigg]^{-1} \;;\; 1 \leq i \leq n, \; 1 \leq l \leq k}{u_{ij} = \Bigg[\sum\limits_{j=1}^k \Big(\frac{d^2(\vec{x}_i, \vec{v}_j)}{d^2(\vec{x}_i, \vec{v}_l)}\Big)^{1/(m-1)} \Bigg]^{-1} \;;\; 1 \leq i \leq n, \; 1 \leq l \leq k}

The typicality degrees in the possibilistic part of the objective function \eqn{J_{PFCM}}{J_{PFCM}} is calculated as follows:

\eqn{t_{ij} =\Bigg[1 + \Big(\frac{b \; d^2(\vec{x}_i, \vec{v}_j)}{\Omega_j}\Big)^{1/(\eta -1)}\Bigg]^{-1} \;;\; 1 \leq i \leq n, \; 1 \leq j \leq k}{t_{ij} = \Bigg[1 + \Big(\frac{b \; d^2(\vec{x}_i, \vec{v}_j)}{\Omega_j}\Big)^{1/(\eta -1)}\Bigg]^{-1} \;;\; 1 \leq i \leq n, \; 1 \leq j \leq k}

The constraints with PFCM are:

\eqn{0 \leq u_{ij}, t_{ij} \leq 1}{0 \leq u_{ij}, t_{ij} \leq 1}

\eqn{0 \leq \sum\limits_{i=1}^n u_{ij} \leq n \;\;;\; 1 \leq j \leq k}{0 \leq \sum\limits_{i=1}^n u_{ij} \leq n \;\;;\; 1 \leq j \leq k}

\eqn{0 \leq \sum\limits_{j=1}^k t_{ij} \leq k \;\;;\; 1 \leq i \leq n}{0 \leq \sum\limits_{j=1}^k t_{ij} \leq k \;\;;\; 1 \leq i \leq n}

\eqn{\sum\limits_{j=1}^k u_{ij} = 1 \;\;;\; 1 \leq i \leq n}{\sum\limits_{j=1}^k u_{ij} = 1 \;\;;\; 1 \leq i\leq n}

\eqn{a}{a} and \eqn{b}{b} are the coefficients to define the relative importance of fuzzy membership and typicality degrees for weighting the probabilistic and possibilistic terms of the objective function, \eqn{ a > 0; \; b > 0}{ a > 0; \; b > 0}.

\eqn{m}{m} is the fuzzifier to specify the amount of fuzziness for the clustering; \eqn{1\leq m\leq \infty}. It is usually chosen as 2. 

\eqn{\eta}{\eta} is the typicality exponent to specify the amount of typicality for the clustering; \eqn{1\leq \eta\leq \infty}. It is usually chosen as 2. 

\eqn{\Omega}{\Omega} is the possibilistic penalty terms for controlling the variance of the clusters.

The update equation for cluster prototypes:

\eqn{\vec{v}_j =\frac{\sum\limits_{i=1}^n (a \; u_{ij}^m + b \; t_{ij}^\eta) \; \vec{x}_i}{\sum\limits_{i=1}^n (a \; u_{ij}^m + b \; t_{ij}^\eta)} \;;\; 1 \leq j \leq k}{\vec{v}_j =\frac{\sum\limits_{i=1}^n (a \; u_{ij}^m + b \; t_{ij}^\eta) \; \vec{x}_i}{\sum\limits_{i=1}^n (a \; u_{ij}^m + b \; t_{ij}^\eta)} \;;\; 1 \leq j \leq k}
}

\value{an object of class \sQuote{ppclust}, which is a list consists of the following items:
   \item{v}{a numeric matrix containing the final cluster prototypes.}
   \item{t}{a numeric matrix containing the typicality degrees of the data objects.}
   \item{d}{a numeric matrix containing the distances of objects to the final cluster prototypes.}
   \item{x}{a numeric matrix containing the processed data set.}
   \item{cluster}{a numeric vector containing the cluster labels found by defuzzifying the typicality degrees of the objects.}
   \item{csize}{a numeric vector for the number of objects in the clusters.}
   \item{k}{an integer for the number of clusters.}
   \item{m}{a number for the used fuzziness exponent.}
   \item{eta}{a number for the used typicality exponent.}
   \item{a}{a number for the fuzzy part of the objective function.}
   \item{b}{a number for the possibilistic part of the objective function.}
   \item{omega}{a numeric vector of reference distances.}
   \item{iter}{an integer vector for the number of iterations in each start of the algorithm.}
   \item{best.start}{an integer for the index of start that produced the minimum objective functional.}
   \item{func.val}{a numeric vector for the objective function values in each start of the algorithm.}
   \item{comp.time}{a numeric vector for the execution time in each start of the algorithm.}
   \item{stand}{a logical value, \code{TRUE} shows that \code{x} data set contains the standardized values of raw data.}
   \item{wss}{a number for the within-cluster sum of squares for each cluster.}
   \item{bwss}{a number for the between-cluster sum of squares.}
   \item{tss}{a number for the total within-cluster sum of squares.}
   \item{twss}{a number for the total sum of squares.}
   \item{algorithm}{a string for the name of partitioning algorithm. It is \sQuote{PCM} with this function.}
   \item{call}{a string for the matched function call generating this \sQuote{ppclust} object.}
}

\author{
Zeynel Cebeci, Alper Tuna Kavlak & Figen Yildiz
}

\references{
Arthur, D. & Vassilvitskii, S. (2007). K-means++: The advantages of careful seeding, in \emph{Proc. of the 18th Annual ACM-SIAM Symposium on Discrete Algorithms}, p. 1027-1035. <\url{http://ilpubs.stanford.edu:8090/778/1/2006-13.pdf}>

Pal, N. R., Pal, K. & Bezdek, J. C. (2005). A possibilistic fuzzy c-means clustering algorithm. \emph{IEEE
Trans. Fuzzy Systems}, 13 (4): 517-530. <doi:10.1109/TFUZZ.2004.840099>
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
 \code{\link{pca}},
 \code{\link{pcm}},
 \code{\link{pcmr}},
 \code{\link{upfc}}
}

\examples{
# Load the dataset X12
data(x12)

# Set the initial centers of clusters
v0 <- matrix(nrow=2, ncol=2, c(-3.34, 1.67, 1.67, 0.00), byrow=FALSE)

# Run FCM with the initial centers in v0
res.fcm <- fcm(x12, centers=v0, m=2)

# Run PFCM with the final centers and memberhips from FCM
res.pfcm <- pfcm(x12, centers=res.fcm$v, memberships=res.fcm$u, m=2, eta=2)

# Show the typicality and fuzzy membership degrees from PFCM
res.pfcm$t
res.pfcm$u
}

\concept{probabilistic clustering}
\concept{possibilistic clustering}
\concept{prototype-based clustering}
\concept{partitioning clustering}
\concept{cluster analysis}

\keyword{cluster}