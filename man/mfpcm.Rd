\name{mfpcm}
\alias{mfpcm}
\title{
Modified Fuzzy Possibilistic C-Means Clustering
}
\description{
Partitions a numeric data set by using the Modified Fuzzy and Possibilistic C-Means (MFPCM) clustering algorithm (Saad & Alimi, 2009).
}
\usage{
mfpcm(x, centers, memberships, m=2, eta=2,  
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
Modified Fuzzy and Possibilistic C Means (MFPCM) algorithm was proposed by Pal et al (1997) intented to incorporate a weight parameter to the objective function of FPCM as follows:

\eqn{J_{MFPCM}(\mathbf{X}; \mathbf{V}, \mathbf{U}, \mathbf{T})=\sum\limits_{i=1}^n u_{ij}^m w_{ij}^m \; d^{2m}(\vec{x}_i, \vec{v}_j) + t_{ij}^\eta w_{ij}^\eta \; d^{2\eta}(\vec{x}_i, \vec{v}_j)}{J_{MFPCM}(\mathbf{X}; \mathbf{V}, \mathbf{U}, \mathbf{T})=\sum\limits_{i=1}^n u_{ij}^m w_{ij}^m \; d^{2m}(\vec{x}_i, \vec{v}_j) + t_{ij}^\eta w_{ij}^\eta \; d^{2\eta}(\vec{x}_i, \vec{v}_j)}

In the above ojective function, every data object is considered to has its own weight in relation to every cluster. Therefore it is expected that the weight permits to have a better classification especially in the case of noise data (Saad & Alimi, 2009). The weight is calculated with the following equation: 

\eqn{w_{ij} = exp \Bigg[- \frac{d^2(\vec{x}_i, \vec{v}_j)}{\sum\limits_{i=1}^n d^2(\vec{x}_i, \bar{v}) \frac{k}{n}} \Bigg]}{w_{ij} = exp \Bigg[- \frac{d^2(\vec{x}_i, \vec{v}_j)}{\sum\limits_{i=1}^n d^2(\vec{x}_i, \bar{v}) \frac{k}{n}} \Bigg]}

The objective function of MFPCM is minimized by using the following update equations:

\eqn{u_{ij} =\Bigg[\sum\limits_{j=1}^k \Big(\frac{d^2(\vec{x}_i, \vec{v}_j)}{d^2(\vec{x}_i, \vec{v}_l)}\Big)^{2m/(m-1)} \Bigg]^{-1} \;\;; 1 \leq i \leq n,\; 1 \leq l \leq k}{u_{ij} =\Bigg[\sum\limits_{j=1}^k \Big(\frac{d^2(\vec{x}_i, \vec{v}_j)}{d^2(\vec{x}_i, \vec{v}_l)}\Big)^{2m/(m-1)} \Bigg]^{-1} \;\;; 1 \leq i \leq n,\; 1 \leq l \leq k}

\eqn{t_{ij} =\Bigg[\sum\limits_{l=1}^n \Big(\frac{d^2(\vec{x}_i, \vec{v}_j)}{d^2(\vec{x}_i, \vec{v}_l)}\Big)^{2\eta/(\eta-1)} \Bigg]^{-1} \;\;; 1 \leq i \leq n, \; 1 \leq j \leq k}{t_{ij} =\Bigg[\sum\limits_{l=1}^n \Big(\frac{d^2(\vec{x}_i, \vec{v}_j)}{d^2(\vec{x}_i, \vec{v}_l)}\Big)^{2\eta/(\eta-1)} \Bigg]^{-1} \;\;; 1\leq i \leq n, \; 1 \leq j \leq k}

\eqn{\vec{v}_{j} =\frac{\sum\limits_{i=1}^n (u_{ij}^m w_{ij}^m + t_{ij}^\eta w_{ij}^\eta) \vec{x}_i}{\sum\limits_{i=1}^n (u_{ij}^m w_{ij}^m + t_{ij}^\eta) w_{ij}^\eta} \;\;; {1\leq j\leq k}}{\vec{v}_{j} =\frac{\sum\limits_{i=1}^n (u_{ij}^m w_{ij}^m + t_{ij}^\eta w_{ij}^\eta) \vec{x}_i}{\sum\limits_{i=1}^n (u_{ij}^m w_{ij}^m + t_{ij}^\eta) w_{ij}^\eta} \;\;; {1\leq j\leq k}}
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

Saad, M. F. & Alimi, A. M. (2009). Modified fuzzy possibilistic c-means. In \emph{Proc. of the Int.  Multiconference of Engineers and Computer Scientists}, 1: 18-20. <ISBN:978-988-17012-2-0>
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

# Run FCM with the initial prototypes and memberships
mfpcm.res <- mfpcm(x, centers=v, memberships=u, m=2, eta=2)

# Show the fuzzy membership degrees for the top 5 objects
head(mfpcm.res$u, 5)

# Show the possibilistic membership degrees for the top 5 objects
head(mfpcm.res$t, 5)
}

\concept{fuzzy c-means clustering}
\concept{possibilistic c-means clustering}
\concept{mixed c-means clustering}
\concept{prototype-based clustering}
\concept{partitioning clustering}
\concept{cluster analysis}

\keyword{cluster}