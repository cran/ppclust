\name{fcm}
\alias{fcm}
\title{
Fuzzy C-Means Clustering
}
\description{
Partitions a numeric data set by using the Fuzzy C-Means (FCM) clustering algorithm (Bezdek, 1974;1981).
}
\usage{
fcm(x, centers, memberships, m=2, dmetric="sqeuclidean", pw = 2, 
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
  \item{alginitv}{a string for the initialization of cluster prototypes matrix. The default is \option{kmpp} for K-means++ initialization method (Arthur & Vassilvitskii, 2007). For the list of alternative options see \code{\link[inaparc]{get.algorithms}}.}
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
Fuzzy C-Means (FCM) clustering algorithm was firstly studied by Dunn (1973) and generalized by Bezdek in 1974 (Bezdek, 1981). Unlike K-means algorithm, each data object is not the member of only one cluster but is the member of all clusters with varying degrees of memberhip between 0 and 1. It is an iterative clustering algorithm that partitions the data set into a predefined \emph{k} partitions by minimizing the weighted within group sum of squared errors. The objective function of FCM is:

\eqn{J_{FCM}(\mathbf{X}; \mathbf{V}, \mathbf{U})=\sum\limits_{i=1}^n u_{ij}^m d^2(\vec{x}_i, \vec{v}_j)}{J_{FCM}(\mathbf{X}; \mathbf{V}, \mathbf{U})=\sum\limits_{i=1}^n u_{ij}^m d^2(\vec{x}_i, \vec{v}_j)}

In the objective function, \eqn{m}{m} is the fuzzifier to specify the amount of 'fuzziness' of the clustering result; \eqn{1 \leq m \leq \infty}. It is usually chosen as 2. The higher values of \eqn{m}{m} result with the more fuzzy clusters while the lower values give harder clusters. If it is 1, FCM becomes an hard algorithm and produces the same results with K-means.

FCM must satisfy the following constraints:

\eqn{u_{ij}=[0,1] \;\;;\; 1 \leq i\leq n \;, 1 \leq j\leq k}{u_{ij}=[0,1] \;\;;\; 1 \leq i\leq n \;, 1 \leq j\leq k} 

\eqn{0 \leq \sum\limits_{i=1}^n u_{ij} \leq n \;\;;\; 1 \leq j\leq k}{0 \leq \sum\limits_{i=1}^n u_{ij} \leq n \;\;;\; 1 \leq j\leq k}

\eqn{\sum\limits_{j=1}^k u_{ij} = 1 \;\;;\; 1 \leq i\leq n}{\sum\limits_{j=1}^k u_{ij} = 1 \;\;;\; 1 \leq i\leq n}

The objective function of FCM is minimized by using the following update equations:

\eqn{u_{ij} =\Bigg[\sum\limits_{j=1}^k \Big(\frac{d^2(\vec{x}_i, \vec{v}_j)}{d^2(\vec{x}_i, \vec{v}_l)}\Big)^{1/(m-1)} \Bigg]^{-1} \;\;; {1\leq i\leq n},\; {1\leq l \leq k}}{u_{ij} =\Bigg[\sum\limits_{j=1}^k \Big(\frac{d^2(\vec{x}_i, \vec{v}_j)}{d^2(\vec{x}_i, \vec{v}_k)}\Bil)^{1/(m-1)} \Bigg]^{-1} \;\;; {1\leq i\leq n},\; {1\leq l \leq k}} 

\eqn{\vec{v}_{j} =\frac{\sum\limits_{i=1}^n u_{ij}^m \vec{x}_i}{\sum\limits_{i=1}^n u_{ij}^m} \;\;; {1\leq j\leq k}}{\vec{v}_{j} =\frac{\sum\limits_{i=1}^n u_{ij}^m \vec{x}_i}{\sum\limits_{i=1}^n u_{ij}^m} \;\;; {1\leq j\leq k}}
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
Zeynel Cebeci, Figen Yildiz & Alper Tuna Kavlak
}

\references{
Arthur, D. & Vassilvitskii, S. (2007). K-means++: The advantages of careful seeding, in \emph{Proc. of the 18th Annual ACM-SIAM Symposium on Discrete Algorithms}, pp. 1027-1035. <{http://ilpubs.stanford.edu:8090/778/1/2006-13.pdf}>

Dunn, J.C. (1973). A fuzzy relative of the ISODATA process and its use in detecting compact well-separated clusters. \emph{J. Cybernetics}, 3(3):32-57. <doi:10.1080/01969727308546046>

Bezdek, J.C. (1974). Cluster validity with fuzzy sets. \emph{J. Cybernetics}, 3: 58-73. <doi:10.1080/01969727308546047>

Bezdek J.C. (1981). \emph{Pattern recognition with fuzzy objective function algorithms}. Plenum, NY. <\href{https://isbnsearch.org/isbn/0306406713}{ISBN:0306406713}>
}

\seealso{
 \code{\link{ekm}},
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

# Initialize the prototype matrix using K-means++ algorithm
v <- inaparc::kmpp(x, k=3)$v

# Initialize the memberships degrees matrix 
u <- inaparc::imembrand(nrow(x), k=3)$u

# Run FCM with the initial prototypes and memberships
fcm.res <- fcm(x, centers=v, memberships=u, m=2)

# Show the fuzzy membership degrees for the top 5 objects
head(fcm.res$u, 5)
}

\concept{fuzzy c-means clustering}
\concept{prototype-based clustering}
\concept{partitioning clustering}
\concept{cluster analysis}

\keyword{cluster}