\name{fpppcm}
\alias{fpppcm}
\title{
Fuzzy Possibilistic Product Partition C-Means Clustering
}
\description{
Partitions a numeric data set by using the Fuzzy Possibilistic Product Partition C-Means (FPPPCM) clustering algorithm which has been proposed by Szilagyi & Szilagyi (2014).
}
\usage{
fpppcm(x, centers, memberships, m=2, eta=2, K=1, omega, 
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
Fuzzy Possibilistic Product Partition C-Means (FPPPCM) clustering algorithm aimed to eliminate the effect of outliers in the other fuzzy and possibilistic clustering algorithms. The algorithm includes a probabilistic and a possibilistic term via multiplicative way instead of additive combination (Gosztolya & Szilagyi, 2015). The objective function of the algorithm as follows:

\eqn{J_{FPPPCM}(\mathbf{X}; \mathbf{V}, \mathbf{U}, \mathbf{T})=\sum\limits_{j=1}^k \sum\limits_{i=1}^n u_{ij}^m \big[ t_{ij}^\eta \; d^2(\vec{x}_i, \vec{v}_j) + \Omega_j (1-t_{ij})^\eta \big]}{J_{FPPPCM}(\mathbf{X}; \mathbf{V}, \mathbf{U}, \mathbf{T})=\sum\limits_{j=1}^k \sum\limits_{i=1}^n u_{ij}^m \big[ t_{ij}^\eta \; d^2(\vec{x}_i, \vec{v}_j) + \Omega_j (1-t_{ij})^\eta \big]}

The fuzzy membership degrees in the probabilistic part of the objective function \eqn{J_{FPPPCM}}{J_{FPPPCM}} is updated as follows:

\eqn{u_{ij} = \frac{\Big[t_{ij}^\eta \; d^2(\vec{x}_i, \vec{v}_j) \; + \; \Omega_j (1-t_{ij})^\eta \Big]^{-1/(m-1)}}{\Big[ \sum\limits_{l=1}^k t_{il}^\eta \; d^2(\vec{x}_i, \vec{v}_l) \; + \; \Omega_l (1-t_{il})^\eta \Big]^{-1/(m-1)}} \;;\; 1 \leq i \leq n, \; 1 \leq j \leq k}{u_{ij} = \frac{\Big[t_{ij}^\eta \; d^2(\vec{x}_i, \vec{v}_j) + \Omega_j (1-t_{ij})^\eta \Big]^{-1/(m-1)}}{\Big[ \sum\limits_{l=1}^k t_{il}^\eta \; d^2(\vec{x}_i, \vec{v}_l) + \Omega_l (1-t_{il})^\eta \Big]^{-1/(m-1)}} \;;\; 1 \leq i \leq n, \; 1 \leq j \leq k}

The typicality degrees in the possibilistic part of the objective function \eqn{J_{FPPPCM}}{J_{FPPPCM}} is calculated as follows:

\eqn{t_{ij} =\Bigg[1 + \Big(\frac{d^2(\vec{x}_i, \vec{v}_j)}{\Omega_j}\Big)^{1/(\eta -1)}\Bigg]^{-1} \;;\; 1 \leq i \leq n, \; 1 \leq j \leq k}{t_{ij} = \Bigg[1 + \Big(\frac{d^2(\vec{x}_i, \vec{v}_j)}{\Omega_j}\Big)^{1/(\eta -1)}\Bigg]^{-1} \;;\; 1 \leq i \leq n, \; 1 \leq j \leq k}

\eqn{m}{m} is the fuzzifier to specify the amount of fuzziness for the clustering; \eqn{1\leq m\leq \infty}. It is usually chosen as 2. 

\eqn{\eta}{\eta} is the typicality exponent to specify the amount of typicality for the clustering; \eqn{1\leq \eta\leq \infty}. It is usually chosen as 2. 

\eqn{\Omega}{\Omega} is the possibilistic penalty term to control the variance of the clusters.

The update equation for cluster prototypes:

\eqn{\vec{v}_j =\frac{\sum\limits_{i=1}^n u_{ij}^m \; t_{ij}^\eta \; \vec{x}_i}{\sum\limits_{i=1}^n u_{ij}^m \; t_{ij}^\eta} \;;\; 1 \leq j \leq k}{\vec{v}_j =\frac{\sum\limits_{i=1}^n u_{ij}^m \; t_{ij}^\eta \; \vec{x}_i}{\sum\limits_{i=1}^n u_{ij}^m \; t_{ij}^\eta} \;;\; 1 \leq j \leq k}
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

Szilagyi, L. & Szilagyi, S. M. (2014). Generalization rules for the suppressed fuzzy c-means clustering algorithm. \emph{Neurocomputing}, 139:298-309. <doi:10.1016/j.neucom.2014.02.027>

Gosztolya, G. & Szilagyi, L. (2015). Application of fuzzy and possibilistic c-means clustering models in blind speaker clustering. \emph{Acta Polytechnica Hungarica}, 12(7):41-56. <\url{http://publicatio.bibl.u-szeged.hu/6151/1/2015-acta-polytechnica.pdf}>
}

\seealso{
 \code{\link{ekm}},
 \code{\link{fcm}},
 \code{\link{fcm2}},
 \code{\link{fpcm}},
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
# Load dataset X16
data(x16)
x <- x16[,-3]
# Initialize the prototype matrix using K-means++
v <- inaparc::kmpp(x, k=2)$v
# Initialize the memberships degrees matrix 
u <- inaparc::imembrand(nrow(x), k=2)$u

# Run FPPPCM 
res.fpppcm <- fpppcm(x, centers=v, memberships=u, m=2, eta=2)

# Display typicality degrees 
res.fpppcm$t

# Run FPPPCM for eta=3
res.fpppcm <- fpppcm(x, centers=v, memberships=u, m=2, eta=3)

# Display typicality degrees 
res.fpppcm$t
}

\concept{probabilistic clustering}
\concept{possibilistic clustering}
\concept{prototype-based clustering}
\concept{partitioning clustering}
\concept{cluster analysis}

\keyword{cluster}