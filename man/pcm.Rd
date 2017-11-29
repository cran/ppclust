\name{pcm}
\alias{pcm}
\title{
Possibilistic C-Means Clustering
}
\description{
Partitions a numeric multidimensional data set by using the Possibilistic C-Means (PCM) clustering algorithm which has been proposed by Krishnapuram & Keller (1993, 1996).
}
\usage{
pcm(x, centers, memberships, eta=2, K=1, omega, oftype = 1,
    dmetric="sqeuclidean", pw=2, fcmrun=TRUE, 
    alginitv="kmpp", alginitu="imembrand", 
    nstart=1, iter.max=1000, con.val=1e-09, 
    fixcent=FALSE, fixmemb=FALSE, stand=FALSE, numseed)
}

\arguments{
  \item{x}{a numeric vector, data frame or matrix.}
  \item{centers}{an integer specifying the number of clusters or a numeric matrix containing the initial cluster centers.}
  \item{memberships}{a numeric matrix containing the initial membership degrees. If missing, it is internally generated.}
  \item{eta}{a number greater than 1 to be used as the typicality exponent. The default is 2.}
  \item{K}{a number greater than 0 to be used as the weight of penalty term. The default is 1.}
  \item{omega}{a numeric vector of reference distances. If missing, it is internally generated.}
  \item{oftype}{an integer for the type of objective function. The default option is \option{1} for the PCM objective function in Krishnapuram & Keller (1993). Use \option{2} for PCM objective function in Krishnapuram & Keller (1996).}
  \item{dmetric}{a string for the distance metric. The default is \option{sqeuclidean} for the squared Euclidean distances. See \code{\link{get.dmetrics}} for the alternative options.}
  \item{pw}{a number for the power of Minkowski distance calculation. The default is \option{2} if the \code{dmetric} is \option{minkowski}.}
  \item{alginitv}{a string for the initialization of cluster prototypes matrix. The default is \option{kmpp} for K-means++ initialization method (Arthur & Vassilvitskii, 2007). For the list of alternative options see \code{\link[inaparc]{get.algorithms}}.}
  \item{alginitu}{a string for the initialization of memberships degrees matrix. The default is \option{imembrand} for random sampling of initial membership degrees.}
  \item{fcmrun}{a logical value for using the results from a previous FCM run. The default is \option{TRUE} for starting the algorithm with the initial centers and memberships from FCM run. Set the value to \option{FALSE} to cancel to use the results of FCM run.}
  \item{nstart}{an integer for the number of starts for clustering. The default is 1.}
  \item{iter.max}{an integer for the maximum number of iterations allowed. The default is 1000.}
  \item{con.val}{a number for the convergence value between the iterations. The default is 1e-09.}
  \item{fixcent}{a logical value to fix the initial cluster centers. The default is \code{FALSE}. If it is \code{TRUE}, the initial centers are not changed in the successive starts of the algorithm when the \code{nstart} is greater than 1.}
  \item{fixmemb}{a logical value to fix the initial membership degrees. The default is \code{FALSE}. If it is \code{TRUE}, the initial memberships are not changed in the successive starts of the algorithm when the \code{nstart} is greater than 1.}
  \item{stand}{a logical value to standardize data. The default value is \code{FALSE}. If its value is \code{TRUE}, the data matrix \code{x} is standardized.}
  \item{numseed}{a seeding number to set the seed of R's random number generator.}
}

\details{
Possibilistic C-Means (PCM) clustering algorithm was proposed by Krishnapuram and Keller (1993) to overcome the noise problem of FCM (Nasraoui, 1999). The algorithm differs from the other clustering algorithms in that the resulting partition of the data can be interpreted as a possibilistic partition, and the membership values can be interpreted as degrees of possibility of the points belonging to the classes, i.e., the compatibilities of the points with the class prototypes (Krishnapuram & Keller, 1993). PCM relaxes the probabilistic constraint on the sum of the memberships of an object to all clusters in FCM. However, PCM must run on the fuzzy clustering results of FCM in order to calculate the parameter \eqn{\Omega}{\Omega}. Although it solves the noise sensitivity problem of FCM, the performance of PCM depends heavily on the initialization and often deteriorates due to the coincident clustering problem (Filippone et al, 2007). 

The objective function of PCM is:

\eqn{J_{PCM}(\mathbf{X}; \mathbf{V}, \mathbf{T})=\sum\limits_{i=1}^n t_{ij}^\eta \; d^2(\vec{x}_i, \vec{v}_j) + \sum\limits_{j=1}^k \Omega_j \sum\limits_{i=1}^n (1 - t_{ij})^\eta}{J_{PCM}(\mathbf{X}; \mathbf{V}, \mathbf{T})=\sum\limits_{i=1}^n t_{ij}^\eta \; d^2(\vec{x}_i, \vec{v}_j) + \sum\limits_{j=1}^k \Omega_j \sum\limits_{i=1}^n (1 - t_{ij})^\eta}

In the objective function above, the first term leads to a minimization of the weighted distances while the second term suppresses the trivial solution (Timm et al, 2004). 

Krishnapuram and Keller later proposed an alternative objective function for PCM (Krishnapuram & Keller, 1996):

\eqn{J_{PCM_2}(\mathbf{X}; \mathbf{V}, \mathbf{T})=\sum\limits_{i=1}^n t_{ij}^\eta \; d^2(\vec{x}_i, \vec{v}_j) + \sum\limits_{j=1}^k \Omega_j \sum\limits_{i=1}^n (t_{ij}^\eta \; \log{t_{ij}^\eta}  - t_{ij}^\eta)}{J_{PCM_2}(\mathbf{X}; \mathbf{V}, \mathbf{T})=\sum\limits_{i=1}^n t_{ij}^\eta \; d^2(\vec{x}_i, \vec{v}_j) + \sum\limits_{j=1}^k \Omega_j \sum\limits_{i=1}^n (t_{ij}^\eta \; \log{t_{ij}^\eta}  - t_{ij}^\eta)}

Where:

\eqn{\vec{\Omega} = K \sum\limits_{i=1}^n u_{ij}^m \; d^2(\vec{x}_i, \vec{v}_j) / \sum\limits_{i=1}^n u_{ij}^m}{\vec{\Omega} = K \sum\limits_{i=1}^n u_{ij}^m \; d^2(\vec{x}_i, \vec{v}_j) / \sum\limits_{i=1}^n u_{ij}^m}

Since the membership calculation in PCM is possibilistic, the membership degrees can be treated as typicality values which measure how typical a data object is for a given cluster, regardless of all other clusters. The update equation of typicality degrees remains identical to those of FCM, and is derived from the objective function of PCM is:

\eqn{t_{ij} =\Bigg[1 + \Big(\frac{d^2(\vec{x}_i, \vec{v}_j)}{\Omega_j}\Big)^{(1/(m-1)}\Bigg]^{-1} \;\;; 1 \leq i \leq  n,\; 1 \leq j \leq k}{t_{ij} =\Bigg[1 + \Big(\frac{d^2(\vec{x}_i, \vec{v}_j)}{\Omega_j}\Big)^{(1/(m-1)}\Bigg]^{-1} \;\;; 1 \leq i \leq  n,\; 1 \leq j \leq k}

The update equation for cluster prototypes is the same with those of FCM:

\eqn{\vec{v}_{j} =\frac{\sum\limits_{i=1}^n t_{ij}^m \vec{x}_i}{\sum\limits_{i=1}^n t_{ij}^m} \;\;; 1 \leq j \leq k}{\vec{v}_{j} =\frac{\sum\limits_{i=1}^n t_{ij}^m \vec{x}_i}{\sum\limits_{i=1}^n t_{ij}^m} \;\;; 1 \leq j \leq k}
}

\value{an object of class \sQuote{ppclust}, which is a list consists of the following items:
   \item{v}{a numeric matrix containing the final cluster prototypes.}
   \item{t}{a numeric matrix containing the typicality degrees of the data objects.}
   \item{d}{a numeric matrix containing the distances of objects to the final cluster prototypes.}
   \item{x}{a numeric matrix containing the processed data set.}
   \item{cluster}{a numeric vector containing the cluster labels found by defuzzifying the typicality degrees of the objects.}
   \item{csize}{a numeric vector containing the number of objects in the clusters.}
   \item{k}{an integer for the number of clusters.}
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

\references{
Arthur, D. & Vassilvitskii, S. (2007). K-means++: The advantages of careful seeding, in \emph{Proc. of the 18th Annual ACM-SIAM Symposium on Discrete Algorithms}, p. 1027-1035. <\url{http://ilpubs.stanford.edu:8090/778/1/2006-13.pdf}>

Krishnapuram, R. & Keller, J.M. (1993). A possibilistic approach to clustering. \emph{IEEE Transactions on Fuzzy Systems}, 1(2):98-110. <doi:10.1109/91.227387>

Krishnapuram, R. & Keller, J.M. (1996). The possibilistic c-means algorithm: insights and recommendations. \emph{IEEE Transactions on Fuzzy Systems}, 4(3):385-393. <doi:10.1109/91.531779>

Timm, H., Borgelt, C., Doring, C. & Kruse, R. (2004). An extension to possibilistic fuzzy cluster analysis. \emph{Fuzzy Sets and Systems}, 147 (1): 3-16. <doi:10.1016/j.fss.2003.11.009>

Filippone, M., Masulli, F., & Rovetta, S. (2007). Possibilistic clustering in feature space. In \emph{Int. Workshop on Fuzzy Logic and Applications}, pp. 219-226. Springer, Berlin, Heidelberg. <doi:10.1007/978-3-540-73400-0_27>
}

\note{PCM usually leads to acceptable results, although it suffers from stability problems if it is not initialized with the corresponding probabilistic algorithm. Therefore, it needs fuzzy partitioning results of a probabilistic algorithm, i.e. FCM, in order to compute some input arguments such as the cluster prototypes and mobilization scale parameter (Timm et al, 2004).
}

\author{
Zeynel Cebeci, Alper Tuna Kavlak
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
 \code{\link{pcmr}},
 \code{\link{pfcm}},
 \code{\link{upfc}}
}

\examples{
# Load the dataset X12
data(x12)

# Initialize the prototype matrix using K-means++
v <- inaparc::kmpp(x12, k=2)$v
# Initialize the memberships degrees matrix 
u <- inaparc::imembrand(nrow(x12), k=2)$u

# Run FCM with the initial prototypes and memberships
fcm.res <- fcm(x12, centers=v, memberships=u, m=2)

# Run PCM with the prototypes and memberships from FCM run
pcm.res <- pcm(x12, centers=fcm.res$v, memberships=fcm.res$u, eta=2)

# Display the typicality degrees
print(pcm.res$t)

# Plot the crisp memberships using maximum typicality degrees
plotcluster(pcm.res, mt="t", trans=TRUE )

# Plot the crisp memberships using the typicality degrees > 0.5
plotcluster(pcm.res, mt="t", cm="threshold", tv=0.5)
}

\concept{possibilistic c-means clustering}
\concept{prototype-based clustering}
\concept{partitioning clustering}
\concept{cluster analysis}

\keyword{cluster}
