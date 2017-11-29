\name{pcmr}
\alias{pcmr}
\title{
Possibilistic C-Means Clustering with Repulsion
}
\description{
Partitions a numeric data set by using the Possibilistic C-Means with Repulsion (PCMR) clustering algorithm which has been proposed by Wachs et al (2006).
}
\usage{
pcmr(x, centers, memberships, eta=2, K=1, omega, gamma=15,
    dmetric="sqeuclidean", pw=2, alginitv="kmpp", alginitu="imembrand",
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
  \item{gamma}{a number for normalization. Gamma value can be in the range of 0.1 and 200, but generally 10 is used. In Shapira & Wachs(2004)
gamma = 15 gave the best accuracy for PCMR.}
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
Possibilistic C-Means with Repulsion (PCMR) aims to  minimize  the  intracluster  distances while maximizing the intercluster distances without using implicitly the constraints of FCM, but by adding a cluster repulsion term to the objective function of PCM (Wachs et al, 2006).  

\eqn{J_{PCMR}(\mathbf{X}; \mathbf{V}, \mathbf{T})=\sum\limits_{i=1}^n t_{ij}^\eta \; d^2(\vec{x}_i, \vec{v}_j) + \sum\limits_{j=1}^k \Omega_j \sum\limits_{i=1}^n (1-t_{ij})^\eta + \gamma \sum\limits_{j=1}^k \sum\limits_{l=1, l \neq j}^k (1/d^2(\vec{v}_j, \vec{v}_l))}{J_{PCMR}(\mathbf{X}; \mathbf{V}, \mathbf{T})=\sum\limits_{i=1}^n t_{ij}^\eta \; d^2(\vec{x}_i, \vec{v}_j) + \sum\limits_{j=1}^k \Omega_j \sum\limits_{i=1}^n (1-t_{ij})^\eta + \gamma \sum\limits_{j=1}^k \sum\limits_{l=1, l \neq j}^k (1/d^2(\vec{v}_j, \vec{v}_l))}

Where \eqn{\gamma}{\gamma} is a weighting factor, and \eqn{t_{ij}}{t_{ij}} satisfies:

\eqn{t_{ij} \in [0,1], \forall j}{t_{ij} \in [0,1], \forall j}

The repulsion term is relevant if the clusters are close enough. When the distance increases it becomes smaller until it  is compensated by the  attraction of  the clusters. On the other  hand, if the clusters are sufficiently spread out, and the intercluster distance decreases (due to the first two terms), the attraction of the cluster can be compensated only by the repulsion term. 

The update equation for the cluster prototypes:

\eqn{\vec{v}_j =\frac{\sum\limits_{i=1}^n t_{ij} \vec{x}_i - \gamma \sum\limits_{j=1}^k v_j \; (1/ d^2(\vec{v}_j, \vec{v}_l))}{\sum\limits_{i=1}^n t_{ij} - \gamma \sum\limits_{j=1}^k v_j \; (1/ d^2(\vec{v}_j, \vec{v}_l))} \;;\; 1 \leq l \leq k}{\vec{v}_j =\frac{\sum\limits_{i=1}^n t_{ij} \vec{x}_i - \gamma \sum\limits_{j=1}^k v_j \; (1/ d^2(\vec{v}_j, \vec{v}_l))}{\sum\limits_{i=1}^n t_{ij} - \gamma \sum\limits_{j=1}^k v_j \; (1/ d^2(\vec{v}_j, \vec{v}_l))} \;;\; 1 \leq l \leq k}
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
   \item{gamma}{a number for normalization.}
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
Zeynel Cebeci, A. Tuna Kavlak, Figen Yildiz
}
\references{
Arthur, D. & Vassilvitskii, S. (2007). K-means++: The advantages of careful seeding, in \emph{Proc. of the 18th Annual ACM-SIAM Symposium on Discrete Algorithms}, p. 1027-1035. <\url{http://ilpubs.stanford.edu:8090/778/1/2006-13.pdf}>

Wachs, J., Shapira, O. & Stern, H. (2006). A method to enhance the 'Possibilistic C-Means with Repulsion' algorithm based on cluster validity index. In \emph{Applied Soft Computing Technologies: The Challenge of Complexity}, pp. 77-87. Springer, Berlin, Heidelberg. <doi:10.1007/3-540-31662-0_6>
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
 \code{\link{pfcm}},
 \code{\link{upfc}}
}

\examples{
# Load data set X12
data(x12)

# Initialize the prototype matrix using K-means++
v <- inaparc::kmpp(x12, k=2)$v
# Initialize the memberships degrees matrix 
u <- inaparc::imembrand(nrow(x12), k=2)$u

# Run FCM with the initial prototypes and memberships
fcm.res <- fcm(x12, centers=v, memberships=u, m=2)

# Run PCMR with the prototypes and memberships from FCM run
pcmr.res <- pcmr(x12, centers=fcm.res$v, memberships=fcm.res$u, eta=2)

# Show the typicality degrees for the top 5 objects
head(pcmr.res$t, 5)

# Plot the crisp memberships using maximum typicality degrees
plotcluster(pcmr.res, mt="t", cm="max")

# Plot the crisp memberships using the typicality degrees > 0.5
plotcluster(pcmr.res, mt="t", cm="threshold", tv=0.5)
}

\concept{possibilistic c-means clustering}
\concept{prototype-based clustering}
\concept{partitioning clustering}
\concept{cluster analysis}

\keyword{cluster}