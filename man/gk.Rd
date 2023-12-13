\name{gk}
\alias{gk}
\title{
Gustafson-Kessel Clustering
}
\description{
Partitions a numeric data set by using the Gustafson-Kessel (GK) clustering algorithm (Gustafson & Kessel, 1979). Unlike FCM using the Euclidean distance, GK uses cluster specific Mahalanobis distance.
}
\usage{
gk(x, centers, memberships, m=2,  
   dmetric="sqeuclidean", pw = 2, alginitv="kmpp", 
   alginitu="imembrand", nstart=1, iter.max=1e03, con.val=1e-09, 
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
  \item{numseed}{a seeding number to set the seed of R's random number generator.}
}

\details{
As an extension of basic FCM algorithm, Gustafson and Kessel (GK) clustering algorithm employs an adaptive distance norm in order to detect clusters with different geometrical shapes (Babuska, 2001; Correa et al, 2011). 

The objective function of GK is:

\eqn{J_{GK}(\mathbf{X}; \mathbf{V}, \mathbf{A}, \mathbf{U}) = \sum\limits_{i=1}^n \sum\limits_{j=1}^k  u_{ij}^m d_{A_j}(\vec{x}_i, \vec{v}_j)}{J_{GK}(\mathbf{X}; \mathbf{V}, \mathbf{A}, \mathbf{U}) = \sum\limits_{i=1}^n \sum\limits_{j=1}^k  u_{ij}^m d_{A_j}(\vec{x}_i, \vec{v}_j)}

In the above equation \eqn{d_{A_j}(\vec{x}_i, \vec{v}_j)}{d_{A_j}(\vec{x}_i, \vec{v}_j)} is the Mahalanobis distance between the data object \eqn{\vec{x}_j}{\vec{x}_j} and cluster prototype \eqn{\vec{v}_i}{\vec{v}_i}.

\eqn{d_{A_j}(\vec{x}_i, \vec{v}_j) = (\vec{x}_i - \vec{v}_j)^T \mathbf{A}_j (\vec{x}_i - \vec{v}_j)}{d_{A_j}(\vec{x}_i, \vec{v}_j) = (\vec{x}_i - \vec{v}_j)^T \mathbf{A}_j (\vec{x}_i - \vec{v}_j)} 

As it is understood from the above equation, each cluster has its own norm-inducing matrix in GK, so that the norm matrix \eqn{\mathbf{A}}{\mathbf{A}} is a \emph{k}-length tuples of the cluster-specific norm-inducing matrices:

\eqn{\mathbf{A}=\{\mathbf{A}_1, \mathbf{A}_2, \dots, \mathbf{A}_k\}}{\mathbf{A}=\{\mathbf{A}_1, \mathbf{A}_2, \dots, \mathbf{A}_k\}}. 

The objective function of GK cannot be directly minimized with respect to \eqn{\mathbf{A}_j}{\mathbf{A}_j} since it is
linear in \eqn{\mathbf{A}_j}{\mathbf{A}_j}. For obtaining a feasible solution, \eqn{\mathbf{A}_j}{\mathbf{A}_j} must be constrained in some way. One of the usual ways of accomplishing this is to constrain the determinant of \eqn{\mathbf{A}_j}{\mathbf{A}_j} (Babuska, 2001):

\eqn{\mid\mathbf{A}_j\mid=\rho_j \;,\; \rho_j > 0 \;,\; 1 \leq j \leq k}{\mid\mathbf{A}_j\mid=\rho_j \;,\; \rho_j > 0 \;,\; 1 \leq j \leq k}

Allowing the matrix \eqn{\mathbf{A}_j}{\mathbf{A}_j} to vary with its determinant fixed corresponds to optimizing the shape of cluster while its volume remains constant. By using the Lagrange-multiplier method, the norm-inducing matrix for the cluster \eqn{j}{j} is defined as follows (Babuska, 2001):

\eqn{\mathbf{A}_j = [\rho_i \; det(\mathbf{F}_j)]^{-1/n}{\mathbf{F}_j}^{-1}}{\mathbf{A}_j = [\rho_i \; det(\mathbf{F}_j)]^{-1/n}{\mathbf{F}_j}^{-1}},

where:

\eqn{n}{n} is the number of data objects,

\eqn{\rho_j}{\rho_j} represents the volume of cluster \eqn{j}{j},

\eqn{\mathbf{F}_j}{\mathbf{F}_j} is the fuzzy covariance matrix for the cluster \eqn{j}{j}, and is calculated as follows:

\eqn{\mathbf{F}_j = \frac{\sum\limits_{i=1}^n u_{ij}^m \; d^2(\vec{x}_i, \vec{v}_j)}{\sum\limits_{i=1}^n u_{ij}^m}}{\mathbf{F}_j=\frac{\sum\limits_{i=1}^n u_{ij}^m \; d^2(\vec{x}_i, \vec{v}_j)}{\sum\limits_{i=1}^n u_{ij}^m}}

\eqn{m}{m} is the fuzzifier to specify the amount of fuzziness for the clustering; \eqn{1\leq m\leq \infty}. It is usually chosen as 2. 

GK must satisfy the following constraints:

\eqn{\sum\limits_{j=1}^k u_{ij} = 1 \;\;;\; 1 \leq i \leq n}{\sum\limits_{j=1}^k u_{ij} = 1 \;\;;\; 1 \leq i\leq n}

\eqn{\sum\limits_{i=1}^n u_{ij} > 0 \;\;;\; 1 \leq j \leq k}{\sum\limits_{i=1}^n u_{ij} > 0 \;\;;\; 1 \leq j \leq k}

The objective function of GK is minimized by using the following update equations:

\eqn{u_{ij} =\Bigg[\sum\limits_{j=1}^k \Big(\frac{d_{A_j}(\vec{x}_i, \vec{v}_j)}{d_{A_l}(\vec{x}_i, \vec{v}_l)}\Big)^{1/(m-1)} \Bigg]^{-1} \;\;; 1 \leq i \leq n,\; 1 \leq l \leq k}{u_{ij} =\Bigg[\sum\limits_{j=1}^k \Big(\frac{d_{A_j}(\vec{x}_i, \vec{v}_j)}{d_{A_l}(\vec{x}_i, \vec{v}_l)}\Big)^{1/(m-1)} \Bigg]^{-1} \;\;; 1 \leq i \leq n,\; 1 \leq l \leq k}

\eqn{\vec{v}_{j} =\frac{\sum\limits_{i=1}^n u_{ij}^m \vec{x}_i}{\sum\limits_{i=1}^n u_{ij}^m} \;\;; 1 \leq j \leq k}{\vec{v}_{j} =\frac{\sum\limits_{i=1}^n u_{ij}^m \vec{x}_i}{\sum\limits_{i=1}^n u_{ij}^m} \;\;; 1 \leq j \leq k}
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

\references{
Arthur, D. & Vassilvitskii, S. (2007). K-means++: The advantages of careful seeding, in \emph{Proc. of the 18th Annual ACM-SIAM Symposium on Discrete Algorithms}, pp. 1027-1035. \url{http://ilpubs.stanford.edu:8090/778/1/2006-13.pdf}

Gustafson, D. E. & Kessel, W. C. (1979). Fuzzy clustering with a fuzzy covariance matrix. In \emph{Proc. of IEEE Conf. on Decision and Control including the 17th Symposium on Adaptive Processes}, San Diego. pp. 761-766. \doi{10.1109/CDC.1978.268028}

Babuska, R. (2001). Fuzzy and neural control. DISC Course Lecture Notes. Delft University of Technology. Delft, the Netherlands. \href{https://tr.scribd.com/document/209211977/Fuzzy-and-Neural-Control}{https://tr.scribd.com/document/209211977/Fuzzy-and-Neural-Control}.

Correa, C., Valero, C., Barreiro, P., Diago, M. P., & Tardáguila, J. (2011). A comparison of fuzzy clustering algorithms applied to feature extraction on vineyard. In \emph{Proc. of the 14th Conf. of the Spanish Assoc. for Artificial Intelligence}. \url{http://oa.upm.es/9246/}
}

\author{
Zeynel Cebeci & Cagatay Cebeci
}

\seealso{
 \code{\link{ekm}},
 \code{\link{fcm}},
 \code{\link{fcm2}},
 \code{\link{fpcm}},
 \code{\link{fpppcm}},
 \code{\link{gg}},
 \code{\link{gkpfcm}},
 \code{\link{hcm}},
 \code{\link{pca}},
 \code{\link{pcm}},
 \code{\link{pcmr}},
 \code{\link{pfcm}},
 \code{\link{upfc}}
}

\examples{
\dontrun{
# Load dataset iris 
data(iris)
x <- iris[,-5]

# Initialize the prototype matrix using K-means++
v <- inaparc::kmpp(x, k=3)$v

# Initialize the membership degrees matrix 
u <- inaparc::imembrand(nrow(x), k=3)$u

# Run FCM with the initial prototypes and memberships
gk.res <- gk(x, centers=v, memberships=u, m=2)

# Show the fuzzy membership degrees for the top 5 objects
head(gk.res$u, 5)
}
}

\concept{fuzzy clustering}
\concept{gustafson-kessel clustering}
\concept{prototype-based clustering}
\concept{partitioning clustering}
\concept{cluster analysis}

\keyword{cluster}
