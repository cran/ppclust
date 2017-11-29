\name{gg}
\alias{gg}
\title{Gath-Geva Clustering Algorithm}
\description{
Partitions a numeric data set by using the Gath-Geva (GG) clustering algorithm (Gath & Geva, 1989). The function \code{gg} is based on the code in \code{\link[advclust]{fuzzy.GG}} function of the package \pkg{advclust} by Bagus and Pramana (2016) with some more extended initialization methods and distance metrics.
}
\usage{
gg(x, centers, memberships, m=2, ggversion="simple", 
   dmetric="sqeuclidean", pw = 2, alginitv="kmpp", 
   alginitu="imembrand", nstart=1, iter.max=1e03, con.val=1e-09, 
   fixcent=FALSE, fixmemb=FALSE, stand=FALSE, numseed)
}

\arguments{
  \item{x}{a numeric vector, data frame or matrix.}
  \item{centers}{an integer specifying the number of clusters or a numeric matrix containing the initial cluster centers.}
  \item{memberships}{a numeric matrix containing the initial membership degrees. If missing, it is internally generated.}
  \item{ggversion}{a string for the version of Gath-Geva algorithm. The default is \option{simple} for the simplified version (Hoepner et al (1999). Use \option{original} for the original Gath-Geva algorithm (Gath & Geva, 1989).}
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
Gath and Geva (1989) proposed that the fuzzy maximum likelihood estimates (FMLE) clustering algorithm can be used to detect clusters of varying shapes, sizes and densities. Instead of using Euclidean distance as in FCM, Gath-Geva (GG) algorithm uses a distance norm based on fuzzy maximum likelihood estimates as described by Bezdek & Dunn (1975). These are the Gauss distances, and hence, GG is also so-called Gaussian Mixture Decomposition.  

The objective function of GG is:

\eqn{J_{GG}(\mathbf{X}; \mathbf{V}, \mathbf{A}, \mathbf{U}) = \sum\limits_{i=1}^n \sum\limits_{j=1}^k  u_{ij}^m d_{A_j}(\vec{x}_i, \vec{v}_j)}{J_{GG}(\mathbf{X}; \mathbf{V}, \mathbf{A}, \mathbf{U}) = \sum\limits_{i=1}^n \sum\limits_{j=1}^k  u_{ij}^m d_{A_j}(\vec{x}_i, \vec{v}_j)}

In the above equation, \eqn{d_{A_j}(\vec{x}_i, \vec{v}_j)}{d_{A_j}(\vec{x}_i, \vec{v}_j)} is the Gauss distance between the object \eqn{\vec{x}_j}{\vec{x}_j} and cluster prototype \eqn{\vec{v}_i}{\vec{v}_i}. 

\eqn{d_{A_j}(\vec{x}_i, \vec{v}_j) = \frac{(2 \pi)^{\frac{n}{2}}\sqrt{\det{\mathbf{A}_j}}}{\alpha_j} \; exp\big(\frac{1}{2} (\vec{x}_i - \vec{v}_j)^T \mathbf{A}_j^{-1}(\vec{x}_i - \vec{v}_j)\big)}{d_{A_j}(\vec{x}_i, \vec{v}_j) = \frac{(2 \pi)^{\frac{n}{2}}\sqrt{\det{\mathbf{A}_j}}}{\alpha_j} \; exp\big(\frac{1}{2} (\vec{x}_i - \vec{v}_j)^T \mathbf{A}_j^{-1}(\vec{x}_i - \vec{v}_j)\big)}

\eqn{\mathbf{A}_j = \frac{\sum\limits_{i=1}^n u_{ij}^m (\vec{x}_i - \vec{v}_j)^T (\vec{x}_i - \vec{v}_j)}{\sum\limits_{i=1}^n u_{ij}^m} \;;\; 1 \leq j \leq k}{\mathbf{A}_j = \frac{\sum\limits_{i=1}^n u_{ij}^m (\vec{x}_i - \vec{v}_j)^T (\vec{x}_i - \vec{v}_j)}{\sum\limits_{i=1}^n u_{ij}^m} \;;\; 1 \leq j \leq k}

The argument \eqn{\alpha_j}{\alpha_j} is the prior probability for belonging \eqn{\vec{x}_i}{\vec{x}_i} to the cluster \eqn{j}{j}:

\eqn{{\alpha_j} = \frac{\sum\limits_{i=1}^n u_{ij}^m}{n}}{{\alpha_j} = \ frac{\sum\limits_{i=1}^n u_{ij}^m}{n}}

The argument \eqn{\alpha_j}{\alpha_j} takes into account the influence of the size of a cluster, and bigger clusters attract more elements (Cruz, 2013).

\eqn{m}{m} is the fuzzifier to specify the amount of fuzziness for the clustering. Although it is 1 in the original FMLE algorithm (Bezdek & Dunn, 1975) a higher value of it (\eqn{1\leq m\leq \infty}{1\leq m\leq \infty}) can be used to make the partition more fuzzy compensating the exponential term of the distance norm (Balasko et al, 2005).  

The objective function of GG is minimized by using the following update equations:

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
Arthur, D. & Vassilvitskii, S. (2007). K-means++: The advantages of careful seeding, in \emph{Proc. of the 18th Annual ACM-SIAM Symposium on Discrete Algorithms}, p. 1027-1035. <\url{http://ilpubs.stanford.edu:8090/778/1/2006-13.pdf}>

Bezdek, J.C. & Dunn J.C. (1975). Optimal fuzzy partitions: A heuristic for estimating the parameters in a mixture of normal dustrubutions. \emph{IEEE Transactions on Computers}, C-24(8):835-838. <doi: 10.1109/T-C.1975.224317>

Gath, I. & Geva, A.B. (1989). Unsupervised optimal fuzzy clustering. \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence}, 11 (7): 773-781. <doi:10.1109/34.192473>

Hoeppner, F., Klawonn, F., Kruse, R. & Runkler, T. (1999). \emph{Fuzzy cluster analysis}. New York: John Wiley and Sons. <ISBN:0471988642>

Balasko, B., Abonyi, J. & Feil, B. (2005). Fuzzy clustering and data analysis toolbox. Department of Process Eng., Univ. of Veszprem, Veszprem.

Cruz, A. (2013). Cluster Algorithms. <\url{http://equipe.nce.ufrj.br/adriano/fuzzy/transparencias/clusteralgorithms/clusteralgorithms.pdf}>

Bagus, A. F. & Pramana, S. (2016). advclust: Object Oriented Advanced Clustering. R package version 0.4. \url{https://CRAN.R-project.org/package=advclust}
}

\note{
Due to the exponential distance norm, this algorithm needs a good initialization since it converges to a near local optimum. So, usually another clustering algorithm, i.e. FCM, is used to initialize the partition matrix \eqn{\mathbf{U}}{\mathbf{U}} (Balasko et al, 2005; Cruz, 2013).
}

\author{
Zeynel Cebeci
}

\seealso{
 \code{\link{ekm}},
 \code{\link{fcm}},
 \code{\link{fcm2}},
 \code{\link{fpcm}},
 \code{\link{fpppcm}},
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
\dontrun{
# Load dataset iris 
data(iris)
x <- iris[,-5]

# Initialize the prototype matrix using Inofrep algorithm
v <- inaparc::inofrep(x, k=3)$v
# Initialize the memberships degrees matrix 
u <- inaparc::imembrand(nrow(x), k=3)$u

# Run Gath & Geva with the initial prototypes and memberships
gg.res <- gg(x, centers=v, memberships=u, m=2)

# Show the fuzzy memberships degrees for the top 5 objects
head(gg.res$u, 5)
}
}

\concept{fuzzy clustering}
\concept{gath-geva clustering}
\concept{prototype-based clustering}
\concept{partitioning clustering}
\concept{cluster analysis}

\keyword{cluster}