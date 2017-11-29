\name{gkpfcm}
\alias{gkpfcm}
\title{
Gustafson-Kessel Clustering Using PFCM
}
\description{
Partitions a numeric data set by using the Gustafson-Kessel (GK) algorithm within the PFCM (Possibilistic Fuzzy C-Means) clustering algorithm (Ojeda-Magaina et al, 2006).
}
\usage{
gkpfcm(x, centers, memberships, m=2, eta=2, K=1, omega, a, b, 
    dmetric="sqeuclidean", pw = 2, alginitv="kmpp", 
    alginitu="imembrand", nstart=1, iter.max=1000, con.val=1e-09, 
    fixcent=FALSE, fixmemb=FALSE, stand=FALSE, numseed)
}

\arguments{
  \item{x}{a numeric vector, data frame or matrix.}
  \item{centers}{an integer specifying the number of clusters or a numeric matrix containing the initial cluster centers.}
  \item{memberships}{a numeric matrix containing the initial membership degrees. If missing, it is internally generated.}
  \item{m}{a number greater than 1 to be used as the fuzziness exponent or fuzzifier. The default is 2.}
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
  \item{fixcent}{a logical flag to make the initial cluster centers not changed along the different starts of the algorithm. The default is \code{FALSE}. If it is \code{TRUE}, the initial centers are not changed in the successive starts of the algorithm when the \code{nstart} is greater than 1.}
  \item{fixmemb}{a logical flag to make the initial membership degrees not changed along the different starts of the algorithm. The default is \code{FALSE}. If it is \code{TRUE}, the initial memberships are not changed in the successive starts of the algorithm when the \code{nstart} is greater than 1.}
  \item{stand}{a logical flag to standardize data. Its default value is \code{FALSE}. If its value is \code{TRUE}, the data matrix \code{x} is standardized.}
  \item{numseed}{a seeding number to set the seed of R's random number generator.}
}

\details{
Gustafson-Kessel clustering within Possibilistic Fuzzy C-Means (GKPFCM) algorithm is an improvement for PFCM algorithm that consists of the modification of the distance metric for \eqn{d_{ij_A}}{d_{ij_A}}. The original PFCM uses the Euclidean distance whereas GKPFCM uses the Mahalanobis distance with GK algorithm. Babuska et al (2002) have proposed an improvement for calculating the covariance matrix \eqn{\mathbf{F}_j}{\mathbf{F}_j} as follows:

\eqn{\mathbf{F}_j := (1 - \gamma) \mathbf{F}_j + \gamma (\mathbf{F}_0)^{1/n} \mathbf{I}}{\mathbf{F}_j := (1 - \gamma) \mathbf{F}_j + \gamma (\mathbf{F}_0)^{1/n} \mathbf{I}}

In the above equation, \eqn{\mathbf{F}_j}{\mathbf{F}_j} involves a weighted identity matrix. The eigenvalues \eqn{\lambda_{ij}}{\lambda_{ij}} and the eigenvectors \eqn{\Phi_{ij}}{\Phi_{ij}} of the resulting matrix are calculated, and the maximum eigenvalue \eqn{\lambda_{i,max} = max_{j}/ \lambda_{ij}}{\lambda_{i,max} = max_{j}/ \lambda_{ij}} is determined. With the obtained results, \eqn{\lambda_{i,max} = \lambda_{ij}/\beta, \forall j}{\lambda_{i,max} = \lambda_{ij}/\beta, \forall j}, which satisfies \eqn{\lambda_{i,max} / \lambda_{ij} \geq \beta} is calculated. Finally, \eqn{\mathbf{F}_j}{\mathbf{F}_j} is recomputed with the following equation:

\eqn{\mathbf{F}_j = [\Phi_{j,1},\dots, \Phi_{j,n}] diag(\lambda_{j,1}, \dots, \lambda_{j,n}) [\Phi_{j,1},\dots, \Phi_{j,n}]^{-1} \;\; \forall j}

The objective function of GKPFCM is:

\eqn{J_{GKPFCM}(\mathbf{X}; \mathbf{V}, \mathbf{A}, \mathbf{U}) = \sum\limits_{i=1}^n \sum\limits_{j=1}^k  u_{ij}^m d_{A_j}(\vec{x}_i, \vec{v}_j)}{J_{GK}(\mathbf{X}; \mathbf{V}, \mathbf{A}, \mathbf{U}) = \sum\limits_{i=1}^n \sum\limits_{j=1}^k  u_{ij}^m d_{A_j}(\vec{x}_i, \vec{v}_j)}

\eqn{m}{m} is the fuzzifier to specify the amount of fuzziness for the clustering; \eqn{1\leq m\leq \infty}. It is usually chosen as 2. 

\eqn{\eta}{\eta} is the typicality exponent to specify the amount of typicality for the clustering; \eqn{1\leq \eta\leq \infty}. It is usually chosen as 2. 

The objective function \eqn{J_{GKPFCM}}{J_{GKPFCM}} is minimized by using the following update equations:

\eqn{u_{ij} =\Bigg[\sum\limits_{j=1}^k \Big(\frac{d_{A_j}(\vec{x}_i, \vec{v}_j)}{d_{A_j}(\vec{x}_i, \vec{v}_l)}\Big)^{2/(m-1)} \Bigg]^{-1} \;\;; 1\leq i \leq n \;,\; 1 \leq l \leq k}{u_{ij} =\Bigg[\sum\limits_{j=1}^k \Big(\frac{d_{A_j}(\vec{x}_i, \vec{v}_j)}{d_{A_j}(\vec{x}_i, \vec{v}_l)}\Big)^{2/(m-1)} \Bigg]^{-1} \;\;; 1\leq i \leq n \;,\; 1 \leq l \leq k}

\eqn{t_{ij} =\Bigg[\sum\limits_{j=1}^k \Big(\frac{d_{A_j}(\vec{x}_i, \vec{v}_j))}{d_{A_j}(\vec{x}_i, \vec{v}_l))}\Big)^{2/(\eta-1)} \Bigg]^{-1} \;;\; 1 \leq i \leq n \;;\, 1 \leq l \leq k}{t_{ij} =\Bigg[\sum\limits_{j=1}^k \Big(\frac{d_{A_j}(\vec{x}_i, \vec{v}_j))}{d_{A_j}(\vec{x}_i, \vec{v}_l))}\Big)^{2/(\eta-1)} \Bigg]^{-1} \;;\; 1 \leq i \leq n \;,\; 1 \leq l \leq k}

\eqn{\vec{v}_{j} =\frac{\sum\limits_{i=1}^n (u_{ij}^m + t_{ij}^\eta) \vec{x}_i}{\sum\limits_{i=1}^n (u_{ij}^m + t_{ij}^\eta)} \;\;; {1\leq j\leq k}}{\vec{v}_{j} =\frac{\sum\limits_{i=1}^n (u_{ij}^m + t_{ij}^\eta) \vec{x}_i}{\sum\limits_{i=1}^n (u_{ij}^m + t_{ij}^\eta)} \;\;; {1\leq j\leq k}}
}

\value{an object of class \sQuote{ppclust}, which is a list consists of the following items:
   \item{x}{a numeric matrix containing the processed data set.}
   \item{v}{a numeric matrix containing the final cluster prototypes (centers of clusters).}
   \item{u}{a numeric matrix containing the fuzzy memberships degrees of the data objects.}
   \item{d}{a numeric matrix containing the distances of objects to the final cluster prototypes.}
   \item{k}{an integer for the number of clusters.}
   \item{m}{a number for the fuzzifier.}
   \item{eta}{a number greater than 1 to be used as the typicality exponent.}
   \item{a}{a number for the relative importance of the fuzzy part of the objective function.}
   \item{b}{a number for the relative importance of the possibilistic part of the objective function.}
   \item{omega}{a numeric vector of reference distances.}
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
Zeynel Cebeci & Cagatay Cebeci
}
\references{
Arthur, D. & Vassilvitskii, S. (2007). K-means++: The advantages of careful seeding, in \emph{Proc. of the 18th Annual ACM-SIAM Symposium on Discrete Algorithms}, p. 1027-1035. <\url{http://ilpubs.stanford.edu:8090/778/1/2006-13.pdf}>

Gustafson, D. E. & Kessel, W. C. (1979). Fuzzy clustering with a fuzzy covariance matrix. In \emph{Proc. of IEEE Conf. on Decision and Control including the 17th Symposium on Adaptive Processes}, San Diego. pp. 761-766. <doi:10.1109/CDC.1978.268028>

Babuska, R., van der Veen, P. J. &  Kaymak, U. (2002). Improved covariance estimation for Gustafson-Kessel clustering. In \emph{Proc. of Int. Conf. on Fuzzy Systems}, Hawaii, 2002, pp. 1081-1085. <\href{https://tr.scribd.com/document/209211977/Fuzzy-and-Neural-Control}{https://tr.scribd.com/document/209211977/Fuzzy-and-Neural-Control}>.

Ojeda-Magaina, B., Ruelas, R., Corona-Nakamura, M. A. & Andina, D. (2006). An improvement to the possibilistic fuzzy c-means clustering algorithm. In \emph{Proc. of IEEE World Automation Congress, 2006 (WAC'06)}. pp. 1-8. <doi:10.1109/WAC.2006.376056>
}

\seealso{
 \code{\link{ekm}},
 \code{\link{fcm}},
 \code{\link{fcm2}},
 \code{\link{fpcm}},
 \code{\link{fpppcm}},
 \code{\link{gg}},
 \code{\link{gk}},
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

# Initialize the memberships degrees matrix 
u <- inaparc::imembrand(nrow(x), k=3)$u

# Run FCM with the initial prototypes and memberships
gkpfcm.res <- gkpfcm(x, centers=v, memberships=u, m=2)

# Show the fuzzy membership degrees for the top 5 objects
head(gkpfcm.res$u, 5)
}
}

\concept{fuzzy c-means clustering}
\concept{gustafson-kessel clustering}
\concept{mixed c-means clustering}
\concept{prototype-based clustering}
\concept{partitioning clustering}
\concept{cluster analysis}

\keyword{cluster}