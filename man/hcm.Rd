\name{hcm}
\alias{hcm}
\title{
Hard C-Means Clustering
}
\description{
Partitions a numeric data set by using Hard C-Means (HCM) clustering algorithm (or K-Means) which has been proposed by MacQueen(1967). The function \code{hcm} is an extension of the basic \code{\link[stats]{kmeans}} with more input arguments and output values in order to make the clustering results comparable with those of other fuzzy and possibilistic algorithms. For instance, not only the Euclidean distance metric but also a number of distance metrics such as the squared Euclidean distance, the squared Chord distance etc. can be employed with the function \code{hcm}.
}
\usage{
hcm(x, centers, dmetric="euclidean", pw=2, alginitv="kmpp",  
   nstart=1, iter.max=1000, con.val=1e-9, stand=FALSE, numseed)
}

\arguments{
  \item{x}{a numeric vector, data frame or matrix.}
  \item{centers}{an integer specifying the number of clusters or a numeric matrix containing the initial cluster centers.}
  \item{dmetric}{a string for the distance metric. The default is \option{euclidean} for the squared Euclidean distances. See \code{\link{get.dmetrics}} for the alternative options.}
  \item{pw}{a number for the power of Minkowski distance calculation. The default is 2 if the \code{dmetric} is \option{minkowski}.}
  \item{alginitv}{a string for the initialization of cluster prototypes matrix. The default is \option{kmpp} for K-means++ initialization method (Arthur & Vassilvitskii, 2007). For the list of alternative options see \code{\link[inaparc]{get.algorithms}}.}
  \item{nstart}{an integer for the number of starts for clustering. The default is 1.}
  \item{iter.max}{an integer for the maximum number of iterations allowed. The default is 1000.}
  \item{con.val}{a number for the convergence value between the iterations. The default is 1e-09.}
  \item{stand}{a logical flag to standardize data. Its default value is \code{FALSE}. If its value is \code{TRUE}, the data matrix \code{x} is standardized.}
  \item{numseed}{a seeding number to set the seed of R's random number generator.}
}

\details{
Hard C-Means (HCM) clustering algorithm (or K-means) partitions a data set into \emph{k} groups, so-called clusters. The objective function of HCM is:

\eqn{J_{HCM}(\mathbf{X}; \mathbf{V})=\sum\limits_{i=1}^n d^2(\vec{x}_i, \vec{v}_j)}{J_{HCM}(\mathbf{X}; \mathbf{V})=\sum\limits_{i=1}^n d^2(\vec{x}_i, \vec{v}_j)}

See \code{\link{ppclust-package}} for the details about the terms in the above equation of \eqn{J_{HCM}}{J_{HCM}}.

The update equation for membership degrees is:

\eqn{u_{ij} = \left\{
                 \begin{array}{rl}
                   1 & if \; d^2(\vec{x}_i, \vec{v}_j) = min_{1\leq l\leq k} \; (d^2(\vec{x}_i, \vec{v}_l)) \\
                   0 & otherwise
                 \end{array}
               \right. }{u_{ij} = \left\{
                 \begin{array}{rl}
                   1 & if \; d^2(\vec{x}_i, \vec{v}_j) = min_{1\leq l\leq k} \; (d^2(\vec{x}_i, \vec{v}_l)) \\
                   0 & otherwise
                 \end{array}
               \right. }


The update equation for cluster prototypes is:
  
\eqn{\vec{v}_{j} =\frac{\sum\limits_{i=1}^n u_{ij} \vec{x}_i}{\sum\limits_{i=1}^n u_{ij}} \;\;; {1\leq j\leq k}}{\vec{v}_{j} =\frac{\sum\limits_{i=1}^n u_{ij} \vec{x}_i}{\sum\limits_{i=1}^n u_{ij}} \;\;; {1\leq j\leq k}}
}

\value{an object of class \sQuote{ppclust}, which is a list consists of the following items:
   \item{x}{a numeric matrix containing the processed data set.}
   \item{v}{a numeric matrix containing the final cluster prototypes (centers of clusters).}
   \item{u}{a numeric matrix containing the hard membership degrees of the data objects.}
   \item{d}{a numeric matrix containing the distances of objects to the final cluster prototypes.}
   \item{k}{an integer for the number of clusters.}
   \item{cluster}{a numeric vector containing the cluster labels of the data objects.}
   \item{csize}{a numeric vector containing the number of objects in the clusters.}
   \item{best.start}{an integer for the index of start with the minimum objective functional.}
   \item{iter}{an integer vector for the number of iterations in each start of the algorithm.}
   \item{func.val}{a numeric vector for the objective function values of each start of the algorithm.}
   \item{comp.time}{a numeric vector for the execution time of each start of the algorithm.}
   \item{wss}{a numeric vector containing the within-cluster sum of squares for each cluster.}
   \item{bwss}{a number for the between-cluster sum of squares.}
   \item{tss}{a number for the total within-cluster sum of squares.}
   \item{twss}{a number for the total sum of squares.}
   \item{stand}{a logical value, \code{TRUE} shows that \code{x} data set contains the standardized values of raw data.}
   \item{algorithm}{a string for the name of partitioning algorithm. It is \sQuote{HCM} with this function.}
   \item{call}{a string for the matched function call generating this \sQuote{ppclust} object.}
}

\author{
Zeynel Cebeci & Figen Yildiz
}

\references{
Arthur, D. & Vassilvitskii, S. (2007). K-means++: The advantages of careful seeding, in \emph{Proc. of the 18th Annual ACM-SIAM Symposium on Discrete Algorithms}, p. 1027-1035. <\url{http://ilpubs.stanford.edu:8090/778/1/2006-13.pdf}>

MacQueen, J.B. (1967). Some methods for classification and analysis of multivariate observations. In \emph{Proc. of 5th Berkeley Symp. on Mathematical Statistics and Probability}, Berkeley, Univ. of California Press, 1: 281-297. <\url{http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.308.8619&rep=rep1&type=pdf}>
}

\seealso{
 \code{\link[stats]{kmeans}},
 \code{\link{ekm}},
 \code{\link{fcm}},
 \code{\link{fcm2}},
 \code{\link{fpcm}},
 \code{\link{fpppcm}},
 \code{\link{gg}},
 \code{\link{gk}},
 \code{\link{gkpfcm}},
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

# Run HCM with the initial prototypes
res.hcm <- hcm(x, centers=v)

# Print, summarize and plot the clustering result
res.hcm$cluster
summary(res.hcm$cluster)
plot(x, col=res.hcm$cluster, pch=16)
}
}

\concept{hard c-means clustering}
\concept{k-means clustering}
\concept{prototype-based clustering}
\concept{partitioning clustering}
\concept{cluster analysis}

\keyword{cluster}
