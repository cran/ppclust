\name{comp.omega}
\alias{comp.omega}
\title{
Compute the possibilistic penalty argument for PCM
}
\description{
Computes the vector \eqn{\vec{\Omega}}{\vec{\Omega}}, the possibilistic penalty argument for Possibilistic C-Means clustering analysis.
}
\usage{
comp.omega(d, u, m=2, pco=NULL, K=1)
}

\arguments{
  \item{d}{a numeric matrix containing the distances of objects to the centers of clusters.}
  \item{u}{a numeric matrix containing the fuzzy or possibilistic membership degrees of the data objects.}
  \item{pco}{an object of \sQuote{ppclust} class that contains the results from a clustering algorithm. If it is supplied, there is no need to input the arguments \code{d}, \code{u} and \code{m} separately because they are parsed from the object \code{pco} itself.}
  \item{m}{a number for the fuzzy exponent. The default is 2.}
  \item{K}{a number greater than 0 to be used as the weight of penalty term. The default is 1.}
}

\value{
  \item{omega}{a numeric vector containing the mobilization scale values for each cluster.}
}

\details{
The vector \eqn{\vec{\Omega}}{\vec{\Omega}} is called possibilistic penalty term that controls the variance of clusters (Gosztolya & Szilagyi, 2015). In general, it is computed by using the fuzzy intra cluster distances from a previous run of FCM. This vector of mobilization scale parameters related with the spread of clusters contains the distance values at where membership degrees becomes 0.5 for each cluster (Timm et al, 2004; Wachs et al, 2006; Alamelumangai, 2014). 

\deqn{\vec{\Omega} = K \frac{\sum\limits_{i=1}^n u_{ij}^m \; d^2(\vec{x}_i,\vec{v}_j)}{\sum\limits_{i=1}^n u_{ij}^m}\; ;\; 1\leq j\leq k}{\vec{\Omega} = K \frac{\sum\limits_{i=1}^n u_{ij}^m \; d^2(\vec{x}_i,\vec{v}_j)}{\sum\limits_{i=1}^n u_{ij}^m}\; ;\; 1\leq j\leq k}

Where:
\eqn{K}{K} is a coefficent, \eqn{K \in (0,\infty)}. It is generally chosen to be 1.
}

\seealso{
 \code{\link{pcm}},
 \code{\link{pcmr}},
 \code{\link{fpcm}},
 \code{\link{gkpfcm}},
 \code{\link{pfcm}}
}

\author{
Zeynel Cebeci
}

\references{
Alamelumangai N. (2014). Computer aided segmentation of mammary carcinoma on ultrasound images using soft computing techniques. PhD Thesis, Anna Univ., IN. <\url{http://hdl.handle.net/10603/50590}>

Wachs, J., Shapira, O., & Stern, H. (2006). A method to enhance the \sQuote{Possibilistic C-Means with Repulsion} algorithm based on cluster validity index. In \emph{Applied Soft Computing Technologies: The Challenge of Complexity}, pp. 77-87. Springer, Berlin, Heidelberg. <doi: 10.1007/3-540-31662-0_6>

Gosztolya, G. & Szilagyi, L. (2015). Application of fuzzy and possibilistic c-means clustering models in blind speaker clustering. \emph{Acta Polytechnica Hungarica}, 12(7):41-56. <\url{http://publicatio.bibl.u-szeged.hu/6151/1/2015-acta-polytechnica.pdf}>
}

\examples{
data(iris)
x <- iris[,-5]

# Run FCM 
res.fcm <- fcm(x=x, centers=3)

# Compute the mobilization scale values using the results from FCM
vomg1 <- comp.omega(pco=res.fcm)
vomg1

# Compute the mobilization scale values using the distances and memberships from FCM
vomg2 <- comp.omega(res.fcm$d, res.fcm$u, m=3)
vomg2

# Compute the mobilization scale values with the K value of 10
vomg3 <- comp.omega(res.fcm$d, res.fcm$u, m=2, K=10)
vomg3
}

\concept{possibilistic c-means clustering}
\concept{prototype-based clustering}
\concept{partitioning clustering}
\concept{cluster analysis}

\keyword{cluster}
