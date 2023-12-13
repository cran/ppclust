\name{ppclust-package}
\alias{ppclust-package}
\docType{package}
\title{Probabilistic and Possibilistic Cluster Analysis}

\description{
Partitioning clustering simply divides the objects in a data set into non-overlapping subsets or clusters by using the prototype-based probabilistic and possibilistic clustering algorithms. The package covers various functions for K-Means (MacQueen, 1967), Fuzzy C-Means (Bezdek, 1974), Possibilitic C-Means (Krishnapuram & Keller, 1993;1996), Possibilistic and Fuzzy C-Means (Pal et al, 2005), Possibilistic Clustering Algorithm (Yang et al, 2006), Possibilistic C-Means with Repulsion (Wachs et al, 2006), Unsupervised Possibilistic Fuzzy C-Means (Wu et al, 2010) and the other variant algorithms which produce hard, fuzzy and possibilistic partitions of numeric data sets. The cluster prototypes and membership matrices required by the partitioning algorithms can be initialized with many initialization techniques that are available in the package \sQuote{\pkg{inaparc}}. As the distance metrics, not only the Euclidean distance but also a set of the commonly used distance metrics are available to use with some of the algorithms in the package.
}

\details{
The goal of prototype-based algorithms as the most widely-used group of partitioning clustering algorithms is to partition a data set of \emph{n} objects with \emph{p} features into \emph{k}, a pre-defined number of clusters which are the non-hierchical subsets of data. On the clustering context, a prototype is a data item that represents or characterizes a cluster. Usually, a prototype can be regarded as the most central point in a data subspace (Tan et al. 2006). 

Among the partitioning-based algorithms, the hard clustering algorithms, i.e. K-means, assume that each data point belongs to one cluster;  however, in practice clusters may overlap and data points may belong to more than one cluster.  In this case the membership degrees of a data point to clusters should be a value between zero and one. This idea has been modelled with the fuzzy clustering algorithms. The Fuzzy C-means proposed by Bezdek (FCM) (Bezdek, 1981), is the well-known partitioning based fuzzy algorithm. It assigns a fuzzy membership degree to each data point based on its distances to the cluster centers. If a data point is closer to a cluster center, its membership to this cluster be higher than its memberships to the other clusters. 

As an extension of the basic FCM algorithm, Gustafson and Kessel (GK) clustering algorithm employs an adaptive distance norm in order to detect clusters with different geometrical shapes (Babuska, 2001; Correa et al, 2011). Gath and Geva (1989) revealed that the fuzzy maximum likelihood estimates (FMLE) clustering algorithm can be used to detect clusters of varying shapes, sizes and densities. 

In the related literature it has been revealed that FCM is sensitive to noise and outliers in the data sets. Krishnapuram and Keller proposed and later improved the Possibilistic C-means (PCM) algorithm in order to avoid the effect of outliers or noises on the clustering performance (Krishnapuram & Keller, 1993;1996). Althouh PCM solves the problems due to noise and outliers by relaxing the probabilistic constraint of FCM, it has the disadvantage to  generate coincident clusters with poor initializations. In order to overcome this problem, several variants of FCM and PCM have been proposed. Fuzzy Possibilistic C-means (FPCM) (Pal et al, 1997) is one of the mixed algorithms for simultaneously constructing memberships and typicalities. However FPCM has some problems because of the row sum constraint with the typicality values that produces unrealistic typicality values for large data sets. By adding a repulsion term forcing clusters far away from each other, an extension of PCM with repulsion has been introduced by Timm et al (2004). Pal et al (2005) proposed the Possibilistic Fuzzy C-means (PFCM) to overcome the noise sensitivity defect of FCM, to eliminate the coincident clusters problem of PCM and to eliminate the row sum constraints of FPCM.

Possibilistic Clustering Algorithm (PCA) proposed by Yang and Wu (2006) was in the front lines of another direction of algorithms improving FCM and PCM. The authors argued that the resulting membership of their proposed algorithm becomes an exponential function, so that it is robust to noise and outliers. Recently, Wu et al (2010) introduced the Unsupervised Possibilistic Fuzzy Clustering (UPFC) algorithm. UPFC is an extension of PCA by combining FCM and PCA algorithms in order to overcome the noise sensitivity problem of FCM and the coincident clusters problem of PCM. When compared to previous algorithms, UPFC seems promising algorithm for clustering in fuzzy and possibilistic environments since it needs not a probabilistic initialization.
}

\section{Basic Notations}{
Given a data set \eqn{\mathbf{X}}{\mathbf{X}} describing \eqn{n}{n} data objects, the probabilistic and possibilistic clustering algorithms partition data into \eqn{k}{k}, a predefined number of clusters through the minimization of their related objective functions with some probabilistic or possibilistic constraints. 

\eqn{\mathbf{X} = \{\vec{x}_1, \vec{x}_2,\dots, \vec{x}_n\} \subseteq \Re^p}{\vec{X} = \{\vec{x}_1, \vec{x}_2,\dots, \vec{x}_n\} \subseteq\Re^p} is the data set for \eqn{n}{n} objects in the \emph{p}-dimensional data space \eqn{\Re}{\Re}, 
where:

\itemize{
 \item \eqn{n}{n} is the number of objects in the data set, \eqn{1\leq n\leq \infty}{1\leq n\leq \infty} 
 \item \eqn{p}{p} is the number of features or variables which describes the data objects,
 \item \eqn{\vec{x}_i}{\vec{x}_i} is \emph{p}-length data vector for the i\emph{th} data object. 
}

On the clustering context, clusters are mostly represented by their prototypes. The prototypes are generally the centers of clusters which can be either centroids or medoids. The probabilistic and possibilistic partitioning clustering algorithms start with initialization of a cluster prototype matrix \eqn{\mathbf{V}}{\mathbf{V}}, and updates it through the iteration steps until it is stabilized. 

\eqn{\mathbf{V} = \{\vec{v}_1, \vec{v}_2, \dots, \vec{v}_k\} \subseteq\Re^n}{\vec{V} = \{\vec{v}_1, \vec{v}_2, \dots, \vec{v}_k\} \subseteq\Re^n} is the protoype matrix of the clusters, where:

\itemize{
  \item \eqn{k}{k} is the number of clusters, \eqn{1\leq k\leq n}{1\leq k\leq n} 
  \item \eqn{\vec{v}_j}{\vec{v}_j} is the \emph{p}-length prototype vector for the j\emph{th} cluster
} 

The clustering algorithms compute the membership degrees of data objects by using some distance metrics for calculation of their proximities to the cluster centers. 

\eqn{d(\vec{x}_i, \vec{v}_j)}{d(\vec{x}_i, \vec{v}_j)} is the distance measure between the data object \eqn{\vec{x}_i}{\vec{x}_i} and cluster prototype \eqn{\vec{v}_j}{\vec{v}_j}. In general, the squared Euclidean distance metric are used in most of the applications: 

\eqn{d_{sq.euclidean}(\vec{x}_i, \vec{v}_j) = d^2(\vec{x}_i, \vec{v}_j) = \mid\mid \vec{x}_i - \vec{v}_j\mid \mid^2 \; = \; (\vec{x}_i - \vec{v}_j)^T \cdot (\vec{x}_i - \vec{v}_j)}{d_{sq.euclidean}(\vec{x}_i, \vec{v}_j) = d^2(\vec{x}_i, \vec{v}_j) = \mid\mid \vec{x}_i - \vec{v}_j\mid \mid^2 \; = \; (\vec{x}_i - \vec{v}_j)^T \cdot (\vec{x}_i - \vec{v}_j)} 

The clustering algorithms usually are run with the standard and squared Euclidean distance norms, which induce hyperspherical clusters. Therefore they are able find the clusters with the same shape and orientation because the norm inducing matrix is an identity matrix \eqn{\mathbf{A} = \mathbf{I}}{\mathbf{A} = \mathbf{I}}. On the other hand, the distance metrics can be employed with a \eqn{n \times n}{n \times n} diagonal norm inducing matrix \eqn{\mathbf{A} = \mathbf{I} \; 1/\sigma_j^2}{\mathbf{A} = \mathbf{I} \; 1/\sigma_j^2} which modifies the distances depending on the direction in which the distance is measured (Timm et al, 2004; Balasko et al 2005). In this case, the squared Euclidean distance with the norm matrix \eqn{\mathbf{A}}{\mathbf{A}} becomes:

\eqn{d_{sq.euclidean}(\vec{x}_i, \vec{v}_j) = d_{A}^2(\vec{x}_i, \vec{v}_j) = \mid\mid \vec{x}_i - \vec{v}_j\mid \mid_A^2 \; = \; (\vec{x}_i - \vec{v}_j)^T \mathbf{A} (\vec{x}_i - \vec{v}_j)}{d_{sq.euclidean}(\vec{x}_i, \vec{v}_j) = d_{A}^2(\vec{x}_i, \vec{v}_j) = \mid\mid \vec{x}_i - \vec{v}_j\mid \mid_A^2 \; = \; (\vec{x}_i - \vec{v}_j)^T \mathbf{A} (\vec{x}_i - \vec{v}_j)} 

\eqn{\mathbf{A}}{\mathbf{A}} can be also formed as the inverse of the \eqn{n \times n}{n \times n} covariance matrix \eqn{\mathbf{F}}{\mathbf{F}}: 

\eqn{\mathbf{A} = \mathbf{F}^{-1}}{\mathbf{A} = \mathbf{F}^{-1}}, where:

\eqn{\mathbf{F} = \frac{1}{n} \sum\limits_{i=1}^n (\vec{x}_i - \bar{x})^T (\vec{x}_i - \bar{x})}{\mathbf{F} = \frac{1}{n} \sum\limits_{i=1}^n (\vec{x}_i - \bar{x})^T (\vec{x}_i - \bar{x})}.  

Where \eqn{\bar{x}}{\bar{x}} stands for the sample mean of the data. When the distance is induced with the norm matrix \eqn{\mathbf{A}}{\mathbf{A}} the Mahalanobis norm on \eqn{\Re^n}{\Re^n} can be written as follows: 

\eqn{d_{mahalanobis}(\vec{x}_i, \vec{v}_j)_A = \; (\vec{x}_i - \vec{v}_j)^T \mathbf{A}_j (\vec{x}_i - \vec{v}_j)}{d_{mahalanobis}(\vec{x}_i, \vec{v}_j)_A = \; (\vec{x}_i - \vec{v}_j)^T \mathbf{A}_j (\vec{x}_i - \vec{v}_j)} 

The membership degrees are the measures specifying the amount of belongingness of the data objects to the different clusters. A data point nearer to the center of a cluster has a higher degree of membership to this cluster. 

\eqn{\mathbf{U} = [u_{ij}]}{\vec{U} = [u_{ij}]} is the matrix for an hard, fuzzy or possibilistic partition of \eqn{\mathbf{X}}{\mathbf{X}}, where:

\itemize{
 \item \eqn{u_{ij} = u_{i}(\vec{x}_{j})}{u_{ij} = u_{i}(\vec{x}_{j})} is the membership degree of \eqn{\vec{x}_i}{\vec{x}_i} to the j\emph{th} cluster 
 }
}

\references{

Babuska, R. (2001). Fuzzy and neural control. DISC Course Lecture Notes. Delft University of Technology. Delft, the Netherlands. <\href{https://tr.scribd.com/document/209211977/Fuzzy-and-Neural-Control}{https://tr.scribd.com/document/209211977/Fuzzy-and-Neural-Control}>.

Balasko, B., Abonyi, J. & Feil, B. (2005). Fuzzy clustering and data analysis toolbox. Department of Process Eng., Univ. of Veszprem, Veszprem. 

Bezdek, J.C. (1981). \emph{Pattern recognition with fuzzy objective function algorithms}. Plenum, NY. <\href{https://isbnsearch.org/isbn/0306406713}{ISBN:0306406713}>

Cebeci, Z. & Yildiz, F. (2015). Bulanik C-Ortalamalar algoritmasýnýn farklý kume buyuklukleri icin hesaplama performansi ve kumeleme gecerliliginin karsilastirilmasi, In Proc.pf \emph{9. Ulusal Zootekni Bilim Kongresi}, Sep. 2015, Konya. pp. 227-239. \doi{10.13140/RG.2.1.2909.9288}
 
Cebeci, Z., Kavlak, A.T. & Yildiz, F. (2017). Validation of fuzzy and possibilistic clustering results, in Proc. of \emph{2017 Int. Artificial Intelligence & Data Processing Symposium}, IEEE. pp. 1-7. \doi{10.1109/IDAP.2017.8090183}
 
Cebeci, Z. (2018). Initialization of membership degree matrix for fast convergence of Fuzzy C-Means clustering, in Proc. of \emph{2018 Int.Conf.on Artificial Intelligence & Data Processing}, pp. 1-5. IEEE. \doi{10.1109/IDAP.2018.8620920}
 
Cebeci, Z. & Cebeci, C. (2018) kpeaks: an R package for quick selection of k for cluster analysis, in Proc. of \emph{2018 Int. Conf. on Artificial Intelligence & Data Processing}, pp. 1-7, IEEE. \doi{10.1109/IDAP.2018.8620896}

Cebeci, Z. (2019). Comparison of internal validity indices for fuzzy clustering. \emph{Journal of Agricultural Informatics}, 10(2):1-14. \doi{10.17700/jai.2019.10.2.537}
 
Correa, C., Valero, C., Barreiro, P., Diago, M. P., & Tardáguila, J. (2011). A comparison of fuzzy clustering algorithms applied to feature extraction on vineyard. In \emph{Proc. of the 14th Conf. of the Spanish Assoc. for Artificial Intelligence}. <\url{http://oa.upm.es/9246/}>.

Gath, I. & Geva, A.B. (1989). Unsupervised optimal fuzzy clustering. \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence}, 11 (7): 773-781. <doi:10.1109/34.192473>

Gustafson, D. E. & Kessel, W. C. (1979). Fuzzy clustering with a fuzzy covariance matrix. In \emph{Proc. of IEEE Conf. on Decision and Control including the 17th Symposium on Adaptive Processes}, San Diego. pp. 761-766. <doi:10.1109/CDC.1978.268028>

Pal, N.R., Pal, K., & Bezdek, J.C. (1997). A mixed c-means clustering model. In \emph{Proc. of the 6th IEEE Int. Conf. on Fuzzy Systems}, 1, pp. 11-21. <doi:10.1109/FUZZY.1997.616338>

Pal, N.R., Pal, S.K., Keller,J.M. & Bezdek, J.C. (2005). A possibilistic fuzzy c-means clustering algorithm. \emph{IEEE Transactions on Fuzzy Systems}, 13 (4): 517-530. <doi: 10.1109/TFUZZ.2004.840099>

Tan, P. N., Steinbach, M., & Kumar, V. (2006). Cluster analysis: Basic concepts and algorithms. In \emph{Introduction to Data Mining}. Pearson Addison Wesley. <\url{http://www-users.cs.umn.edu/~kumar/dmbook/ch8.pdf}>

Timm, H., Borgelt, C., Doring, C. & Kruse, R. (2004). An extension to possibilistic fuzzy cluster analysis. \emph{Fuzzy Sets and Systems}, 147 (1): 3-16. <doi:10.1016/j.fss.2003.11.009>

Yang, M. S. & Wu, K. L. (2006). Unsupervised possibilistic clustering. \emph{Pattern Recognition}, 39(1): 5-21. <doi:10.1016/j.patcog.2005.07.005>

Wu, X., Wu, B., Sun, J. & Fu, H. (2010). Unsupervised possibilistic fuzzy clustering. \emph{J. of Information & Computational Sci.}, 7 (5): 1075-1080.
}

\author{Zeynel Cebeci, Figen Yildiz, A. Tuna Kavlak, Cagatay Cebeci & Hasan Onder}

\seealso{
 \code{\link{as.ppclust}},
 \code{\link{comp.omega}},
 \code{\link{fcm}},
 \code{\link{fcm2}},
 \code{\link{fpcm}},
 \code{\link{fpppcm}},
 \code{\link{gg}},
 \code{\link{gk}},
 \code{\link{gkpfcm}},
 \code{\link{hcm}},
 \code{\link{is.ppclust}}
 \code{\link{pca}},
 \code{\link{ppclust2}}
 \code{\link{pcm}},
 \code{\link{pcmr}},
 \code{\link{pfcm}},
 \code{\link{summary.ppclust}},
 \code{\link{upfc}}
}

\concept{fuzzy clustering}
\concept{crisp clustering}
\concept{hard clustering}
\concept{probabilistic clustering}
\concept{possibilistic clustering}
\concept{prototype-based clustering}
\concept{partitioning clustering}
\concept{partitional clustering}
\concept{cluster analysis}
\concept{unsupervised learning}
\concept{flat clustering}
\concept{non-hiearchical clustering}

\keyword{cluster}
\keyword{multivariete}
