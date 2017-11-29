\name{plotcluster}
\alias{plotcluster}
\title{
Plot Clustering Results
}
\description{
Plots clustering results from a cluster analysis with \sQuote{ppclust}.
}
\usage{
plotcluster(objx, mt, cm, tv, cp=1, pt=19, trans=FALSE)
}

\arguments{
  \item{objx}{an object of \sQuote{ppclust} class.}
  \item{mt}{a character to specify the membership type. The default is \option{u} for fuzzy membership degrees. The alternative option is  \option{t} for typicality degrees. The default is \option{t} for the algorithms which produce both types of membership matrices.}
  \item{cm}{a character to specify the crisping method. The default is \option{max} for using maximum degree of memberships for each data obejct. The alternative option is  \option{threshold} for the membership degrees exceeding a user-specified threshold value with \code{tv}.}
  \item{tv}{a number specifying the threshold membership degree when the value of \code{cm} is selected as \option{threshold}. The value of \code{tv} should be between 0 and 1. The default is \option{max} for using maximum degree of memberships for each data obejct. The alternative option is \option{threshold} for the 0.5.}
  \item{cp}{an integer for the index of available color palettes. The default is 1. The options are 2, 3, 4 and 5 for different color themes.}
 \item{pt}{an integer specifying the ASCII code of a point character to be used in plotting. The default is 19 for solid circles. Use \option{#} for displaying the cluster with their cluster labels.}
  \item{trans}{a logical value for the type of plots. The default is \option{FALSE} for solid point colors. The alternative option is \option{TRUE} for transparent point colors.}
}

\author{
Zeynel Cebeci
}

\examples{
# Run FCM for 3 clusters on the data set Iris
res.fcm <- fcm(x=iris[,-5], centers=3)

par(ask=TRUE)
# Plot the clustering results with solid colors
plotcluster(res.fcm, cp=1)

# Plot the same clustering results with transparent colors
plotcluster(res.fcm, cp=1, trans=TRUE)

# Plot the same clustering results for the memberships > 0.75
plotcluster(res.fcm, cp=1, cm="threshold", tv=0.75, trans=TRUE)
par(ask=FALSE)
}

\concept{possibilistic c-means clustering}
\concept{prototype-based clustering}
\concept{partitioning clustering}
\concept{cluster analysis}

\keyword{cluster}
