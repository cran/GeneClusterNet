\name{GeneCluster}
\alias{GeneCluster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Clustering  dynamic gene expression
%%  ~~function to do ... ~~
}
\description{ Providing functional clustering of time course gene expressions by using the Legendre orthogonal polynomials (LOP) to model cluster-specific curves for each cluster 
}
\usage{
GeneCluster(mExpression, times, NumberOfCluster, orderLOP)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mExpression}{ a gene expression matrix with \eqn{p} columns (length of time vector) and \eqn{n} rows ( number of genes).
%%     ~~Describe \code{mExpression} here~~
}
  \item{times}{ time vector specifies the time points of measurements.
%%     ~~Describe \code{times} here~~
}
  \item{NumberOfCluster}{ number of cluster (J)
%%     ~~Describe \code{NumberOfCluster} here~~
}
  \item{orderLOP}{ order of Legendre Polynomials (r)
%%     ~~Describe \code{orderLOP} here~~
}
}
\value{
        A list that contains $MeanExpression is the matrix with J rows ( J is number of Cluster) and n columns ( length of time vector), each rows is the mean expression of a cluster. $LOPCoefficient is the coefficient matrix of LOP with J rows and r+1 columns.$Classifications indicates the cluster label for each of genes. 
        A list of Plots, first plot is the mean expression plot for every clusters and each of the rest plots displays the trajectories of gene expressions in each cluster. 
}

\author{Yaqun Wang \email{yw505@sph.rutgers.edu}, Zhengyang Shi
}

\references{
Wang, Y., Xu, M., Wang, Z., Tao, M., Zhu, J., Wang, L., et al. (2012). \emph{How to cluster gene expression dynamics in response to environmental signals}. Briefings in bioinformatics, 13(2), 162-174.

Wang, Y., Berceli, S. A., Garbey, M. and Wu, R. (2016). \emph{Inference of gene regulatory network through adaptive dynamic Beyesian networm modeling}. Technical Report.
}

\examples{
 # load the package 
 library(GeneClusterNet)

 # Set the number of Cluster is 3 and order of Legendre Polynomials is 5.
 set.seed(1234)
 data(mExpression)
 
 Sample=mExpression[sample(1:nrow(mExpression),50,replace=FALSE),]

 GeneCluster(Sample, times=c(1:18), NumberOfCluster=3,orderLOP=5)


}
