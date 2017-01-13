\name{GeneClusterInterp}
\alias{GeneClusterInterp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Interpolating gene expression measurements 
%%  ~~function to do ... ~~
}
\description{To perform the dynamic Bayesian network analysis, the time points of gene expression measurements have to be evenly spaced.  If the original time points in a dataset are not even, this function can not only interpolate measurements to have evenly spaced time course gene expressions, but also allow users to specify the number of time points.
}
\usage{
GeneClusterInterp(LOPCoefficient, OriginalTime, outLen = 20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{LOPCoefficient}{ coefficient matrix of LOP.
%%     ~~Describe \code{LOPCoefficient} here~~
}
  \item{OriginalTime}{ time vector 
%%     ~~Describe \code{OriginalTime} here~~
}
  \item{outLen}{ number of new time points.
%%     ~~Describe \code{outLen} here~~
}
}

\value{
A matrix with J + 1 rows ( J is number of Cluster) and outLen columns ( number of equal time space). The first row of it provides the new time vector.
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
 set.seed(1234)
 data(mExpression)
 Sample=mExpression[sample(1:nrow(mExpression),50,replace=FALSE),]
 LOPCoefficient =GeneCluster (Sample, times=c(1:18), NumberOfCluster=3,orderLOP=5)$ LOPCoefficient
 GeneClusterInterp (LOPCoefficient, OriginalTime=c(1:18),outLen=20)


}

