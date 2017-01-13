\name{GeneClusterBIC}
\alias{GeneClusterBIC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Optimal number of Gene Clusters 
%%  ~~function to do ... ~~
}
\description{Given time cousre expressions of \eqn{n} genes, time vector, order of Legendre Polynomials and a range of cluster numbers, e.g. from 1 to 15, the function can identify the optimal number of clusters, which has the smallest BIC value.
}
\usage{
GeneClusterBIC(mExpression, times, G = c(1:15), orderLOP)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mExpression}{ a gene expression matrix with \eqn{p} columns (length of time vector) and \eqn{n} rows ( number of genes).
%%     ~~Describe \code{mExpression} here~~
}
  \item{times}{ time vector specifies the time points of measurements.
%%     ~~Describe \code{times} here~~
}
  \item{G}{ range of number of clusters
%%     ~~Describe \code{G} here~~
}
  \item{orderLOP}{ order of Legendre Polynomials 
%%     ~~Describe \code{orderLOP} here~~
}
}

\value{ A list of BIC corresponds to every number of clusters and the optimal BIC. 
        A plot shows the smallest BIC.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
 GeneClusterBIC(Sample, times=c(1:18), G=c(1:5), orderLOP=5)

}

