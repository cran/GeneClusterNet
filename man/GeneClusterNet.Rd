\name{GeneClusterNet}
\alias{GeneClusterNet}
\title{Inference of gene regulatory network}
\description{ This function clusters gene expressions and makes inference of gene regulatory network with dynamic Bayesian network method based on G1DBN. It integrates the analyses of GeneCluster, GeneClusterBIC, GeneClusterInterp and G1DBN.  

}
\usage{
GeneClusterNet(mExpression, times, orderLOP, alpha1 = 0.5, alpha2 = 0.05, 
	realign = F, cutoff = c(lowCut = -0.35, upCut = 0.2), 
	NumberOfCluster = 0, sLabels = NULL)
}
\arguments{
  \item{mExpression}{ a gene expression matrix with \eqn{p} columns (length of time vector) and \eqn{n} rows ( number of genes).
}
  \item{times}{ time vector 
}
  \item{orderLOP}{ order of Legendre Polynomials
}
  \item{alpha1}{threshold that use for edge selection in the 1st order dependence score matrix S1 and its default setting is 0.5. See https://cran.r-project.org/web/packages/G1DBN/G1DBN.pdf.
}
  \item{alpha2}{ threshold that use for edge selection in the score matrix S2 and its default setting is 0.05. See https://cran.r-project.org/web/packages/G1DBN/G1DBN.pdf
}
  \item{realign}{ defalut is FALSE, realign is true, the cutoff variable indicates a cutoff range 
}
  \item{cutoff}{ cutoff range for determining the time of initial up or down regulaiotn.
}
  \item{NumberOfCluster}{ number of cluster. if it is set to zero, the function will determine the optimal number of gene expression clusters automatically.
}
  \item{sLabels}{ defalut is NULL, the user could specify the cluster labels when number of cluster is specified. 
}
}

\value{
        A score matrix and an adjacency matrix. It also creates a list of plots same as the function GeneCluster. In addtition, it creates a plot of gene regulatory network.
}
\author{Yaqun Wang \email{yw505@sph.rutgers.edu}, Zhengyang Shi
}

\references{
Wang, Y., Xu, M., Wang, Z., Tao, M., Zhu, J., Wang, L., et al. (2012). \emph{How to cluster gene expression dynamics in response to environmental signals}. Briefings in bioinformatics, 13(2), 162-174.

Wang, Y., Berceli, S. A., Garbey, M. and Wu, R. (2016). \emph{Inference of gene regulatory network through adaptive dynamic Beyesian networm modeling}. Technical Report.

R package G1DBN available at \url{https://cran.r-project.org/package=G1DBN}
}


\examples{
 # load the package 
 library(GeneClusterNet)
 set.seed(1234)
 data(mExpression)
 Sample=mExpression[sample(1:nrow(mExpression),50,replace=FALSE),]
 GeneClusterNet (Sample, times=c(1:18), orderLOP=5 ,alpha1=0.5, alpha2=0.05,NumberOfCluster = 3)

}

