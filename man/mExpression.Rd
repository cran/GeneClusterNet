\name{mExpression}
\alias{mExpression}
\docType{data}
\title{ Gene expression data set of Yeast.
%%   ~~ data name/kind ... ~~
}
\description{
 Sample data with 789 genes and 18 evenly spaced gene expression measurements. 
}
\references{
Spellman, P. T., Sherlock, G., Zhang, M. Q., Iyer, V. R., Anders, K., Eisen, M. B., ... & Futcher, B. (1998). \emph{Comprehensive identification of cell cycle-regulated genes of the yeast Saccharomyces cerevisiae by microarray hybridization}. Molecular biology of the cell, 9(12), 3273-3297.
}
\usage{data("mExpression")}

\examples{

# load the package 
library(GeneClusterNet)

data(mExpression)

# first 5 rows

Sample=mExpression[1:5,]

plot(1:18,Sample[1,],type="l", ylim=c(min(Sample),max(Sample)),xlab="Time",ylab="Expression",lwd=2)
color=1
for (i in 2:5){
  color=color+1
  lines(1:18,Sample[i,],col=color,lwd=2)
}
}
\keyword{datasets}
