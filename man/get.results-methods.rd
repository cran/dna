\name{get.results-methods}
\docType{methods}
\alias{get.results}
\alias{get.results-methods}
\alias{get.results,resultsIndTest}
\alias{get.results,resultsClassTest}
\alias{get.results,resultsModTest}

\title{Method for Function \code{get.results}}
\description{get.results-methods}
\section{Methods}{
\describe{

\item{\code{signature(object = "resultsIndTest")}}{
returns the p-values and test statistics for tests for
differential connectivity of individual genes. 
}

\item{\code{signature(object = "resultsClassTest")}}{
returns the p-value, test statistic, and the class of genes
for a test for differential connectivity of the class of genes. 
}

\item{\code{signature(object = "resultsModTest")}}{
returns the p-value, test statistic, and the modules for each network
for a test for overall modular structure. 
}
}}
\references{
Gill, R., Datta, S., and Datta, S. (2010) A statistical framework for differential network analysis from microarray data. \emph{BMC Bioinformatics}, \bold{11}, 95.
}
\examples{
# small example illustrating test procedures
X1=rbind(
c(2.5,6.7,4.5,2.3,8.4,3.1),
c(1.2,0.7,4.0,9.1,6.6,7.1),
c(4.3,-1.2,7.5,3.8,1.0,9.3),
c(9.5,7.6,5.4,2.3,1.1,0.2))
colnames(X1)=paste("G",1:6,sep="")

X2=rbind(
c(4.5,2.4,6.8,5.6,4.5,1.2,4.5),
c(7.6,9.0,0.1,3.4,5.6,5.5,1.2),
c(8.3,4.5,7.0,1.2,4.3,3.7,6.8),
c(3.4,1.1,6.9,7.2,3.1,0.9,6.6),
c(3.4,2.2,1.3,5.5,9.8,6.7,0.6))
colnames(X2)=paste("G",8:2,sep="")

# perform a test for differential connectivity of individual genes 
# with PLS connectivity scores and squared distances
## Not run: tig=test.individual.genes(X1,X2)
## Not run: results.tig=get.results(tig)
## Not run: results.tig

# perform a test for differential connectivity of all genes 
# with PLS connectivity scores and squared distances
## Not run: tcg=test.class.genes(X1,X2)
## Not run: results.tcg=get.results(tcg)
## Not run: results.tcg

# perform a test for modular structure using a minimum module size of 2
# and threshold of .5 with PLS connectivity scores
## Not run: test.modular.structure(X1,X2,min.module.size=2)
## Not run: results.tms=get.results(tms)
## Not run: results.tms
}
\keyword{methods}
\keyword{tests}
