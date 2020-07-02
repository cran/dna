\name{resultsClassTest-class}
\docType{class}
\alias{resultsClassTest-class}
\alias{show,resultsClassTest-method}
\alias{get.results,resultsClassTest-method}

\title{Class \code{"resultsClassTest"}}
\description{
Tests whether the connectivity scores for a set of important genes differ between two networks.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("resultsClassTest", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{p.value}:}{Object of class \code{"numeric"}; the p-value for the significance test.}
    \item{\code{delta}:}{Object of class \code{"numeric"}; the test statistic for the significance test.}
    \item{\code{class.genes}:}{Object of class \code{"character"}; the list of important genes.}
  }
}
\section{Methods}{
  \describe{
    \item{get.results}{\code{signature(object = "resultsClassTest")}: returns the p-value, test statistic, and the class of genes for a test for differential connectivity of the class of genes. }
    \item{show}{\code{signature(object = "resultsClassTest")}: summarizes the test by outputing the class of important genes, the value of the test statistic, and its p-value.}
	}
}

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

# perform a test for differential connectivity of all genes 
# with PLS connectivity scores and squared distances
## Not run: tcg=test.class.genes(X1,X2)
## Not run: results.tcg=get.results(tcg)
## Not run: results.tcg
}
\keyword{classes}
