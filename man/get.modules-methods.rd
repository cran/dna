\name{get.modules-methods}
\docType{methods}
\alias{get.modules}
\alias{get.modules-methods}
\alias{get.modules,modules}
\title{Methods for Function \code{get.modules}}
\description{get.modules-methods}
\section{Methods}{
\describe{

\item{\code{signature(object = "modules")}}{
returns the module number for each gene in a network.  For the genes
that are in none of the modules, the module number is listed as 0. 
}
}}

\references{
Gill, R., Datta, S., and Datta, S. (2010) A statistical framework for differential network analysis from microarray data. \emph{BMC Bioinformatics}, \bold{11}, 95.
}
\examples{
# artificial example to show how to obtain modules from a matrix of
# connectivity scores
set.seed(26)
s=matrix(runif(100,-1,1),10,10);diag(s)=1;s=round((s+t(s))/2,1)
the.modules=network.modules(s,m=3,epsilon=.7)
the.modules

# method for extracting the modules
get.modules(the.modules)
}
\keyword{methods}
\keyword{modules}
