\name{modules-class}
\docType{class}
\alias{modules-class}
\alias{show,modules-method}
\alias{summary,modules-method}
\alias{get.modules,modules-method}

\title{Class \code{"modules"}}
\description{The class of modules of a network.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("modules", ...)}.}
\section{Slots}{
  \describe{
    \item{\code{module}:}{Object of class \code{"factor"} which represents the
module number of the genes in a network.  For genes that are in none of the
modules, the module number is listed as 0.}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "modules")}: lists the genes in each of the modules.}
    \item{summary}{\code{signature(x = "modules")}: summarizes a network by listing the number of genes in each module.}
    \item{get.modules}{\code{signature(object = "modules")}: returns the module number for each gene in a network.  For the genes that are in none of the modules, the module number is listed as 0. 
}
	 }
}
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

# summary method useful for large networks
summary(the.modules)

# method for extracting the modules
get.modules(the.modules)

# plot a graph of the modules
## Not run: network.modules(s,m=3,epsilon=.7,plot=TRUE)
## Not run: tkplot.close('1')
}
\keyword{classes}
