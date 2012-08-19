\name{network.modules}
\alias{network.modules}
\title{Determine modules for network}
\description{
   Determine the modular structure for a network.
}
\usage{
network.modules(s,m,epsilon,plot=FALSE,...)
}
\arguments{
   \item{s}{scores for a network.}
   \item{m}{minimum cluster size parameter.}
   \item{epsilon}{threshold parameter.}
   \item{plot}{indicates whether to create a graph for the network using the
tkplot function in the igraph package.}
   \item{...}{additional arguments passed to the tkplot function in the igraph package.}
}
\value{
   \item{modules}{an object of S4-class "modules" for the network}   
}
\author{
The authors are Ryan Gill, Somnath Datta, and Susmita Datta.
The software is maintained by Ryan Gill \email{rsgill01@louisville.edu}.
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
}

