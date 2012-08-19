\name{gennet}
\alias{gennet}
\title{General Regression network}
\description{
   Computes the connectivity scores for a network based on a specified regression method.
}
\usage{
gennet(data, f, recenter.data=FALSE, rescale.data=FALSE, 
symmetrize.scores=FALSE, rescale.scores = FALSE, ...)
}
\arguments{
   \item{data}{microarray dataset with genes in columns and samples in rows.}
   \item{f}{regression method.}
   \item{recenter.data}{indicates whether data should be recentered.}
   \item{rescale.data}{indicates whether data should be rescaled.}
   \item{symmetrize.scores}{indicates whether PLS scores should be made to be symmetric.}
   \item{rescale.scores}{indicates whether PLS scores should be rescaled so that the largest score for each gene should be 1 in magnitude.}
   \item{...}{Any additional arguments for f.}
}
\value{
   \item{gennet}{a matrix of interactions between gene pairs based on the regression method supplied by the user.}   
}
\author{
The authors are Ryan Gill, Somnath Datta, and Susmita Datta.
The software is maintained by Ryan Gill \email{rsgill01@louisville.edu}.
}
\references{
Gill, R., Datta, S., and Datta, S. (2010) A statistical framework for differential network analysis from microarray data. \emph{BMC Bioinformatics}, \bold{11}, 95.

Hastie, T., Tibshirani, R., and Friedman, J. (2009) \emph{The Elements of Statistical Learning}. Springer: New York.
}
\examples{
# small example using gennet with a user-defined ridge regression
X1=rbind(
c(2.5,6.7,4.5,2.3,8.4,3.1),
c(1.2,0.7,4.0,9.1,6.6,7.1),
c(4.3,-1.2,7.5,3.8,1.0,9.3),
c(9.5,7.6,5.4,2.3,1.1,0.2))

## Not run: ourRR=function(X,y,lambda=1){solve(t(X)%*%X+diag(ncol(X)))%*%t(X)%*%y}
## Not run: gennet(X1,f=ourRR,recenter.data=TRUE,rescale.data=TRUE,symmetrize.scores=TRUE,rescale.scores=FALSE)

# compare results with RRnet
RRnet(X1)
}
