\name{PLSnet}
\alias{PLSnet}
\title{Partial Least Squares network}
\description{
   Computes the connectivity scores for a network based on partial least squares (PLS).
}
\usage{
PLSnet(data,ncom=3,rescale.data=TRUE,symmetrize.scores=TRUE,rescale.scores=FALSE)
}
\arguments{
   \item{data}{microarray dataset with genes in columns and samples in rows.}
   \item{ncom}{the number of PLS components (latent variables) in PLS models.}
   \item{rescale.data}{indicates whether data should be rescaled.}
   \item{symmetrize.scores}{indicates whether PLS scores should be made to be symmetric.}
   \item{rescale.scores}{indicates whether PLS scores should be rescaled so that the largest score for each gene should be 1 in magnitude.}
}
\value{
   \item{PLSnet}{a matrix of interactions between gene pairs based on partial least squares.}   
}
\author{
The authors are Ryan Gill, Somnath Datta, and Susmita Datta.
The software is maintained by Ryan Gill \email{rsgill01@louisville.edu}.
}
\references{
Gill, R., Datta, S., and Datta, S. (2010) A statistical framework for differential network analysis from microarray data. \emph{BMC Bioinformatics}, \bold{11}, 95.

Pihur, V., Datta, S., and Datta, S. (2008) Reconstruction of genetic association networks from microarray data: a partial least squares approach. \emph{Bioinformatics}, \bold{24}(4), 561--568.
}
\examples{
# small example using PLSnet with 3 latent PLS components,
# data rescaled, and scores symmetrized but not rescaled
X1=rbind(
c(2.5,6.7,4.5,2.3,8.4,3.1),
c(1.2,0.7,4.0,9.1,6.6,7.1),
c(4.3,-1.2,7.5,3.8,1.0,9.3),
c(9.5,7.6,5.4,2.3,1.1,0.2))
s=PLSnet(X1)
print(round(s,4))

# small example using PLSnet with 2 latent PLS components,
# data rescaled, and scores symmetrized and rescaled
s2=PLSnet(X1,ncom=2,rescale.data=TRUE,symmetrize.scores=TRUE,rescale.scores=TRUE)
print(round(s2,4))
}
