\name{cornet}
\alias{cornet}
\title{Correlation network}
\description{
   Computes the connectivity scores for a network based on correlation.
}
\usage{
cornet(data,rescale.scores=FALSE)
}
\arguments{
   \item{data}{microarray dataset with genes in columns and samples in rows.}
   \item{rescale.scores}{indicates whether PLS scores should be rescaled so that the largest score for each gene should be 1 in magnitude.}
}
\value{
   \item{cornet}{a correlation matrix measuring the interactions between gene pairs.}   
}
\author{
The authors are Ryan Gill, Somnath Datta, and Susmita Datta.
The software is maintained by Ryan Gill \email{rsgill01@louisville.edu}.
}
\references{
Gill, R., Datta, S., and Datta, S. (2010) A statistical framework for differential network analysis from microarray data. \emph{BMC Bioinformatics}, \bold{11}, 95.
}
\examples{
# small example using cornet without rescaled scores
X1=rbind(
c(2.5,6.7,4.5,2.3,8.4,3.1),
c(1.2,0.7,4.0,9.1,6.6,7.1),
c(4.3,-1.2,7.5,3.8,1.0,9.3),
c(9.5,7.6,5.4,2.3,1.1,0.2))
s=cornet(X1)
print(round(s,4))

# small example using cornet with rescaled scores
s2=cornet(X1,rescale.scores=TRUE)
print(round(s2,4))
}
