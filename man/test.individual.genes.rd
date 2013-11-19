\name{test.individual.genes}
\alias{test.individual.genes}
\title{Test for differential connectivity of individual genes}
\description{
   Tests for differential connectivity of individual genes between 
two networks using PLS scores.
}
\usage{
test.individual.genes(X1,X2,scores,distance="abs",
num.permutations=1000,check.networks=TRUE,...)
}
\arguments{
   \item{X1}{network 1 with genes in columns and samples in rows.}
   \item{X2}{network 2 with genes in columns and samples in rows.}
   \item{scores}{type of connectivity score to be used.  Either one of the built-in methods ("PLS", "PC", "RR", or "cor") can be used or a user-defined method can be supplied.}
   \item{distance}{distance function to be used. Either one of the built-in functions ("abs" or "sqr") can be used or a user-defined distance function can be supplied.}
   \item{num.permutations}{the number of random permutations.}
   \item{check.networks}{indicates whether get.common.networks is
used to preprocess the networks before the test is performed.}
   \item{...}{additional arguments for scores or distance.}
}
\value{
   \item{results}{result of test (of class resultsIndTest).}
}
\author{
The authors are Ryan Gill, Somnath Datta, and Susmita Datta.
The software is maintained by Ryan Gill \email{rsgill01@louisville.edu}.
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

# perform a test for differential connectivity of individual genes 
# with PLS connectivity scores and squared distances
## Not run: tig=test.individual.genes(X1,X2)
## Not run: summary(tig)

# extract results for a test for differential connectivity of individual genes
## Not run: results.tig=get.results(tig)
## Not run: results.tig

# perform a test for differential connectivity of individual genes 
# with PLS connectivity scores without rescaling the data,
# symmetrizing the scores, or rescaling the scores and with squared distances
# based on 10000 permutations
## Not run: tig2=test.individual.genes(X1,X2,scores="PLS",distance="sqr",
## num.permutations=10000,rescale.data=FALSE,symmetrize.scores=FALSE,
## rescale.scores=FALSE)
## Not run: summary(tig2)

# perform a test for differential connectivity of individual genes 
# with PLS connectivity scores and with custom distances
## Not run: our.dist=function(score1,score2){pmin(abs(score1-score2),1)}
## Not run: tig3=test.individual.genes(X1,X2,scores=PLSnet,distance=our.dist)
## Not run: summary(tig3)

# perform a test for differential connectivity of individual genes 
# with correlation connectivity scores
## Not run: tig4=test.individual.genes(X1,X2,scores="cor")
## Not run: summary(tig4)

# perform a test for differential connectivity of individual genes 
# with principal components regression connectivity scores
## Not run: tig5=test.individual.genes(X1,X2,scores="PC")
## Not run: summary(tig5)

# perform a test for differential connectivity of individual genes 
# with ridge regression connectivity scores with rescaled data,
# symmetrized and rescaled scores and a penalty parameter equal to 3
## Not run: tig6=test.individual.genes(X1,X2,scores="RR",
## rescale.scores=TRUE,lambda=3)
## Not run: summary(tig6)

# perform a test for differential connectivity of individual genes  
# with custom ridge regression connectivity scores with 
# centered and rescaled data and symmetrized and rescaled scores
## Not run: ourRR=function(X,y,lambda=3){
## solve(t(X)%*%X+lambda*diag(ncol(X)))%*%t(X)%*%y}
## Not run: ourRRnet=function(X){gennet(X,f=ourRR,recenter.data=TRUE,
## rescale.data=TRUE,symmetrize.scores=TRUE,rescale.scores=TRUE)}
## Not run: tig7=test.individual.genes(X1,X2,scores=ourRRnet)
## Not run: summary(tig7)
}
