\name{pairOfNetworks-class}
\docType{class}
\alias{pairOfNetworks-class}
\alias{get.common.networks,pairOfNetworks-method}
\alias{get.network1,pairOfNetworks-method}
\alias{get.network2,pairOfNetworks-method}
\alias{show,pairOfNetworks-method}

\title{Class \code{"pairOfNetworks"}}
\description{The class of a pair of networks that can be used for differential network analysis.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("pairOfNetworks", ...)}.}
\section{Slots}{
  \describe{
    \item{\code{network1}:}{Object of class \code{"matrix"}; the first network of expression values with genes in columns and subjects in rows.}
    \item{\code{network2}:}{Object of class \code{"matrix"}; the second network of expression values with genes in columns and subjects in rows.}
  }
}
\section{Methods}{
  \describe{
    \item{get.common.networks}{\code{signature(object = "pairOfNetworks")}:
 }
    \item{get.network1}{\code{signature(object = "pairOfNetworks")}: returns the first network from an object of class "pairOfNetworks".}
    \item{get.network2}{\code{signature(object = "pairOfNetworks")}: returns the second network from an object of class "pairOfNetworks".}
    \item{show}{\code{signature(object = "pairOfNetworks")}: lists the number of subjects and genes in each network as well as the number of genes common to both networks.}
	 }
}
\examples{
# small example illustrating how a pair of networks is
# preprocessed to obtain and align a set of common genes
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
networks=new("pairOfNetworks",network1=X1,network2=X2)
get.common.networks(networks)

# extract network 1
get.network1(networks)
# extract network 2
get.network2(networks)
}
\keyword{classes}
