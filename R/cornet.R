# file dna/R/cornet.R
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

cornet = function(data,rescale.scores=FALSE){
 data=as.matrix(data)
 n=as.integer(nrow(data))
 p=as.integer(ncol(data))
 pp=as.integer(p*p) 
 gene.names=colnames(data)
 if (is.null(gene.names))
  gene.names=paste("Gene",1:p)
 out=.C("rcor",as.double(data),s=double(pp),n,p)
 s=matrix(out$s,p,p,byrow=FALSE) 
 rownames(s)=gene.names
 colnames(s)=gene.names
 if (rescale.scores==TRUE){
  for (i in 1:p) 
   s[i,i]=0
  ss=max(abs(s))
  diag(s)=rep(ss,length(diag(s)))
  s=s/ss
  for (i in 1:p)
   s[i,i]=1
 }
 s
}

