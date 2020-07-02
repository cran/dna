# file dna/R/PCnet.R
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

PCnet = function(data, ncom=3, rescale.data=TRUE, symmetrize.scores=TRUE,
rescale.scores=FALSE){
 data=as.matrix(data)
 n=as.integer(nrow(data))
 p=as.integer(ncol(data))
 pp=as.integer(p*p) 
 gene.names=colnames(data)
 if (is.null(gene.names))
  gene.names=paste("Gene",1:p)
 out=.C("rpcnet",as.double(data),s=double(pp),as.integer(ncom),n,p,
as.integer(rescale.data),as.integer(symmetrize.scores),
as.integer(rescale.scores))
 s=matrix(out$s,p,p,byrow=FALSE) 
 rownames(s)=gene.names
 colnames(s)=gene.names
 s
}

