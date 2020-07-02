# file dna/R/gennet.R
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

gennet=function (data, f, recenter.data=FALSE, rescale.data=FALSE, 
symmetrize.scores=FALSE, rescale.scores=FALSE, ...) 
{
    gene.names = colnames(data)
    data=t(as.matrix(data))
    data=scale(t(data),center=recenter.data,scale=rescale.data)  
    n = ncol(data)
    if (is.null(gene.names)) 
        gene.names = paste("Gene", 1:n)
    s = matrix(0, n, n)
    rownames(s) = gene.names
    colnames(s) = gene.names
    for (i in 1:n) {
        X = data[, -i]
        y = data[, i]
        s[i, i] = 1 - rescale.scores
        s[i, -i] = f(X, y, ...)
    }
   if (symmetrize.scores==TRUE)
    s=0.5*(s+t(s))
   if (rescale.scores==TRUE){
    ss=max(abs(s))
    diag(s)=rep(ss,length(diag(s)))
    s=s/ss
   }
   s
}

