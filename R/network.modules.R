# file dna/R/network.modules.R
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

network.modules=function(s,m,epsilon,plot=FALSE,interactive=FALSE,...){
 p=as.integer(nrow(s))
 if (is.null(colnames(s))){
  colnames(s)=paste("Gene",1:p)
  rownames(s)=paste("Gene",1:p)
 }
 out=.C("rgmd",as.double(s),module=integer(p),as.integer(m),as.double(epsilon),p) 
 out$module=as.factor(out$module)
 names(out$module)=rownames(s)
 if (plot==TRUE){
  if (requireNamespace("igraph",quietly=TRUE)){
   if (sum(out$module!=0)>0){
    graph.genenames=names(out$module)[out$module!=0]
    graph.s=s[out$module!=0,out$module!=0]
    Is=abs(graph.s)>=epsilon
    rs=row(graph.s)
    cs=col(graph.s)
    vIs=c(Is)
    vrs=c(rs)
    vcs=c(cs)
    orc=rs<cs
    ex=vrs[vIs&orc]
    ey=vcs[vIs&orc]
    edges=rbind(ex,ey)
    g=igraph::graph.empty(directed=FALSE)
    g=igraph::add.vertices(g,length(graph.genenames),names=graph.genenames)
    g=igraph::add.edges(g,edges)
    if (interactive==TRUE)
     igraph::tkplot(g,vertex.label=igraph::V(g)$names,...)
    else
     plot(g,vertex.label=igraph::V(g)$names,...)
   }
  }
  else{
   cat("No plot created since there are no modules in this network.\n")
  }
 }
 new("modules",module=out$module)
}

