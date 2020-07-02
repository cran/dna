# file dna/R/AllClass.R
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

setClass("modules",representation(module="factor"),
validity=function(object){
 if (sum(is.na(as.integer(levels(object@module))))>0)
  stop("[modules: validation] module is not integer-valued")
 a=sort(as.integer(levels(object@module)))
 n.a=length(a)
 if (sum(((a[1]!=0)&(a[1]!=1))|(a[-1]-a[-n.a]!=1))!=0)
  stop("[modules: validation] module is not properly defined")
 return(TRUE)
})

setClass("pairOfNetworks",representation(network1="matrix",network2="matrix"),
validity=function(object){
 if (is.null(colnames(object@network1))|is.null(colnames(object@network2))){
  if (ncol(object@network1)!=ncol(object@network2))
   stop("[pairOfNetworks: validation] No labels are given for the genes \n and the number of genes in each network do not match.")
 }
 else{
  if (sum(duplicated(colnames(object@network1)))>0)
   stop(paste("[pairOfNetworks: validation] There are duplicated gene names in Network 1:",colnames(object@network1)[duplicated(colnames(object@network1))]))
  if (sum(duplicated(colnames(object@network2)))>0)
   stop(paste("[pairOfNetworks: validation] There are duplicated gene names in Network 2:",colnames(object@network2)[duplicated(colnames(object@network2))]))
  if (length(intersect(colnames(object@network1),colnames(object@network2)))==0)
   stop("[pairOfNetworks: validation] The networks have no common genes.")
 }
 return(TRUE)
})

setClass("resultsClassTest",representation(p.value="numeric",
delta="numeric",class.genes="character"))

setClass("resultsIndTest",representation(p.values="numeric",
d="numeric"))

setClass("resultsModTest",representation(p.value="numeric",
N="numeric",modules1="modules",modules2="modules"))

