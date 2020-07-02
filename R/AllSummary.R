# file dna/R/AllSummary.R
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

setMethod("summary","modules",
function(object){
 cat("Class: modules\n")
 F=slot(object,"module")
 mF=length(levels(F))-as.numeric("0"%in%levels(F))
 if (mF>0){
  for (i in 1:mF)
   cat(sum(F==i),"genes in Module",i,"\n")
 }
 else{
  cat("No modules\n")
 }
})

setMethod("summary","resultsIndTest",
function(object){
 cat("Tests for differential connectivity of individual genes\n\n")
 p.val=slot(object,"p.values")
 cat(sum(p.val<.001),"genes are significant at level 0.001\n")
 cat(sum(p.val<.005),"genes are significant at level 0.005\n")
 cat(sum(p.val<.01),"genes are significant at level 0.01\n")
 cat(sum(p.val<.05),"genes are significant at level 0.05\n")
})

setMethod("summary","resultsModTest",
function(object){
 cat("Tests for differential modular structure in two networks of genes\n\n")
 p.val=slot(object,"p.value")
 N=slot(object,"N")
 mod1=slot(object,"modules1")
 mod2=slot(object,"modules2")
 cat("Network 1:\n")
 summary(mod1)
 cat("\n")
 cat("Network 2:\n")
 summary(mod2)
 cat("\n")
 cat("Test statistic: N=",N,"\n")
 cat("P-value=",p.val,"\n") 
})


