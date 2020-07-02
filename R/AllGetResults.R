# file dna/R/AllGetResults.R
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

setMethod("get.results","resultsClassTest",function(object)
{
 p.val=slot(object,"p.value")
 delta=slot(object,"delta")
 genelist=slot(object,"class.genes")
 return(list(p.value=p.val,delta=delta,class.genes=genelist))
})

setMethod("get.results","resultsIndTest",function(object)
{
 p.val=slot(object,"p.values")
 test.stat=slot(object,"d")
 gene.names=names(p.val)
 p=length(p.val)
 out=data.frame(d=test.stat,p.value=p.val)
 return(data.frame(out[order(p.val),]))
})

setMethod("get.results","resultsModTest",function(object)
{
 p.val=slot(object,"p.value")
 N=slot(object,"N")
 mod1=slot(object,"modules1")
 mod2=slot(object,"modules2")
 return(list(p.value=p.val,N=N,modules1=mod1,modules2=mod2))
})

