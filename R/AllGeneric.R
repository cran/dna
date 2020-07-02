# file dna/R/AllGeneric.R
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

#if (!isGeneric("summary"))
# setGeneric("summary",function(object){standardGeneric("summary")})

if (!isGeneric("get.modules"))
 setGeneric("get.modules",function(object){standardGeneric("get.modules")})

if (!isGeneric("get.network1"))
 setGeneric("get.network1",function(object){standardGeneric("get.network1")})

if (!isGeneric("get.network2"))
 setGeneric("get.network2",function(object){standardGeneric("get.network2")})

if (!isGeneric("get.common.networks"))
 setGeneric("get.common.networks",function(object){standardGeneric("get.common.networks")})

if (!isGeneric("get.results"))
 setGeneric("get.results",function(object){standardGeneric("get.results")})


