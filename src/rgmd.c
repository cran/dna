/*  file dna/src/rgmd.c
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 or 3 of the License
 *  (at your option).
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  A copy of the GNU General Public License is available at
 *  http://www.r-project.org/Licenses/
 *
 *
 * Exports
 *	rgmd(...)
 *
 * to be called as  .C(.)  in ../R/network.modules.R
 */

#include "plslib.h"

void rgmd(double *s, int *module, int *m, double *ep, int *p){
 int i,j,k,group,gcount,isdone;
 group=1;
 for (i=0;i<*p;i++)
  module[i]=0;
 for (i=0;i<*p-1;i++)
  if (module[i]==0){
   module[i]=group;
   gcount=1;
   isdone=0;
   while (isdone==0){
    isdone=1;
    for (j=i;j<*p;j++)
     if (module[j]==group){
      for (k=i+1;k<*p;k++)
       if (module[k]==0)
	if (fabs(s[j*(*p)+k])>=*ep){
	 module[k]=group;
	 gcount++;
	 isdone=0;
	}
     }
   }
   if (gcount>=*m)
    group++;
   else
    for (j=0;j<*p;j++)
     if (module[j]==group)
      module[j]=-1;
  }
 for (j=0;j<*p;j++)
  if (module[j]==-1)
   module[j]=0;
}

