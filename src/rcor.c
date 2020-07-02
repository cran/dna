/*  file dna/src/rcor.c
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
 *	rcor(...)
 *
 * to be called as  .C(.)  in ../R/cornet.R
 */

#include "plslib.h"
#include <R_ext/Lapack.h>

void rcor(double *origdata, double *s, int *n, int *p){
 double *data;
 int count;
 int i,j,k,jip;
 double tempsum, tempsum2;
 double tempmean, tempsd, tempcor;

 int np=(*n)*(*p); 
 int incx=1;
 int incy=1;

 data=Calloc(np,double); 
 F77_CALL(dcopy)(&np,origdata,&incx,data,&incy);

 count=0;
 for (j=0;j<*p;j++){
  tempsum=0;
  tempsum2=0;
  for (i=0;i<*n;i++){
   tempsum+=data[count];
   tempsum2+=data[count]*data[count];
   count++;
  }
  tempmean=tempsum/(*n);
  tempsd=sqrt((tempsum2-(*n)*tempmean*tempmean)/((*n)-1));
  count-=*n;
  for (i=0;i<*n;i++){
   data[count]=(origdata[count]-tempmean)/tempsd;    
   count++;
  }
 }
 for (j=0;j<*p;j++){
  jip=j*(*p)+j;
  s[jip]=1;
 }
 for (j=0;j<*p-1;j++)
  for (i=j+1;i<*p;i++){
   count=0;
   tempsum=0;
   for (k=0;k<*n;k++)
    tempsum+=data[j*(*n)+k]*data[i*(*n)+k];
   tempcor=tempsum/(*n-1.0);
   jip=i*(*p)+j;
   s[jip]=tempcor;
   jip=j*(*p)+i;
   s[jip]=tempcor;
  }
 Free(data);
}

