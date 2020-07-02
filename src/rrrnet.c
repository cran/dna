/*  file dna/src/rrrnet.c
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
 *	rrrnet(...)
 *
 * to be called as  .C(.)  in ../R/RRnet.R
 */

#include "plslib.h"
#include <R_ext/Lapack.h>

void rrrnet(double *origdata, double *s, double *lambda, int *n, int *p, int *rescaleData, int *symmetrizeScores, int *rescaleScores){
 double *data;
 double *X;
 double *XTXpI;
 double *y;
 double *b;

 char uplo='U';
 char trans;
 int info;
 int intone=1;
 double doubleone=1.0;
 double doublezero=0.0;

 int i, j, ii, jip;
 int count, xcount;
 double tempsum, tempsum2, tempmean, tempdouble;
 double tempsd=0.0;

 int pm1=*p-1;
 int np=(*n)*(*p);
 
 data=Calloc(np,double);
 X=Calloc((*n)*pm1,double);
 XTXpI=Calloc(pm1*pm1,double);
 y=Calloc(*n,double);
 b=Calloc(pm1,double);

 F77_CALL(dcopy)(&np,origdata,&intone,data,&intone);

 count=0;
 for (j=0;j<*p;j++){
  tempsum=0;
  tempsum2=0;
  for (i=0;i<*n;i++){
   tempsum+=data[count];
   if (*rescaleData==1)
    tempsum2+=data[count]*data[count];
   count++;
  }
  tempmean=tempsum/(*n);
  if (*rescaleData==1)
   tempsd=sqrt((tempsum2-(*n)*tempmean*tempmean)/((*n)-1));
  count-=*n;
  for (i=0;i<*n;i++){
   if (*rescaleData==1)
    data[count]=(origdata[count]-tempmean)/tempsd;
   else
    data[count]=origdata[count]-tempmean;
   count++;
  }
 }

// put first column of data into y
 count=0;
 for (i=0;i<*n;i++){
  y[i]=data[count];
  count++;
 }

// put remaining columns of data into X
 xcount=0;
 for (i=0;i<*n;i++){
  for (j=0; j<pm1;j++){
   X[xcount]=data[count];
   xcount++;
   count++;
  }
 }

 for (j=0;j<*p;j++){
  jip=(j-1)*(*n);
  if (j>0)
   for (i=0;i<*n;i++){
    tempdouble=X[jip];
    X[jip]=y[i];
    y[i]=tempdouble;
    jip++;
   }
  count=0;
  for (i=0;i<pm1;i++)
   for (ii=0;ii<pm1;ii++){
    if (i==ii)
     XTXpI[count]=1;
    else
     XTXpI[count]=0;
    count++;
   }
  trans='t';
  F77_CALL(dsyrk)(&uplo,&trans,&pm1,n,&doubleone,X,n,lambda,XTXpI,&pm1 FCLEN FCLEN);
  F77_CALL(dgemv)(&trans,n,&pm1,&doubleone,X,n,y,&intone,&doublezero,b,&intone FCLEN);
  F77_CALL(dposv)(&uplo,&pm1,&intone,XTXpI,&pm1,b,&pm1,&info FCLEN);
  count=0;
  for (i=0;i<*p;i++)
   if (i!=j){
    s[i*(*p)+j]=b[count];
    count++;
   }
 }

 if (*symmetrizeScores==1){
  for (i=0;i<*p-1;i++)
   for (j=i+1;j<*p;j++){
    s[i*(*p)+j]=(s[i*(*p)+j]+s[j*(*p)+i])/2;
    s[j*(*p)+i]=s[i*(*p)+j];
   }
 }
 if (*rescaleScores==1){
  tempdouble=0;
  int pp=(*p)*(*p);
  for (i=0;i<*p;i++)
   for (j=0;j<*p;j++)
    if (i!=j)
     if (fabs(s[i*(*p)+j])>tempdouble)
      tempdouble=fabs(s[i*(*p)+j]);
  tempdouble=1/tempdouble;
  F77_CALL(dscal)(&pp,&tempdouble,s,&intone); 
 }
 for (j=0;j<*p;j++)
  s[j*(*p)+j]=1;

 Free(b);
 Free(XTXpI);
 Free(X);
 Free(y);
 Free(data);
}

