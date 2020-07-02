/*  file dna/src/rplsnet.c
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
 *	rplsnet(...)
 *
 * to be called as  .C(.)  in ../R/PLSnet.R
 */

#include "plslib.h"
#include <R_ext/Lapack.h>

void rplsnet(double *origdata, double *s, int *ncom, int *n, int *p, 
int *rescaleData, int *symmetrizeScores, int *rescaleScores){
 double *data;
 double *X;
 double *y;
 double *b;
 double *tempb;
 double *cX;
 double *tc;
 double *cc;
 double *tTT;
 double *TTTT;
 double *tX;
 double *c;
 double *tsk;
 
 char uplo = 'U';
 char side = 'L';
 char trans; 
 double alpha = 1.0;
 double beta = 0.0;
 int incx = 1;
 int incy = 1;

 int np=(*n)*(*p);
 int pm1=*p-1;
 int npm1=(*n)*pm1;
 int nn=(*n)*(*n);
 double m1=-1.0;
 int ncompm1=(*ncom)*pm1;

 int i, j, k, l, count, xcount, jip;

 double tempdouble, tempsum, tempsum2, tempmean;
 double tempsd=0.0;

 data=Calloc(np,double);
 X=Calloc((*n)*pm1,double);
 y=Calloc(*n,double);
 b=Calloc(*ncom,double);
 cX=Calloc((*n)*pm1,double);
 tc=Calloc(pm1*(*ncom),double);
 tTT=Calloc((*ncom)*(*n),double);
 TTTT=Calloc((*n)*(*n),double);
 tX=Calloc((*n)*pm1,double);
 c=Calloc(pm1*(*ncom),double);
 tsk=Calloc(pm1,double);

 F77_CALL(dcopy)(&np,origdata,&incx,data,&incy);

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
  tempb=b;
  jip=(j-1)*(*n);
  if (j>0)
   for (i=0;i<*n;i++){
    tempdouble=X[jip];
    X[jip]=y[i];
    y[i]=tempdouble;
    jip++;
   }
  F77_CALL(dcopy)(&npm1,X,&incx,cX,&incy);
  for (k=0;k<*ncom;k++){
   cc=tc+k*pm1;
   trans='t';
// column k+1 of cc <- t(cX) * y 
   F77_CALL(dgemv)(&trans,n,&pm1,&alpha,cX,n,y,&incx,&beta,cc,&incy FCLEN);
// standardize column k+1 of cc
   tempdouble=1/F77_CALL(dnrm2)(&pm1,cc,&incx);
   F77_CALL(dscal)(&pm1,&tempdouble,cc,&incx);
// column k+1 of tTT <- cX * cc
   trans='n';
   F77_CALL(dgemv)(&trans,n,&pm1,&alpha,cX,n,cc,&incx,&beta,tTT,&incy FCLEN);
// component k+1 of tempb <- (t(column k+1 of tTT) * y)/(t(column k+1 of tTT)*column k+1 of tTT)
   *tempb=F77_CALL(ddot)(n,tTT,&incx,y,&incy);
   tempdouble=1/F77_CALL(ddot)(n,tTT,&incx,tTT,&incy);
   *tempb*=tempdouble;
   tempb++;
// TTTT <- column k+1 of tTT * t(column k+1 of tTT)
   for (i=0;i<(*n)*(*n);i++)
    TTTT[i]=0;
   F77_CALL(dsyr)(&uplo,n,&alpha,tTT,&incx,TTTT,n FCLEN);
   for (i=0;i<(*n)-1;i++)
    for (l=i+1;l<(*n);l++)
     TTTT[i*(*n)+l]=TTTT[l*(*n)+i];
// divide TTTT by t(column k+1 of tTT) * column k+1 of tTT
   F77_CALL(dscal)(&nn,&tempdouble,TTTT,&incx);
// subtract TTTT by I
   for (i=0;i<*n;i++)
    TTTT[i*(*n+1)]-=1;
   F77_CALL(dscal)(&nn,&m1,TTTT,&incx);
// tX <- TTTT*cX
   F77_CALL(dsymm)(&side,&uplo,n,&pm1,&alpha,TTTT,n,cX,n,&beta,tX,n FCLEN FCLEN);
   F77_CALL(dcopy)(&npm1,tX,&incx,cX,&incy);
  }
  F77_CALL(dcopy)(&ncompm1,tc,&incx,c,&incy);

// tsk <- c * b
  trans='n';
  F77_CALL(dgemv)(&trans,&pm1,ncom,&alpha,c,&pm1,b,&incx,&beta,tsk,&incy FCLEN);
  count=0;
// column j+1 of s <- tsk
  for (i=0;i<*p;i++)
   if (i!=j){
    s[i*(*p)+j]=tsk[count];
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
  F77_CALL(dscal)(&pp,&tempdouble,s,&incx); 
 }
 for (j=0;j<*p;j++)
  s[j*(*p)+j]=1;

 Free(tsk);
 Free(c);
 Free(tX);
 Free(TTTT);
 Free(tTT);
 Free(tc);
 Free(cX);
 Free(b);
 Free(X);
 Free(y);
 Free(data); 
}
