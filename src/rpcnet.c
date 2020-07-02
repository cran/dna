/*  file dna/src/rpcnet.c
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
 *	rpcnet(...)
 *
 * to be called as  .C(.)  in ../R/PCnet.R
 */

#include "plslib.h"
#include <R_ext/Lapack.h>

void rpcnet(double *origdata, double *s, int *ncom, int *n, int *p, int *rescaleData, int *symmetricScores, int *rescaleScores){
 double *data;
 double *X;
 double *y;
 double *b;
 double *G;
 double *T;
 double *eigenvalues;
 double *eigenvectors;
 int *isuppz;

 char uplo = 'U';
 char trans; 
 char jobz = 'V';
 char range = 'I';

 int pm1=*p-1;
 int np=(*n)*(*p);
 int npm1=(*n)*pm1;
 int i, j, k, count, xcount, jip, m;
 int incx=1;
 int incy=1;
 int lwork;
 int liwork;
 int il=*p-*ncom;
 int iu=pm1;
 
 int info;
 int iworkin;
 int *iwork;
 double tempdouble, tempsum, tempsum2, tempmean;
 double tempsd=0.0;
 double alpha=1.0;
 double beta=0.0;
 double doubleone=1.0;
 double doublezero=0.0;
 double doublem1=-1.0; 
 double vl, vu;
 double workin;
 double *work;
 double *tsk;

 data=Calloc(np,double);
 X=Calloc(npm1,double);
 y=Calloc(*n,double);
 b=Calloc(*ncom,double);
 G=Calloc(pm1*pm1,double);
 T=Calloc((*n)*(*ncom),double);
 eigenvalues=Calloc(*ncom+1,double);
 eigenvectors=Calloc(pm1*(*ncom)+1,double);
 isuppz=Calloc(pm1,int);
 tsk=Calloc(pm1,double);
 work=Calloc(1,double);
 iwork=Calloc(1,int);

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

 count=0;
 for (i=0;i<*n;i++){
  y[i]=data[count];
  count++;
 }

 xcount=0;
 for (i=0;i<*n;i++)
  for (j=0;j<pm1;j++){
   X[xcount]=data[count];
   xcount++;
   count++;
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
  trans='t';
  F77_CALL(dsyrk)(&uplo,&trans,&pm1,n,&doubleone,X,n,&doublezero,G,&pm1 FCLEN FCLEN);
  lwork=-1;
  liwork=-1;
  F77_CALL(dsyevr)(&jobz,&range,&uplo,&pm1,G,&pm1,&vl,&vu,&il,&iu,&doublem1,&m,eigenvalues,eigenvectors,&pm1,isuppz,&workin,&lwork,&iworkin,&liwork,&info FCLEN FCLEN FCLEN);
  lwork = (int)workin;
  work=Realloc(work,lwork,double);
  liwork = (int)iworkin;
  iwork=Realloc(iwork,liwork,int);
  F77_CALL(dsyevr)(&jobz,&range,&uplo,&pm1,G,&pm1,&vl,&vu,&il,&iu,&doublem1,&m,eigenvalues,eigenvectors,&pm1,isuppz,work,&lwork,iwork,&liwork,&info FCLEN FCLEN FCLEN);
  trans='n';
  F77_CALL(dgemm)(&trans,&trans,n,ncom,&pm1,&alpha,X,n,eigenvectors,&pm1,&beta,T,n FCLEN FCLEN);
  count=0;
  for (k=0;k<*ncom;k++){
   tempsum=0;
   tempsum2=0;
   xcount=0;
   for (i=0;i<*n;i++){
    tempsum+=T[count]*y[xcount];
    tempsum2+=T[count]*T[count];
    count++;
    xcount++;
   }
   b[k]=tempsum/tempsum2;
  }
  F77_CALL(dgemv)(&trans,&pm1,ncom,&alpha,eigenvectors,&pm1,b,&incx,&beta,tsk,&incy FCLEN);
  count=0;
  for (i=0;i<*p;i++)
   if (i!=j){
    s[i*(*p)+j]=tsk[count];
    count++;
   }
 }

 if (*symmetricScores==1){
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
 Free(iwork);
 Free(work);
 Free(tsk);
 Free(isuppz);
 Free(eigenvectors);
 Free(eigenvalues);
 Free(T);
 Free(G);
 Free(b);
 Free(y);
 Free(X);
 Free(data);
}
