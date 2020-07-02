/*  file dna/src/dna.c
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
 */

#include "plslib.h"
#include <Rmath.h>

void rrunif(int *x, int *n)
{
    double a = 0;
    double b;
    double tx;
    b=*n;
    GetRNGstate();
    tx=0.0+runif(a, b);
    x[0]=(int) tx;
    PutRNGstate();
}

void absDISTclassgenes(double *s1,double *s2,int *f,int *nf,double *dlt,int *p){
 int i,j,mi;
 *dlt=0.0;
 for (i=0;i<*nf;i++)
  for (j=0;j<*nf;j++)
   if (i!=j){
    mi=(f[i]-1)*(*p)+(f[j]-1);
    (*dlt)+=fabs(s1[mi]-s2[mi]);
   }
}

void absDISTsinglegene(double *s1,double *s2,double *d,int *p){
 int i,j,ii;
 for (i=0;i<*p;i++)
  d[i]=0.0;
 for (i=0;i<*p;i++)
  for (j=0;j<*p;j++)
   if (i!=j){
    ii=i*(*p)+j;
    d[i]+=fabs(s1[ii]-s2[ii]);
   }
}

void sqrDISTclassgenes(double *s1,double *s2,int *f,int *nf,double *dlt,int *p){
 int i,j,mi;
 *dlt=0.0;
 for (i=0;i<*nf;i++)
  for (j=0;j<*nf;j++)
   if (i!=j){
    mi=(f[i]-1)*(*p)+(f[j]-1);
    (*dlt)+=(s1[mi]-s2[mi])*(s1[mi]-s2[mi]);
   }
}

void sqrDISTsinglegene(double *s1,double *s2,double *d,int *p){
 int i,j,ii;
 for (i=0;i<*p;i++)
  d[i]=0.0;
 for (i=0;i<*p;i++)
  for (j=0;j<*p;j++)
   if (i!=j){
    ii=i*(*p)+j;
    d[i]+=(s1[ii]-s2[ii])*(s1[ii]-s2[ii]);
   }
}

void perm(double *x1, double *x2, double *px1, double *px2, int *n1, int *n2, int *p){
 int i,j,k,kk,v,m,w,pr,r;
 int n;
 int *urn;
 n=(*n1)+(*n2);
 urn=Calloc(n,int);
 for (i=0;i<n;i++)
  urn[i]=i;
 for (j=0;j<*n1;j++){
  pr=j;
  m=n-j;
  rrunif(&v,&m);
  w=urn[v];
  if (w<*n1){
   r=w;
   for (k=0;k<*p;k++){
    px1[pr]=x1[r];
    pr+=*n1;
    r+=*n1;
   }
  }
  else{
   r=w-*n1;
   for (k=0;k<*p;k++){
    px1[pr]=x2[r];
    pr+=*n1;
    r+=*n2;    
   }
  }
  for (k=v;k<m-1;k++)
   urn[k]=urn[k+1];
 }
 for (j=0;j<*n2;j++)
  for (k=0;k<*p;k++){
   pr=j;
   w=urn[j];
   if (w<*n1){
    r=w;
    for (kk=0;kk<*p;kk++){
     px2[pr]=x1[r];
     pr+=*n2;
     r+=*n1;
    }
   }
   else{
    r=w-*n1;
    for (kk=0;kk<*p;kk++){
     px2[pr]=x2[r];
     pr+=*n2;
     r+=*n2;    
    }
   }
  }
 Free(urn);
}

void UnionIntersectionStat(int *module1, int *module2, double *sN, int *p){
 int i,j,gnum,gden,G0;
 G0=0;
 *sN=0;
 for (i=0;i<*p;i++){
  if ((module1[i]>0)||(module2[i]>0)){
   G0++;
   gnum=0;
   gden=0;
   for (j=0;j<*p;j++){
    if ((module1[i]>0)&&(module1[i]==module1[j])){
     gden++;
     if ((module2[i]>0)&&(module2[i]==module2[j]))
      gnum++;
    }
    else if ((module2[i]>0)&&(module2[i]==module2[j]))
     gden++;
   }
   (*sN)+=gnum/(gden+0.); 
  }
 }
 if (G0>0){
  (*sN)=1-(*sN)/(G0+0.);
 }
}

#include <R_ext/Rdynload.h>

static const R_CMethodDef CEntries[] = {
    {"rcor", (DL_FUNC) &rcor, 4},
    {"rpcnet", (DL_FUNC) &rpcnet, 8},
    {"rplsnet", (DL_FUNC) &rplsnet, 8},
    {"rrrnet", (DL_FUNC) &rrrnet, 8}, 
    {"rgmd", (DL_FUNC) &rgmd, 5},
    {"tdcclassPC", (DL_FUNC) &tdcclassPC, 15},
    {"tdcclassPLS", (DL_FUNC) &tdcclassPLS, 15},
    {"tdcclassRR", (DL_FUNC) &tdcclassRR, 15},
    {"tdcindPC", (DL_FUNC) &tdcindPC, 13},
    {"tdcindPLS", (DL_FUNC) &tdcindPLS, 13},
    {"tdcindRR", (DL_FUNC) &tdcindRR, 13},
    {"tdmsPC", (DL_FUNC) &tdmsPC, 16},
    {"tdmsPLS", (DL_FUNC) &tdmsPLS, 16},
    {"tdmsRR", (DL_FUNC) &tdmsRR, 16},
    {NULL, NULL, 0}
};

void R_init_dna(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
