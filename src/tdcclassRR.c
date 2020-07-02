/*  file dna/src/tdcclassRR.c
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
 *	tdcclassRR(...)
 *
 * to be called as  .C(.)  in ../R/test.class.genes.R
 */

#include "plslib.h"

void tdcclassRR(double *x1, double *x2, int *f, int *nf, double *pval, double *dlt, int *n1, int *n2, int *p, double *lambda, int *nperm, int *rescaleData, int *symmetrizeScores, int *rescaleScores, int *distancetype){
 int i,pvalnum;
 double pdlt;
 double *px1;
 double *px2;
 double *s1;
 double *s2;
 double *ps1;
 double *ps2;
 px1=Calloc((*n1)*(*p),double);
 px2=Calloc((*n2)*(*p),double);
 s1=Calloc((*p)*(*p),double);
 s2=Calloc((*p)*(*p),double);
 ps1=Calloc((*p)*(*p),double);
 ps2=Calloc((*p)*(*p),double);
 rrrnet(x1,s1,lambda,n1,p,rescaleData,symmetrizeScores,rescaleScores);
 rrrnet(x2,s2,lambda,n2,p,rescaleData,symmetrizeScores,rescaleScores);
 if (*distancetype==1)
  absDISTclassgenes(s1,s2,f,nf,dlt,p);
 else if (*distancetype==2)
  sqrDISTclassgenes(s1,s2,f,nf,dlt,p);
 (*dlt)/=*nf;
 (*dlt)/=(*nf)-1.0;
 pvalnum=0;
 Rprintf("Starting permutation test:\n");
 for (i=0;i<*nperm;i++){
  Rprintf("permutation %i out of %i\n",i+1,*nperm);
  perm(x1,x2,px1,px2,n1,n2,p);
  rrrnet(px1,ps1,lambda,n1,p,rescaleData,symmetrizeScores,rescaleScores);
  rrrnet(px2,ps2,lambda,n2,p,rescaleData,symmetrizeScores,rescaleScores);
  if (*distancetype==1)
   absDISTclassgenes(ps1,ps2,f,nf,&pdlt,p);
  else if (*distancetype==2)
   sqrDISTclassgenes(ps1,ps2,f,nf,&pdlt,p);
  pdlt/=*nf;
  pdlt/=(*nf)-1.0;
  if (pdlt>=*dlt)
   pvalnum++;
 }
 *pval=(0.0+pvalnum)/(*nperm);
 Free(s1);
 Free(s2);
 Free(ps1);
 Free(ps2);
 Free(px1);
 Free(px2);
}
