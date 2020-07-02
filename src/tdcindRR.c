/*  file dna/src/tdcindRR.c
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
 *	tdcindRR(...)
 *
 * to be called as  .C(.)  in ../R/test.individual.genes.R
 */

#include "plslib.h"

void tdcindRR(double *x1, double *x2, double *pval, double *d, int *n1, int *n2, int *p, double *lambda, int *nperm, int *rescaleData, int *symmetrizeScores, int *rescaleScores, int *distancetype){
 int i,j;
 int *pvalnum;
 double *pd;
 double *px1;
 double *px2;
 double *s1;
 double *s2;
 double *ps1;
 double *ps2;
 pvalnum=Calloc(*p,int);
 pd=Calloc(*p,double);
 px1=Calloc((*n1)*(*p),double);
 px2=Calloc((*n2)*(*p),double);
 s1=Calloc((*p)*(*p),double);
 s2=Calloc((*p)*(*p),double);
 ps1=Calloc((*p)*(*p),double);
 ps2=Calloc((*p)*(*p),double);
 rrrnet(x1,s1,lambda,n1,p,rescaleData,symmetrizeScores,rescaleScores);
 rrrnet(x2,s2,lambda,n2,p,rescaleData,symmetrizeScores,rescaleScores);
 if (*distancetype==1)
  absDISTsinglegene(s1,s2,d,p);
 else if (*distancetype==2)
  sqrDISTsinglegene(s1,s2,d,p);
 for (j=0;j<*p;j++)
  d[j]/=(*p-1.0);
 for (j=0;j<*p;j++)
  pvalnum[j]=0;
 Rprintf("Starting permutation test:\n");
 for (i=0;i<*nperm;i++){
  Rprintf("permutation %i out of %i\n",i+1,*nperm);
  perm(x1,x2,px1,px2,n1,n2,p);
  rrrnet(px1,ps1,lambda,n1,p,rescaleData,symmetrizeScores,rescaleScores);
  rrrnet(px2,ps2,lambda,n2,p,rescaleData,symmetrizeScores,rescaleScores);
  if (*distancetype==1)
   absDISTsinglegene(ps1,ps2,pd,p);
  else if (*distancetype==2)
   sqrDISTsinglegene(ps1,ps2,pd,p); 
  for (j=0;j<*p;j++)
   pd[j]/=(*p-1.0);
  for (j=0;j<*p;j++)
   if (pd[j]>=d[j])
    pvalnum[j]++;
 }
 for (j=0;j<*p;j++)
  pval[j]=(0.0+pvalnum[j])/(*nperm);
 Free(pvalnum);
 Free(pd);
 Free(s1);
 Free(s2);
 Free(ps1);
 Free(ps2);
 Free(px1);
 Free(px2);
}
