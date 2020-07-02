/*  file dna/src/tdmsPLS.c
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
 *	tdmsPLS(...)
 *
 * to be called as  .C(.)  in ../R/test.modular.structure.R
 */
 
#include "plslib.h"

void tdmsPLS(double *x1, double *x2, int *m, double *epsilon, double *pval, double *sN, int *module1, int *module2, int *n1, int *n2, int *p, int *ncom, int *nperm, int *rescaleData, int *symmetrizeScores, int *rescaleScores){
 int i,pvalnum;
 double psN;
 double *px1;
 double *px2;
 double *s1;
 double *s2;
 double *ps1;
 double *ps2;
 int *smod1;
 int *smod2;
 px1=Calloc((*n1)*(*p),double);
 px2=Calloc((*n2)*(*p),double);
 s1=Calloc((*p)*(*p),double);
 s2=Calloc((*p)*(*p),double);
 ps1=Calloc((*p)*(*p),double);
 ps2=Calloc((*p)*(*p),double);
 smod1=Calloc(*p,int);
 smod2=Calloc(*p,int);
 rplsnet(x1,s1,ncom,n1,p,rescaleData,symmetrizeScores,rescaleScores);
 rplsnet(x2,s2,ncom,n2,p,rescaleData,symmetrizeScores,rescaleScores);
 rgmd(s1,module1,m,epsilon,p);
 rgmd(s2,module2,m,epsilon,p);
 UnionIntersectionStat(module1,module2,sN,p);
 pvalnum=0;
 Rprintf("Starting permutation test:\n");
 for (i=0;i<*nperm;i++){
  Rprintf("permutation %i out of %i\n",i+1,*nperm);
  perm(x1,x2,px1,px2,n1,n2,p);
  rplsnet(px1,ps1,ncom,n1,p,rescaleData,symmetrizeScores,rescaleScores);
  rplsnet(px2,ps2,ncom,n2,p,rescaleData,symmetrizeScores,rescaleScores);
  rgmd(ps1,smod1,m,epsilon,p);
  rgmd(ps2,smod2,m,epsilon,p);
  UnionIntersectionStat(smod1,smod2,&psN,p);  
  if (psN>=*sN)
   pvalnum++;
 }
 *pval=(0.0+pvalnum)/(*nperm);
 Free(smod1);
 Free(smod2);
 Free(s1);
 Free(s2);
 Free(ps1);
 Free(ps2);
 Free(px1);
 Free(px2);
}
