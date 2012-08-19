#include "plslib.h"
#include "R.h"
void tdcindPLS(double *x1, double *x2, double *pval, double *d, int *n1, int *n2, int *p, int *ncom, int *nperm, int *rescaleData, int *symmetrizeScores, int *rescaleScores, int *distancetype){
 int i,j;
 int *pvalnum;
 double *pd;
 double *px1;
 double *px2;
 double *s1;
 double *s2;
 double *ps1;
 double *ps2;
 pvalnum=malloc((*p)*sizeof(int));
 pd=malloc((*p)*sizeof(double));
 px1=malloc((*n1)*(*p)*sizeof(double));
 px2=malloc((*n2)*(*p)*sizeof(double));
 s1=malloc((*p)*(*p)*sizeof(double));
 s2=malloc((*p)*(*p)*sizeof(double));
 ps1=malloc((*p)*(*p)*sizeof(double));
 ps2=malloc((*p)*(*p)*sizeof(double));
 rplsnet(x1,s1,ncom,n1,p,rescaleData,symmetrizeScores,rescaleScores);
 rplsnet(x2,s2,ncom,n2,p,rescaleData,symmetrizeScores,rescaleScores);
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
  rplsnet(px1,ps1,ncom,n1,p,rescaleData,symmetrizeScores,rescaleScores);
  rplsnet(px2,ps2,ncom,n2,p,rescaleData,symmetrizeScores,rescaleScores);
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
 free(pvalnum);
 free(pd);
 free(s1);
 free(s2);
 free(ps1);
 free(ps2);
 free(px1);
 free(px2);
}