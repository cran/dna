#include "plslib.h"
#include "R.h"
void tdcclassPLS(double *x1, double *x2, int *f, int *nf, double *pval, double *dlt, int *n1, int *n2, int *p, int *ncom, int *nperm, int *rescaleData, int *symmetrizeScores, int *rescaleScores, int *distancetype){
 int i,pvalnum;
 double pdlt;
 double *px1;
 double *px2;
 double *s1;
 double *s2;
 double *ps1;
 double *ps2;
 px1=malloc((*n1)*(*p)*sizeof(double));
 px2=malloc((*n2)*(*p)*sizeof(double));
 s1=malloc((*p)*(*p)*sizeof(double));
 s2=malloc((*p)*(*p)*sizeof(double));
 ps1=malloc((*p)*(*p)*sizeof(double));
 ps2=malloc((*p)*(*p)*sizeof(double));
 rplsnet(x1,s1,ncom,n1,p,rescaleData,symmetrizeScores,rescaleScores);
 rplsnet(x2,s2,ncom,n2,p,rescaleData,symmetrizeScores,rescaleScores);
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
  rplsnet(px1,ps1,ncom,n1,p,rescaleData,symmetrizeScores,rescaleScores);
  rplsnet(px2,ps2,ncom,n2,p,rescaleData,symmetrizeScores,rescaleScores);
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
 free(s1);
 free(s2);
 free(ps1);
 free(ps2);
 free(px1);
 free(px2);
}
