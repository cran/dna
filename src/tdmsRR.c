#include "plslib.h"
#include "R.h"
void tdmsRR(double *x1, double *x2, int *m, double *epsilon, double *pval, double *sN, int *module1, int *module2, int *n1, int *n2, int *p, double *lambda, int *nperm, int *rescaleData, int *symmetrizeScores, int *rescaleScores){
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
 px1=malloc((*n1)*(*p)*sizeof(double));
 px2=malloc((*n2)*(*p)*sizeof(double));
 s1=malloc((*p)*(*p)*sizeof(double));
 s2=malloc((*p)*(*p)*sizeof(double));
 ps1=malloc((*p)*(*p)*sizeof(double));
 ps2=malloc((*p)*(*p)*sizeof(double));
 smod1=malloc((*p)*sizeof(int));
 smod2=malloc((*p)*sizeof(int));
 rrrnet(x1,s1,lambda,n1,p,rescaleData,symmetrizeScores,rescaleScores);
 rrrnet(x2,s2,lambda,n2,p,rescaleData,symmetrizeScores,rescaleScores);
 rgmd(s1,module1,m,epsilon,p);
 rgmd(s2,module2,m,epsilon,p);
 UnionIntersectionStat(module1,module2,sN,p);
 pvalnum=0;
 Rprintf("Starting permutation test:\n");
 for (i=0;i<*nperm;i++){
  Rprintf("permutation %i out of %i\n",i+1,*nperm);
  perm(x1,x2,px1,px2,n1,n2,p);
  rrrnet(px1,ps1,lambda,n1,p,rescaleData,symmetrizeScores,rescaleScores);
  rrrnet(px2,ps2,lambda,n2,p,rescaleData,symmetrizeScores,rescaleScores);
  rgmd(ps1,smod1,m,epsilon,p);
  rgmd(ps2,smod2,m,epsilon,p);
  UnionIntersectionStat(smod1,smod2,&psN,p);  
  if (psN>=*sN)
   pvalnum++;
 }
 *pval=(0.0+pvalnum)/(*nperm);
 free(smod1);
 free(smod2);
 free(s1);
 free(s2);
 free(ps1);
 free(ps2);
 free(px1);
 free(px2);
}
