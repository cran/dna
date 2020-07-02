#include <R.h>

void rplsnet(double *origdata, double *s, int *ncom, int *n, int *p, int *rescaleData, int *symmetrizeScores, int *rescaleScores);

void rrrnet(double *origdata, double *s, double *lambda, int *n, int *p, int *rescaleData, int *symmetrizeScores, int *rescaleScores);

void rpcnet(double *origdata, double *s, int *ncom, int *n, int *p, int *rescaleData, int *symmetricScores, int *rescaleScores);

void rcor(double *origdata, double *s, int *n, int *p);

void rgmd(double *s, int *module, int *m, double *ep, int *p);

void rrunif(int *x, int *n);

void perm(double *x1, double *x2, double *px1, double *px2, int *n1, int *n2, int *p);

void absDISTsinglegene(double *s1,double *s2,double *d,int *p);

void sqrDISTsinglegene(double *s1,double *s2,double *d,int *p);

void absDISTclassgenes(double *s1,double *s2,int *f,int *nf,double *dlt,int *p);

void sqrDISTclassgenes(double *s1,double *s2,int *f,int *nf,double *dlt,int *p);

void tdcindPLS(double *x1, double *x2, double *pval, double *d, int *n1, int *n2, int *p, int *ncom, int *nperm, int *rescaleData, int *symmetrizeScores, int *rescaleScores, int *distancetype);

void tdcindPC(double *x1, double *x2, double *pval, double *d, int *n1, int *n2, int *p, int *ncom, int *nperm, int *rescaleData, int *symmetrizeScores, int *rescaleScores, int *distancetype);

void tdcindRR(double *x1, double *x2, double *pval, double *d, int *n1, int *n2, int *p, double *lambda, int *nperm, int *rescaleData, int *symmetrizeScores, int *rescaleScores, int *distancetype);

void tdcclassPLS(double *x1, double *x2, int *f, int *nf, double *pval, double *dlt, int *n1, int *n2, int *p, int *ncom, int *nperm, int *rescaleData, int *symmetrizeScores, int *rescaleScores, int *distancetype);

void tdcclassPC(double *x1, double *x2, int *f, int *nf, double *pval, double *dlt, int *n1, int *n2, int *p, int *ncom, int *nperm, int *rescaleData, int *symmetrizeScores, int *rescaleScores, int *distancetype);

void tdcclassRR(double *x1, double *x2, int *f, int *nf, double *pval, double *dlt, int *n1, int *n2, int *p, double *lambda, int *nperm, int *rescaleData, int *symmetrizeScores, int *rescaleScores, int *distancetype);

void UnionIntersectionStat(int *module1, int *module2, double *sN, int *p);

void tdmsPLS(double *x1, double *x2, int *m, double *epsilon, double *pval, double *sN, int *module1, int *module2, int *n1, int *n2, int *p, int *ncom, int *nperm, int *rescaleData, int *symmetrizeScores, int *rescaleScores);

void tdmsPC(double *x1, double *x2, int *m, double *epsilon, double *pval, double *sN, int *module1, int *module2, int *n1, int *n2, int *p, int *ncom, int *nperm, int *rescaleData, int *symmetrizeScores, int *rescaleScores);

void tdmsRR(double *x1, double *x2, int *m, double *epsilon, double *pval, double *sN, int *module1, int *module2, int *n1, int *n2, int *p, double *lambda, int *nperm, int *rescaleData, int *symmetrizeScores, int *rescaleScores);
