#include "plslib.h"
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
