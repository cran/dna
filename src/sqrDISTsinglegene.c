#include "plslib.h"
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