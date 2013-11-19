#include "plslib.h"
void perm(double *x1, double *x2, double *px1, double *px2, int *n1, int *n2, int *p){
 int i,j,k,kk,v,m,w,pr,r;
 int n;
 int *urn;
 n=(*n1)+(*n2);
 urn=malloc(n*sizeof(int));
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
 free(urn);
}
