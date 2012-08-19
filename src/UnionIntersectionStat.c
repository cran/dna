#include "plslib.h"
void UnionIntersectionStat(int *module1, int *module2, double *sN, int *p){
 int i,j,gnum,gden,G0;
 G0=0;
 *sN=0;
 for (i=0;i<*p;i++)
  if ((module1[i]>0)||(module2[i]>0)){
   G0++;
   gnum=0;
   gden=0;
   for (j=0;j<*p;j++)
    if ((module1[i]>0)&&(module1[i]==module1[j])){
     gden++;
     if ((module2[i]>0)&&(module2[i]==module2[j]))
      gnum++;
    }
    else if ((module2[i]>0)&&(module2[i]==module2[j]))
     gden++;
   (*sN)+=gnum/(gden+0.); 
  }
  if (G0>0)
   (*sN)=1-(*sN)/(G0+0.);
}
