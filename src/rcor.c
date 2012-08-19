#include "plslib.h"
extern int dcopy_(int *n, double *x, int *incx, double *y, int *incy);

void rcor(double *origdata, double *s, int *n, int *p){
 double *data;
 int count;
 int i,j,k,jip;
 double tempsum, tempsum2;
 double tempmean, tempsd, tempcor;

 int np=(*n)*(*p); 
 int incx=1;
 int incy=1;

 data=malloc(np*sizeof(double)); 
 dcopy_(&np,origdata,&incx,data,&incy);

 count=0;
 for (j=0;j<*p;j++){
  tempsum=0;
  tempsum2=0;
  for (i=0;i<*n;i++){
   tempsum+=data[count];
   tempsum2+=data[count]*data[count];
   count++;
  }
  tempmean=tempsum/(*n);
  tempsd=sqrt((tempsum2-(*n)*tempmean*tempmean)/((*n)-1));
  count-=*n;
  for (i=0;i<*n;i++){
   data[count]=(origdata[count]-tempmean)/tempsd;    
   count++;
  }
 }
 for (j=0;j<*p;j++){
  jip=j*(*p)+j;
  s[jip]=1;
 }
 for (j=0;j<*p-1;j++)
  for (i=j+1;i<*p;i++){
   count=0;
   tempsum=0;
   for (k=0;k<*n;k++)
    tempsum+=data[j*(*n)+k]*data[i*(*n)+k];
   tempcor=tempsum/(*n-1.0);
   jip=i*(*p)+j;
   s[jip]=tempcor;
   jip=j*(*p)+i;
   s[jip]=tempcor;
  }
 free(data);
}

