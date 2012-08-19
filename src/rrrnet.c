#include "plslib.h"
extern int dscal_(int *n, double *alpha, double *x, int *incx);

extern int dcopy_(int *n, double *x, int *incx, double *y, int *incy);

extern int dgemv_(char *trans, int *m, int *n, double *alpha, double *a, int *lda, double *x, int *incx, double *beta, double *y, int *incy);

extern int dposv_(char *uplo, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, int *info);

extern int dsyrk_(char *uplo, char *trans, int *n, int *k, double *alpha, double *a, int *lda, double *beta, double *c, int *ldc);

void rrrnet(double *origdata, double *s, double *lambda, int *n, int *p, int *rescaleData, int *symmetrizeScores, int *rescaleScores){
 double *data;
 double *X;
 double *XTXpI;
 double *y;
 double *b;

 char uplo='U';
 char trans;
 int info;
 int intone=1;
 double doubleone=1.0;
 double doublezero=0.0;

 int i, j, ii, jip;
 int count, xcount;
 double tempsum, tempsum2, tempmean, tempdouble;
 double tempsd=0.0;

 int pm1=*p-1;
 int np=(*n)*(*p);
 
 data=malloc(np*sizeof(double));
 X=malloc((*n)*pm1*sizeof(double));
 XTXpI=malloc(pm1*pm1*sizeof(double));
 y=malloc((*n)*sizeof(double));
 b=malloc(pm1*sizeof(double));

 dcopy_(&np,origdata,&intone,data,&intone);

 count=0;
 for (j=0;j<*p;j++){
  tempsum=0;
  tempsum2=0;
  for (i=0;i<*n;i++){
   tempsum+=data[count];
   if (*rescaleData==1)
    tempsum2+=data[count]*data[count];
   count++;
  }
  tempmean=tempsum/(*n);
  if (*rescaleData==1)
   tempsd=sqrt((tempsum2-(*n)*tempmean*tempmean)/((*n)-1));
  count-=*n;
  for (i=0;i<*n;i++){
   if (*rescaleData==1)
    data[count]=(origdata[count]-tempmean)/tempsd;
   else
    data[count]=origdata[count]-tempmean;
   count++;
  }
 }

// put first column of data into y
 count=0;
 for (i=0;i<*n;i++){
  y[i]=data[count];
  count++;
 }

// put remaining columns of data into X
 xcount=0;
 for (i=0;i<*n;i++){
  for (j=0; j<pm1;j++){
   X[xcount]=data[count];
   xcount++;
   count++;
  }
 }

 for (j=0;j<*p;j++){
  jip=(j-1)*(*n);
  if (j>0)
   for (i=0;i<*n;i++){
    tempdouble=X[jip];
    X[jip]=y[i];
    y[i]=tempdouble;
    jip++;
   }
  count=0;
  for (i=0;i<pm1;i++)
   for (ii=0;ii<pm1;ii++){
    if (i==ii)
     XTXpI[count]=1;
    else
     XTXpI[count]=0;
    count++;
   }
  trans='t';
  dsyrk_(&uplo,&trans,&pm1,n,&doubleone,X,n,lambda,XTXpI,&pm1);
  dgemv_(&trans,n,&pm1,&doubleone,X,n,y,&intone,&doublezero,b,&intone);
  dposv_(&uplo,&pm1,&intone,XTXpI,&pm1,b,&pm1,&info);
  count=0;
  for (i=0;i<*p;i++)
   if (i!=j){
    s[i*(*p)+j]=b[count];
    count++;
   }
 }

 if (*symmetrizeScores==1){
  for (i=0;i<*p-1;i++)
   for (j=i+1;j<*p;j++){
    s[i*(*p)+j]=(s[i*(*p)+j]+s[j*(*p)+i])/2;
    s[j*(*p)+i]=s[i*(*p)+j];
   }
 }
 if (*rescaleScores==1){
  tempdouble=0;
  int pp=(*p)*(*p);
  for (i=0;i<*p;i++)
   for (j=0;j<*p;j++)
    if (i!=j)
     if (fabs(s[i*(*p)+j])>tempdouble)
      tempdouble=fabs(s[i*(*p)+j]);
  tempdouble=1/tempdouble;
  dscal_(&pp,&tempdouble,s,&intone); 
 }
 for (j=0;j<*p;j++)
  s[j*(*p)+j]=1;

 free(b);
 free(XTXpI);
 free(X);
 free(y);
 free(data);
}

