#include "plslib.h"
extern double ddot_(int *n, double *dx, int *incx, double *dy, int *incy);

extern double dnrm2_(int *n, double *x, int *incx);

extern int dscal_(int *n, double *alpha, double *x, int *incx);

extern int dcopy_(int *n, double *x, int *incx, double *y, int *incy);

extern int dsyr_(char *uplo, int *n, double *alpha, double *x, int *incx, double *a, int *lda);

extern int dgemv_(char *trans, int *m, int *n, double *alpha, double *a, int *lda, double *x, int *incx, double *beta, double *y, int *incy);

extern int dsymm_(char *side, char *uplo, int *m, int *n, double *alpha, double *a,int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);

void rplsnet(double *origdata, double *s, int *ncom, int *n, int *p, 
int *rescaleData, int *symmetrizeScores, int *rescaleScores){
 double *data;
 double *X;
 double *y;
 double *b;
 double *tempb;
 double *cX;
 double *tc;
 double *cc;
 double *tTT;
 double *TTTT;
 double *tX;
 double *c;
 double *tsk;
 
 char uplo = 'U';
 char side = 'L';
 char trans; 
 double alpha = 1.0;
 double beta = 0.0;
 int incx = 1;
 int incy = 1;

 int np=(*n)*(*p);
 int pm1=*p-1;
 int npm1=(*n)*pm1;
 int nn=(*n)*(*n);
 double m1=-1.0;
 int ncompm1=(*ncom)*pm1;

 int i, j, k, l, count, xcount, jip;

 double tempdouble, tempsum, tempsum2, tempmean;
 double tempsd=0.0;

 data=malloc(np*sizeof(double));
 X=malloc((*n)*pm1*sizeof(double));
 y=malloc((*n)*sizeof(double));
 b=malloc((*ncom)*sizeof(double));
 cX=malloc((*n)*pm1*sizeof(double));
 tc=malloc(pm1*(*ncom)*sizeof(double));
 tTT=malloc((*ncom)*(*n)*sizeof(double));
 TTTT=malloc((*n)*(*n)*sizeof(double));
 tX=malloc((*n)*pm1*sizeof(double));
 c=malloc(pm1*(*ncom)*sizeof(double));
 tsk=malloc(pm1*sizeof(double));

 dcopy_(&np,origdata,&incx,data,&incy);

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
  tempb=b;
  jip=(j-1)*(*n);
  if (j>0)
   for (i=0;i<*n;i++){
    tempdouble=X[jip];
    X[jip]=y[i];
    y[i]=tempdouble;
    jip++;
   }
  dcopy_(&npm1,X,&incx,cX,&incy);
  for (k=0;k<*ncom;k++){
   cc=tc+k*pm1;
   trans='t';
// column k+1 of cc <- t(cX) * y 
   dgemv_(&trans,n,&pm1,&alpha,cX,n,y,&incx,&beta,cc,&incy);
// standardize column k+1 of cc
   tempdouble=1/dnrm2_(&pm1,cc,&incx);
   dscal_(&pm1,&tempdouble,cc,&incx);
// column k+1 of tTT <- cX * cc
   trans='n';
   dgemv_(&trans,n,&pm1,&alpha,cX,n,cc,&incx,&beta,tTT,&incy);
// component k+1 of tempb <- (t(column k+1 of tTT) * y)/(t(column k+1 of tTT)*column k+1 of tTT)
   *tempb=ddot_(n,tTT,&incx,y,&incy);
   tempdouble=1/ddot_(n,tTT,&incx,tTT,&incy);
   *tempb*=tempdouble;
   tempb++;
// TTTT <- column k+1 of tTT * t(column k+1 of tTT)
   for (i=0;i<(*n)*(*n);i++)
    TTTT[i]=0;
   dsyr_(&uplo,n,&alpha,tTT,&incx,TTTT,n);
   for (i=0;i<(*n)-1;i++)
    for (l=i+1;l<(*n);l++)
     TTTT[i*(*n)+l]=TTTT[l*(*n)+i];
// divide TTTT by t(column k+1 of tTT) * column k+1 of tTT
   dscal_(&nn,&tempdouble,TTTT,&incx);
// subtract TTTT by I
   for (i=0;i<*n;i++)
    TTTT[i*(*n+1)]-=1;
   dscal_(&nn,&m1,TTTT,&incx);
// tX <- TTTT*cX
   dsymm_(&side,&uplo,n,&pm1,&alpha,TTTT,n,cX,n,&beta,tX,n);
   dcopy_(&npm1,tX,&incx,cX,&incy);
  }
  dcopy_(&ncompm1,tc,&incx,c,&incy);

// tsk <- c * b
  trans='n';
  dgemv_(&trans,&pm1,ncom,&alpha,c,&pm1,b,&incx,&beta,tsk,&incy);
  count=0;
// column j+1 of s <- tsk
  for (i=0;i<*p;i++)
   if (i!=j){
    s[i*(*p)+j]=tsk[count];
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
  dscal_(&pp,&tempdouble,s,&incx); 
 }
 for (j=0;j<*p;j++)
  s[j*(*p)+j]=1;

 free(tsk);
 free(c);
 free(tX);
 free(TTTT);
 free(tTT);
 free(tc);
 free(cX);
 free(b);
 free(X);
 free(y);
 free(data); 
}
