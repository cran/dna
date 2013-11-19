#include "plslib.h"

extern int dscal_(int *n, double *alpha, double *x, int *incx);

extern int dcopy_(int *n, double *x, int *incx, double *y, int *incy);

extern int dsyrk_(char *uplo, char *trans, int *n, int *k, double *alpha, double *a, int *lda, double *beta, double *c, int *ldc);

extern int dsyevr_(char *jobz, char *range, char *uplo, int *n, double *A, int *lda, double *vl, double *vu, int *il, int *iu, double *abstol, int *m, double *w, double *z, int *ldz, int *isuppz, double *work, int *lwork, int *iwork, int *liwork, int *info);

extern int dgemv_(char *trans, int *m, int *n, double *alpha, double *a, int *lda, double *x, int *incx, double *beta, double *y, int *incy);

extern int dgemm_(char *transa, char *transb, int* m, int* n, int* k, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);

void rpcnet(double *origdata, double *s, int *ncom, int *n, int *p, int *rescaleData, int *symmetricScores, int *rescaleScores){
 double *data;
 double *X;
 double *y;
 double *b;
 double *G;
 double *T;
 double *eigenvalues;
 double *eigenvectors;
 int *isuppz;

 char uplo = 'U';
 char trans; 
 char jobz = 'V';
 char range = 'I';

 int pm1=*p-1;
 int np=(*n)*(*p);
 int npm1=(*n)*pm1;
 int i, j, k, count, xcount, jip, m;
 int incx=1;
 int incy=1;
 int lwork;
 int liwork;
 int il=*p-*ncom;
 int iu=pm1;
 
 int info;
 int iworkin;
 int *iwork;
 double tempdouble, tempsum, tempsum2, tempmean;
 double tempsd=0.0;
 double alpha=1.0;
 double beta=0.0;
 double doubleone=1.0;
 double doublezero=0.0;
 double doublem1=-1.0; 
 double vl, vu;
 double workin;
 double *work;
 double *tsk;

 data=malloc(np*sizeof(double));
 X=malloc(npm1*sizeof(double));
 y=malloc((*n)*sizeof(double));
 b=malloc((*ncom)*sizeof(double));
 G=malloc(pm1*pm1*sizeof(double));
 T=malloc((*n)*(*ncom)*sizeof(double));
 eigenvalues=malloc((*ncom+1)*sizeof(double));
 eigenvectors=malloc((pm1*(*ncom)+1)*sizeof(double));
 isuppz=malloc(pm1*sizeof(int));
 tsk=malloc(pm1*sizeof(double));
 work=malloc(sizeof(double));
 iwork=malloc(sizeof(double));

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

 count=0;
 for (i=0;i<*n;i++){
  y[i]=data[count];
  count++;
 }

 xcount=0;
 for (i=0;i<*n;i++)
  for (j=0;j<pm1;j++){
   X[xcount]=data[count];
   xcount++;
   count++;
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
  trans='t';
  dsyrk_(&uplo,&trans,&pm1,n,&doubleone,X,n,&doublezero,G,&pm1);
  lwork=-1;
  liwork=-1;
  dsyevr_(&jobz,&range,&uplo,&pm1,G,&pm1,&vl,&vu,&il,&iu,&doublem1,&m,eigenvalues,eigenvectors,&pm1,isuppz,&workin,&lwork,&iworkin,&liwork,&info);
  lwork = (int)workin;
  work=realloc(work,lwork*sizeof(double));
  liwork = (int)iworkin;
  iwork=realloc(iwork,liwork*sizeof(int));
  dsyevr_(&jobz,&range,&uplo,&pm1,G,&pm1,&vl,&vu,&il,&iu,&doublem1,&m,eigenvalues,eigenvectors,&pm1,isuppz,work,&lwork,iwork,&liwork,&info);
  trans='n';
  dgemm_(&trans,&trans,n,ncom,&pm1,&alpha,X,n,eigenvectors,&pm1,&beta,T,n);
  count=0;
  for (k=0;k<*ncom;k++){
   tempsum=0;
   tempsum2=0;
   xcount=0;
   for (i=0;i<*n;i++){
    tempsum+=T[count]*y[xcount];
    tempsum2+=T[count]*T[count];
    count++;
    xcount++;
   }
   b[k]=tempsum/tempsum2;
  }
  dgemv_(&trans,&pm1,ncom,&alpha,eigenvectors,&pm1,b,&incx,&beta,tsk,&incy);
  count=0;
  for (i=0;i<*p;i++)
   if (i!=j){
    s[i*(*p)+j]=tsk[count];
    count++;
   }
 }

 if (*symmetricScores==1){
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
 free(iwork);
 free(work);
 free(tsk);
 free(isuppz);
 free(eigenvectors);
 free(eigenvalues);
 free(T);
 free(G);
 free(b);
 free(y);
 free(X);
 free(data);
}
