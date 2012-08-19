#include <R.h>
#include <Rmath.h>

void rrunif(int *x, int *n)
{
    double a = 0;
    double b;
    double tx;
    b=*n;
    GetRNGstate();
    tx=0.0+runif(a, b);
    x[0]=(int) tx;
    PutRNGstate();
}
