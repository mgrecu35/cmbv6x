#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;


extern "C" void interp_arm_(float *x, float *y, int *n,
                            float *xi, float *yi, int *ni)
{
  int i;
  vec xa=vec(*n);
  vec ya=vec(*n);
  vec xia=vec(*ni);
  vec yia=vec(*ni);
  for(i=0;i<*n;i++)
    {
      xa(i)=x[i];
      ya(i)=y[i];
    }
  for(i=0;i<*ni;i++)
    xia(i)=xi[i];
  
  
  interp1(xa, ya, xia, yia);
  for(i=0;i<*ni;i++)
    {
      yi[i]=yia(i);
      printf("%i %g \n",i,yia(i));
    }
}

extern "C" void kgainc_(float *dtb, int *n, float *s, float *kgain)
{
  mat A(*n,*n);
  int i,j;
  for(i=0;i<*n;i++)
    {
      for(j=0;j<*n;j++)
        A(i,j)=dtb[i]*dtb[j];
      A(i,i)=A(i,i)+(*s);
    }
  //A.print("A=:");
  mat B=pinv(A,0.0001);
  //B.print("B=:");
  
  for(i=0;i<*n;i++)
    {
      kgain[i]=0;
      for(j=0;j<*n;j++)
        kgain[i]+=dtb[j]*B(j,i);
    }
}
//pinv=linalg.pinv(dot(dtb.T,dtb)+eye(6)*4)
//  kgain=dot(dtb,pinv)

extern "C" void interp_armi_(int *x, float *y, int *n,
                             int *xi, float *yi, int *ni)
{
  int i;
  vec xa=vec(*n);
  vec ya=vec(*n);
  vec xia=vec(*ni);
  vec yia;
  for(i=0;i<*n;i++)
    {
      xa(i)=x[i];
      ya(i)=y[i];
    }
  for(i=0;i<*ni;i++)
    {
      xia(i)=xi[i];
      //printf("%g \n",xi[i]);
    }
  //  xia.print("xi:");
  
  interp1(xa, ya, xia, yia);
  //printf(" %i %i \n",*n,*ni);
  for(i=0;i<*ni;i++)
    {
      yi[i]=yia(i);
      //  printf("%i %lg %lg \n",i,xia(i),yia(i));
    }
}
