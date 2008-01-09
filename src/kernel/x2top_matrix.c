#include <stdio.h>
#include "typedefs.h"
#include "smalloc.h"
#include "x2top_matrix.h"
#include "gmx_lapack.h"

/*extern void gaussj(double **a,int n,double **b,int m);*/
	
double **alloc_matrix(int n,int m)
{
  double **ptr,*ptr0;
  int i;
  
  /* There's always time for more pointer arithmetic! */
  snew(ptr,n);
  snew(ptr0,n*n);
  for(i=0; (i<n); i++) 
    ptr[i] = ptr0 + i*n;
  return ptr;
}

void free_matrix(double **a,int n)
{
  sfree(a[0]);
  sfree(a);
}

void mat_mult(FILE *fp,int n,int m,double **x,double **y,double **z)
{
  int i,j,k;
  
  if (fp)
    fprintf(fp,"Multiplying %d x %d matrix with a %d x %d matrix\n",n,m,m,n);
  if (fp)
    for(i=0; (i<n); i++) {
      for(j=0; (j<m); j++) 
	fprintf(fp,"%8.1e",x[i][j]);
      fprintf(fp,"\n");
    }
  for(i=0; (i<m); i++) {
    for(j=0; (j<m); j++) {
      z[i][j] = 0;
      for(k=0; (k<n); k++)
	z[i][j] += x[k][i]*y[j][k];
    }
  }
}

void mat_inv(FILE *fp,int n,double **a)
{
  int i,j,m,lda,*ipiv,lwork,info;
  double **test,**id,**b,*work;

  test = alloc_matrix(n,n);  
  b = alloc_matrix(n,n);
  for(i=0; (i<n); i++) 
    b[i][i] = 1.0;
  if (fp)
    fprintf(fp,"Inverting %d square matrix\n",n);
  for(i=0; (i<n); i++) {
    b[i][i] = 1;
    for(j=0; (j<n); j++)  {
      if (fp) 
	fprintf(fp,"%8.1f",a[i][j]);
      test[i][j] = a[i][j];
    }
    if (fp)
      fprintf(fp,"\n");
  }
  /*  gaussj(a,n,b,n);*/
  snew(ipiv,n);
  lwork = n*n;
  snew(work,lwork);
  m = lda   = n;
  info  = 0;
  F77_FUNC(dgetrf,DGETRF)(&n,&m,a[0],&lda,ipiv,&info);
  if (info != 0)
    gmx_fatal(FARGS,"LU decomposition failed. Info = %d",info);
  F77_FUNC(dgetri,DGETRI)(&n,a[0],&lda,ipiv,work,&lwork,&info);
  if (info != 0)
    gmx_fatal(FARGS,"Matrix inversion failed. Info = %d",info);
    
  if (fp) {
    fprintf(fp,"And here is the result:\n");
    for(i=0; (i<n); i++) {
      for(j=0; (j<n); j++) 
	fprintf(fp,"  %8.1e",a[i][j]);
      fprintf(fp,"\n");
    }
  }
  id = alloc_matrix(n,n);
  mat_mult(fp,n,n,test,a,id);
  if (fp) {
    fprintf(fp,"And here is the product of A and Ainv\n");
    for(i=0; (i<n); i++) {
      for(j=0; (j<n); j++) 
	fprintf(fp,"  %8.1e",id[i][j]);
      fprintf(fp,"\n");
    }
  }
}

