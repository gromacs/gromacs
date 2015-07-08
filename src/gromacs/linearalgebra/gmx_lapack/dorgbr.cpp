#include "../gmx_lapack.h"
#include "lapack_limits.h"

void
F77_FUNC(dorgbr,DORGBR)(const char *vect,
	int *m,
	int *n,
	int *k,
	double *a,
	int *lda,
	double *tau,
	double *work,
	int *lwork,
	int *info)
{
  int wantq,iinfo,j,i,i1,wrksz;
  int mn = (*m < *n) ? *m : *n;

  wantq = (*vect=='Q' || *vect=='q');

  *info = 0;
  wrksz = mn*DORGBR_BLOCKSIZE;
  if(*lwork==-1) {
    work[0] = wrksz;
    return;
  }
  
  if(*m==0 || *n==0)
    return;

  if(wantq) {
    if(*m>=*k)
      F77_FUNC(dorgqr,DORGQR)(m,n,k,a,lda,tau,work,lwork,&iinfo);
    else {
      for(j=*m;j>=2;j--) {
	a[(j-1)*(*lda)+0] = 0.0;
	for(i=j+1;i<=*m;i++)
	  a[(j-1)*(*lda)+(i-1)] = a[(j-2)*(*lda)+(i-1)]; 
      }
      a[0] = 1.0;
      for(i=2;i<=*m;i++)
	a[i-1] = 0.0;
      if(*m>1) {
	i1 = *m-1;
	F77_FUNC(dorgqr,DORGQR)(&i1,&i1,&i1,&(a[*lda+1]),lda,tau,work,lwork,&iinfo);
      }
    }
  } else {
    if(*k<*n)
      F77_FUNC(dorglq,DORGLQ)(m,n,k,a,lda,tau,work,lwork,&iinfo);
    else {
      a[0] = 1.0;
      for(i=2;i<=*m;i++)
	a[i-1] = 0.0;
      for(j=2;j<=*n;j++) {
	for(i=j-1;i>=2;i--)
	  a[(j-1)*(*lda)+(i-1)] = a[(j-1)*(*lda)+(i-2)]; 
	a[(j-1)*(*lda)+0] = 0.0;
      }
      if(*n>1) {
	i1 = *n-1;
	F77_FUNC(dorglq,DORGLQ)(&i1,&i1,&i1,&(a[*lda+1]),lda,tau,work,lwork,&iinfo);
      }
    }
  }
  work[0] = wrksz;
  return;
}
 
