#include<ctype.h>
#include "../gmx_lapack.h"

/* LAPACK */
void
F77_FUNC(dlacpy,DLACPY)(const char *uplo,
	int *m,
	int *n,
	double *a,
	int *lda,
	double *b,
	int *ldb)
{
  int i,j,minjm;
  const char ch=toupper(*uplo);

  if(ch=='U') {
    for(j=0;j<*n;j++) {
      minjm = (j < (*m-1)) ? j : (*m-1);
      for(i=0;i<=minjm;i++)
	b[j*(*ldb)+i] = a[j*(*lda)+i];
    }
  } else if(ch=='L') {
    for(j=0;j<*n;j++) {
      for(i=j;i<*m;i++)
	b[j*(*ldb)+i] = a[j*(*lda)+i];
    }
  } else {
    for(j=0;j<*n;j++) {
      for(i=0;i<*m;i++)
	b[j*(*ldb)+i] = a[j*(*lda)+i];
    }    
  }
}
