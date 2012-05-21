#include "../gmx_lapack.h"

void
F77_FUNC(dgebd2,DGEBD2)(int *m,
	int *n,
	double *a,
	int *lda,
	double *d,
	double *e,
	double *tauq,
	double *taup,
	double *work,
	int *info)
{
  int i,i1,i2,i3;
    
    *info = 0;

  if(*m>=*n) {
    /* reduce to upper bidiag. form */
    for(i=0;i<*n;i++) {
      i1 = *m - i;
      i2 = ( (i+1) < (*m-1)) ? (i+1) : (*m-1);
      i3 = 1;
      F77_FUNC(dlarfg,DLARFG)(&i1,&(a[i*(*lda)+i]),&(a[i*(*lda)+i2]),&i3,&(tauq[i]));
      d[i] = a[i*(*lda)+i];
      a[i*(*lda)+i] = 1.0;
      i2 = *n - i - 1;
      F77_FUNC(dlarf,DLARF)("L",&i1,&i2,&(a[i*(*lda)+i]),&i3,&(tauq[i]),&(a[(i+1)*(*lda)+i]),lda,work);
      a[i*(*lda)+i] = d[i];

      if(i<(*n-1)) {

	i1 = *n - i -1;
	i2 = ( (i+2) < (*n-1)) ? (i+2) : (*n-1); 
	F77_FUNC(dlarfg,DLARFG)(&i1,&(a[(i+1)*(*lda)+i]),&(a[i2*(*lda)+i]),lda,&(taup[i]));

	e[i] = a[(i+1)*(*lda)+i];
	a[(i+1)*(*lda)+i] = 1.0;

	i1 = *m - i - 1;
	i2 = *n - i - 1;
	F77_FUNC(dlarf,DLARF)("R",&i1,&i2,&(a[(i+1)*(*lda)+i]),lda,&(taup[i]),&(a[(i+1)*(*lda)+i+1]),lda,work);
	a[(i+1)*(*lda)+i] = e[i];
      } else
	taup[i] = 0.0;
    }
  } else {
    /* reduce to lower bidiag. form */
    for(i=0;i<*m;i++) {
      i1 = *n - i;
      i2 = ( (i+1) < (*n-1)) ? (i+1) : (*n-1);
      i3 = 1;
      F77_FUNC(dlarfg,DLARFG)(&i1,&(a[i*(*lda)+i]),&(a[i2*(*lda)+i]),lda,&(taup[i]));
      d[i] = a[i*(*lda)+i];
      a[i*(*lda)+i] = 1.0;

      i2 = *m - i - 1;
      i3 = ( (i+1) < (*m-1)) ? (i+1) : (*m-1);
      F77_FUNC(dlarf,DLARF)("R",&i2,&i1,&(a[i*(*lda)+i]),lda,&(taup[i]),&(a[(i)*(*lda)+i3]),lda,work);
      a[i*(*lda)+i] = d[i];

      if(i<(*m-1)) {

	i1 = *m - i - 1;
	i2 = ( (i+2) < (*m-1)) ? (i+2) : (*m-1);
	i3 = 1;
	F77_FUNC(dlarfg,DLARFG)(&i1,&(a[(i)*(*lda)+i+1]),&(a[i*(*lda)+i2]),&i3,&(tauq[i]));

	e[i] = a[(i)*(*lda)+i+1];
	a[(i)*(*lda)+i+1] = 1.0;

	i1 = *m - i - 1;
	i2 = *n - i - 1;
	i3 = 1;
	F77_FUNC(dlarf,DLARF)("L",&i1,&i2,&(a[(i)*(*lda)+i+1]),&i3,&(tauq[i]),&(a[(i+1)*(*lda)+i+1]),lda,work);
	a[(i)*(*lda)+i+1] = e[i];
      } else
	tauq[i] = 0.0;
    }
  }
  return;
}
