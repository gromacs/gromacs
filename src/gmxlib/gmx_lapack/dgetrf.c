#include "gmx_blas.h"
#include "gmx_lapack.h"
#include "lapack_limits.h"

void
F77_FUNC(dgetrf,DGETRF)(int *m,
	int *n,
	double *a,
	int *lda,
	int *ipiv,
	int *info)
{
  int mindim,jb;
  int i,j,k,l;
  int iinfo;
  double minusone = -1.0;
  double one = 1.0;

  if(*m<=0 || *n<=0)
    return;

  *info = 0;

  mindim = (*m < *n) ? *m : *n;

  if(DGETRF_BLOCKSIZE>=mindim) {

    /* unblocked code */
    F77_FUNC(dgetf2,DGETF2)(m,n,a,lda,ipiv,info);

  } else {

    /* blocked case */

    for(j=1;j<=mindim;j+=DGETRF_BLOCKSIZE) {
      jb = ( DGETRF_BLOCKSIZE < (mindim-j+1)) ? DGETRF_BLOCKSIZE : (mindim-j+1);
      /* factor diag. and subdiag blocks and test for singularity */
      k = *m-j+1;
      F77_FUNC(dgetf2,DGETF2)(&k,&jb,&(a[(j-1)*(*lda)+(j-1)]),lda,&(ipiv[j-1]),&iinfo);
      
      if(*info==0 && iinfo>0)
	*info = iinfo + j - 1;

      /* adjust pivot indices */
      k = (*m < (j+jb-1)) ? *m : (j+jb-1);
      for(i=j;i<=k;i++)
	ipiv[i-1] += j - 1;

      /* Apply to columns 1 throughj j-1 */
      k = j - 1;
      i = j + jb - 1;
      l = 1;
      F77_FUNC(dlaswp,DLASWP)(&k,a,lda,&j,&i,ipiv,&l);
      if((j+jb)<=*n) {
	/* Apply to cols. j+jb through n */
	k = *n-j-jb+1;
	i = j+jb-1;
	l = 1;
	F77_FUNC(dlaswp,DLASWP)(&k,&(a[(j+jb-1)*(*lda)+0]),lda,&j,&i,ipiv,&l);
	/* Compute block row of U */
	k = *n-j-jb+1;
	F77_FUNC(dtrsm,DTRSM)("Left","Lower","No transpose","Unit",&jb,&k,&one,
	       &(a[(j-1)*(*lda)+(j-1)]),lda,&(a[(j+jb-1)*(*lda)+(j-1)]),lda);

	if((j+jb)<=*m) {
	  /* Update trailing submatrix */
	  k = *m-j-jb+1;
	  i = *n-j-jb+1;
	  F77_FUNC(dgemm,DGEMM)("No transpose","No transpose",&k,&i,&jb,&minusone,
		 &(a[(j-1)*(*lda)+(j+jb-1)]),lda,
		 &(a[(j+jb-1)*(*lda)+(j-1)]),lda,&one,
		 &(a[(j+jb-1)*(*lda)+(j+jb-1)]),lda);
	}

      }
    }
  }
}
