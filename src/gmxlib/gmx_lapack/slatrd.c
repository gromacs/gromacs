#include <ctype.h>
#include "gmx_blas.h"
#include "gmx_lapack.h"
#include "lapack_limits.h"


void
F77_FUNC(slatrd,SLATRD)(char *  uplo,
       int  *   n,
       int  *   nb,
       float * a,
       int *    lda,
       float * e,
       float * tau,
       float * w,
       int *    ldw)
{
  int i,iw;
  int ti1,ti2,ti3;
  float one,zero,minusone,alpha;
  char ch=toupper(*uplo);

  one=1.0;
  minusone=-1.0;
  zero=0.0;

  if(*n<=0)
    return;

  if(ch=='U') {
    for(i=*n;i>=(*n-*nb+1);i--) {
      iw = i -*n + *nb;
      
      if(i<*n) {
	ti1 = *n-i;
	ti2 = 1;
	/* BLAS */
	F77_FUNC(sgemv,SGEMV)("N",&i,&ti1,&minusone, &(a[ i*(*lda) + 0]),lda,&(w[iw*(*ldw)+(i-1)]),
	       ldw,&one, &(a[ (i-1)*(*lda) + 0]), &ti2);
	/* BLAS */
	F77_FUNC(sgemv,SGEMV)("N",&i,&ti1,&minusone, &(w[ iw*(*ldw) + 0]),ldw,&(a[i*(*lda)+(i-1)]),
	       lda,&one, &(a[ (i-1)*(*lda) + 0]), &ti2);
      }

      if(i>1) {
	/*  Generate elementary reflector H(i) to annihilate
	 *              A(1:i-2,i) 
	 */
	ti1 = i-1;
	ti2 = 1;

	/* LAPACK */
	F77_FUNC(slarfg,SLARFG)(&ti1,&(a[(i-1)*(*lda)+(i-2)]),&(a[(i-1)*(*lda)+0]),&ti2,&(tau[i-2]));
      
	e[i-2] = a[(i-1)*(*lda)+(i-2)];
	a[(i-1)*(*lda)+(i-2)] = 1.0;

	/* Compute W(1:i-1,i) */
	ti1 = i-1;
	ti2 = 1;

	/* BLAS */
	F77_FUNC(ssymv,SSYMV)("U",&ti1,&one,a,lda,&(a[(i-1)*(*lda)+0]),&ti2,&zero,
	       &(w[(iw-1)*(*ldw)+0]),&ti2);
	if(i<*n) {
	  ti1 = i-1;
	  ti2 = *n-i;
	  ti3 = 1;
	  /* BLAS */
	  F77_FUNC(sgemv,SGEMV)("T",&ti1,&ti2,&one,&(w[iw*(*ldw)+0]),ldw,&(a[(i-1)*(*lda)+0]),&ti3,
		 &zero,&(w[(iw-1)*(*ldw)+i]),&ti3);
	
	  /* BLAS */
	  F77_FUNC(sgemv,SGEMV)("N",&ti1,&ti2,&minusone,&(a[i*(*lda)+0]),lda,&(w[(iw-1)*(*ldw)+i]),&ti3,
		 &one,&(w[(iw-1)*(*ldw)+0]),&ti3);
	
	  /* BLAS */
	  F77_FUNC(sgemv,SGEMV)("T",&ti1,&ti2,&one,&(a[i*(*lda)+0]),lda,&(a[(i-1)*(*lda)+0]),&ti3,
		 &zero,&(w[(iw-1)*(*ldw)+i]),&ti3);
	
	  /* BLAS */
	  F77_FUNC(sgemv,SGEMV)("N",&ti1,&ti2,&minusone,&(w[iw*(*ldw)+0]),ldw,&(w[(iw-1)*(*ldw)+i]),&ti3,
		 &one,&(w[(iw-1)*(*ldw)+0]),&ti3);
	}
      
	ti1 = i-1;
	ti2 = 1;
	/* BLAS */
	F77_FUNC(sscal,SSCAL)(&ti1,&(tau[i-2]),&(w[(iw-1)*(*ldw)+0]),&ti2);
      
	alpha = -0.5*tau[i-2]*F77_FUNC(sdot,SDOT)(&ti1,&(w[(iw-1)*(*ldw)+0]),&ti2,
				    &(a[(i-1)*(*lda)+0]),&ti2);
      
	/* BLAS */
	F77_FUNC(saxpy,SAXPY)(&ti1,&alpha,&(a[(i-1)*(*lda)+0]),&ti2,&(w[(iw-1)*(*ldw)+0]),&ti2);

      }
    }
  } else {
    /* lower */
    for(i=1;i<=*nb;i++) {

      ti1 = *n-i+1;
      ti2 = i-1;
      ti3 = 1;
      /* BLAS */
      F77_FUNC(sgemv,SGEMV)("N",&ti1,&ti2,&minusone, &(a[ i-1 ]),lda,&(w[ i-1 ]),
	       ldw,&one, &(a[ (i-1)*(*lda) + (i-1)]), &ti3);
      /* BLAS */
      F77_FUNC(sgemv,SGEMV)("N",&ti1,&ti2,&minusone, &(w[ i-1 ]),ldw,&(a[ i-1 ]),
	       lda,&one, &(a[ (i-1)*(*lda) + (i-1)]), &ti3);

      if(i<*n) {
	ti1 = *n - i;
	ti2 = (*n < i+2 ) ? *n : (i+2);
	ti3 = 1;
	/* LAPACK */
	F77_FUNC(slarfg,SLARFG)(&ti1,&(a[(i-1)*(*lda)+(i)]),&(a[(i-1)*(*lda)+(ti2-1)]),&ti3,&(tau[i-1]));
	e[i-1] = a[(i-1)*(*lda)+(i)];
	a[(i-1)*(*lda)+(i)] = 1.0;
	
	ti1 = *n - i;
	ti2 = 1;
	F77_FUNC(ssymv,SSYMV)("L",&ti1,&one,&(a[i*(*lda)+i]),lda,&(a[(i-1)*(*lda)+i]),&ti2,
	       &zero,&(w[(i-1)*(*ldw)+i]),&ti2);
	ti1 = *n - i;
	ti2 = i-1;
	ti3 = 1;
	/* BLAS */
	F77_FUNC(sgemv,SGEMV)("T",&ti1,&ti2,&one,&(w[ i ]),ldw,&(a[(i-1)*(*lda)+i]),&ti3,
	       &zero,&(w[(i-1)*(*ldw)+0]),&ti3);
	
	/* BLAS */
	F77_FUNC(sgemv,SGEMV)("N",&ti1,&ti2,&minusone,&(a[ i ]),lda,&(w[(i-1)*(*ldw)+0]),&ti3,
	       &one,&(w[(i-1)*(*ldw)+i]),&ti3);
	
	/* BLAS */
	F77_FUNC(sgemv,SGEMV)("T",&ti1,&ti2,&one,&(a[ i ]),lda,&(a[(i-1)*(*lda)+i]),&ti3,
	       &zero,&(w[(i-1)*(*ldw)+0]),&ti3);
	
	/* BLAS */
	F77_FUNC(sgemv,SGEMV)("N",&ti1,&ti2,&minusone,&(w[ i ]),ldw,&(w[(i-1)*(*ldw)+0]),&ti3,
	       &one,&(w[(i-1)*(*ldw)+i]),&ti3);

	F77_FUNC(sscal,SSCAL)(&ti1,&(tau[i-1]),&(w[(i-1)*(*ldw)+i]),&ti3);
	alpha = -0.5*tau[i-1]*F77_FUNC(sdot,SDOT)(&ti1,&(w[(i-1)*(*ldw)+i]),&ti3,
				   &(a[(i-1)*(*lda)+i]),&ti3);
	
	F77_FUNC(saxpy,SAXPY)(&ti1,&alpha,&(a[(i-1)*(*lda)+i]),&ti3,&(w[(i-1)*(*ldw)+i]),&ti3);
      }
    }
  }
  return;
}
	


  
