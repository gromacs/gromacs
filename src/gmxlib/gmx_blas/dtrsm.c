#include <ctype.h>
#include <math.h>

#include <types/simple.h>
#include "gmx_blas.h"

void
F77_FUNC(dtrsm,DTRSM)(char * side,
       char * uplo,
       char * transa,
       char * diag,
       int *  m__,
       int *  n__,
       double *alpha__,
       double *a,
       int *  lda__,
       double *b,
       int *  ldb__)
{
  char xside  = toupper(*side);
  char xuplo  = toupper(*uplo);
  char xtrans = toupper(*transa);
  char xdiag  = toupper(*diag);
  int i,j,k;
  double temp;

  
  int m = *m__;
  int n = *n__;
  int lda = *lda__;
  int ldb = *ldb__;
  double alpha = *alpha__;

  if(n<=0)
    return;
  
  if(fabs(alpha)<GMX_DOUBLE_MIN) { 
    for(j=0;j<n;j++)
      for(i=0;i<m;i++)
	b[j*(ldb)+i] = 0.0;
    return;
  }

  if(xside=='L') {
    /* left side */
    if(xtrans=='N') {
      /* No transpose */
      if(xuplo=='U') {
	/* upper */
	for(j=0;j<n;j++) {
	  if(fabs(alpha-1.0)>GMX_DOUBLE_EPS) {
	    for(i=0;i<m;i++)
	      b[j*(ldb)+i] *= alpha;
	  }
	  for(k=m-1;k>=0;k--) {
	    if(fabs(b[j*(ldb)+k])>GMX_DOUBLE_MIN) {
	      if(xdiag=='N')
		b[j*(ldb)+k] /= a[k*(lda)+k];
	      for(i=0;i<k;i++)
		b[j*(ldb)+i] -= b[j*(ldb)+k]*a[k*(lda)+i];
	    }
	  }
	}
      } else {
	/* lower */
	for(j=0;j<n;j++) {
	  if(fabs(alpha-1.0)>GMX_DOUBLE_EPS)
	    for(i=0;i<m;i++)
	      b[j*(ldb)+i] *= alpha;
	  for(k=0;k<m;k++) {
	    if(fabs(b[j*(ldb)+k])>GMX_DOUBLE_MIN) {
	      if(xdiag=='N')
		b[j*(ldb)+k] /= a[k*(lda)+k];
	      for(i=k+1;i<m;i++)
		b[j*(ldb)+i] -= b[j*(ldb)+k]*a[k*(lda)+i];
	    }
	  }
	}
      }
    } else {
      /* Transpose */
      if(xuplo=='U') {
	/* upper */
	for(j=0;j<n;j++) {
	  for(i=0;i<m;i++) {
	    temp = alpha * b[j*(ldb)+i];
	    for(k=0;k<i;k++)
	      temp -= a[i*(lda)+k] * b[j*(ldb)+k];
	    if(xdiag=='N')
		temp /= a[i*(lda)+i];
	    b[j*(ldb)+i] = temp;
	  }
	}
      } else {
	/* lower */
	for(j=0;j<n;j++) {
	  for(i=m-1;i>=0;i--) {
	    temp = alpha * b[j*(ldb)+i];
	    for(k=i+1;k<m;k++)
	      temp -= a[i*(lda)+k] * b[j*(ldb)+k];
	    if(xdiag=='N')
		temp /= a[i*(lda)+i];
	    b[j*(ldb)+i] = temp;
	  }
	}
      }
    }
  } else {
    /* right side */
    if(xtrans=='N') {
      /* No transpose */
      if(xuplo=='U') {
	/* upper */
	for(j=0;j<n;j++) {
	  if(fabs(alpha-1.0)>GMX_DOUBLE_EPS)
	    for(i=0;i<m;i++)
	      b[j*(ldb)+i] *= alpha;
	  for(k=0;k<j;k++) {
	    if(fabs(a[j*(lda)+k])>GMX_DOUBLE_MIN) {
	      for(i=0;i<m;i++)
		b[j*(ldb)+i] -= a[j*(lda)+k]*b[k*(ldb)+i];
	    }
	  }
	  if(xdiag=='N') {
	    temp = 1.0/a[j*(lda)+j];
	    for(i=0;i<m;i++)
	      b[j*(ldb)+i] *= temp;
	  }
	}
      } else {
	/* lower */
	for(j=n-1;j>=0;j--) {
	  if(fabs(alpha)>GMX_DOUBLE_MIN)
	    for(i=0;i<m;i++)
	      b[j*(ldb)+i] *= alpha;
	  for(k=j+1;k<n;k++) {
	    if(fabs(a[j*(lda)+k])>GMX_DOUBLE_MIN) {
	      for(i=0;i<m;i++)
		b[j*(ldb)+i] -= a[j*(lda)+k]*b[k*(ldb)+i];
	    }
	  }
	  if(xdiag=='N') {
	    temp = 1.0/a[j*(lda)+j];
	    for(i=0;i<m;i++)
	      b[j*(ldb)+i] *= temp;
	  }
	}
      }
    } else {
      /* Transpose */
      if(xuplo=='U') {
	/* upper */
	for(k=n-1;k>=0;k--) {
	  if(xdiag=='N') {
	    temp = 1.0/a[k*(lda)+k];
	    for(i=0;i<m;i++)
	      b[k*(ldb)+i] *= temp;
	  }
	  for(j=0;j<k;j++) {
	    if(fabs(a[k*(lda)+j])>GMX_DOUBLE_MIN) {
	      temp = a[k*(lda)+j];
	      for(i=0;i<m;i++)
		b[j*(ldb)+i] -= temp * b[k*(ldb)+i];
	    }
	  }
	  if(fabs(alpha-1.0)>GMX_DOUBLE_EPS)
	    for(i=0;i<m;i++)
	      b[k*(ldb)+i] *= alpha;
	}
      } else {
	/* lower */
	for(k=0;k<n;k++) {
	  if(xdiag=='N') {
	    temp = 1.0/a[k*(lda)+k];
	    for(i=0;i<m;i++)
	      b[k*(ldb)+i] *= temp;
	  }
	  for(j=k+1;j<n;j++) {
	    if(fabs(a[k*(lda)+j])>GMX_DOUBLE_MIN) {
	      temp = a[k*(lda)+j];
	      for(i=0;i<m;i++)
		b[j*(ldb)+i] -= temp * b[k*(ldb)+i];
	    }
	  }
	  if(fabs(alpha-1.0)>GMX_DOUBLE_EPS)
	    for(i=0;i<m;i++)
	      b[k*(ldb)+i] *= alpha;
	}
      }      
    }
  }    
}
