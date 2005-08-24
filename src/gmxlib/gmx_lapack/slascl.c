#include <ctype.h>
#include <math.h>
#include <types/simple.h>

#include "gmx_lapack.h"
#include "lapack_limits.h"


void
F77_FUNC(slascl,SLASCL)(char *type,
                        int *kl,
                        int *ku,
                        float *cfrom,
                        float *cto,
                        int *m,
                        int *n,
                        float *a,
                        int *lda,
                        int *info)
{
  char ch=toupper(*type);
  int i,j,k,l,k1,k2,k3,k4;
  int done=0;
  float minval,smlnum,bignum;
  float cfromc, ctoc, cfrom1, cto1, mul;

  if(*n<=0 || *m<=0)
    return;

  minval = GMX_FLOAT_MIN;
  smlnum = minval / GMX_FLOAT_EPS;
  bignum = 1.0 / smlnum;

  cfromc = *cfrom;
  ctoc   = *cto;

  while(!done) {
    
    cfrom1 = cfromc * smlnum;
    cto1   = ctoc / bignum;

    if(fabs(cfrom1)>fabs(ctoc) && fabs(ctoc)>GMX_FLOAT_MIN) {
      mul = smlnum;
      done = 0;
      cfromc = cfrom1;
    } else if(fabs(cto1)>fabs(cfromc)) {
      mul = bignum;
      done = 0;
      ctoc = cto1;
    } else {
      mul = ctoc / cfromc;
      done = 1;
    }

    switch(ch) {
    case 'G': 
      /* Full matrix */
      for(j=0;j<*n;j++)
	for(i=0;i<*m;i++)
	  a[j*(*lda)+i] *= mul;
      break;

    case 'L': 
      /* Lower triangular matrix */
      for(j=0;j<*n;j++)
	for(i=j;i<*m;i++)
	  a[j*(*lda)+i] *= mul;
      break;

    case 'U': 
      /* Upper triangular matrix */
      for(j=0;j<*n;j++) {
	k = (j < (*m-1)) ? j : (*m-1);
	for(i=0;i<=k;i++)
	  a[j*(*lda)+i] *= mul;
      }
      break;

    case 'H': 
      /* Upper Hessenberg matrix */
      for(j=0;j<*n;j++) {
	k = ((j+1) < (*m-1)) ? (j+1) : (*m-1);
	for(i=0;i<=k;i++)
	  a[j*(*lda)+i] *= mul;
      }
      break;

    case 'B': 
      /* Symmetric band matrix, lower bandwidth KL, upper KU,
       * only the lower half stored.
       */
      k3 = *kl;
      k4 = *n - 1;
      for(j=0;j<*n;j++) {
	k = (k3 < (k4-j)) ? k3 : (k4-j);
	for(i=0;i<=k;i++)
	  a[j*(*lda)+i] *= mul;
      }
      break;

    case 'Q': 
      /* Symmetric band matrix, lower bandwidth KL, upper KU,
       * only the upper half stored.
       */
      k1 = *ku;
      k3 = *ku;
      for(j=0;j<*n;j++) {
	k = ((k1-j) > 0) ? (k1-j) : 0;
	for(i=k;i<=k3;i++)
	  a[j*(*lda)+i] *= mul;
      }
      break;

    case 'Z': 
      /* Band matrix, lower bandwidth KL, upper KU. */

      k1 = *kl + *ku;
      k2 = *kl;
      k3 = 2*(*kl) + *ku;
      k4 = *kl + *ku - 1 + *m;
      for(j=0;j<*n;j++) {
	k = ((k1-j) > k2) ? (k1-j) : k2;
	l = (k3 < (k4-j)) ? k3 : (k4-j);
	for(i=k;i<=l;i++)
	  a[j*(*lda)+i] *= mul;
      }
      break;

    default:
      *info = -1;
      return;
    }
  } /* finished */

  *info = 0;
  return;
}
