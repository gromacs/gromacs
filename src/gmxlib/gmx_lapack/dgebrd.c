#include "gmx_lapack.h"
#include "gmx_blas.h"
#include "lapack_limits.h"


void
F77_FUNC(dgebrd,DGEBRD)(int *m, 
	int *n, 
	double *a, 
	int *lda, 
	double *d__, 
	double *e,
	double *tauq, 
	double *taup,
	double *work, 
	int *lwork,
	int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, i_1, i_2, i_3, i_4;

    /* Local variables */
    int i_, j, nx,nb;
    double ws;
    int nbmin, iinfo, minmn;
    int ldwrkx, ldwrky;
    double one = 1.0;
    double minusone = -1.0;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --d__;
    --e;
    --tauq;
    --taup;
    --work;

    nb = DGEBRD_BLOCKSIZE;
    *info = 0;
    if (*lwork==-1) {
      work[1] = (double) ( (*m + *n) * nb);
      return;
    }
    minmn = (*m < *n) ? *m : *n;
    if (minmn == 0) {
      work[1] = 1.;
      return;
    }

    ws = (*m > *n) ? *m : *n;
    ldwrkx = *m;
    ldwrky = *n;

    if (nb > 1 && nb < minmn) {
	nx = DGEBRD_CROSSOVER;
	if (nx < minmn) {
	    ws = (double) ((*m + *n) * nb);
	    if ((double) (*lwork) < ws) {
	      nbmin = DGEBRD_MINBLOCKSIZE;
		if (*lwork >= (*m + *n) * nbmin) {
		    nb = *lwork / (*m + *n);
		} else {
		    nb = 1;
		    nx = minmn;
		}
	    }
	}
    } else {
	nx = minmn;
    }

    i_1 = minmn - nx;
    i_2 = nb;
    for (i_ = 1; i_2 < 0 ? i_ >= i_1 : i_ <= i_1; i_ += i_2) {

	i_3 = *m - i_ + 1;
	i_4 = *n - i_ + 1;
	F77_FUNC(dlabrd,DLABRD)(&i_3, &i_4, &nb, &a[i_ + i_ * a_dim1], lda, &d__[i_], 
		&e[i_], &tauq[i_], &taup[i_], &work[1], &ldwrkx, 
		&work[ldwrkx * nb + 1], &ldwrky);

	i_3 = *m - i_ - nb + 1;
	i_4 = *n - i_ - nb + 1;
	F77_FUNC(dgemm,DGEMM)("N", "T", &i_3, &i_4, &nb, &minusone, 
	       &a[i_ + nb + i_ * a_dim1], lda, &work[ldwrkx * nb + nb + 1],
	       &ldwrky, &one, &a[i_ + nb + (i_ + nb) * a_dim1], lda);
	i_3 = *m - i_ - nb + 1;
	i_4 = *n - i_ - nb + 1;
	F77_FUNC(dgemm,DGEMM)("N", "N", &i_3, &i_4, &nb, &minusone, &work[nb + 1], &ldwrkx,
	       &a[i_ + (i_ + nb) * a_dim1], lda, &one, 
	       &a[i_ + nb + (i_ + nb) * a_dim1], lda);

	if (*m >= *n) {
	    i_3 = i_ + nb - 1;
	    for (j = i_; j <= i_3; ++j) {
		a[j + j * a_dim1] = d__[j];
		a[j + (j + 1) * a_dim1] = e[j];
	    }
	} else {
	    i_3 = i_ + nb - 1;
	    for (j = i_; j <= i_3; ++j) {
		a[j + j * a_dim1] = d__[j];
		a[j + 1 + j * a_dim1] = e[j];
	    }
	}
    }

    i_2 = *m - i_ + 1;
    i_1 = *n - i_ + 1;
    F77_FUNC(dgebd2,DGEBD2)(&i_2, &i_1, &a[i_ + i_ * a_dim1], lda, &d__[i_], &e[i_], &
	    tauq[i_], &taup[i_], &work[1], &iinfo);
    work[1] = ws;
    return;

}
