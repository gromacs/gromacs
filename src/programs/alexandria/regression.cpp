#include "regression.h"
    
#include <cmath>
#include <cstdio>
#include <cstring>

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/smalloc.h"

extern "C" void dgelsd_(int* m, int* n, int* nrhs, double* a, int* lda,
                        double* b, int* ldb, double* s, double* rcond, int* rank,
                        double* work, int* lwork, int* iwork, int* info );

void multi_regression2(int nrow, double y[], int ncol,
                       double **a, double x[])
{
    /* Query and allocate the optimal workspace */
    int     lwork = -1;
    int     lda   = nrow;
    int     ldb   = nrow;
    int     nrhs  = 1;
    int     rank;
    double  rcond = -1.0;
    double  wkopt;
    double *s;
    snew(s, nrow);
    // Compute length of integer array iwork according to
    // https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_lapack_examples/dgelsd_ex.c.htm
    int  smlsiz = 25;
    int  nlvl   = std::max(0L, std::lround(std::log2(std::min(nrow, ncol)/smlsiz + 1) ) + 1);
    int  liwork = 3*std::min(ncol,nrow)*nlvl + 11*std::min(nrow, ncol);
    int *iwork;
    snew(iwork, liwork);
    int info;
    dgelsd_ (&nrow, &ncol, &nrhs, a[0], &lda, y, &ldb, s, &rcond, &rank, &wkopt, &lwork,
             iwork, &info );
    lwork = (int)wkopt;
    double *work;
    snew(work, lwork);
    /* Solve the equations A*X = B */
    dgelsd_ (&nrow, &ncol, &nrhs, a[0], &lda, y, &ldb, s, &rcond, &rank, work, &lwork,
             iwork, &info );
    /* Check for convergence */
    if( info > 0 ) 
    {
        char buf[256];
        snprintf(buf, sizeof(buf), "The algorithm computing SVD failed to converge and the least squares solution could not be computed. Info = %d.", info);
        GMX_THROW(gmx::InvalidInputError(buf));
    }
    for(int i = 0; i < ncol; i++)
    {
        x[i] = y[i];
    }
 
    sfree(s);
    sfree(work);
    sfree(iwork);
}

