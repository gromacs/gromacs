/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "types/simple.h"

static void nrerror(const char error_text[], gmx_bool bExit)
{
    fprintf(stderr, "Numerical Recipes run-time error...\n");
    fprintf(stderr, "%s\n", error_text);
    if (bExit)
    {
        fprintf(stderr, "...now exiting to system...\n");
        exit(1);
    }
}

/* dont use the keyword vector - it will clash with the
 * altivec extensions used for powerpc processors.
 */

static real *rvector(int nl, int nh)
{
    real *v;

    v = (real *)malloc((unsigned) (nh-nl+1)*sizeof(real));
    if (!v)
    {
        nrerror("allocation failure in rvector()", TRUE);
    }
    return v-nl;
}

static int *ivector(int nl, int nh)
{
    int *v;

    v = (int *)malloc((unsigned) (nh-nl+1)*sizeof(int));
    if (!v)
    {
        nrerror("allocation failure in ivector()", TRUE);
    }
    return v-nl;
}


static real **matrix1(int nrl, int nrh, int ncl, int nch)
{
    int    i;
    real **m;

    m = (real **) malloc((unsigned) (nrh-nrl+1)*sizeof(real*));
    if (!m)
    {
        nrerror("allocation failure 1 in matrix1()", TRUE);
    }
    m -= nrl;

    for (i = nrl; i <= nrh; i++)
    {
        m[i] = (real *) malloc((unsigned) (nch-ncl+1)*sizeof(real));
        if (!m[i])
        {
            nrerror("allocation failure 2 in matrix1()", TRUE);
        }
        m[i] -= ncl;
    }
    return m;
}

static double **dmatrix(int nrl, int nrh, int ncl, int nch)
{
    int      i;
    double **m;

    m = (double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
    if (!m)
    {
        nrerror("allocation failure 1 in dmatrix()", TRUE);
    }
    m -= nrl;

    for (i = nrl; i <= nrh; i++)
    {
        m[i] = (double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
        if (!m[i])
        {
            nrerror("allocation failure 2 in dmatrix()", TRUE);
        }
        m[i] -= ncl;
    }
    return m;
}

static int **imatrix1(int nrl, int nrh, int ncl, int nch)
{
    int i, **m;

    m = (int **)malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
    if (!m)
    {
        nrerror("allocation failure 1 in imatrix1()", TRUE);
    }
    m -= nrl;

    for (i = nrl; i <= nrh; i++)
    {
        m[i] = (int *)malloc((unsigned) (nch-ncl+1)*sizeof(int));
        if (!m[i])
        {
            nrerror("allocation failure 2 in imatrix1()", TRUE);
        }
        m[i] -= ncl;
    }
    return m;
}



static real **submatrix(real **a, int oldrl, int oldrh, int oldcl,
                        int newrl, int newcl)
{
    int    i, j;
    real **m;

    m = (real **) malloc((unsigned) (oldrh-oldrl+1)*sizeof(real*));
    if (!m)
    {
        nrerror("allocation failure in submatrix()", TRUE);
    }
    m -= newrl;

    for (i = oldrl, j = newrl; i <= oldrh; i++, j++)
    {
        m[j] = a[i]+oldcl-newcl;
    }

    return m;
}



static void free_vector(real *v, int nl)
{
    free((char*) (v+nl));
}

static void free_ivector(int *v, int nl)
{
    free((char*) (v+nl));
}

static void free_dvector(int *v, int nl)
{
    free((char*) (v+nl));
}



static void free_matrix(real **m, int nrl, int nrh, int ncl)
{
    int i;

    for (i = nrh; i >= nrl; i--)
    {
        free((char*) (m[i]+ncl));
    }
    free((char*) (m+nrl));
}

static real **convert_matrix(real *a, int nrl, int nrh, int ncl, int nch)
{
    int    i, j, nrow, ncol;
    real **m;

    nrow = nrh-nrl+1;
    ncol = nch-ncl+1;
    m    = (real **) malloc((unsigned) (nrow)*sizeof(real*));
    if (!m)
    {
        nrerror("allocation failure in convert_matrix()", TRUE);
    }
    m -= nrl;
    for (i = 0, j = nrl; i <= nrow-1; i++, j++)
    {
        m[j] = a+ncol*i-ncl;
    }
    return m;
}



static void free_convert_matrix(real **b, int nrl)
{
    free((char*) (b+nrl));
}

#define SWAP(a, b) {real temp = (a); (a) = (b); (b) = temp; }

static void dump_mat(int n, real **a)
{
    int i, j;

    for (i = 1; (i <= n); i++)
    {
        for (j = 1; (j <= n); j++)
        {
            fprintf(stderr, "  %10.3f", a[i][j]);
        }
        fprintf(stderr, "\n");
    }
}

gmx_bool gaussj(real **a, int n, real **b, int m)
{
    int *indxc, *indxr, *ipiv;
    int  i, icol = 0, irow = 0, j, k, l, ll;
    real big, dum, pivinv;

    indxc = ivector(1, n);
    indxr = ivector(1, n);
    ipiv  = ivector(1, n);
    for (j = 1; j <= n; j++)
    {
        ipiv[j] = 0;
    }
    for (i = 1; i <= n; i++)
    {
        big = 0.0;
        for (j = 1; j <= n; j++)
        {
            if (ipiv[j] != 1)
            {
                for (k = 1; k <= n; k++)
                {
                    if (ipiv[k] == 0)
                    {
                        if (fabs(a[j][k]) >= big)
                        {
                            big  = fabs(a[j][k]);
                            irow = j;
                            icol = k;
                        }
                    }
                    else if (ipiv[k] > 1)
                    {
                        nrerror("GAUSSJ: Singular Matrix-1", FALSE);
                        return FALSE;
                    }
                }
            }
        }
        ++(ipiv[icol]);
        if (irow != icol)
        {
            for (l = 1; l <= n; l++)
            {
                SWAP(a[irow][l], a[icol][l]);
            }
            for (l = 1; l <= m; l++)
            {
                SWAP(b[irow][l], b[icol][l]);
            }
        }
        indxr[i] = irow;
        indxc[i] = icol;
        if (a[icol][icol] == 0.0)
        {
            fprintf(stderr, "irow = %d, icol = %d\n", irow, icol);
            dump_mat(n, a);
            nrerror("GAUSSJ: Singular Matrix-2", FALSE);
            return FALSE;
        }
        pivinv        = 1.0/a[icol][icol];
        a[icol][icol] = 1.0;
        for (l = 1; l <= n; l++)
        {
            a[icol][l] *= pivinv;
        }
        for (l = 1; l <= m; l++)
        {
            b[icol][l] *= pivinv;
        }
        for (ll = 1; ll <= n; ll++)
        {
            if (ll != icol)
            {
                dum         = a[ll][icol];
                a[ll][icol] = 0.0;
                for (l = 1; l <= n; l++)
                {
                    a[ll][l] -= a[icol][l]*dum;
                }
                for (l = 1; l <= m; l++)
                {
                    b[ll][l] -= b[icol][l]*dum;
                }
            }
        }
    }
    for (l = n; l >= 1; l--)
    {
        if (indxr[l] != indxc[l])
        {
            for (k = 1; k <= n; k++)
            {
                SWAP(a[k][indxr[l]], a[k][indxc[l]]);
            }
        }
    }
    free_ivector(ipiv, 1);
    free_ivector(indxr, 1);
    free_ivector(indxc, 1);

    return TRUE;
}

#undef SWAP


static void covsrt(real **covar, int ma, int lista[], int mfit)
{
    int  i, j;
    real swap;

    for (j = 1; j < ma; j++)
    {
        for (i = j+1; i <= ma; i++)
        {
            covar[i][j] = 0.0;
        }
    }
    for (i = 1; i < mfit; i++)
    {
        for (j = i+1; j <= mfit; j++)
        {
            if (lista[j] > lista[i])
            {
                covar[lista[j]][lista[i]] = covar[i][j];
            }
            else
            {
                covar[lista[i]][lista[j]] = covar[i][j];
            }
        }
    }
    swap = covar[1][1];
    for (j = 1; j <= ma; j++)
    {
        covar[1][j] = covar[j][j];
        covar[j][j] = 0.0;
    }
    covar[lista[1]][lista[1]] = swap;
    for (j = 2; j <= mfit; j++)
    {
        covar[lista[j]][lista[j]] = covar[1][j];
    }
    for (j = 2; j <= ma; j++)
    {
        for (i = 1; i <= j-1; i++)
        {
            covar[i][j] = covar[j][i];
        }
    }
}

#define SWAP(a, b) {swap = (a); (a) = (b); (b) = swap; }

static void covsrt_new(real **covar, int ma, int ia[], int mfit)
/* Expand in storage the covariance matrix covar, so as to take
 * into account parameters that are being held fixed. (For the
 * latter, return zero covariances.)
 */
{
    int  i, j, k;
    real swap;
    for (i = mfit+1; i <= ma; i++)
    {
        for (j = 1; j <= i; j++)
        {
            covar[i][j] = covar[j][i] = 0.0;
        }
    }
    k = mfit;
    for (j = ma; j >= 1; j--)
    {
        if (ia[j])
        {
            for (i = 1; i <= ma; i++)
            {
                SWAP(covar[i][k], covar[i][j]);
            }
            for (i = 1; i <= ma; i++)
            {
                SWAP(covar[k][i], covar[j][i]);
            }
            k--;
        }
    }
}
#undef SWAP

static void mrqcof(real x[], real y[], real sig[], int ndata, real a[],
                   int ma, int lista[], int mfit,
                   real **alpha, real beta[], real *chisq,
                   void (*funcs)(real, real *, real *, real *))
{
    int  k, j, i;
    real ymod, wt, sig2i, dy, *dyda;

    dyda = rvector(1, ma);
    for (j = 1; j <= mfit; j++)
    {
        for (k = 1; k <= j; k++)
        {
            alpha[j][k] = 0.0;
        }
        beta[j] = 0.0;
    }
    *chisq = 0.0;
    for (i = 1; i <= ndata; i++)
    {
        (*funcs)(x[i], a, &ymod, dyda);
        sig2i = 1.0/(sig[i]*sig[i]);
        dy    = y[i]-ymod;
        for (j = 1; j <= mfit; j++)
        {
            wt = dyda[lista[j]]*sig2i;
            for (k = 1; k <= j; k++)
            {
                alpha[j][k] += wt*dyda[lista[k]];
            }
            beta[j] += dy*wt;
        }
        (*chisq) += dy*dy*sig2i;
    }
    for (j = 2; j <= mfit; j++)
    {
        for (k = 1; k <= j-1; k++)
        {
            alpha[k][j] = alpha[j][k];
        }
    }
    free_vector(dyda, 1);
}


gmx_bool mrqmin(real x[], real y[], real sig[], int ndata, real a[],
                int ma, int lista[], int mfit,
                real **covar, real **alpha, real *chisq,
                void (*funcs)(real, real *, real *, real *),
                real *alamda)
{
    int          k, kk, j, ihit;
    static real *da, *atry, **oneda, *beta, ochisq;

    if (*alamda < 0.0)
    {
        oneda = matrix1(1, mfit, 1, 1);
        atry  = rvector(1, ma);
        da    = rvector(1, ma);
        beta  = rvector(1, ma);
        kk    = mfit+1;
        for (j = 1; j <= ma; j++)
        {
            ihit = 0;
            for (k = 1; k <= mfit; k++)
            {
                if (lista[k] == j)
                {
                    ihit++;
                }
            }
            if (ihit == 0)
            {
                lista[kk++] = j;
            }
            else if (ihit > 1)
            {
                nrerror("Bad LISTA permutation in MRQMIN-1", FALSE);
                return FALSE;
            }
        }
        if (kk != ma+1)
        {
            nrerror("Bad LISTA permutation in MRQMIN-2", FALSE);
            return FALSE;
        }
        *alamda = 0.001;
        mrqcof(x, y, sig, ndata, a, ma, lista, mfit, alpha, beta, chisq, funcs);
        ochisq = (*chisq);
    }
    for (j = 1; j <= mfit; j++)
    {
        for (k = 1; k <= mfit; k++)
        {
            covar[j][k] = alpha[j][k];
        }
        covar[j][j] = alpha[j][j]*(1.0+(*alamda));
        oneda[j][1] = beta[j];
    }
    if (!gaussj(covar, mfit, oneda, 1))
    {
        return FALSE;
    }
    for (j = 1; j <= mfit; j++)
    {
        da[j] = oneda[j][1];
    }
    if (*alamda == 0.0)
    {
        covsrt(covar, ma, lista, mfit);
        free_vector(beta, 1);
        free_vector(da, 1);
        free_vector(atry, 1);
        free_matrix(oneda, 1, mfit, 1);
        return TRUE;
    }
    for (j = 1; j <= ma; j++)
    {
        atry[j] = a[j];
    }
    for (j = 1; j <= mfit; j++)
    {
        atry[lista[j]] = a[lista[j]]+da[j];
    }
    mrqcof(x, y, sig, ndata, atry, ma, lista, mfit, covar, da, chisq, funcs);
    if (*chisq < ochisq)
    {
        *alamda *= 0.1;
        ochisq   = (*chisq);
        for (j = 1; j <= mfit; j++)
        {
            for (k = 1; k <= mfit; k++)
            {
                alpha[j][k] = covar[j][k];
            }
            beta[j]     = da[j];
            a[lista[j]] = atry[lista[j]];
        }
    }
    else
    {
        *alamda *= 10.0;
        *chisq   = ochisq;
    }
    return TRUE;
}


gmx_bool mrqmin_new(real x[], real y[], real sig[], int ndata, real a[],
                    int ia[], int ma, real **covar, real **alpha, real *chisq,
                    void (*funcs)(real, real [], real *, real []),
                    real *alamda)
/* Levenberg-Marquardt method, attempting to reduce the value Chi^2
 * of a fit between a set of data points x[1..ndata], y[1..ndata]
 * with individual standard deviations sig[1..ndata], and a nonlinear
 * function dependent on ma coefficients a[1..ma]. The input array
 * ia[1..ma] indicates by nonzero entries those components of a that
 * should be fitted for, and by zero entries those components that
 * should be held fixed at their input values. The program returns
 * current best-fit values for the parameters a[1..ma], and
 * Chi^2 = chisq. The arrays covar[1..ma][1..ma], alpha[1..ma][1..ma]
 * are used as working space during most iterations. Supply a routine
 * funcs(x,a,yfit,dyda,ma) that evaluates the fitting function yfit,
 * and its derivatives dyda[1..ma] with respect to the fitting
 * parameters a at x. On the first call provide an initial guess for
 * the parameters a, and set alamda < 0 for initialization (which then
 * sets alamda=.001). If a step succeeds chisq becomes smaller and
 * alamda de-creases by a factor of 10. If a step fails alamda grows by
 * a factor of 10. You must call this routine repeatedly until
 * convergence is achieved. Then, make one final call with alamda=0,
 * so that covar[1..ma][1..ma] returns the covariance matrix, and alpha
 * the curvature matrix.
 * (Parameters held fixed will return zero covariances.)
 */
{
    void covsrt(real **covar, int ma, int ia[], int mfit);
    gmx_bool gaussj(real **a, int n, real **b, int m);
    void mrqcof_new(real x[], real y[], real sig[], int ndata, real a[],
                    int ia[], int ma, real **alpha, real beta[], real *chisq,
                    void (*funcs)(real, real [], real *, real []));
    int         j, k, l;
    static int  mfit;
    static real ochisq, *atry, *beta, *da, **oneda;

    if (*alamda < 0.0)                    /* Initialization. */
    {
        atry = rvector(1, ma);
        beta = rvector(1, ma);
        da   = rvector(1, ma);
        for (mfit = 0, j = 1; j <= ma; j++)
        {
            if (ia[j])
            {
                mfit++;
            }
        }
        oneda   = matrix1(1, mfit, 1, 1);
        *alamda = 0.001;
        mrqcof_new(x, y, sig, ndata, a, ia, ma, alpha, beta, chisq, funcs);
        ochisq = (*chisq);
        for (j = 1; j <= ma; j++)
        {
            atry[j] = a[j];
        }
    }
    for (j = 1; j <= mfit; j++) /* Alter linearized fitting matrix, by augmenting. */
    {
        for (k = 1; k <= mfit; k++)
        {
            covar[j][k] = alpha[j][k]; /* diagonal elements. */
        }
        covar[j][j] = alpha[j][j]*(1.0+(*alamda));
        oneda[j][1] = beta[j];
    }
    if (!gaussj(covar, mfit, oneda, 1)) /* Matrix solution. */
    {
        return FALSE;
    }
    for (j = 1; j <= mfit; j++)
    {
        da[j] = oneda[j][1];
    }
    if (*alamda == 0.0) /* Once converged, evaluate covariance matrix. */
    {
        covsrt_new(covar, ma, ia, mfit);
        free_matrix(oneda, 1, mfit, 1);
        free_vector(da, 1);
        free_vector(beta, 1);
        free_vector(atry, 1);
        return TRUE;
    }
    for (j = 0, l = 1; l <= ma; l++) /* Did the trial succeed? */
    {
        if (ia[l])
        {
            atry[l] = a[l]+da[++j];
        }
    }
    mrqcof_new(x, y, sig, ndata, atry, ia, ma, covar, da, chisq, funcs);
    if (*chisq < ochisq)
    {
        /* Success, accept the new solution. */
        *alamda *= 0.1;
        ochisq   = (*chisq);
        for (j = 1; j <= mfit; j++)
        {
            for (k = 1; k <= mfit; k++)
            {
                alpha[j][k] = covar[j][k];
            }
            beta[j] = da[j];
        }
        for (l = 1; l <= ma; l++)
        {
            a[l] = atry[l];
        }
    }
    else                 /* Failure, increase alamda and return. */
    {
        *alamda *= 10.0;
        *chisq   = ochisq;
    }
    return TRUE;
}

void mrqcof_new(real x[], real y[], real sig[], int ndata, real a[],
                int ia[], int ma, real **alpha, real beta[], real *chisq,
                void (*funcs)(real, real [], real *, real[]))
/* Used by mrqmin to evaluate the linearized fitting matrix alpha, and
 * vector beta as in (15.5.8), and calculate Chi^2.
 */
{
    int  i, j, k, l, m, mfit = 0;
    real ymod, wt, sig2i, dy, *dyda;

    dyda = rvector(1, ma);
    for (j = 1; j <= ma; j++)
    {
        if (ia[j])
        {
            mfit++;
        }
    }
    for (j = 1; j <= mfit; j++) /* Initialize (symmetric) alpha), beta. */
    {
        for (k = 1; k <= j; k++)
        {
            alpha[j][k] = 0.0;
        }
        beta[j] = 0.0;
    }
    *chisq = 0.0;
    for (i = 1; i <= ndata; i++) /* Summation loop over all data. */
    {
        (*funcs)(x[i], a, &ymod, dyda);
        sig2i = 1.0/(sig[i]*sig[i]);
        dy    = y[i]-ymod;
        for (j = 0, l = 1; l <= ma; l++)
        {
            if (ia[l])
            {
                wt = dyda[l]*sig2i;
                for (j++, k = 0, m = 1; m <= l; m++)
                {
                    if (ia[m])
                    {
                        alpha[j][++k] += wt*dyda[m];
                    }
                }
                beta[j] += dy*wt;
            }
        }
        *chisq += dy*dy*sig2i;  /* And find Chi^2. */
    }
    for (j = 2; j <= mfit; j++) /* Fill in the symmetric side. */
    {
        for (k = 1; k < j; k++)
        {
            alpha[k][j] = alpha[j][k];
        }
    }
    free_vector(dyda, 1);
}
