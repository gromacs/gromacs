/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2014-2018 
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour, 
 *             Paul J. van Maaren, 
 *             David van der Spoel (Project leader)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, 
 * Boston, MA  02110-1301, USA.
 */
 
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>

#include "gromacs/math/vec.h"
#include "gromacs/linearalgebra/matrix.h"
#include "gromacs/utility/exceptions.h"

#include "regression.h"

extern "C"
{
void dgelsd_(int* m, int* n, int* nrhs, double* a, int* lda,
             double* b, int* ldb, double* s, double* rcond, int* rank,
             double* work, int* lwork, int* iwork, int* info );
             
void dgels_(const char* trans, int* m, int* n, int* nrhs, double* a, int* lda,
            double* b, int* ldb, double* work, int* lwork, int* info );
            
void dgesvd_(const char* jobu, const char* jobvt, int* m, int* n, double* a,
             int* lda, double* s, double* u, int* ldu, double* vt, int* ldvt,
             double* work, int* lwork, int* info );
}

static void multi_regression2(int nrow, double y[], int ncol,
                              double **a, double x[])
{
    /* Query and allocate the optimal workspace */
    int                 lwork = -1;
    int                 lda   = nrow;
    int                 ldb   = nrow;
    int                 nrhs  = 1;
    int                 rank;
    double              rcond = -1.0;
    double              wkopt;
    std::vector<double> s;
    s.resize(nrow);
    // Compute length of integer array iwork according to
    // https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_lapack_examples/dgelsd_ex.c.htm
    int              smlsiz = 25;
    int              nlvl   = std::max(0L, std::lround(std::log2(std::min(nrow, ncol)/smlsiz + 1) ) + 1);
    int              liwork = 3*std::min(ncol, nrow)*nlvl + 11*std::min(nrow, ncol);
    std::vector<int> iwork;
    iwork.resize(liwork);
    int              info;
    bool             bDgelsd = true;
    if (bDgelsd)
    {
        dgelsd_ (&nrow, &ncol, &nrhs, a[0], &lda, y, &ldb, s.data(),
                 &rcond, &rank, &wkopt, &lwork,
                 iwork.data(), &info );
    }
    else
    {
        dgels_ ("No transpose", &nrow, &ncol, &nrhs, a[0], &lda, y, &ldb,
                &wkopt, &lwork, &info );
    }
    lwork = (int)wkopt;
    std::vector<double> work;
    work.resize(lwork);
    /* Solve the equations A*X = B */
    if (bDgelsd)
    {
        dgelsd_ (&nrow, &ncol, &nrhs, a[0], &lda, y, &ldb, s.data(),
                 &rcond, &rank, work.data(), &lwork,
                 iwork.data(), &info );
    }
    else
    {
        dgels_ ("No transpose", &nrow, &ncol, &nrhs, a[0], &lda, y, &ldb,
                work.data(), &lwork, &info );
    }
    /* Check for convergence */
    if (info > 0)
    {
        char buf[256];
        snprintf(buf, sizeof(buf), "The algorithm computing SVD failed to converge and the least squares solution could not be computed. Info = %d.", info);
        GMX_THROW(gmx::InvalidInputError(buf));
    }
    for (int i = 0; i < ncol; i++)
    {
        x[i] = y[i];
    }
}

static void tensor2matrix(tensor c, double **a)
{
    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            a[i][j] = c[i][j];
        }
    }

}

static void array2tensor(double **a, tensor c)
{
    tensor tmp;    
    tmp[XX][XX] = a[XX][ZZ];
    tmp[XX][YY] = a[XX][YY];
    tmp[XX][ZZ] = a[XX][XX];   
    tmp[YY][XX] = a[ZZ][ZZ];
    tmp[YY][YY] = a[ZZ][YY];
    tmp[YY][ZZ] = a[ZZ][XX];    
    tmp[ZZ][XX] = a[YY][ZZ];
    tmp[ZZ][YY] = a[YY][YY];
    tmp[ZZ][ZZ] = a[YY][XX];   
    transpose(tmp, c);
}

static void unitTensor(tensor I, int ncol, int nrow)
{
    for (int i = 0; i < ncol; i++)
    {
        for (int j = 0; j < nrow; j++)
        {
            I[i][j]  = 1 ? (ncol == nrow) : 0;
        }
    }
}

static void SVD(int*    m,    int*     n,     double** a,
                int*    lda,  double** s,     double** u, 
                int*    ldu,  double** vt,    int*     ldvt)
{    
    int       lwork;   
    int       info; 
    double    wkopt;    
    std::vector<double> work;
    
    lwork = -1;
    dgesvd_("All", "All", m, n, a[0], lda, s[0], u[0], ldu, vt[0], ldvt, &wkopt, &lwork, &info);
    lwork = (int)wkopt;
    work.resize(lwork);
    dgesvd_("All", "All", m, n, a[0], lda, s[0], u[0], ldu, vt[0], ldvt, work.data(), &lwork, &info);
    if (info > 0)
    {
        char buf[256];
        snprintf(buf, sizeof(buf), "SVD algorithm failed to converge. Info = %d.", info);
        GMX_THROW(gmx::InvalidInputError(buf));
    }   
}

/*
 * Kabsch algorithm to find the optimal rotation matrix to 
 * rotate matrix P (N, 3) onto matrix Q (N, 3). 
 * https://en.wikipedia.org/wiki/Kabsch_algorithm
 *
 * tensor here is 3 by 3 matrix! 
 */
void kabsch_rotation(tensor P, tensor Q, tensor rotated_P)
{
    int       m      = DIM;
    int       n      = DIM;
    int       lda    = m;
    int       ldu    = m;
    int       ldvt   = n;        
    double  **A      = alloc_matrix(m, n);;        
    double  **s      = alloc_matrix(m, n);;
    double  **u      = alloc_matrix(m, n);;
    double  **vt     = alloc_matrix(m, n);;   
    tensor    C, U, Vt, UVt;
    tensor    Pt, I, UI, R;
    
        
    /*transpose of p (Pt)*/
    transpose(P, Pt);
    
    /*Covariance matrix (A)*/
    mmul(Pt, Q, C);    
    tensor2matrix(C, A);
    
    /* SVD of the covariance matrix (A)
     * u : left singular vectors
     * vt: right singular vectors, trasnposed
     * s : singular values
     * A = u*s*vt
     */
    SVD(&m, &n, A, &lda, s, u, &ldu, vt, &ldvt);       
    array2tensor(u, U);
    array2tensor(vt, Vt);
    
    /*Check for correction*/
    mmul(U, Vt, UVt);
    int d;
    if (det(UVt) > 0)
    {
        d = 1;
    }
    else
    {
        d = -1;
    }
    
    /*Unit tensor (I)*/
    unitTensor(I, m, n);
    I[ZZ][ZZ] = d;
    
    /*Rotation matrix (R)*/
    mmul(U, I, UI);
    mmul(UI, Vt, R);
    
    /*Rotate*/
    mmul(R, P, rotated_P);
  
    free_matrix(A);
    free_matrix(s);
    free_matrix(u);
    free_matrix(vt);
}

MatrixWrapper::MatrixWrapper(int ncolumn, int nrow)
{
    ncolumn_ = ncolumn;
    nrow_    = nrow;
    a_       = alloc_matrix(ncolumn_, nrow_);
}

MatrixWrapper::~MatrixWrapper()
{
    free_matrix(a_);
}

void MatrixWrapper::set(int col, int row, double value)
{
    a_[col][row] = value;
}

void MatrixWrapper::setRow(int row, const double value[])
{
    for (int i = 0; i < ncolumn_; i++)
    {
        set(i, row, value[i]);
    }
}

double MatrixWrapper::get(int col, int row) const
{
    return a_[col][row];
}

void MatrixWrapper::solve(std::vector<double> rhs, std::vector<double> *solution)
{
    multi_regression2(nrow_, rhs.data(), ncolumn_, a_, solution->data());
}
