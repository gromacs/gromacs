/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Implements test of regression code
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_alexandria
 */


#include <cmath>
#include <cstdlib>

#include <gtest/gtest.h>

#include "programs/alexandria/regression.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

extern "C" void dgelsd_(int* m, int* n, int* nrhs, double* a, int* lda,
                        double* b, int* ldb, double* s, double* rcond, int* rank,
                        double* work, int* lwork, int* iwork, int* info );


class RegressionTest : public gmx::test::CommandLineTestBase
{
    protected:
        gmx::test::TestReferenceChecker checker_;

        //init set tolecrance
        RegressionTest () : checker_(this->rootChecker())
        {
            auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 5e-2);
            checker_.setDefaultTolerance(tolerance);
        }

        void testRegression(int m, int n, double a[], double b[])
        {
            // Code and data from
            // https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_lapack_examples/dgelsd_ex.c.htm
            int    lwork = -1;
            int    lda   = m;
            int    ldb   = n;
            int    nrhs  = 1;
            int    rank  = 0;
            int    info;
            double rcond = -1.0;
            double wkopt = 0.0;
            std::vector<double> s;
            s.resize(m);
            int  smlsiz = 25;
            int  nlvl   = std::max(0L, std::lround(std::log2(std::min(n, m)/smlsiz + 1) ) + 1);
            int  liwork = 3*std::min(m, n)*nlvl + 11*std::min(m, n);
            std::vector<int> iwork;
            iwork.resize(liwork);
 
            dgelsd_( &m, &n, &nrhs, a, &lda, b, &ldb, s.data(), 
                     &rcond, &rank, &wkopt, &lwork,
                     iwork.data(), &info );
            lwork = (int)wkopt;
            std::vector<double> work;
            work.resize(lwork);
            /* Solve the equations A*X = B */
            dgelsd_( &m, &n, &nrhs, a, &lda, b, &ldb, s.data(),
                     &rcond, &rank, work.data(), &lwork,
                     iwork.data(), &info );
            /* Check for convergence */
            if ( info > 0 ) {
                printf( "The algorithm computing SVD failed to converge;\n" );
                printf( "the least squares solution could not be computed.\n" );
                exit( 1 );
            }
            std::vector<double> x;
            for(int i = 0; i < n; i++)
            {
                x.push_back(b[i]);
            }
            checker_.checkSequence(x.begin(), x.end(), "solution");
        }
};

TEST_F (RegressionTest, Solve_A_x_is_B_1)
{
#define M 4
#define N 5
#define NRHS 1
#define LDA M
#define LDB N
    double a[LDA*N] = {
        0.12, -6.91, -3.33,  3.97,
        -8.19,  2.22, -8.94,  3.33,
        7.69, -5.12, -6.72, -2.74,
        -2.26, -9.08, -4.40, -7.92,
        -4.71,  9.96, -9.98, -3.20
    };
    double b[LDB*NRHS] = {
        7.30,  1.33,  2.68, -9.62,  0.00,
    };

    testRegression(M, N, a, b);
#undef M
#undef N
}

TEST_F (RegressionTest, Solve_A_x_is_B_2)
{
#define M 3
#define N 3
    double a[M*N] = {
        3.0,  0.0, 0.0,
        0.0,  2.0, 0.0,
        0.0,  0.0, 1.0 
    };
    double b[N] = {
        3.0, 4.0, 5.0,
    };

    testRegression(M, N, a, b);
#undef M
#undef N
}

