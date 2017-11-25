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
#include "gromacs/linearalgebra/matrix.h"

#include "programs/alexandria/regression.h"

#include <cmath>
#include <cstdlib>

#include <gtest/gtest.h>

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

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

        void testRegression(int ncol, int nrow, double **a, double b[])
        {
            // Code and data from
            std::vector<double> x;
            x.resize(ncol);
            
            multi_regression2(nrow, b, ncol, a, x.data());

            checker_.checkSequence(x.begin(), x.end(), "solution");
        }
};

TEST_F (RegressionTest, Solve_A_x_is_B_1)
{
#define NCOL 4
#define NROW 5
    double **a = alloc_matrix(NROW, NCOL);
    a[0][0] =  0.12; a[0][1] = -6.91; a[0][2] = -3.33; a[0][3] =  3.97;
    a[1][0] = -8.19; a[1][1] =  2.22; a[1][2] = -8.94; a[1][3] =  3.33;
    a[2][0] =  7.69; a[2][1] = -5.12; a[2][2] = -6.72; a[2][3] = -2.74;
    a[3][0] = -2.26; a[3][1] = -9.08; a[3][2] = -4.40; a[3][4] = -7.92;
    a[4][0] = -4.71; a[4][1] =  9.96; a[4][2] = -9.98; a[4][3] = -3.20;
    
    double b[NROW] = {
        7.30,  1.33,  2.68, -9.62,  0.00,
    };

    testRegression(NCOL, NROW, a, b);
    free_matrix(a);
#undef NCOL
#undef NROW
}

TEST_F (RegressionTest, Solve_A_x_is_B_2)
{
#define NCOL 3
#define NROW 3
    double **a = alloc_matrix(NROW, NCOL);
    a[0][0] = 3.0;
    a[1][1] = 2.0;
    a[2][2] = 1.0;
    double b[NROW] = {
        3.0, 4.0, 5.0,
    };
    // Diagonal matrix, should give as answer ( 1, 2, 5 )
    testRegression(NCOL, NROW, a, b);
    free_matrix(a);
#undef NCOL
#undef NROW
}

TEST_F (RegressionTest, Solve_A_x_is_B_3)
{
#define NCOL 2
#define NROW 4
    double **a = alloc_matrix(NROW, NCOL);
    a[0][0] = 3.0;
    a[1][1] = 2.0; 
    a[2][0] = 1.0;
    a[3][0] = 1.0;
    a[3][1] = 2.0;
    double b[NROW] = {
        3.0, 4.0, 1.0, 5.0
    };
    // Answer should be ( 1, 2 )
    testRegression(NCOL, NROW, a, b);
    free_matrix(a);
#undef NCOL
#undef NROW
}

TEST_F (RegressionTest, Solve_A_x_is_B_4)
{
#define NCOL 1
#define NROW 2
    double **a = alloc_matrix(NROW, NCOL);
    a[0][0] = 3.0; a[0][1] = 2.0;
    double b[NROW] = {
        6.0, 4.0 
    };
    // Answer should be ( 2 )
    testRegression(NCOL, NROW, a, b);
    free_matrix(a);
#undef NCOL
#undef NROW
}

TEST_F (RegressionTest, Solve_A_x_is_B_5)
{
#define NCOL 2
#define NROW 2
    double **a = alloc_matrix(NROW, NCOL);
    a[0][0] = 3.0; a[1][1] = 2.0;
    double b[NROW] = {
        5.0, -5.0
    };
    // Answer should be ( 1, -1)
    testRegression(NCOL, NROW, a, b);
    free_matrix(a);
#undef NCOL
#undef NROW
}

