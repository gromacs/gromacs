/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
 * Tests various corners of splinetable implementations.
 *
 * The main point of these tests is to check that different functions
 * using table routines compile and work as intended.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_tables
 */
#include "gmxpre.h"

#include "gromacs/tables/splinetable.h"

#include <cmath>

#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/utilities.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/snprintf.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace gmx
{

namespace
{
//! Function returning the square of the argument.
static double my_function(double x)
{
    return x*x;
}

class QuadraticSplineTableTest : public ::testing::Test
{
    protected:
        test::TestReferenceData        refData_;
        test::TestReferenceChecker     checker_;
        QuadraticSplineTableTest( )
            : checker_(refData_.rootChecker())
        {
            checker_.setDefaultTolerance(test::relativeToleranceAsFloatingPoint(1, 1e-3));
        }

        void testSplineTableV(double                     dx,
                              const std::vector<double> &table_in_v)
        {
            std::vector<double>  table_in_f;
            QuadraticSplineTable st(dx, table_in_v, table_in_f);
            checker_.checkInteger(st.size(), "tableSize");
            checker_.checkDouble(st.spacing(), "tableSpacing");
            const std::vector<real> f = st.F();
            checker_.checkSequence(f.begin(), f.end(), "F");
            const std::vector<real> v = st.V();
            checker_.checkSequence(v.begin(), v.end(), "V");
            const std::vector<real> fvd0 = st.FDV0();
            checker_.checkSequence(fvd0.begin(), fvd0.end(), "FDV0");
        }
        void testAnalytically()
        {
            const int           ntable = 8;
            double              dx     = 0.1;
            std::vector<double> v, f;
            for (int i = 0; (i < ntable); i++)
            {
                v.push_back(my_function(i*dx));
            }
            QuadraticSplineTable         st(dx, v, f);
            double                       rtab  = st.size()*st.spacing();
            double                       rr[6] = { 0.123, 0.371, 0.478, 0.511, 0.739, 0.911 };
            test::FloatingPointTolerance toler =
                test::relativeToleranceAsFloatingPoint(1, 1e-3);
            for (unsigned int i = 0; (i < 6); i++)
            {
                double r     = rr[i] * rtab;
                double etab  = st.energy(r);
                double eref  = my_function(r);
                EXPECT_DOUBLE_EQ_TOL(etab, eref, toler);
                double ftab  = st.force(r);
                double fref  = -2*r;
                EXPECT_DOUBLE_EQ_TOL(ftab, fref, toler);
            }
        }
        void testSplineTableVF(double                     dx,
                               const std::vector<double> &table_in_v,
                               const std::vector<double> &table_in_f)
        {
            QuadraticSplineTable st(dx, table_in_v, table_in_f);
            checker_.checkInteger(st.size(), "tableSize");
            checker_.checkDouble(st.spacing(), "tableSpacing");
            const std::vector<real> f = st.F();
            checker_.checkSequence(f.begin(), f.end(), "F");
            const std::vector<real> v = st.V();
            checker_.checkSequence(v.begin(), v.end(), "V");
            const std::vector<real> fdv0 = st.FDV0();
            checker_.checkSequence(fdv0.begin(), fdv0.end(), "FDV0");
        }
};

TEST_F(QuadraticSplineTableTest, VHarmonic)
{
    /* Harmonic potential */
    std::vector<double> v;
    double              dx = 0.1;
    for (unsigned int i = 0; (i < 8); i++)
    {
        v.push_back(pow(dx*i, 2));
    }
    testSplineTableV(dx, v);
    testAnalytically();
}

TEST_F(QuadraticSplineTableTest, VFHarmonic)
{
    /* Harmonic potential */
    std::vector<double> v, f;
    double              dx = 0.03;
    for (unsigned int i = 0; (i < 8); i++)
    {
        v.push_back(pow(dx*i, 2));
        f.push_back(-2*dx*i);
    }
    testSplineTableVF(dx, v, f);
}

TEST_F(QuadraticSplineTableTest, EwaldQ)
{
    /* Long range Coulomb potential */
    std::vector<double> v;
    double              b  = 4.5;
    double              dx = 0.001;
    v.push_back(2*b/sqrt(M_PI));
    for (int i = 1; (i < 8); i++)
    {
        double r = dx*i;
        v.push_back(std::erf(b*r)/r);
    }
    testSplineTableV(dx, v);
}

TEST_F(QuadraticSplineTableTest, EwaldLJ)
{
    /* Long range LJ potential */
    std::vector<double> v;
    double              b  = 4.5;
    double              dx = 0.001;
    v.push_back(pow(b, 6.0)/6.0);
    for (int i = 1; (i < 8); i++)
    {
        double r      = dx*i;
        double br     = b*r;
        double br2    = br*br;
        double br4    = br2*br2;
        double r6     = pow(r, 6.0);
        double expbr2 = exp(-br2);
        v.push_back((-std::expm1(-br2) - expbr2*(br2 + 0.5*br4))/r6);
    }
    testSplineTableV(dx, v);
}

} // namespace

} // namespace gmx
