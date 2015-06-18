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

#include "gromacs/legacyheaders/names.h"
#include "gromacs/legacyheaders/types/enums.h"
#include "gromacs/math/utilities.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/snprintf.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace gmx
{

namespace
{

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

        void testSplineTableV(int           ntable,
                              double        dx,
                              const double *table_in_v)
        {
            QuadraticSplineTable st(ntable, dx, table_in_v, NULL);
            checker_.checkInteger(st.size(), "tableSize");
            checker_.checkDouble(st.spacing(), "tableSpacing");
            double r = 0.735*st.size()*st.spacing();
            checker_.checkDouble(r, "distance");
            checker_.checkDouble(st.energy(r), "energy");
            checker_.checkDouble(st.force(r), "force");
            const std::vector<double> f = st.F();
            checker_.checkSequence(f.begin(), f.end(), "F");
            const std::vector<double> v = st.V();
            checker_.checkSequence(v.begin(), v.end(), "V");
            const std::vector<double> fvd0 = st.FDV0();
            checker_.checkSequence(fvd0.begin(), fvd0.end(), "FDV0");
        }
        void testAnalytically()
        {
            const int ntable = 8;
            double    dx     = 0.1, v[ntable];
            for (int i = 0; (i < ntable); i++)
            {
                v[i] = my_function(i*dx);
            }
            QuadraticSplineTable st(ntable, dx, v, NULL);
            double               rtab  = st.size()*st.spacing();
            double               rr[6] = { 0.123, 0.371, 0.478, 0.511, 0.739, 0.911 };
            double               toler = 1e-4;
            for (unsigned int i = 0; (i < 6); i++)
            {
                char   buf[256];
                double r     = rr[i] * rtab;
                double etab  = st.energy(r);
                double eref  = my_function(r);
                //EXPECT_DOUBLE_EQ_TOL(etab, eref, toler);
                double ftab  = st.force(r);
                double fref  = -2*r;
                //EXPECT_DOUBLE_EQ_TOL(ftab, fref, toler);
            }
        }
        void testSplineTableVF(int           ntable,
                               double        dx,
                               const double *table_in_v,
                               const double *table_in_f)
        {
            QuadraticSplineTable st(ntable, dx, table_in_v, table_in_f);
            checker_.checkInteger(st.size(), "tableSize");
            checker_.checkDouble(st.spacing(), "tableSpacing");
            double r = 0.735*st.size()*st.spacing();
            checker_.checkDouble(r, "distance");
            checker_.checkDouble(st.energy(r), "energy");
            checker_.checkDouble(st.force(r), "force");
            const std::vector<double> f = st.F();
            checker_.checkSequence(f.begin(), f.end(), "F");
            const std::vector<double> v = st.V();
            checker_.checkSequence(v.begin(), v.end(), "V");
            const std::vector<double> fdv0 = st.FDV0();
            checker_.checkSequence(fdv0.begin(), fdv0.end(), "FDV0");
        }
};

TEST_F(QuadraticSplineTableTest, VHarmonic)
{
    /* Harmonic potential */
    double v[8];
    double dx = 0.1;
    for (unsigned int i = 0; (i < 8); i++)
    {
        v[i] = pow(dx*i, 2);
    }
    testSplineTableV(8, dx, v);
}

TEST_F(QuadraticSplineTableTest, VFHarmonic)
{
    /* Harmonic potential */
    double v[8], f[8];
    double dx = 0.03;
    for (unsigned int i = 0; (i < 8); i++)
    {
        v[i] = pow(dx*i, 2);
        f[i] = -2*dx*i;
    }
    testSplineTableVF(8, dx, v, f);
}

TEST_F(QuadraticSplineTableTest, EwaldQ)
{
    /* Long range Coulomb potential */
    double v[8];
    double b  = 4.5;
    double dx = 0.001;
    v[0] = 2*b/sqrt(M_PI);
    for (int i = 1; (i < 8); i++)
    {
        double r = dx*i;
        v[i] = gmx_erfd(b*r)/r;
    }
    testSplineTableV(8, dx, v);
}

TEST_F(QuadraticSplineTableTest, EwaldLJ)
{
    /* Long range LJ potential */
    double v[8];
    double b  = 4.5;
    double dx = 0.001;
    v[0] = pow(b, 6.0)/6.0;
    for (int i = 1; (i < 8); i++)
    {
        double r      = dx*i;
        double br     = b*r;
        double br2    = br*br;
        double br4    = br2*br2;
        double r6     = pow(r, 6.0);
        double expbr2 = exp(-br2);
        v[i]          = (-gmx_expm1(-br2) - expbr2*(br2 + 0.5*br4))/r6;
        /* v[i]          = (1 - expbr2*(1 + br2 + 0.5*br4))/r6; */
    }
    testSplineTableV(8, dx, v);
}

} // namespace

} // namespace gmx
