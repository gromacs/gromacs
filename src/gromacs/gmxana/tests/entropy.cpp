/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2017, by the GROMACS development team, led by
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
 * Tests for entropy code
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "gmxpre.h"

#include "gromacs/gmxana/entropy.h"

#include "config.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

//#include "testutils/cmdlinetest.h"
//#include "testutils/stdiohelper.h"
//#include "testutils/textblockmatchers.h"

namespace gmx
{

namespace
{

class Entropy : public ::testing::Test
{
    protected:
        test::TestReferenceData        refData_;
        test::TestReferenceChecker     checker_;
        Entropy( ) : checker_(refData_.rootChecker())
        {
        }
    public:
        void runTest(EntropyAlgorithm ea, real temperature, int nskip)
        {
            real   eigenvalue[] =
            {
                -1.94784e-07, 2.82555e-07, 4.22414e-07,
                398.335, 725.784, 763.472, 7000.48, 7572.49,
                8500.42, 8956.78, 12872.3, 15920.4, 22010.5,
                22695.3, 38346, 43599.4, 76653.8, 105882,
                116206, 119411, 345267, 349635, 370323, 372737
            };
            int    neigval = sizeof(eigenvalue)/sizeof(eigenvalue[0]);
            real   S       = calc_entropy(ea, neigval, eigenvalue,
                                          temperature, nskip);
            double tol   = 1e-3;
            checker_.setDefaultTolerance(test::relativeToleranceAsFloatingPoint(1, tol));
            checker_.checkReal(S, "entropy");

        }
};

TEST_F(Entropy, Schliter_300_Skip_0)
{
    runTest(eaSchlitter, 300.0, 0);
}

TEST_F(Entropy, Schliter_300_Skip_3)
{
    runTest(eaSchlitter, 300.0, 3);
}

TEST_F(Entropy, QuasiHarmonic_200_Skip_0)
{
    runTest(eaQuasiHarmonic, 200.0, 0);
}

TEST_F(Entropy, QuasiHarmonic_200_Skip_3)
{
    runTest(eaQuasiHarmonic, 300.0, 3);
}

} // namespace

} // namespace gmx
