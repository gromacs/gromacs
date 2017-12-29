/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2017,2018, by the GROMACS development team, led by
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

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

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
            double            ev[] = {
                -1.94784025325e-07, 2.8255525230e-07, 4.2241402525e-07,
                398.33505325, 725.7835240, 763.472250, 7000.4825320, 7572.4925230,
                8500.42024567, 8956.7832450, 12872.305233, 15920.23540, 22010.50235,
                22695.303579, 38346.53520, 43599.402352, 76653.823520, 105882.52230,
                116206.03578, 119411.25470, 345267.43220, 349635.02352, 370323.32550, 372737.0253
            };
            int               neigval = sizeof(ev)/sizeof(ev[0]);
            std::vector<real> eigenvalue;
            for (int i = 0; i < neigval; i++)
            {
                eigenvalue.push_back(ev[i]);
            }
            real   S       = calc_entropy(ea, neigval, eigenvalue.data(),
                                          temperature, nskip);
            checker_.setDefaultTolerance(test::relativeToleranceAsFloatingPoint(1, 1e-7));
            checker_.checkReal(S, "entropy");
        }
};

TEST_F(Entropy, Schlitter_300_Skip_0)
{
    runTest(eaSchlitter, 300.0, 0);
}

TEST_F(Entropy, Schlitter_300_Skip_3)
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
