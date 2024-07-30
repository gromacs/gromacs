/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2017- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Tests for entropy code
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "gmxpre.h"

#include <algorithm>
#include <iterator>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/gmxana/thermochemistry.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{

class Entropy : public ::testing::Test
{
protected:
    test::TestReferenceData    refData_;
    test::TestReferenceChecker checker_;
    Entropy() : checker_(refData_.rootChecker()) {}

public:
    void runSchlitter(real temperature, gmx_bool bLinear)
    {
        std::vector<double> ev = { 0.00197861,  0.000389439, 0.000316043,  0.000150392, 0.000110254,
                                   8.99659e-05, 8.06572e-05, 5.14339e-05,  4.34268e-05, 2.16063e-05,
                                   1.65182e-05, 1.3965e-05,  1.01937e-05,  5.83113e-06, 4.59067e-06,
                                   3.39577e-06, 2.14837e-06, 1.20468e-06,  9.60817e-18, 1.48941e-19,
                                   1.13939e-19, 5.02332e-20, -6.90708e-21, -1.91354e-18 };
        std::vector<real>   eigenvalue;
        std::copy(ev.begin(), ev.end(), std::back_inserter(eigenvalue));

        real S = calcSchlitterEntropy(eigenvalue, temperature, bLinear);
        checker_.setDefaultTolerance(test::relativeToleranceAsFloatingPoint(1, 1e-7));
        checker_.checkReal(S, "entropy");
    }

    void runQuasiHarmonic(real temperature, gmx_bool bLinear)
    {
        std::vector<double> ev = { -31.403, -7.73169, -3.80315, -2.15659e-06, -1.70991e-07, 236,
                                   4609.83, 6718.07,  6966.27,  8587.85,      10736.3,      13543.7,
                                   17721.3, 22868,    35517.8,  44118.1,      75827.9,      106277,
                                   115132,  120782,   445118,   451149,       481058,       484576 };
        std::vector<real>   eigenvalue;
        std::copy(ev.begin(), ev.end(), std::back_inserter(eigenvalue));

        real S = calcQuasiHarmonicEntropy(eigenvalue, temperature, bLinear, 1.0);

        checker_.setDefaultTolerance(test::relativeToleranceAsFloatingPoint(1, 1e-7));
        checker_.checkReal(S, "entropy");
    }
};

TEST_F(Entropy, Schlitter_300_NoLinear)
{
    runSchlitter(300.0, FALSE);
}

TEST_F(Entropy, Schlitter_300_Linear)
{
    runSchlitter(300.0, TRUE);
}

TEST_F(Entropy, QuasiHarmonic_300_NoLinear)
{
    runQuasiHarmonic(300.0, FALSE);
}

TEST_F(Entropy, QuasiHarmonic_200_NoLinear)
{
    runQuasiHarmonic(200.0, FALSE);
}

TEST_F(Entropy, QuasiHarmonic_200_Linear)
{
    runQuasiHarmonic(200.0, TRUE);
}

} // namespace
} // namespace test
} // namespace gmx
