/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2014- The GROMACS Authors
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
 * Implements low level test of manyautocorrelation routines
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_correlationfunctions
 */
#include "gmxpre.h"

#include "gromacs/correlationfunctions/manyautocorrelation.h"

#include <cmath>

#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/real.h"

#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{
namespace
{

class ManyAutocorrelationTest : public ::testing::Test
{
};

TEST_F(ManyAutocorrelationTest, Empty)
{
    std::vector<std::vector<real>> c;
    EXPECT_THROW_GMX(many_auto_correl(&c), gmx::InconsistentInputError);
}

#ifndef NDEBUG
TEST_F(ManyAutocorrelationTest, DifferentLength)
{
    std::vector<std::vector<real>> c;
    c.resize(3);
    c[0].resize(10);
    c[1].resize(10);
    c[2].resize(8);
    EXPECT_THROW_GMX(many_auto_correl(&c), gmx::InconsistentInputError);
}
#endif

TEST_F(ManyAutocorrelationTest, MultipleSeries)
{
    // Regression test: a thread that processes more than one series (the
    // usual case with GMX_OPENMP=OFF, and whenever the series count
    // exceeds the thread count) must not leave power-spectrum values from
    // the previous series in the zero-padding area of its work array.
    const int ndata   = 8;
    const int nseries = 16;
    const int nfft    = (3 * ndata / 2) + 1;

    std::vector<std::vector<real>> original(nseries);
    for (int i = 0; i < nseries; ++i)
    {
        original[i].resize(ndata);
        for (int j = 0; j < ndata; ++j)
        {
            original[i][j] = std::sin(0.7 * (i + 1) + 1.1 * j) + 0.3 * std::cos(2.3 * j);
        }
    }

    std::vector<std::vector<real>> data = original;
    ASSERT_NO_THROW_GMX(many_auto_correl(&data));

    for (int i = 0; i < nseries; ++i)
    {
        real scale = 0.0;
        for (int j = 0; j < ndata; ++j)
        {
            scale += original[i][j] * original[i][j];
        }
        for (int k = 0; k < ndata; ++k)
        {
            // Correlation of the sequence zero-padded to nfft entries.
            real expected = 0.0;
            for (int m = 0; m < ndata; ++m)
            {
                const int mshifted = (m + k) % nfft;
                if (mshifted < ndata)
                {
                    expected += original[i][m] * original[i][mshifted];
                }
            }
            EXPECT_REAL_EQ_TOL(expected,
                               data[i][k],
                               relativeToleranceAsPrecisionDependentFloatingPoint(scale, 1e-3, 1e-6))
                    << "series " << i << ", lag " << k;
        }
    }
}

} // namespace
} // namespace test
} // namespace gmx
