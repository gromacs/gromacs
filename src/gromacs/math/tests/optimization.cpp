/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
 *
 * \brief Tests routines in neldermead.h .
 *
 * \author Christian Blau <blau@kth.se>
 */

#include "gmxpre.h"

#include "gromacs/math/optimization.h"

#include <cmath>

#include <functional>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/functions.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"

#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{

real mcCormick(ArrayRef<const real> x)
{
    return std::sin(x[0] + x[1]) + square(x[0] - x[1]) - 1.5_real * x[0] + 2.5_real * x[1] + 1._real;
}

struct RosenBrock3d
{
    real operator()(ArrayRef<const real> x)
    {
        return 100 * square(x[1] - square(x[0])) + square(1 - x[0])
               + 100 * square(x[2] - square(x[1])) + square(1 - x[1]);
    }
};

TEST(NelderMead, Optimizes2DFunctionCorrectly)
{
    std::vector<real> initalPoint = { 1, 1 };
    auto              result      = nelderMead(mcCormick, initalPoint);
    EXPECT_REAL_EQ_TOL(result.functionValue_, -1.91329, relativeToleranceAsFloatingPoint(1, 5e-5));

    initalPoint = { 0, 0 };
    result      = nelderMead(mcCormick, initalPoint);
    EXPECT_REAL_EQ_TOL(result.functionValue_, -1.91329, relativeToleranceAsFloatingPoint(1, 5e-5));
}

TEST(NelderMead, Optimizes3DFunctorCorrectly)
{
    std::vector<real> initalPoint = { 0, 0, 0 };
    auto              result      = nelderMead(RosenBrock3d(), initalPoint);
    EXPECT_REAL_EQ_TOL(result.coordinates_[0], 1.00, relativeToleranceAsFloatingPoint(1, 1e-6));
    EXPECT_REAL_EQ_TOL(result.coordinates_[1], 1.00, relativeToleranceAsFloatingPoint(1, 2e-6));
    EXPECT_REAL_EQ_TOL(result.coordinates_[2], 1.00, relativeToleranceAsFloatingPoint(1, 5e-6));
    EXPECT_REAL_EQ_TOL(result.functionValue_, 0, relativeToleranceAsFloatingPoint(1, 1e-7));
}

} // namespace
} // namespace test
} // namespace gmx
