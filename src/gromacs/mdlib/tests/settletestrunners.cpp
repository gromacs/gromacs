/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
 * \brief Defines the test runner for CPU version of SETTLE.
 *
 * Also adds stub for the GPU version to keep the compiler happy.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "settletestrunners.h"

#include "config.h"

#include <array>
#include <memory>
#include <utility>

#include <gtest/gtest.h>

#include "gromacs/math/arrayrefwithpadding.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/settle.h"
#include "gromacs/mdlib/tests/settletestdata.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/arrayref.h"

#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{

void SettleHostTestRunner::applySettle(SettleTestData*    testData,
                                       const t_pbc        pbc,
                                       const bool         updateVelocities,
                                       const bool         calcVirial,
                                       const std::string& testDescription)
{
    SettleData settled(testData->mtop_);

    settled.setConstraints(
            testData->idef_->il[F_SETTLE], testData->numAtoms_, testData->masses_, testData->inverseMasses_);

    bool errorOccured;
    int  numThreads  = 1;
    int  threadIndex = 0;
    csettle(settled,
            numThreads,
            threadIndex,
            &pbc,
            testData->x_.arrayRefWithPadding(),
            testData->xPrime_.arrayRefWithPadding(),
            testData->reciprocalTimeStep_,
            updateVelocities ? testData->v_.arrayRefWithPadding() : ArrayRefWithPadding<RVec>(),
            calcVirial,
            testData->virial_,
            &errorOccured);
    EXPECT_FALSE(errorOccured) << testDescription;
}

} // namespace test
} // namespace gmx
