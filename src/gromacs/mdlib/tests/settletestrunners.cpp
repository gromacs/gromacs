/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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

#include <gtest/gtest.h>

#include "gromacs/mdlib/settle.h"

#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{

void applySettle(SettleTestData*    testData,
                 const t_pbc        pbc,
                 const bool         updateVelocities,
                 const bool         calcVirial,
                 const std::string& testDescription)
{
    settledata* settled = settle_init(testData->mtop_);

    settle_set_constraints(settled, &testData->ilist_, testData->mdatoms_);

    bool errorOccured;
    int  numThreads  = 1;
    int  threadIndex = 0;
    csettle(settled, numThreads, threadIndex, &pbc,
            reinterpret_cast<real*>(as_rvec_array(testData->x_.data())),
            reinterpret_cast<real*>(as_rvec_array(testData->xPrime_.data())), testData->reciprocalTimeStep_,
            updateVelocities ? reinterpret_cast<real*>(as_rvec_array(testData->v_.data())) : nullptr,
            calcVirial, testData->virial_, &errorOccured);
    settle_free(settled);
    EXPECT_FALSE(errorOccured) << testDescription;
}

#if GMX_GPU != GMX_GPU_CUDA

void applySettleGpu(gmx_unused SettleTestData* testData,
                    gmx_unused const t_pbc pbc,
                    gmx_unused const bool  updateVelocities,
                    gmx_unused const bool  calcVirial,
                    gmx_unused const std::string& testDescription)
{
    FAIL() << "Dummy SETTLE GPU function was called instead of the real one in the SETTLE test.";
}

#endif
} // namespace test
} // namespace gmx
