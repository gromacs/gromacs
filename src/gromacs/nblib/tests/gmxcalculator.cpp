/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
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
 * This implements basic nblib utility tests
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#include "gmxpre.h"

#include <gtest/gtest.h>

#include "gromacs/nblib/gmxsetup.h"
#include "gromacs/nblib/simulationstate.h"

#include "testhelpers.h"
#include "testsystems.h"

namespace nblib
{
namespace test
{
namespace
{

TEST(NBlibTest, CanConstructGmxForceCalculator)
{
    ArgonSimulationStateBuilder      argonSystemBuilder;
    SimulationState                  simState = argonSystemBuilder.setupSimulationState();
    std::shared_ptr<NBKernelOptions> options = std::make_shared<NBKernelOptions>(NBKernelOptions());
    EXPECT_NO_THROW(GmxForceCalculator(simState, options));
}

TEST(NBlibTest, GmxForceCalculatorCanCompute)
{
    ArgonSimulationStateBuilder argonSystemBuilder;
    SimulationState             simState       = argonSystemBuilder.setupSimulationState();
    NBKernelOptions             options        = NBKernelOptions();
    options.nbnxmSimd                          = BenchMarkKernels::SimdNo;
    std::unique_ptr<NbvSetupUtil> nbvSetupUtil = std::make_unique<NbvSetupUtil>(simState, options);
    std::unique_ptr<GmxForceCalculator> gmxForceCalculator = nbvSetupUtil->setupGmxForceCalculator();
    EXPECT_NO_THROW(gmxForceCalculator->compute());
}

TEST(NBlibTest, CanSetupStepWorkload)
{
    std::shared_ptr<NBKernelOptions> options = std::make_shared<NBKernelOptions>(NBKernelOptions());
    EXPECT_NO_THROW(setupStepWorkload(options));
}

TEST(NBlibTest, GmxForceCalculatorCanSetupInteractionConst)
{
    std::shared_ptr<NBKernelOptions> options = std::make_shared<NBKernelOptions>(NBKernelOptions());
    EXPECT_NO_THROW(setupInteractionConst(options));
}
} // namespace
} // namespace test
} // namespace nblib
