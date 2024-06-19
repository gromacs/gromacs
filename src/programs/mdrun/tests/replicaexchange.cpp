/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2013- The GROMACS Authors
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
 * Tests for the mdrun replica-exchange functionality
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "config.h"

#include <memory>
#include <regex>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/filestream.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testfilemanager.h"

#include "energycomparison.h"
#include "multisimtest.h"
#include "trajectorycomparison.h"

namespace gmx
{
namespace test
{

//! Convenience typedef
typedef MultiSimTest ReplicaExchangeEnsembleTest;

TEST_P(ReplicaExchangeEnsembleTest, ExitsNormally)
{
    mdrunCaller_->addOption("-replex", 1);
    runExitsNormallyTest();
}

/* Note, not all preprocessor implementations nest macro expansions
   the same way / at all, if we would try to duplicate less code. */

#if GMX_LIB_MPI
INSTANTIATE_TEST_SUITE_P(
        WithDifferentControlVariables,
        ReplicaExchangeEnsembleTest,
        ::testing::Combine(::testing::Values(NumRanksPerSimulation(1), NumRanksPerSimulation(2)),
                           ::testing::Values(IntegrationAlgorithm::MD),
                           ::testing::Values(TemperatureCoupling::VRescale),
                           ::testing::Values(PressureCoupling::No, PressureCoupling::Berendsen)));
#else
INSTANTIATE_TEST_SUITE_P(
        DISABLED_WithDifferentControlVariables,
        ReplicaExchangeEnsembleTest,
        ::testing::Combine(::testing::Values(NumRanksPerSimulation(1), NumRanksPerSimulation(2)),
                           ::testing::Values(IntegrationAlgorithm::MD),
                           ::testing::Values(TemperatureCoupling::VRescale),
                           ::testing::Values(PressureCoupling::No, PressureCoupling::Berendsen)));
#endif

//! Convenience typedef
typedef MultiSimTest ReplicaExchangeTerminationTest;

TEST_P(ReplicaExchangeTerminationTest, WritesCheckpointAfterMaxhTerminationAndThenRestarts)
{
    mdrunCaller_->addOption("-replex", 1);
    runMaxhTest();
}

INSTANTIATE_TEST_SUITE_P(InNvt,
                         ReplicaExchangeTerminationTest,
                         ::testing::Combine(::testing::Values(NumRanksPerSimulation(1)),
                                            ::testing::Values(IntegrationAlgorithm::MD),
                                            ::testing::Values(TemperatureCoupling::VRescale),
                                            ::testing::Values(PressureCoupling::No)));

} // namespace test
} // namespace gmx
