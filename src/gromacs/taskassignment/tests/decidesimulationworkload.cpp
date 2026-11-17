/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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
 * Tests for setupStepWorkload(), in particular where virtual-site forces are spread
 * and the force copies between host and GPU that this requires.
 *
 * setupStepWorkload() also checks the consistency of the flags it sets with release
 * assertions, so every configuration tested here also exercises those checks.
 *
 * \ingroup module_taskassignment
 */
#include "gmxpre.h"

#include "gromacs/taskassignment/decidesimulationworkload.h"

#include <cstdint>

#include <optional>
#include <string>
#include <tuple>

#include <gtest/gtest.h>

#include "gromacs/mdlib/force_flags.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/multipletimestepping.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{
namespace test
{
namespace
{

//! The kinds of MD step that affect where virtual-site forces are spread
enum class StepKind : int
{
    Plain,
    Virial,
    ForceOutput,
    VirialAndForceOutput,
    NeighborSearch,
    Count
};

//! Names of the step kinds, for test names
const char* enumValueToString(StepKind enumValue)
{
    const EnumerationArray<StepKind, const char*> stepKindName = {
        "Plain", "Virial", "ForceOutput", "VirialAndForceOutput", "NeighborSearch"
    };
    return stepKindName[enumValue];
}

//! How PME-PP communication is done
enum class PmeSetup : int
{
    //! No separate PME rank
    NoSeparatePmeRank,
    //! GPU PME on a separate rank, with PME-PP communication staged through the host
    StagedGpuPmeOnSeparateRank,
    Count
};

//! Names of the PME setups, for test names
const char* pmeSetupName(PmeSetup pmeSetup)
{
    switch (pmeSetup)
    {
        case PmeSetup::NoSeparatePmeRank: return "NoSeparatePme";
        case PmeSetup::StagedGpuPmeOnSeparateRank: return "StagedGpuPme";
        default: return "Unknown";
    }
}

//! The output interval for forces
constexpr int c_nstfout = 10;

//! Parameters: have vsites, use GPU update, allow GPU buffer ops, PME setup, step kind
using StepWorkloadTestParameters = std::tuple<bool, bool, bool, PmeSetup, StepKind>;

//! Test fixture for setupStepWorkload()
class StepWorkloadTest : public ::testing::TestWithParam<StepWorkloadTestParameters>
{
public:
    StepWorkloadTest()
    {
        std::tie(haveVirtualSites_, useGpuUpdate_, allowGpuBufferOps_, pmeSetup_, stepKind_) = GetParam();

        simulationWork_.computeNonbonded            = true;
        simulationWork_.useGpuNonbonded             = true;
        simulationWork_.useGpuXBufferOpsWhenAllowed = allowGpuBufferOps_;
        simulationWork_.useGpuFBufferOpsWhenAllowed = allowGpuBufferOps_;
        simulationWork_.useGpuUpdate                = useGpuUpdate_;
        simulationWork_.haveVirtualSites            = haveVirtualSites_;
        if (pmeSetup_ == PmeSetup::StagedGpuPmeOnSeparateRank)
        {
            simulationWork_.useGpuPme                = true;
            simulationWork_.haveSeparatePmeRank      = true;
            simulationWork_.useCpuPmePpCommunication = true;
        }

        inputrec_.outputControl.nstfout = c_nstfout;
    }

    //! Returns the step workload for the step kind of this test
    StepWorkload setupStepWorkloadForStepKind(const bool computeForces) const
    {
        int legacyFlags = GMX_FORCE_STATECHANGED | GMX_FORCE_LISTED | GMX_FORCE_NONBONDED;
        if (computeForces)
        {
            legacyFlags |= GMX_FORCE_FORCES;
        }
        if (stepKind_ == StepKind::Virial || stepKind_ == StepKind::VirialAndForceOutput)
        {
            legacyFlags |= GMX_FORCE_VIRIAL | GMX_FORCE_ENERGY;
        }
        if (stepKind_ == StepKind::NeighborSearch)
        {
            legacyFlags |= GMX_FORCE_NS;
        }
        const bool isForceOutputStep =
                (stepKind_ == StepKind::ForceOutput || stepKind_ == StepKind::VirialAndForceOutput);
        const int64_t  step            = isForceOutputStep ? 2 * c_nstfout : c_nstfout + 3;
        constexpr bool writeCheckpoint = false;
        return setupStepWorkload(legacyFlags,
                                 ArrayRef<const MtsLevel>{},
                                 step,
                                 std::optional<bool>(writeCheckpoint),
                                 domainWork_,
                                 simulationWork_,
                                 inputrec_);
    }

protected:
    //! Whether the system has virtual sites
    bool haveVirtualSites_ = false;
    //! Whether update runs on the GPU
    bool useGpuUpdate_ = false;
    //! Whether GPU buffer ops may be used
    bool allowGpuBufferOps_ = false;
    //! How PME-PP communication is done
    PmeSetup pmeSetup_ = PmeSetup::NoSeparatePmeRank;
    //! The kind of step
    StepKind stepKind_ = StepKind::Plain;
    //! The simulation workload
    SimulationWorkload simulationWork_;
    //! The domain workload, without CPU force work
    DomainLifetimeWorkload domainWork_;
    //! The input record, only used for the output intervals
    t_inputrec inputrec_;
};

TEST_P(StepWorkloadTest, VsiteForcesAreSpreadOnceInTheRightPlace)
{
    const StepWorkload stepWork = setupStepWorkloadForStepKind(true);

    // With GPU update, vsite forces should be spread on the GPU, except on virial
    // and force-output steps, where the spread forces are needed on the host.
    const bool isVirialStep =
            (stepKind_ == StepKind::Virial || stepKind_ == StepKind::VirialAndForceOutput);
    const bool isForceOutputStep =
            (stepKind_ == StepKind::ForceOutput || stepKind_ == StepKind::VirialAndForceOutput);
    const bool expectSpreadOnGpu = haveVirtualSites_ && useGpuUpdate_ && !isVirialStep && !isForceOutputStep;
    EXPECT_EQ(stepWork.spreadVsiteForcesOnGpu, expectSpreadOnGpu);
    EXPECT_EQ(stepWork.spreadVsiteForcesOnHost, haveVirtualSites_ && !expectSpreadOnGpu);

    if (useGpuUpdate_ && stepWork.spreadVsiteForcesOnHost)
    {
        EXPECT_TRUE(stepWork.copyFinalFToGpu)
                << "Forces spread on the host must be uploaded for the GPU update";
        if (stepWork.useGpuFBufferOps)
        {
            EXPECT_TRUE(stepWork.copyReducedFFromGpu)
                    << "Forces reduced on the GPU must be downloaded for spreading on the host";
        }
    }
    if (stepWork.spreadVsiteForcesOnGpu && stepWork.useGpuFBufferOps && pmeSetup_ == PmeSetup::NoSeparatePmeRank)
    {
        EXPECT_FALSE(stepWork.copyReducedFFromGpu)
                << "Forces reduced on the GPU should stay there for spreading on the GPU";
        EXPECT_FALSE(stepWork.copyFinalFToGpu);
    }
}

TEST_P(StepWorkloadTest, VsiteForcesAreNotSpreadWithoutForces)
{
    const StepWorkload stepWork = setupStepWorkloadForStepKind(false);

    EXPECT_FALSE(stepWork.spreadVsiteForcesOnHost);
    EXPECT_FALSE(stepWork.spreadVsiteForcesOnGpu);
}

//! Returns a readable name for the test parameters
std::string nameOfTest(const ::testing::TestParamInfo<StepWorkloadTestParameters>& info)
{
    const auto [haveVirtualSites, useGpuUpdate, allowGpuBufferOps, pmeSetup, stepKind] = info.param;
    return formatString("%s_%s_%s_%s_%s",
                        haveVirtualSites ? "Vsites" : "NoVsites",
                        useGpuUpdate ? "GpuUpdate" : "CpuUpdate",
                        allowGpuBufferOps ? "GpuBufferOps" : "CpuBufferOps",
                        pmeSetupName(pmeSetup),
                        enumValueToString(stepKind));
}

INSTANTIATE_TEST_SUITE_P(WithAndWithoutVsites,
                         StepWorkloadTest,
                         ::testing::Combine(::testing::Bool(),
                                            ::testing::Bool(),
                                            ::testing::Bool(),
                                            ::testing::Values(PmeSetup::NoSeparatePmeRank,
                                                              PmeSetup::StagedGpuPmeOnSeparateRank),
                                            ::testing::Values(StepKind::Plain,
                                                              StepKind::Virial,
                                                              StepKind::ForceOutput,
                                                              StepKind::VirialAndForceOutput,
                                                              StepKind::NeighborSearch)),
                         nameOfTest);

} // namespace
} // namespace test
} // namespace gmx
