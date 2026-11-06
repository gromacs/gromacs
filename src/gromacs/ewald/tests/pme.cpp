/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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
 * Tests for functionality of the pme related classes:
 * class SeparatePmeRanksPermitted
 *
 * \author Dmitry Morozov <dmitry.morozov@jyu.fi>
 * \ingroup module_ewald
 */
#include "gmxpre.h"

#include "gromacs/ewald/pme.h"

#include "config.h"

#include <array>
#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/domdec/domdec.h"
#include "gromacs/ewald/pme_grid.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringcompare.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/testasserts.h"
#include "testutils/testmatchers.h"

#include "pmetestcommon.h"

namespace gmx
{

namespace test
{

class SeparatePmeRanksPermittedTest : public ::testing::Test
{
public:
    void disableFirstReason() { separatePmeRanksPermitted_.disablePmeRanks("First reason"); }

    void disableSecondReason() { separatePmeRanksPermitted_.disablePmeRanks("Second reason"); }

    void disableEmptyReason() { separatePmeRanksPermitted_.disablePmeRanks(""); }

protected:
    SeparatePmeRanksPermitted separatePmeRanksPermitted_;
};

TEST_F(SeparatePmeRanksPermittedTest, ZeroPmeDisableReasons)
{
    // Expect that SeparatePmeRanksPermitted is enabled by default
    EXPECT_TRUE(separatePmeRanksPermitted_.permitSeparatePmeRanks());
}

TEST_F(SeparatePmeRanksPermittedTest, CanBeDisabled)
{
    // Test if disablePmeRanks works
    EXPECT_NO_THROW(disableFirstReason(););
}

TEST_F(SeparatePmeRanksPermittedTest, OneDisableReasonFlag)
{
    disableFirstReason();

    // Expect that SeparatePmeRanksPermitted is disabled now
    EXPECT_FALSE(separatePmeRanksPermitted_.permitSeparatePmeRanks());
}

TEST_F(SeparatePmeRanksPermittedTest, OneDisableReasonText)
{
    disableFirstReason();

    // Expect that reasonsWhyDisabled works with one reason
    EXPECT_TRUE(separatePmeRanksPermitted_.reasonsWhyDisabled() == "First reason");
}

TEST_F(SeparatePmeRanksPermittedTest, TwoDisableReasonText)
{
    disableFirstReason();
    disableSecondReason();

    // Expect that reasonsWhyDisabled works with two reasons
    EXPECT_TRUE(separatePmeRanksPermitted_.reasonsWhyDisabled() == "First reason; Second reason");
}

TEST_F(SeparatePmeRanksPermittedTest, EmptyDisableReasonText)
{
    disableEmptyReason();

    // Expect that reasonsWhyDisabled works with empty reason
    EXPECT_TRUE(separatePmeRanksPermitted_.reasonsWhyDisabled().empty());
}
TEST(PmeStepWorkloadTest, ConvertsRelevantStepWork)
{
    StepWorkload stepWork;
    EXPECT_FALSE(PmeStepWorkload{ stepWork }.computeEnergyAndVirial);

    stepWork.computeEnergy = true;
    EXPECT_TRUE(PmeStepWorkload{ stepWork }.computeEnergyAndVirial);

    stepWork.computeEnergy       = false;
    stepWork.computeVirial       = true;
    stepWork.computeForces       = true;
    stepWork.useGpuPmeFReduction = true;
    const PmeStepWorkload pmeWork{ stepWork };

    EXPECT_TRUE(pmeWork.computeEnergyAndVirial);
    EXPECT_TRUE(pmeWork.computeForces);
    EXPECT_TRUE(pmeWork.useGpuPmeFReduction);
}

TEST(PmeOnlySimulationWorkloadTest, ConvertsRelevantSimulationWork)
{
    SimulationWorkload simulationWork;
    simulationWork.useGpuPmePpCommunication = true;
    simulationWork.useNvshmem               = true;
    simulationWork.useGpuHaloExchange       = true;
    simulationWork.useMdGpuGraph            = true;
    simulationWork.useGpuPme                = true;
    simulationWork.useGpuPmeFft             = false;

    const PmeOnlySimulationWorkload pmeWork{ simulationWork };

    EXPECT_TRUE(pmeWork.useGpu);
    EXPECT_TRUE(pmeWork.useGpuPmePpCommunication);
    EXPECT_TRUE(pmeWork.useNvshmem);
    EXPECT_TRUE(pmeWork.useGpuHaloExchange);
    EXPECT_TRUE(pmeWork.useMdGpuGraph);

    simulationWork.useGpuPme = false;
    const PmeOnlySimulationWorkload cpuPmeWork{ simulationWork };
    EXPECT_FALSE(cpuPmeWork.useGpu);
}


TEST(PmeGridOverlapTest, UsesAlignedZStrideForMajorOverlapBuffer)
{
    constexpr int pmeOrder  = 5;
    constexpr int gridY     = 120;
    constexpr int gridZ     = 121;
    constexpr int minorRank = 1;

    const int unalignedGridZ = gridZ + pmeOrder - 1;
    const int alignedGridZ   = pmeGridAlignedZSize(gridZ, pmeOrder);
    if (alignedGridZ == unalignedGridZ)
    {
        GTEST_SKIP() << "PME z grid alignment is not active in this build";
    }

    const int gridYBlockSize = (gridY + minorRank - 1) / minorRank + pmeOrder;
    const int maxLocalGridY  = gridYBlockSize - 1;

    const int oldBufferSize      = pmeOrder * gridYBlockSize * unalignedGridZ;
    const int requiredBufferSize = pmeOrder * maxLocalGridY * alignedGridZ;
    const int fixedBufferSize    = pmeOrder * gridYBlockSize * alignedGridZ;

    EXPECT_GT(requiredBufferSize, oldBufferSize);
    EXPECT_LE(requiredBufferSize, fixedBufferSize);
}

namespace
{
constexpr matrix        c_fepTestBox = { { 3.0_real, 0.0_real, 0.0_real },
                                         { 0.0_real, 3.0_real, 0.0_real },
                                         { 0.0_real, 0.0_real, 3.0_real } };
const CoordinatesVector c_fepTestX   = { { 0.30_real, 0.70_real, 1.10_real },
                                         { 1.40_real, 0.50_real, 2.20_real },
                                         { 2.10_real, 1.90_real, 0.40_real },
                                         { 0.80_real, 2.50_real, 1.70_real } };
const ChargesVector     c_fepTestQA  = { 0.85_real, -0.85_real, 0.40_real, -0.40_real };
const ChargesVector     c_fepTestQB  = { 0.00_real, 0.00_real, 0.40_real, -0.40_real };
} // namespace

TEST(PmeFepEndpointGridSkipTest, CpuSkipsZeroWeightGridAcrossCoulombAndLjPme)
{
    const ChargesVector c6A    = { 0.0024_real, 0.0018_real, 0.0031_real, 0.0015_real };
    const ChargesVector c6B    = { 0.0000_real, 0.0000_real, 0.0031_real, 0.0015_real };
    const ChargesVector sigmaA = { 0.34_real, 0.31_real, 0.36_real, 0.30_real };
    const ChargesVector sigmaB = { 0.34_real, 0.31_real, 0.36_real, 0.30_real };

    struct CpuCase
    {
        CoulombInteractionType coulombType;
        VanDerWaalsType        vdwType;
        LongRangeVdW           ljComb;
        bool                   fepQ;
        bool                   fepLj;
        real                   lambdaQ;
        real                   lambdaLj;
        bool                   expectHalvedSpreadAndFft;
    };
    const std::array<CpuCase, 8> cases = { {
            { CoulombInteractionType::Pme, VanDerWaalsType::Cut, LongRangeVdW::Geom, true, false, 0.0_real, 0.0_real, true },
            { CoulombInteractionType::Pme, VanDerWaalsType::Cut, LongRangeVdW::Geom, true, false, 1.0_real, 0.0_real, true },
            { CoulombInteractionType::Pme, VanDerWaalsType::Cut, LongRangeVdW::Geom, true, false, 1e-7_real, 0.0_real, false },
            { CoulombInteractionType::Cut, VanDerWaalsType::Pme, LongRangeVdW::Geom, false, true, 0.0_real, 0.0_real, true },
            { CoulombInteractionType::Cut, VanDerWaalsType::Pme, LongRangeVdW::Geom, false, true, 0.0_real, 1.0_real, true },
            { CoulombInteractionType::Cut, VanDerWaalsType::Pme, LongRangeVdW::LB, false, true, 0.0_real, 0.0_real, true },
            { CoulombInteractionType::Cut, VanDerWaalsType::Pme, LongRangeVdW::LB, false, true, 0.0_real, 1.0_real, true },
            { CoulombInteractionType::Pme, VanDerWaalsType::Pme, LongRangeVdW::Geom, true, true, 1.0_real, 0.0_real, true },
    } };

    MDLogger                     dummyLogger;
    const NumPmeDomains          numPmeDomains{ 1, 1 };
    const FloatingPointTolerance tol =
            relativeToleranceAsPrecisionDependentFloatingPoint(1.0, 1e-5F, 1e-12);

#if GMX_OPENMP
    const std::vector<int> threadCounts = { 1, 2 };
#else
    const std::vector<int> threadCounts = { 1 };
#endif

    for (const int numThreads : threadCounts)
    {
        for (size_t caseIndex = 0; caseIndex < cases.size(); ++caseIndex)
        {
            const CpuCase& tc = cases[caseIndex];
            SCOPED_TRACE(formatString("case=%zu, numThreads=%d", caseIndex, numThreads));

            t_inputrec inputRec;
            inputRec.nkx                    = 16;
            inputRec.nky                    = 16;
            inputRec.nkz                    = 16;
            inputRec.pme_order              = 4;
            inputRec.coulombtype            = tc.coulombType;
            inputRec.vdwtype                = tc.vdwType;
            inputRec.ljpme_combination_rule = tc.ljComb;
            inputRec.epsilon_r              = 1.0_real;
            inputRec.efep                   = FreeEnergyPerturbationType::Yes;

            std::unique_ptr<gmx_pme_t> pme = gmx_pme_init(nullptr,
                                                          numPmeDomains,
                                                          &inputRec,
                                                          c_fepTestBox,
                                                          1.0_real,
                                                          tc.fepQ,
                                                          tc.fepLj,
                                                          true,
                                                          2.5_real,
                                                          2.0_real,
                                                          numThreads,
                                                          PmeRunMode::CPU,
                                                          nullptr,
                                                          nullptr,
                                                          nullptr,
                                                          dummyLogger);
            gmx_pme_reinit_atoms(pme.get(), gmx::ssize(c_fepTestX), c_fepTestQA, c_fepTestQB);

            auto evalCpu = [&](bool computeEnergyAndVirial, std::vector<RVec>* fOut, t_nrnb* nrnbOut)
            {
                fOut->assign(c_fepTestX.size(), RVec{ 0.0_real, 0.0_real, 0.0_real });
                *nrnbOut     = t_nrnb{};
                matrix virQ  = { { 0 } };
                matrix virLj = { { 0 } };
                real   eQ = 0, eLj = 0, dvdlQ = 0, dvdlLj = 0;

                StepWorkload sw;
                sw.computeForces = true;
                sw.computeEnergy = computeEnergyAndVirial;
                sw.computeVirial = computeEnergyAndVirial;
                gmx_pme_do(pme.get(),
                           c_fepTestX,
                           *fOut,
                           c_fepTestQA,
                           c_fepTestQB,
                           c6A,
                           c6B,
                           sigmaA,
                           sigmaB,
                           c_fepTestBox,
                           0,
                           0,
                           nrnbOut,
                           nullptr,
                           virQ,
                           virLj,
                           &eQ,
                           &eLj,
                           tc.lambdaQ,
                           tc.lambdaLj,
                           &dvdlQ,
                           &dvdlLj,
                           PmeStepWorkload{ sw });
            };

            std::vector<RVec> fForce, fEnergy;
            t_nrnb            nrnbForce, nrnbEnergy;
            evalCpu(false, &fForce, &nrnbForce);
            evalCpu(true, &fEnergy, &nrnbEnergy);

            EXPECT_GT(nrnbEnergy.n[eNR_SPREADBSP], 0);
            const int mult = tc.expectHalvedSpreadAndFft ? 2 : 1;
            EXPECT_EQ(nrnbForce.n[eNR_SPREADBSP] * mult, nrnbEnergy.n[eNR_SPREADBSP]);
            EXPECT_EQ(nrnbForce.n[eNR_FFT] * mult, nrnbEnergy.n[eNR_FFT]);
            for (size_t i = 0; i < c_fepTestX.size(); ++i)
            {
                for (int d = 0; d < DIM; ++d)
                {
                    EXPECT_REAL_EQ_TOL(fEnergy[i][d], fForce[i][d], tol);
                }
            }
        }
    }
}

} // namespace test

} // namespace gmx
