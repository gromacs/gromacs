/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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
 * Tests for the mdrun -rerun functionality
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "config.h"

#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/options/filenameoption.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/trajectory/energyframe.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/mpitest.h"
#include "testutils/simulationdatabase.h"
#include "testutils/testasserts.h"

#include "energycomparison.h"
#include "energyreader.h"
#include "mdruncomparison.h"
#include "moduletest.h"
#include "trajectorycomparison.h"
#include "trajectoryreader.h"

namespace gmx
{
namespace test
{
namespace
{

//! Functor for comparing reference and test frames on particular energies to match.
class EnergyComparator
{
    public:
        //! Constructor
        EnergyComparator(const EnergyTolerances &energiesToMatch)
            : energiesToMatch_(energiesToMatch) {}
        //! The functor method.
        void operator()(const EnergyFrame &reference, const EnergyFrame &test)
        {
            compareEnergyFrames(reference, test, energiesToMatch_);
        }
        //! Container of the energies to match and the tolerance required.
        const EnergyTolerances &energiesToMatch_;
};

//! Run grompp for a normal mdrun, the run, and its rerun.
void executeRerunTest(TestFileManager        *fileManager,
                      SimulationRunner       *runner,
                      const std::string      &simulationName,
                      int                     maxWarningsTolerated,
                      const MdpFieldValues   &mdpFieldValues,
                      const EnergyTolerances &energiesToMatch)
{
    // TODO At some point we should also test PME-only ranks.
    int numRanksAvailable = getNumberOfTestMpiRanks();
    if (!isNumberOfPpRanksSupported(simulationName, numRanksAvailable))
    {
        fprintf(stdout, "Test system '%s' cannot run with %d ranks.\n"
                "The supported numbers are: %s\n",
                simulationName.c_str(), numRanksAvailable,
                reportNumbersOfPpRanksSupported(simulationName).c_str());
        return;
    }

    auto normalRunTrajectoryFileName = fileManager->getTemporaryFilePath("normal.trr");
    auto normalRunEdrFileName        = fileManager->getTemporaryFilePath("normal.edr");
    auto rerunTrajectoryFileName     = fileManager->getTemporaryFilePath("rerun.trr");
    auto rerunEdrFileName            = fileManager->getTemporaryFilePath("rerun.edr");

    // prepare the .tpr file
    {
        // TODO evolve grompp to report the number of warnings issued, so
        // tests always expect the right number.
        CommandLine caller;
        caller.append("grompp");
        caller.addOption("-maxwarn", maxWarningsTolerated);
        runner->useTopGroAndNdxFromDatabase(simulationName);
        runner->useStringAsMdpFile(prepareMdpFileContents(mdpFieldValues));
        EXPECT_EQ(0, runner->callGrompp(caller));
    }

    // do the normal mdrun
    {
        runner->fullPrecisionTrajectoryFileName_ = normalRunTrajectoryFileName;
        runner->edrFileName_                     = normalRunEdrFileName;
        CommandLine normalRunCaller;
        normalRunCaller.append("mdrun");
        ASSERT_EQ(0, runner->callMdrun(normalRunCaller));
    }

    // do a rerun on the .trr just produced
    {
        runner->fullPrecisionTrajectoryFileName_ = rerunTrajectoryFileName;
        runner->edrFileName_                     = rerunEdrFileName;
        CommandLine rerunCaller;
        rerunCaller.append("mdrun");
        rerunCaller.addOption("-rerun", normalRunTrajectoryFileName);
        ASSERT_EQ(0, runner->callMdrun(rerunCaller));
    }

    // Build the functor that will compare reference and test
    // energy frames on the chosen energy fields.
    //
    // TODO It would be less code if we used a lambda for this, but either
    // clang 3.4 or libstdc++ 5.2.1 have an issue with capturing a
    // std::unordered_map
    EnergyComparator energyComparator(energiesToMatch);
    // Build the manager that will present matching pairs of frames to compare.
    //
    // TODO Here is an unnecessary copy of keys (ie. the energy field
    // names), for convenience. In the future, use a range.
    auto namesOfEnergiesToMatch = getKeys(energiesToMatch);
    FramePairManager<EnergyFrameReader, EnergyFrame>
         energyManager(openEnergyFileToReadFields(normalRunEdrFileName, namesOfEnergiesToMatch),
                  openEnergyFileToReadFields(rerunEdrFileName, namesOfEnergiesToMatch));
    // Compare the energy frames.
    energyManager.compareAllFramePairs(energyComparator);

    // Specify how trajectory frame matching must work.
    TrajectoryFrameMatchSettings trajectoryMatchSettings {
        true, true, true, true, false, true
    };
    /* Specify the default expected tolerances for trajectory
     * components for all simulation systems. */
    TrajectoryTolerances trajectoryTolerances {
        defaultRealTolerance(),                                               // box
        relativeToleranceAsFloatingPoint(1.0, 1.0e-3),                        // positions
        defaultRealTolerance(),                                               // velocities - unused for rerun
        relativeToleranceAsFloatingPoint(100.0, GMX_DOUBLE ? 1.0e-7 : 1.0e-5) // forces
    };

    // Build the functor that will compare reference and test
    // trajectory frames in the chosen way.
    auto trajectoryComparator = [&trajectoryMatchSettings, &trajectoryTolerances](const TrajectoryFrame &reference, const TrajectoryFrame &test)
        {
            compareTrajectoryFrames(reference, test, trajectoryMatchSettings, trajectoryTolerances);
        };
    // Build the manager that will present matching pairs of frames to compare
    FramePairManager<TrajectoryFrameReader, TrajectoryFrame>
    trajectoryManager(std::make_unique<TrajectoryFrameReader>(normalRunTrajectoryFileName),
                      std::make_unique<TrajectoryFrameReader>(rerunTrajectoryFileName));
    // Compare the trajectory frames.
    trajectoryManager.compareAllFramePairs(trajectoryComparator);
}

/*! \brief Test fixture base for mdrun -rerun
 *
 * This test ensures mdrun can run a simulation, writing a trajectory
 * and matching energies, and reproduce the same energies from a rerun
 * to within a tight tolerance. It says nothing about whether a rerun
 * can reproduce energies from a trajectory generated with older code,
 * since that is not a useful property. Whether mdrun produced correct
 * energies then and now needs different kinds of testing, but if
 * true, this test ensures the rerun has the expected property.
 *
 * The limitations of mdrun and its output means that reproducing the
 * same energies is currently only meaningful for integration without
 * thermostats or barostats, however the present form of the test
 * infrastructure has in-principle support for such, if that is ever
 * needed/useful.
 *
 * We should also not compare pressure, because with constraints the
 * non-search steps need a much larger tolerance, and per Redmine 1868
 * we should stop computing pressure in reruns anyway.
 *
 * Similarly, per 1868, in the present implementation the kinetic
 * energy quantities are not generally reproducible, either.
 *
 * The choices for tolerance are arbitrary but sufficient.  Rerun does
 * pair search every frame, so it cannot in general exactly reproduce
 * quantities from a normal run, because the accumulation order
 * differs. (Nor does it reproduce pair-search frames exactly,
 * either). */
class MdrunRerunTest : public MdrunTestFixture,
                       public ::testing::WithParamInterface <
                       std::tuple < std::string, std::string>>
{
};

TEST_P(MdrunRerunTest, WithinTolerances)
{
    auto params         = GetParam();
    auto simulationName = std::get<0>(params);
    auto integrator     = std::get<1>(params);
    SCOPED_TRACE(formatString("Comparing normal and rerun of simulation '%s' "
                              "with integrator '%s'",
                              simulationName.c_str(), integrator.c_str()));

    auto mdpFieldValues = prepareMdpFieldValues(simulationName.c_str(),
                                                integrator.c_str(),
                                                "no", "no");

    // bd is much less reproducible in a rerun than the other integrators
    const int        toleranceScaleFactor = (integrator == "bd") ? 2 : 1;
    EnergyTolerances energiesToMatch
    {{
         {
             interaction_function[F_EPOT].longname,
             relativeToleranceAsPrecisionDependentUlp(10.0, 24 * toleranceScaleFactor, 40 * toleranceScaleFactor)
         },
     }};

    int numWarningsToTolerate = 0;
    executeRerunTest(&fileManager_, &runner_,
                     simulationName, numWarningsToTolerate, mdpFieldValues,
                     energiesToMatch);
}

// TODO The time for OpenCL kernel compilation means these tests time
// out. Once that compilation is cached for the whole process, these
// tests can run in such configurations.
#if GMX_GPU != GMX_GPU_OPENCL
INSTANTIATE_TEST_CASE_P(NormalMdrunIsReproduced, MdrunRerunTest,
                            ::testing::Combine(::testing::Values("argon12", "tip3p5", "alanine_vsite_vacuo"),
                                                   ::testing::Values("md", "md-vv", "bd", "sd")));
#else
INSTANTIATE_TEST_CASE_P(DISABLED_NormalMdrunIsReproduced, MdrunRerunTest,
                            ::testing::Combine(::testing::Values("argon12", "tip3p5", "alanine_vsite_vacuo"),
                                                   ::testing::Values("md", "md-vv", "bd", "sd")));
#endif

class MdrunRerunFreeEnergyTest : public MdrunTestFixture,
                                 public ::testing::WithParamInterface <
                                 std::tuple < std::string, std::string, int>>
{
};

TEST_P(MdrunRerunFreeEnergyTest, WithinTolerances)
{
    auto params          = GetParam();
    auto simulationName  = std::get<0>(params);
    auto integrator      = std::get<1>(params);
    auto initLambdaState = std::get<2>(params);
    SCOPED_TRACE(formatString("Comparing normal and rerun of simulation '%s' "
                              "with integrator '%s' for initial lambda state %d",
                              simulationName.c_str(), integrator.c_str(), initLambdaState));

    auto mdpFieldValues = prepareMdpFieldValues(simulationName.c_str(),
                                                integrator.c_str(),
                                                "no", "no");
    mdpFieldValues["other"] += formatString("\ninit-lambda-state %d", initLambdaState);

    EnergyTolerances energiesToMatch
    {{
         {
             interaction_function[F_EPOT].longname, relativeToleranceAsPrecisionDependentUlp(10.0, 24, 32)
         },
         {
             interaction_function[F_DVDL_COUL].longname, relativeToleranceAsPrecisionDependentUlp(1.0, 8, 8)
         },
         {
             interaction_function[F_DVDL_VDW].longname, relativeToleranceAsPrecisionDependentUlp(1.0, 8, 8)
         },
         {
             interaction_function[F_DVDL_BONDED].longname, relativeToleranceAsPrecisionDependentUlp(1.0, 8, 8)
         },
         {
             interaction_function[F_DVDL_RESTRAINT].longname, relativeToleranceAsPrecisionDependentUlp(1.0, 8, 8)
         }
     }};

    // The md integrator triggers a warning for nearly decoupled
    // states, which we need to suppress. TODO sometimes?
    int numWarningsToTolerate = (integrator == "md") ? 1 : 0;
    executeRerunTest(&fileManager_, &runner_,
                     simulationName, numWarningsToTolerate, mdpFieldValues,
                     energiesToMatch);
}

// TODO The time for OpenCL kernel compilation means these tests time
// out. Once that compilation is cached for the whole process, these
// tests can run in such configurations.
#if GMX_GPU != GMX_GPU_OPENCL
INSTANTIATE_TEST_CASE_P(MdrunIsReproduced, MdrunRerunFreeEnergyTest,
                            ::testing::Combine(::testing::Values("nonanol_vacuo"),
                                                   ::testing::Values("md", "md-vv", "sd"),
                                                   ::testing::Range(0, 11)));
#else
INSTANTIATE_TEST_CASE_P(DISABLED_MdrunIsReproduced, MdrunRerunFreeEnergyTest,
                            ::testing::Combine(::testing::Values("nonanol_vacuo"),
                                                   ::testing::Values("md", "md-vv", "sd"),
                                                   ::testing::Range(0, 11)));
#endif

} // namespace
} // namespace test
} // namespace gmx
