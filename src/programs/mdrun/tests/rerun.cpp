/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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

#include <functional>
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

#include "testutils/testasserts.h"

#include "energycomparison.h"
#include "energyreader.h"
#include "mdruncomparison.h"
#include "trajectorycomparison.h"
#include "trajectoryreader.h"

namespace gmx
{
namespace test
{
namespace
{

//! Describes parameters for individual test cases.
struct RerunTestParams
{
    //! Name of the simulation in the database.
    const char                     *simulationName;
    //! Name of the integrator to test.
    const char                     *integrator;
    //! Names of energy fields to compare with given tolerancecs
    EnergyTolerances energiesToMatch;
    //! Settings for which fields are required to match.
    TrajectoryFrameMatchSettings    matchSettings;
    //! Max warnings grompp must tolerate.
    int                             maxWarningsTolerated;
};

//! Pretty printer so SCOPED_TRACE(p) works.
::std::ostream &operator<<(::std::ostream &os, const RerunTestParams &p)
{
    return os << "Comparing normal and rerun of simulation '" << p.simulationName << "' with integrator '" << p.integrator << "'";
}

/*! \brief Test fixture for mdrun -rerun
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
 * energy quantities are not generally reproducible, either. */
class MdrunRerunTest : public MdrunTestFixture,
                       public ::testing::WithParamInterface<RerunTestParams>
{
    public:
        MdrunRerunTest() :
            normalRunTrajectoryFileName_(fileManager_.getTemporaryFilePath("normal.trr")),
            normalRunEdrFileName_       (fileManager_.getTemporaryFilePath("normal.edr")),
            rerunTrajectoryFileName_    (fileManager_.getTemporaryFilePath("rerun.trr")),
            rerunEdrFileName_           (fileManager_.getTemporaryFilePath("rerun.edr")),
            simulationsWereRun_(false)
        {
        }
        //! Run a normal mdrun and its rerun.
        void runTest()
        {
            // do a normal mdrun
            {
                runner_.fullPrecisionTrajectoryFileName_ = normalRunTrajectoryFileName_;
                runner_.edrFileName_                     = normalRunEdrFileName_;
                CommandLine normalRunCaller;
                normalRunCaller.append("mdrun");
                ASSERT_EQ(0, runner_.callMdrun(normalRunCaller));
            }

            // do a rerun on the .trr just produced
            {
                runner_.fullPrecisionTrajectoryFileName_ = rerunTrajectoryFileName_;
                runner_.edrFileName_                     = rerunEdrFileName_;
                CommandLine rerunCaller;
                rerunCaller.append("mdrun");
                rerunCaller.addOption("-rerun", normalRunTrajectoryFileName_);
                ASSERT_EQ(0, runner_.callMdrun(rerunCaller));
            }

            simulationsWereRun_ = true;
        }
        //! Check the energies of the rerun agree with the normal run within \c tolerance.
        void checkResults(const EnergyTolerances     &energiesToMatch,
                          const TrajectoryFrameMatchSettings &matchSettings)
        {
            // TODO get rid of this
            auto params = GetParam();
            ASSERT_EQ(true, simulationsWereRun_);
            using namespace std::placeholders;

            auto namesOfEnergiesToMatch = getKeys(energiesToMatch);
            // Prepare object to read the energy files and return
            // pairs of frames to compare
            FramePairManager<EnergyFrameReader, EnergyFrame>
            energyManager(openEnergyFileToReadFields(normalRunEdrFileName_, namesOfEnergiesToMatch),
                          openEnergyFileToReadFields(rerunEdrFileName_, namesOfEnergiesToMatch));
            energyManager.compareAllFramePairs(std::bind(compareEnergyFrames, _1, _2, energiesToMatch));

            // Prepare object to read the trajectory files and return
            // pairs of frames to compare
            FramePairManager<TrajectoryFrameReader, TrajectoryFrame>
            trajectoryManager(TrajectoryFrameReaderPtr(new TrajectoryFrameReader(normalRunTrajectoryFileName_)),
                              TrajectoryFrameReaderPtr(new TrajectoryFrameReader(rerunTrajectoryFileName_)));
            trajectoryManager.compareAllFramePairs(std::bind(compareTrajectoryFrames, _1, _2, matchSettings, s_trajectoryTolerances.at(params.simulationName)));
        }
    private:
        //!@{
        //! Names of files used in testing.
        std::string normalRunTrajectoryFileName_;
        std::string normalRunEdrFileName_;
        std::string rerunTrajectoryFileName_;
        std::string rerunEdrFileName_;
        //!@}
        //! True only when simulations have been run and results can be tested.
        bool        simulationsWereRun_;
    public:
        /*! \brief Specifies the expected tolerances for trajectory components for different simulation systems. */
        static std::map<std::string, TrajectoryTolerances> s_trajectoryTolerances;
};

// TODO improve
//! Specify how trajectory frame matching must work.
TrajectoryFrameMatchSettings g_matchSettings {
    true, true, true, true, true, true
};

FloatingPointTolerance g_boxTolerance      = defaultRealTolerance();
FloatingPointTolerance g_positionTolerance = relativeToleranceAsFloatingPoint(1.0, 1.0e-3);
FloatingPointTolerance g_velocityTolerance = defaultRealTolerance(); // unused for rerun
FloatingPointTolerance g_forceTolerance    = relativeToleranceAsFloatingPoint(100.0, GMX_DOUBLE ? 1.0e-7 : 1.0e-5);

std::map<std::string, TrajectoryTolerances> MdrunRerunTest::s_trajectoryTolerances =
{
    { "argon12", { g_boxTolerance, g_positionTolerance, g_velocityTolerance, g_forceTolerance } },
    { "argon12", { g_boxTolerance, g_positionTolerance, g_velocityTolerance, g_forceTolerance } },
    { "spc5", { g_boxTolerance, g_positionTolerance, g_velocityTolerance, g_forceTolerance } },
    { "alanine_vacuo", { g_boxTolerance, g_positionTolerance, g_velocityTolerance, g_forceTolerance } },
    { "alanine_vsite_vacuo", { g_boxTolerance, g_positionTolerance, g_velocityTolerance, g_forceTolerance } },
    { "nonanol_vacuo", { g_boxTolerance, g_positionTolerance, g_velocityTolerance, g_forceTolerance } },
};

//! Specify the energies to compare with a given tolerance.
EnergyTolerances g_energiesToMatch = {{
        { interaction_function[F_EPOT].longname, defaultRealTolerance() }
    }};

/*! \brief Specifies the cases tested.
 *
 * The choices for tolerance are arbitrary but sufficient.  Rerun does
 * pair search every frame, so it cannot in general exactly reproduce
 * quantities from a normal run, because the accumulation order
 * differs. (Nor does it reproduce pair-search frames exactly,
 * either). */
RerunTestParams g_mdrunParams[] = {
    { "argon12", "md", g_energiesToMatch, g_matchSettings, 0 },
    { "argon12", "md-vv", g_energiesToMatch, g_matchSettings, 0 },
    { "argon12", "bd", g_energiesToMatch, g_matchSettings, 0 },
    { "argon12", "sd", g_energiesToMatch, g_matchSettings, 0 },
    { "spc5", "md", g_energiesToMatch, g_matchSettings, 0 },
    { "spc5", "md-vv", g_energiesToMatch, g_matchSettings, 0 },
    { "spc5", "bd", g_energiesToMatch, g_matchSettings, 0 },
    { "spc5", "sd", g_energiesToMatch, g_matchSettings, 0 },
    { "alanine_vsite_vacuo", "md", g_energiesToMatch, g_matchSettings, 0 },
    { "alanine_vsite_vacuo", "md-vv", g_energiesToMatch, g_matchSettings, 0 },
    { "alanine_vsite_vacuo", "bd", g_energiesToMatch, g_matchSettings, 0 },
    { "alanine_vsite_vacuo", "sd", g_energiesToMatch, g_matchSettings, 0 },
};

TEST_P(MdrunRerunTest, WithinTolerances)
{
    auto        params = GetParam();
    SCOPED_TRACE(params);
    {
        // TODO evolve grompp to report the number of warnings issued, so
        // tests always expect the right number.
        CommandLine caller;
        caller.append("grompp");
        caller.addOption("-maxwarn", params.maxWarningsTolerated);
        runner_.useTopGroAndNdxFromDatabase(params.simulationName);
        auto mdpFieldValues = prepareMdpFieldValues(params.simulationName,
                                                    params.integrator,
                                                    "no", "no");
        runner_.useStringAsMdpFile(prepareMdpFileContents(mdpFieldValues));
        EXPECT_EQ(0, runner_.callGrompp(caller));
    }
    runTest();
    checkResults(params.energiesToMatch, params.matchSettings);
}

INSTANTIATE_TEST_CASE_P(NormalMdrunIsReproduced, MdrunRerunTest, ::testing::ValuesIn(g_mdrunParams));

//! Specify the energies to compare with a given tolerance.
EnergyTolerances g_energiesToMatchForFreeEnergyCase = {{
        { interaction_function[F_EPOT].longname, defaultRealTolerance() },
        { interaction_function[F_DVDL].longname, defaultRealTolerance() },
        { interaction_function[F_DVDL_VDW].longname, defaultRealTolerance() }
    }};

//! Specifies the cases tested.
RerunTestParams g_freeEnergyParams[] = {
    // This setup triggers a warning about not using this integrator
    // for nearly decoupled states, which we need to suppress.
    { "nonanol_vacuo", "md", g_energiesToMatchForFreeEnergyCase, g_matchSettings, 1 },
    { "nonanol_vacuo", "md-vv", g_energiesToMatchForFreeEnergyCase, g_matchSettings, 0 },
    { "nonanol_vacuo", "sd", g_energiesToMatchForFreeEnergyCase, g_matchSettings, 0 },
};

INSTANTIATE_TEST_CASE_P(FreeEnergyMdrunIsReproduced, MdrunRerunTest, ::testing::ValuesIn(g_freeEnergyParams));

} // namespace
} // namespace
} // namespace
