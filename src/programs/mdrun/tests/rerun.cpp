/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014,2015,2016, by the GROMACS development team, led by
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
 * \todo Add other tests for mdrun -rerun, e.g.
 * - TpiExitsNormally (since it uses the -rerun machinery)
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/options/filenameoption.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"

#include "testutils/testasserts.h"

#include "energyreader.h"
#include "mdruncomparisonfixture.h"
#include "trajectoryreader.h"

namespace gmx
{
namespace test
{
namespace
{

/*! \brief Manages returning a pair of frames from a normal run
 * and its rerun that are meaningful to compare.
 *
 * This class has an interface similar to one in
 * exactcontinuation.cpp, but not similar enough to warrant extracting
 * an interface class. */
template <class FrameReader, class Frame>
class FramePairManager
{
    public:
        //! Convenience typedef
        typedef std::pair<Frame, Frame> FramePair;
        //! Convenience typedef
        typedef std::unique_ptr<FrameReader> FrameReaderPtr;
        //! Constructor
        FramePairManager(FrameReaderPtr normalRun,
                         FrameReaderPtr rerun) :
            normalRun_(normalRun.release()),
            rerun_(rerun.release())
        {}
        /*! \brief Probe for a pair of valid frames, and return true if both are found.
         *
         * Give a test failure if exactly one frame is found, because
         * that file is longer than the other one, and this is not
         * expected behaviour. */
        bool shouldContinueComparing()
        {
            if (normalRun_->readNextFrame())
            {
                if (rerun_->readNextFrame())
                {
                    // Two valid next frames exist, so we should continue comparing.
                    return true;
                }
                else
                {
                    ADD_FAILURE() << "normal run energy file had at least one more frame than rerun energy file";
                }
            }
            else
            {
                if (rerun_->readNextFrame())
                {
                    ADD_FAILURE() << "rerun energy file had at least one more frame than normal run energy file";
                }
                else
                {
                    // Both files ran out of frames at the same time, which is the expected behaviour.
                }
            }
            // At least one file is out of frames, so should not continue comparing.
            return false;
        }
        /*! \brief Read and return a pair of frames that are
         * meaningful to compare. */
        FramePair getPair()
        {
            return std::make_pair(normalRun_->frame(), rerun_->frame());
        }

    private:
        FrameReaderPtr normalRun_;
        FrameReaderPtr rerun_;
};

//! Describes parameters for individual test cases.
struct RerunTestParams
{
    //! Name of the simulation in the database.
    const char                     *simulationName;
    //! Name of the integrator to test.
    const char                     *integrator;
    //! Names of energies that must compare the same.
    const std::vector<std::string>  energyNames;
    //! Tolerance within with energies and trajectories must compare.
    FloatingPointTolerance          tolerance;
    //! Max warnings grompp must tolerate.
    int                             maxWarningsTolerated;
};

//! Pretty printer for failing tests.
::std::ostream &operator<<(::std::ostream &os, const RerunTestParams &p)
{
    return os << "simulation '" << p.simulationName << "' with integrator '" << p.integrator << "'";
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
 * energy quantities are not generally reproducible, either.
 *
 * \todo Test for valid restart from non-nstlist steps - needs wider
 * tolerance for the first step after restart, because the pair list
 * will differ. */
class MdrunRerunTestFixtureBase : public MdrunComparisonFixture,
                                  public ::testing::WithParamInterface<RerunTestParams>
{
    public:
        MdrunRerunTestFixtureBase() :
            normalRunTrajectoryFileName_(fileManager_.getTemporaryFilePath("normal.trr")),
            normalRunEdrFileName_       (fileManager_.getTemporaryFilePath("normal.edr")),
            rerunTrajectoryFileName_    (fileManager_.getTemporaryFilePath("rerun.trr")),
            rerunEdrFileName_           (fileManager_.getTemporaryFilePath("rerun.edr")),
            simulationsWereRun_(false)
        {
        }
        // Import convenience wrapper from parent class.
        using MdrunComparisonFixture::runGrompp;
        //! Run a normal mdrun and its rerun.
        void runTest()
        {
            ASSERT_EQ(true, gromppWasRun_);

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
        void checkResults(const std::vector<std::string> &requiredEnergies,
                          FloatingPointTolerance          tolerance)
        {
            ASSERT_EQ(true, simulationsWereRun_);

            // Prepare object to read the energy files and return
            // pairs of frames to compare
            FramePairManager<EnergyFrameReader, EnergyFrame>
            energyManager(openEnergyFileToReadFields(normalRunEdrFileName_, requiredEnergies),
                          openEnergyFileToReadFields(rerunEdrFileName_, requiredEnergies));
            while (energyManager.shouldContinueComparing())
            {
                auto frames = energyManager.getPair();
                compareFrames(frames.first, frames.second, tolerance);
            }

            // Prepare object to read the trajectory files and return
            // pairs of frames to compare
            FramePairManager<TrajectoryFrameReader, TrajectoryFrame>
            trajectoryManager(TrajectoryFrameReaderPtr(new TrajectoryFrameReader(normalRunTrajectoryFileName_)),
                              TrajectoryFrameReaderPtr(new TrajectoryFrameReader(rerunTrajectoryFileName_)));
            while (trajectoryManager.shouldContinueComparing())
            {
                auto frames = trajectoryManager.getPair();
                compareFrames(frames.first, frames.second, tolerance);
            }
        }
    private:
        //! Name of file used in testing.
        std::string normalRunTrajectoryFileName_;
        //! Name of file used in testing.
        std::string normalRunEdrFileName_;
        //! Name of file used in testing.
        std::string rerunTrajectoryFileName_;
        //! Name of file used in testing.
        std::string rerunEdrFileName_;
        //! True only when simulations have been run and results can be tested.
        bool        simulationsWereRun_;
};

//! Specify the energies compared.
std::vector<std::string> potentialEnergy_g {{
                                                interaction_function[F_EPOT].longname
                                            }};

/*! \brief Specifies the cases tested.
 *
 * The choices for tolerance are arbitrary but sufficient.  Rerun does
 * pair search every frame, so it cannot in general exactly reproduce
 * quantities from a normal run, because the accumulation order
 * differs. (Nor does it reproduce pair-search frames exactly,
 * either). */
RerunTestParams mdrunParams_g[] = {
    { "argon12", "md", potentialEnergy_g, relativeToleranceAsPrecisionDependentUlp(1.0, 1, 1), 0 },
    { "argon12", "md-vv", potentialEnergy_g, relativeToleranceAsPrecisionDependentUlp(1.0, 1, 1), 0 },
    { "argon12", "bd", potentialEnergy_g, relativeToleranceAsPrecisionDependentUlp(1.0, 1, 1), 0 },
    { "argon12", "sd", potentialEnergy_g, relativeToleranceAsPrecisionDependentUlp(1.0, 1, 1), 0 },
    { "spc5", "md", potentialEnergy_g, relativeToleranceAsPrecisionDependentUlp(1.0, 40, 30), 0 },
    { "spc5", "md-vv", potentialEnergy_g, relativeToleranceAsPrecisionDependentUlp(1.0, 40, 20), 0 },
    { "spc5", "bd", potentialEnergy_g, relativeToleranceAsPrecisionDependentUlp(1.0, 100, 700), 0 },
    { "spc5", "sd", potentialEnergy_g, relativeToleranceAsPrecisionDependentUlp(1.0, 30, 10), 0 },
    { "alanine_vsite_vacuo", "md", potentialEnergy_g, relativeToleranceAsPrecisionDependentUlp(1.0, 1100, 600), 0 },
    { "alanine_vsite_vacuo", "md-vv", potentialEnergy_g, relativeToleranceAsPrecisionDependentUlp(1.0, 60, 50), 0 },
    { "alanine_vsite_vacuo", "bd", potentialEnergy_g, relativeToleranceAsPrecisionDependentUlp(1.0, 140, 150), 0 },
    { "alanine_vsite_vacuo", "sd", potentialEnergy_g, relativeToleranceAsPrecisionDependentUlp(1.0, 140, 100), 0 },
};

//! Name the test case usefully
typedef MdrunRerunTestFixtureBase WithinTolerances;

TEST_P(WithinTolerances, IsReproducedByRerun)
{
    auto        params = GetParam();
    SCOPED_TRACE(params);
    // TODO evolve grompp to report the number of warnings issued, so
    // tests always expect the right number.
    CommandLine caller;
    caller.addOption("-maxwarn", params.maxWarningsTolerated);
    runGrompp(caller, params.simulationName, params.integrator, "no", "no");
    runTest();
    checkResults(params.energyNames, params.tolerance);
}

INSTANTIATE_TEST_CASE_P(NormalMdrun, WithinTolerances, ::testing::ValuesIn(mdrunParams_g));

//! Specify the energies compared.
std::vector<std::string> freeEnergyTestEnergies_g {
    {
        interaction_function[F_EPOT].longname
    },
    {
        interaction_function[F_DVDL].longname
    },
    {
        interaction_function[F_DVDL_VDW].longname
    }
};

//! Specifies the cases tested.
RerunTestParams freeEnergyParams_g[] = {
    // This setup triggers a warning about not using this integrator
    // for nearly decoupled states, which we need to suppress.
    { "nonanol_vacuo", "md", freeEnergyTestEnergies_g, relativeToleranceAsPrecisionDependentUlp(1.0, 10, 20), 1 },
    { "nonanol_vacuo", "md-vv", freeEnergyTestEnergies_g, relativeToleranceAsPrecisionDependentUlp(1.0, 10, 20), 0 },
    { "nonanol_vacuo", "sd", freeEnergyTestEnergies_g, relativeToleranceAsPrecisionDependentUlp(1.0, 10, 40), 0 },
};

INSTANTIATE_TEST_CASE_P(FreeEnergyMdrun, WithinTolerances, ::testing::ValuesIn(freeEnergyParams_g));

} // namespace
} // namespace
} // namespace
