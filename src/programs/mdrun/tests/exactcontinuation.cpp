/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
 * Tests that mdrun restarts are exact, because a two-part run
 * reproduces a single-part run.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include <array>
#include <string>
#include <vector>

#include <tuple>
#include <gtest/gtest.h>

#include "gromacs/options/filenameoption.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
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

/*! \brief Manages returning a pair of frames that are meaningful to
 * compare, having come from a normal run and two equivalent runs that
 * did an exact continuation.
 *
 * This class has an interface similar to one in
 * rerun.cpp, but not similar enough to warrant extracting
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
        FramePairManager(FrameReaderPtr fullRun,
                         FrameReaderPtr firstHalfRun,
                         FrameReaderPtr secondHalfRun) :
            fullRun_(fullRun.release()),
            firstHalfRun_(firstHalfRun.release()),
            secondHalfRun_(secondHalfRun.release()),
            isFirstHalf_(true)
        {}
        /*! \brief Probe for a pair of valid frames, and return true if both are found.
         *
         * Give a test failure if exactly one frame is found, because
         * that file is longer than the other one, and this is not
         * expected behaviour. */
        bool shouldContinueComparing()
        {
            bool haveNormalRunFrame = fullRun_->readNextFrame();
            bool haveHalfRunFrame = isFirstHalf_ ? firstHalfRun_->readNextFrame() : secondHalfRun_->readNextFrame();
            if (isFirstHalf_ && !haveHalfRunFrame)
            {
                // Presumably we ran out of frames for the first half
                isFirstHalf_         = false;
                // The first frame is a duplicate of the last frame in the first half, so discard it.
                haveHalfRunFrame = secondHalfRun_->readNextFrame();
                if (haveHalfRunFrame)
                {
                    secondHalfRun_->frame();
                    haveHalfRunFrame = secondHalfRun_->readNextFrame();
                }
                else
                {
                    ADD_FAILURE() << "No frames in the second half run, so the normal run had at least one more frame than it";
                }
            }

            if (haveNormalRunFrame)
            {
                if (haveHalfRunFrame)
                {
                    // Two valid next frames exist, so we should continue comparing.
                    return true;
                }
                else
                {
                    ADD_FAILURE() << "normal run energy file had at least one more frame than "
                                  << (isFirstHalf_ ? "first" : "second") << " half run energy file";
                }
            }
            else
            {
                if (haveHalfRunFrame)
                {
                    ADD_FAILURE() << (isFirstHalf_ ? "first" : "second")
                                  << " half run energy file had at least one more frame than normal run energy file";
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
            return std::make_pair(fullRun_->frame(), isFirstHalf_ ? firstHalfRun_->frame() : secondHalfRun_->frame());
        }

    private:
        FrameReaderPtr fullRun_;
        FrameReaderPtr firstHalfRun_;
        FrameReaderPtr secondHalfRun_;
        bool isFirstHalf_;
};

//! Describes parameters for individual test cases.
struct ExactContinuationTestParams
{
    //! Name of the simulation in the database.
    const char                     *simulationName;
    //! Name of the integrator to test.
    const char                     *integrator;
        //! Name of the kind of temperature coupling to test.
        const char *temperatureCoupling;
        //! Name of the kind of pressure coupling to test.
        const char *pressureCoupling;
    //! Names of energies that must compare the same.
    const std::vector<std::string>  energyNames;
    //! Settings for which fields are required to match.
    TrajectoryFrameMatchSettings    matchSettings;
    //! Tolerance within with quantities must compare.
    FloatingPointTolerance          tolerance;
    //! Max warnings grompp must tolerate.
    int                             maxWarningsTolerated;
};

//! Pretty printer for failing tests.
::std::ostream &operator<<(::std::ostream &os, const ExactContinuationTestParams &p)
{
    return os << "Comparing a normal simulation using integrator '" << p.integrator <<
        "', temperature coupling '" << p.temperatureCoupling <<
        "', and pressure coupling '" << p.pressureCoupling << "' on system '" <<
        p.simulationName << "' with one that stopped half way and was exactly continued";
}

/*! \brief Test fixture for mdrun exact continuations
 *
 * This test ensures mdrun can run a simulation, writing a trajectory
 * and matching energies, and reproduce to within a tight tolerance
 * the same energies from runs that stopped half way, and restarted
 * from the checkpoint.
 *
 * \todo Is there value in testing with mdrun -reprod? As well as
 * without it?
 *
 * \todo Add constraints, virtual sites peptide-in-water input
 * cases, so they are also tested. Also free energy? */
class MdrunContinuationTestFixture : public MdrunComparisonFixture,
                                     public ::testing::WithParamInterface<ExactContinuationTestParams>
{
    public:
        //! Constructor
        MdrunContinuationTestFixture() :
            fullRunTrajectoryFileName_(fileManager_.getTemporaryFilePath("full.trr")),
            fullRunEdrFileName_(fileManager_.getTemporaryFilePath("full.edr")),
            firstHalfRunTrajectoryFileName_(fileManager_.getTemporaryFilePath("firsthalf.trr")),
            firstHalfRunEdrFileName_(fileManager_.getTemporaryFilePath("firsthalf.edr")),
            firstHalfRunCheckpointFileName_(fileManager_.getTemporaryFilePath("firsthalf.cpt")),
            secondHalfRunTrajectoryFileName_(fileManager_.getTemporaryFilePath("secondhalf")),
            secondHalfRunEdrFileName_(fileManager_.getTemporaryFilePath("secondhalf"))
        {}
        // Import convenience wrapper from parent class.
        using MdrunComparisonFixture::runGrompp;
        /*! \brief Run a normal mdrun, the same mdrun stopping half
         * way, doing a continuation, and compare the results. */
        void runTest()
        {
            ASSERT_EQ(true, gromppWasRun_);

            // do a normal mdrun
            {
                runner_.fullPrecisionTrajectoryFileName_ = fullRunTrajectoryFileName_;
                runner_.edrFileName_ = fullRunEdrFileName_;
                CommandLine fullRunCaller;
                fullRunCaller.append("mdrun");
                ASSERT_EQ(0, runner_.callMdrun(fullRunCaller));
            }

            // do a repeat of the first half of the same mdrun
            {
                // Compute half the number of steps that the normal mdrun did.
                auto nstepsValue = queryMdpKeyValue("nsteps");
                auto halfNsteps = formatString("%" GMX_PRId64, str_to_int64_t(nstepsValue.c_str(), nullptr) / 2);

                runner_.fullPrecisionTrajectoryFileName_ = firstHalfRunTrajectoryFileName_;
                runner_.edrFileName_ = firstHalfRunEdrFileName_;
                CommandLine firstHalfRunCaller;
                firstHalfRunCaller.append("mdrun");
                firstHalfRunCaller.addOption("-nsteps", halfNsteps);
                firstHalfRunCaller.addOption("-cpo", firstHalfRunCheckpointFileName_);
                ASSERT_EQ(0, runner_.callMdrun(firstHalfRunCaller));
            }

            // do a continuation from the first half of that same mdrun
            {
                runner_.fullPrecisionTrajectoryFileName_ = secondHalfRunTrajectoryFileName_;
                runner_.edrFileName_ = secondHalfRunEdrFileName_;
                CommandLine secondHalfRunCaller;
                secondHalfRunCaller.append("mdrun");
                secondHalfRunCaller.append("-noappend");
                secondHalfRunCaller.addOption("-cpi", firstHalfRunCheckpointFileName_);
                ASSERT_EQ(0, runner_.callMdrun(secondHalfRunCaller));
                // Cope with how -noappend works
                secondHalfRunTrajectoryFileName_ += ".part0002.trr";
                secondHalfRunEdrFileName_ += ".part0002.edr";
            }

            simulationsWereRun_ = true;
        }
        //! Check the trajectories and energies of the rerun agree with the normal run within \c tolerance.
        void checkResults(const std::vector<std::string>     &requiredEnergies,
                          const TrajectoryFrameMatchSettings &matchSettings,
                          FloatingPointTolerance              tolerance)
        {
            ASSERT_TRUE(simulationsWereRun_);

            // Prepare object to read the energy files and return
            // pairs of frames to compare
            FramePairManager<EnergyFrameReader, EnergyFrame>
                energyManager(openEnergyFileToReadFields(fullRunEdrFileName_, requiredEnergies),
                              openEnergyFileToReadFields(firstHalfRunEdrFileName_, requiredEnergies),
                              openEnergyFileToReadFields(secondHalfRunEdrFileName_, requiredEnergies));
            while (energyManager.shouldContinueComparing())
            {
                auto frames = energyManager.getPair();
                SCOPED_TRACE("Comparing frames from two runs '" + frames.first.getFrameName() + "' and '" + frames.second.getFrameName() + "'");
                compareFrames(frames.first, frames.second, tolerance);
            }

            // Prepare object to read the trajectory files and return
            // pairs of frames to compare
            FramePairManager<TrajectoryFrameReader, TrajectoryFrame>
            trajectoryManager(TrajectoryFrameReaderPtr(new TrajectoryFrameReader(fullRunTrajectoryFileName_)),
                              TrajectoryFrameReaderPtr(new TrajectoryFrameReader(firstHalfRunTrajectoryFileName_)),
                              TrajectoryFrameReaderPtr(new TrajectoryFrameReader(secondHalfRunTrajectoryFileName_)));
            while (trajectoryManager.shouldContinueComparing())
            {
                auto frames = trajectoryManager.getPair();
                SCOPED_TRACE("Comparing frames from two runs '" + frames.first.getFrameName() + "' and '" + frames.second.getFrameName() + "'");
                compareFrames(frames.first, frames.second, matchSettings, tolerance);
            }
        }
    private:
        //! Name of file used in testing.
        std::string fullRunTrajectoryFileName_;
        //! Name of file used in testing.
        std::string fullRunEdrFileName_;
        //! Name of file used in testing.
        std::string firstHalfRunTrajectoryFileName_;
        //! Name of file used in testing.
        std::string firstHalfRunEdrFileName_;
        //! Name of file used in testing.
        std::string firstHalfRunCheckpointFileName_;
        //! Name of file used in testing.
        std::string secondHalfRunTrajectoryFileName_;
        //! Name of file used in testing.
        std::string secondHalfRunEdrFileName_;
        //! True only when simulations have been run and results can be tested.
        bool        simulationsWereRun_;
};

//! Name the test case usefully
typedef MdrunContinuationTestFixture MdrunContinuationTest;

TEST_P(MdrunContinuationTest, IsExactWithinTolerances)
{
    auto        params = GetParam();
    SCOPED_TRACE(params);

    CommandLine caller;
    caller.addOption("-maxwarn", params.maxWarningsTolerated);
    runGrompp(caller, params.simulationName, params.integrator, params.temperatureCoupling, params.pressureCoupling);
    runTest();
    checkResults(params.energyNames, params.matchSettings, params.tolerance);
}

/*! \brief Specifies the expected tolerance for this test
 *
 * NOTE This is an arbitrary-but-sufficient choice for the
 * tolerance. */
FloatingPointTolerance testTolerance_g = relativeToleranceAsPrecisionDependentUlp(1.0, 160, 1.2e5);

//! Names the energy fields required to compare equal in order to pass this test
std::vector<std::string> namesOfRequiredEnergyFields_g =
    {{interaction_function[F_EPOT].longname,
      interaction_function[F_EKIN].longname,
      interaction_function[F_ETOT].longname,
      interaction_function[F_TEMP].longname,
      interaction_function[F_PRES].longname }};

//! Specify how trajectory frame matching must work.
TrajectoryFrameMatchSettings matchSettings_g {
    true, true, true, true, true, true
};

/*! \brief Specifies the cases tested.
 *
 * Some combinations are not implemented, these are noted here. MTTK
 * pressure coupling is only implemented for one sub-case.
 *
 * \todo In all cases that are not implemented, grompp currently gives
 * a fatal error, so there's no way to check that behaviour, nor
 * automatically flag (via the future change of behaviour of the test)
 * that we should start testing it, when it becomes implemented. */
ExactContinuationTestParams mdrunParams_g[] = {
    { "argon12", "md",    "no",          "no",                namesOfRequiredEnergyFields_g, matchSettings_g, testTolerance_g, 0 },
    { "argon12", "md",    "no",          "berendsen",         namesOfRequiredEnergyFields_g, matchSettings_g, testTolerance_g, 0 },
    { "argon12", "md",    "no",          "parrinello-rahman", namesOfRequiredEnergyFields_g, matchSettings_g, testTolerance_g, 0 },
    { "argon12", "md",    "berendsen",   "no",                namesOfRequiredEnergyFields_g, matchSettings_g, testTolerance_g, 0 },
    { "argon12", "md",    "berendsen",   "berendsen",         namesOfRequiredEnergyFields_g, matchSettings_g, testTolerance_g, 0 },
    { "argon12", "md",    "berendsen",   "parrinello-rahman", namesOfRequiredEnergyFields_g, matchSettings_g, testTolerance_g, 0 },
    { "argon12", "md",    "vrescale",    "no",                namesOfRequiredEnergyFields_g, matchSettings_g, testTolerance_g, 0 },
    { "argon12", "md",    "vrescale",    "berendsen",         namesOfRequiredEnergyFields_g, matchSettings_g, testTolerance_g, 0 },
    { "argon12", "md",    "vrescale",    "parrinello-rahman", namesOfRequiredEnergyFields_g, matchSettings_g, testTolerance_g, 0 },
    { "argon12", "md",    "nose-hoover", "no",                namesOfRequiredEnergyFields_g, matchSettings_g, testTolerance_g, 0 },
    { "argon12", "md",    "nose-hoover", "berendsen",         namesOfRequiredEnergyFields_g, matchSettings_g, testTolerance_g, 1 },
    { "argon12", "md",    "nose-hoover", "parrinello-rahman", namesOfRequiredEnergyFields_g, matchSettings_g, testTolerance_g, 0 },
    { "argon12", "md-vv", "no",          "no",                namesOfRequiredEnergyFields_g, matchSettings_g, testTolerance_g, 0 },
    { "argon12", "md-vv", "no",          "berendsen",         namesOfRequiredEnergyFields_g, matchSettings_g, testTolerance_g, 0 },
    // not implemented
    // { "argon12", "md-vv", "no",          "parrinello-rahman", namesOfRequiredEnergyFields_g, matchSettings_g, testTolerance_g, 0 },
    { "argon12", "md-vv", "berendsen",   "no",                namesOfRequiredEnergyFields_g, matchSettings_g, testTolerance_g, 0 },
    { "argon12", "md-vv", "berendsen",   "berendsen",         namesOfRequiredEnergyFields_g, matchSettings_g, testTolerance_g, 0 },
    // not implemented
    // { "argon12", "md-vv", "berendsen",   "parrinello-rahman", namesOfRequiredEnergyFields_g, matchSettings_g, testTolerance_g, 0 },
    { "argon12", "md-vv", "vrescale",    "no",                namesOfRequiredEnergyFields_g, matchSettings_g, testTolerance_g, 0 },
    { "argon12", "md-vv", "vrescale",    "berendsen",         namesOfRequiredEnergyFields_g, matchSettings_g, testTolerance_g, 0 },
    // not implemented
    // { "argon12", "md-vv", "vrescale",    "parrinello-rahman", namesOfRequiredEnergyFields_g, matchSettings_g, testTolerance_g, 0 },
    { "argon12", "md-vv", "nose-hoover", "no",                namesOfRequiredEnergyFields_g, matchSettings_g, testTolerance_g, 0 },
    // not implemented
    // { "argon12", "md-vv", "nose-hoover", "berendsen",         namesOfRequiredEnergyFields_g, matchSettings_g, testTolerance_g, 0 },
    { "argon12", "md-vv", "nose-hoover", "mttk",              namesOfRequiredEnergyFields_g, matchSettings_g, testTolerance_g, 0 },
    { "argon12", "bd",    "no",          "no",                namesOfRequiredEnergyFields_g, matchSettings_g, testTolerance_g, 0 },
    { "argon12", "bd",    "no",          "berendsen",         namesOfRequiredEnergyFields_g, matchSettings_g, testTolerance_g, 0 },
    { "argon12", "bd",    "no",          "parrinello-rahman", namesOfRequiredEnergyFields_g, matchSettings_g, testTolerance_g, 0 },
    { "argon12", "sd",    "no",          "no",                namesOfRequiredEnergyFields_g, matchSettings_g, testTolerance_g, 0 },
    { "argon12", "sd",    "no",          "berendsen",         namesOfRequiredEnergyFields_g, matchSettings_g, testTolerance_g, 0 },
    { "argon12", "sd",    "no",          "parrinello-rahman", namesOfRequiredEnergyFields_g, matchSettings_g, testTolerance_g, 0 },
};

INSTANTIATE_TEST_CASE_P(ForArgon, MdrunContinuationTest, ::testing::ValuesIn(mdrunParams_g));

} // namespace
} // namespace
} // namespace
