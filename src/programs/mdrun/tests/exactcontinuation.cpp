/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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

#include <string>
#include <tuple>

#include <gtest/gtest.h>

#include "gromacs/topology/idef.h"
#include "gromacs/trajectory/energyframe.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/mpitest.h"
#include "testutils/testasserts.h"

#include "energycomparison.h"
#include "energyreader.h"
#include "moduletest.h"
#include "simulationdatabase.h"

namespace gmx
{
namespace test
{
namespace
{


/*! \brief Manages returning a pair of energy frames from a normal run
 * and the same run in two halves, such that they are meaningful to
 * compare.
 *
 * This class has a similar interface with one in
 * exactcontinuation.cpp, but not enough to warrant extracting an
 * interface class. */
template <class FrameReader, class Frame>
class ContinuationFramePairManager
{
    public:
        //! Convenience typedef
        typedef std::unique_ptr<FrameReader> FrameReaderPtr;
        //! Constructor
        ContinuationFramePairManager(FrameReaderPtr full,
                                     FrameReaderPtr firstPart,
                                     FrameReaderPtr secondPart) :
            full_(std::move(full)),
            firstPart_(std::move(firstPart)),
            secondPart_(std::move(secondPart)),
            isFirstPart_(true)
        {}
        /*! \brief Probe for a pair of valid frames, and return true if both are found.
         *
         * Give a test failure if exactly one frame is found, because
         * that file is longer than the other one, and this is not
         * expected behaviour.
         *
         * \todo This would be straightforward if velocity Verlet
         * behaved like other integrators. */
        bool shouldContinueComparing()
        {
            if (full_->readNextFrame())
            {
                if (isFirstPart_)
                {
                    if (firstPart_->readNextFrame())
                    {
                        // Two valid next frames exist, so we should continue comparing.
                        return true;
                    }
                    else
                    {
                        // First half ran out of frames, move on to the second half
                        isFirstPart_  = false;
                        if (secondPart_->readNextFrame())
                        {
                            // Skip a second-half frame so the one we will
                            // read can compare with the next full-run
                            // frames.
                            secondPart_->frame();
                            if (secondPart_->readNextFrame())
                            {
                                // Two valid next frames exist, so we should continue comparing.
                                return true;
                            }
                            else
                            {
                                ADD_FAILURE() << "Second-half energy file had no (new) frames";
                            }
                        }
                        else
                        {
                            ADD_FAILURE() << "Second-half energy file had no frames";
                        }
                    }
                }
                else
                {
                    if (secondPart_->readNextFrame())
                    {
                        // Two valid next frames exist, so we should continue comparing.
                        return true;
                    }
                    else
                    {
                        ADD_FAILURE() << "Full run energy file had at least one more frame than two-part run energy file";
                    }
                }
            }
            else
            {
                if (isFirstPart_)
                {
                    ADD_FAILURE() << "Full-run energy file ran out of frames before the first half of the two-part run completed";
                }
                else
                {
                    if (secondPart_->readNextFrame())
                    {
                        ADD_FAILURE() << "Two-part run energy file had at least one more frame than full-run energy file";
                    }
                    else
                    {
                        // Both files ran out of frames at the same time, which is the expected behaviour.
                    }
                }
            }
            // At least one file is out of frames, so should not continue comparing.
            return false;
        }
    public:
        //! Compare all possible pairs of frames using \c compareTwoFrames.
        void compareAllFramePairs(std::function<void(const Frame &, const Frame &)> compareTwoFrames)
        {
            while (shouldContinueComparing())
            {
                auto firstFrame  = full_->frame();
                auto secondFrame = isFirstPart_ ? firstPart_->frame() : secondPart_->frame();
                SCOPED_TRACE("Comparing frames from two runs '" + firstFrame.frameName() + "' and '" + secondFrame.frameName() + "'");
                compareTwoFrames(firstFrame, secondFrame);
            }

        }

    private:
        EnergyFrameReaderPtr full_;
        EnergyFrameReaderPtr firstPart_;
        EnergyFrameReaderPtr secondPart_;
        bool                 isFirstPart_;
};

/*! \brief Run grompp for a normal mdrun, the same mdrun stopping half
 * way, doing a continuation, and compare the results. */
void executeTest(TestFileManager        *fileManager,
                 SimulationRunner       *runner,
                 const std::string      &simulationName,
                 int                     maxWarningsTolerated,
                 const MdpFieldValues   &mdpFieldValues,
                 const EnergyTolerances &energiesToMatch)
{
    int numRanksAvailable = getNumberOfTestMpiRanks();
    if (!isNumberOfPpRanksSupported(simulationName, numRanksAvailable))
    {
        fprintf(stdout, "Test system '%s' cannot run with %d ranks.\n"
                "The supported numbers are: %s\n",
                simulationName.c_str(), numRanksAvailable,
                reportNumbersOfPpRanksSupported(simulationName).c_str());
        return;
    }

    // prepare some names for files to use with the two mdrun calls
    std::string fullRunEdrFileName              = fileManager->getTemporaryFilePath("full.edr");
    std::string firstHalfRunEdrFileName         = fileManager->getTemporaryFilePath("firsthalf.edr");
    std::string firstHalfRunCheckpointFileName  = fileManager->getTemporaryFilePath("firsthalf.cpt");
    std::string secondHalfRunEdrFileName        = fileManager->getTemporaryFilePath("secondhalf");

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

    // do a normal mdrun
    {
        runner->edrFileName_ = fullRunEdrFileName;
        CommandLine fullRunCaller;
        fullRunCaller.append("mdrun");
        ASSERT_EQ(0, runner->callMdrun(fullRunCaller));
    }

    // do a repeat of the first half of the same mdrun
    {
        runner->edrFileName_ = firstHalfRunEdrFileName;
        CommandLine firstHalfRunCaller;
        firstHalfRunCaller.append("mdrun");
        // TODO This happens to be half of the default value for
        // nsteps in simulationdatabase.cpp, but we should find a
        // better way to do this.
        firstHalfRunCaller.addOption("-nsteps", 8);
        firstHalfRunCaller.addOption("-cpo", firstHalfRunCheckpointFileName);
        ASSERT_EQ(0, runner->callMdrun(firstHalfRunCaller));
    }

    // do a continuation from the first half of that same mdrun
    {
        runner->edrFileName_ = secondHalfRunEdrFileName;
        CommandLine secondHalfRunCaller;
        secondHalfRunCaller.append("mdrun");
        // TODO We could test with appending but it would need a
        // different implementation.
        secondHalfRunCaller.append("-noappend");
        secondHalfRunCaller.addOption("-cpi", firstHalfRunCheckpointFileName);
        ASSERT_EQ(0, runner->callMdrun(secondHalfRunCaller));
        // Cope with how -noappend works
        secondHalfRunEdrFileName += ".part0002.edr";
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
    ContinuationFramePairManager<EnergyFrameReader, EnergyFrame>
         energyManager(openEnergyFileToReadFields(fullRunEdrFileName, namesOfEnergiesToMatch),
                  openEnergyFileToReadFields(firstHalfRunEdrFileName, namesOfEnergiesToMatch),
                  openEnergyFileToReadFields(secondHalfRunEdrFileName, namesOfEnergiesToMatch));
    // Compare the energy frames.
    energyManager.compareAllFramePairs(energyComparator);
}

/*! \brief Test fixture for mdrun exact continuations
 *
 * This test ensures mdrun can run a simulation, writing a trajectory
 * and matching energies, and reproduce to within a tolerance the same
 * energies from runs that stopped half way, and restarted from the
 * checkpoint.
 *
 * \todo Is there value in testing with mdrun -reprod? As well as
 * without it?
 *
 * \todo Add FEP case. */
class MdrunContinuationIsExact : public MdrunTestFixture,
                                 public ::testing::WithParamInterface <
                                 std::tuple < std::string, std::string, std::string, std::string >>
{
    public:
        //! Constructor
        MdrunContinuationIsExact() {}
};

/* Listing all of these is tedious, but there's no other way to get a
 * usefully descriptive string into the test-case name, so that when
 * one breaks we can find out which one is broken without referring to
 * this source code file.
 *
 * NOTE The choices for the tolerances are arbitrary but sufficient
 * for comparison of runs to work on different hardware, and kinds and
 * degrees parallelism. */

TEST_P(MdrunContinuationIsExact, WithinTolerances)
{
    auto params              = GetParam();
    auto simulationName      = std::get<0>(params);
    auto integrator          = std::get<1>(params);
    auto temperatureCoupling = std::get<2>(params);
    auto pressureCoupling    = std::get<3>(params);
    SCOPED_TRACE(formatString("Comparing normal and two-part run of simulation '%s' "
                              "with integrator '%s'",
                              simulationName.c_str(), integrator.c_str()));

    auto mdpFieldValues = prepareMdpFieldValues(simulationName.c_str(),
                                                integrator.c_str(),
                                                temperatureCoupling.c_str(),
                                                pressureCoupling.c_str());

    EnergyTolerances energiesToMatch
    {{
         {
             interaction_function[F_EPOT].longname,
             relativeToleranceAsPrecisionDependentUlp(10.0, 1, 1)
         },
     }};

    int numWarningsToTolerate = 1;
    executeTest(&fileManager_, &runner_,
                simulationName,
                numWarningsToTolerate, mdpFieldValues,
                energiesToMatch);
}

INSTANTIATE_TEST_CASE_P(NormalIntegrators, MdrunContinuationIsExact,
                            ::testing::Combine(::testing::Values("argon12", "spc5", "alanine_vsite_vacuo"),
                                                   ::testing::Values("md", "md-vv", "bd", "sd"),
                                                   ::testing::Values("no"),
                                                   ::testing::Values("no")));

// TODO Brownian Dynamics diverges by a few ULP at step 16
INSTANTIATE_TEST_CASE_P(NormalIntegratorsWithFEP, MdrunContinuationIsExact,
                            ::testing::Combine(::testing::Values("nonanol_vacuo"),
                                                   ::testing::Values("md", "md-vv", /*"bd",*/ "sd"),
                                                   ::testing::Values("no"),
                                                   ::testing::Values("no")));

INSTANTIATE_TEST_CASE_P(NormalNVT, MdrunContinuationIsExact,
                            ::testing::Combine(::testing::Values("argon12"),
                                                   ::testing::Values("md", "md-vv"),
                                                   ::testing::Values("berendsen", "v-rescale", "nose-hoover"),
                                                   ::testing::Values("no")));

INSTANTIATE_TEST_CASE_P(LeapfrogNPH, MdrunContinuationIsExact,
                            ::testing::Combine(::testing::Values("argon12"),
                                                   ::testing::Values("md"),
                                                   ::testing::Values("no"),
                                                   ::testing::Values("berendsen", "parrinello-rahman")));

INSTANTIATE_TEST_CASE_P(LeapfrogNPT, MdrunContinuationIsExact,
                            ::testing::Combine(::testing::Values("argon12"),
                                                   ::testing::Values("md"),
                                                   ::testing::Values("berendsen", "v-rescale", "nose-hoover"),
                                                   ::testing::Values("berendsen", "parrinello-rahman")));

INSTANTIATE_TEST_CASE_P(VelocityVerletNPH, MdrunContinuationIsExact,
                            ::testing::Combine(::testing::Values("argon12"),
                                                   ::testing::Values("md-vv"),
                                                   ::testing::Values("no"),
                                                   ::testing::Values("berendsen")));

INSTANTIATE_TEST_CASE_P(VelocityVerletNPT, MdrunContinuationIsExact,
                            ::testing::Combine(::testing::Values("argon12"),
                                                   ::testing::Values("md-vv"),
                                                   ::testing::Values("v-rescale"),
                                                   ::testing::Values("berendsen")));

INSTANTIATE_TEST_CASE_P(MTTK, MdrunContinuationIsExact,
                            ::testing::Combine(::testing::Values("argon12"),
                                                   ::testing::Values("md-vv"),
                                                   ::testing::Values("nose-hoover"),
                                                   ::testing::Values("mttk")));

} // namespace
} // namespace
} // namespace
