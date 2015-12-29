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
#include <vector>

#include <tuple>
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
 * Reproducing the same energies is currently only meaningful for
 * integration without thermostats or barostats, however the present
 * form of the test infrastructure has in-principle support for such,
 * if that is ever needed.
 *
 * \todo Test for valid restart from non-nstlist steps - needs wider
 * tolerance for the first step after restart, because the pair list
 * will differ. */
class MdrunRerunTestFixtureBase : public MdrunComparisonFixture
{
    public:
        /*! \brief Describe the energy-file fields that we wish to
         * compare between the normal run and its rerun.
         *
         * Derived classes must specialize this to return the set of
         * fields appropriate for the functionality that they test.
         *
         * We should not compare pressure, because with constraints
         * the non-search steps need a much larger tolerance, and per
         * Redmine 1868 we should stop computing pressure in reruns
         * anyway.
         *
         * Similarly, per 1868, in the present implementation the kinetic
         * energy quantities are not reproducible, either. */
        virtual std::vector<std::string> namesOfRequiredFields() const = 0;
        // Import convenience wrapper from parent class.
        using MdrunComparisonFixture::runTest;
        //! Run a normal mdrun, its rerun, and compare the results
        virtual void runTest(const CommandLine     &gromppCallerRef,
                             const char            *simulationName,
                             const char            *integrator,
                             const char            *tcoupl,
                             const char            *pcoupl,
                             FloatingPointTolerance tolerance)
        {
            /* Note that the database lookups will throw std::out_of_range if
             * the simulationName key is not found */
            runner_.useTopGroAndNdxFromDatabase(simulationName);
            auto mdpFieldValues = prepareMdpFieldValues(simulationName);
            prepareMdpFile(mdpFieldValues, integrator, tcoupl, pcoupl);
            EXPECT_EQ(0, runner_.callGrompp(gromppCallerRef));

            // prepare some names for files to use with the two mdrun calls
            std::string normalRunTrajectoryFileName = fileManager_.getTemporaryFilePath("normal.trr");
            std::string normalRunEdrFileName        = fileManager_.getTemporaryFilePath("normal.edr");
            std::string rerunTrajectoryFileName     = fileManager_.getTemporaryFilePath("rerun.trr");
            std::string rerunEdrFileName            = fileManager_.getTemporaryFilePath("rerun.edr");

            // do a normal mdrun
            {
                runner_.fullPrecisionTrajectoryFileName_ = normalRunTrajectoryFileName;
                runner_.edrFileName_                     = normalRunEdrFileName;
                CommandLine normalRunCaller;
                normalRunCaller.append("mdrun");
                ASSERT_EQ(0, runner_.callMdrun(normalRunCaller));
            }

            // do a rerun on the .trr just produced
            {
                runner_.fullPrecisionTrajectoryFileName_ = rerunTrajectoryFileName;
                runner_.edrFileName_                     = rerunEdrFileName;
                CommandLine rerunCaller;
                rerunCaller.append("mdrun");
                rerunCaller.addOption("-rerun", normalRunTrajectoryFileName);
                ASSERT_EQ(0, runner_.callMdrun(rerunCaller));
            }

            // Prepare object to read the energy files and return
            // pairs of frames to compare
            auto requiredNames = namesOfRequiredFields();
            FramePairManager<EnergyFrameReader, EnergyFrame>
                 energyManager(openEnergyFileToReadFields(normalRunEdrFileName, requiredNames),
                          openEnergyFileToReadFields(rerunEdrFileName, requiredNames));
            while (energyManager.shouldContinueComparing())
            {
                compareFrames(energyManager.getPair(), tolerance);
            }

            // Prepare object to read the trajectory files and return
            // pairs of frames to compare
            FramePairManager<TrajectoryFrameReader, TrajectoryFrame>
            trajectoryManager(TrajectoryFrameReaderPtr(new TrajectoryFrameReader(normalRunTrajectoryFileName)),
                              TrajectoryFrameReaderPtr(new TrajectoryFrameReader(rerunTrajectoryFileName)));
            while (trajectoryManager.shouldContinueComparing())
            {
                compareFrames(trajectoryManager.getPair(), tolerance);
            }
        }
};

class RerunReproducesNormalMdrunFor : public MdrunRerunTestFixtureBase
{
    public:
        //! Compare potential energy only
        virtual std::vector<std::string> namesOfRequiredFields() const
        {
            return std::vector<std::string> {{interaction_function[F_EPOT].longname}};
        }
};

/* Listing all of these is tedious, but there's no other way to get a
 * usefully descriptive string into the test-case name, so that when
 * one breaks, we can find out which one is broken without referring
 * to this source code file.
 *
 * NOTE The choices for tolerance are arbitrary but sufficient.  Rerun
 * does pair search every frame, so it cannot in general exactly
 * reproduce quantities from a normal run, because the accumulation
 * order differs. (Nor does it reproduce pair-search frames exactly,
 * either). */

TEST_F(RerunReproducesNormalMdrunFor, Argon12_Leapfrog)
{
    runTest("argon12", "md", "no", "no", relativeToleranceAsPrecisionDependentUlp(1.0, 1, 1));
}

TEST_F(RerunReproducesNormalMdrunFor, Argon12_VelocityVerlet)
{
    runTest("argon12", "md-vv", "no", "no", relativeToleranceAsPrecisionDependentUlp(1.0, 1, 1));
}

TEST_F(RerunReproducesNormalMdrunFor, Argon12_BrownianDynamics)
{
    runTest("argon12", "bd", "no", "no", relativeToleranceAsPrecisionDependentUlp(1.0, 1, 1));
}

TEST_F(RerunReproducesNormalMdrunFor, Argon12_StochasticDynamics)
{
    runTest("argon12", "sd", "no", "no", relativeToleranceAsPrecisionDependentUlp(1.0, 1, 1));
}

TEST_F(RerunReproducesNormalMdrunFor, Water5_Leapfrog)
{
    runTest("spc5", "md", "no", "no", relativeToleranceAsPrecisionDependentUlp(1.0, 40, 30));
}

TEST_F(RerunReproducesNormalMdrunFor, Water5_VelocityVerlet)
{
    runTest("spc5", "md-vv", "no", "no", relativeToleranceAsPrecisionDependentUlp(1.0, 20, 20));
}

// TODO Some forces are completely wrong, so this test is disabled
TEST_F(RerunReproducesNormalMdrunFor, DISABLED_Water5_BrownianDynamics)
{
    runTest("spc5", "bd", "no", "no", relativeToleranceAsPrecisionDependentUlp(1.0, 40, 50));
}

TEST_F(RerunReproducesNormalMdrunFor, Water5_StochasticDynamics)
{
    runTest("spc5", "sd", "no", "no", relativeToleranceAsPrecisionDependentUlp(1.0, 30, 10));
}

TEST_F(RerunReproducesNormalMdrunFor, AlanineVsiteVacuo_Leapfrog)
{
    runTest("alanine_vsite_vacuo", "md", "no", "no", relativeToleranceAsPrecisionDependentUlp(1.0, 1100, 600));
}

TEST_F(RerunReproducesNormalMdrunFor, AlanineVsiteVacuo_VelocityVerlet)
{
    runTest("alanine_vsite_vacuo", "md-vv", "no", "no", relativeToleranceAsPrecisionDependentUlp(1.0, 60, 50));
}

TEST_F(RerunReproducesNormalMdrunFor, AlanineVsiteVacuo_BrownianDynamics)
{
    runTest("alanine_vsite_vacuo", "bd", "no", "no", relativeToleranceAsPrecisionDependentUlp(1.0, 140, 150));
}

TEST_F(RerunReproducesNormalMdrunFor, AlanineVsiteVacuo_StochasticDynamics)
{
    runTest("alanine_vsite_vacuo", "sd", "no", "no", relativeToleranceAsPrecisionDependentUlp(1.0, 70, 100));
}

class RerunReproducesFreeEnergyMdrunFor : public MdrunRerunTestFixtureBase
{
    public:
        /*! \brief Only compare quantities like the potential energy
         *
         * \todo Perhaps we can compare kinetic energies at some
         * future time (Redmine 1868), but the code is too complicated
         * to fix right now. */
        virtual std::vector<std::string> namesOfRequiredFields() const
        {
            return std::vector<std::string> {{interaction_function[F_EPOT].longname},
                                             {interaction_function[F_DVDL].longname},
                                             {interaction_function[F_DVDL_VDW].longname}};
        }
};

TEST_F(RerunReproducesFreeEnergyMdrunFor, NonanolVacuo_LeapFrog)
{
    CommandLine caller;
    caller.addOption("-maxwarn", 1);
    runTest(caller, "nonanol_vacuo", "md", "no", "no", relativeToleranceAsPrecisionDependentUlp(1.0, 10, 10));
}

TEST_F(RerunReproducesFreeEnergyMdrunFor, NonanolVacuo_VelocityVerlet)
{
    runTest("nonanol_vacuo", "md-vv", "no", "no", relativeToleranceAsPrecisionDependentUlp(1.0, 10, 20));
}

TEST_F(RerunReproducesFreeEnergyMdrunFor, NonanolVacuo_StochasticDynamics)
{
    runTest("nonanol_vacuo", "sd", "no", "no", relativeToleranceAsPrecisionDependentUlp(1.0, 10, 40));
}

} // namespace
} // namespace
} // namespace
