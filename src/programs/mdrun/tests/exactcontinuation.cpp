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
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"

#include "energyreader.h"
#include "moduletest.h"

namespace gmx
{
namespace test
{
namespace
{

/* Do a simulation on argon
 * - writing frames from different kinds of steps: starting, ending, intermediate NS, intermediate non-NS
 * - with other steps between frame-writing steps
 * - with enough buffer that the rerun will compute the same potential energy even though it does NS every frame
 * - without constraints
 * whose restarted run will match with a tight tolerance.
 */
const char *argonSimulationMdpFileFormatString {
    "rcoulomb = 1.0\n"
    "rvdw = 1.0\n"
    "rlist = -1\n"
    "bd-fric = 100\n"
    "cutoff-scheme = Verlet\n"
    "verlet-buffer-tolerance = 0.000001\n"
    "nsteps = 20\n"
    "nstcalcenergy = -1\n"
    "nstenergy = 5\n"
    "nstlist = 10\n"
    "nstxout = 5\n"
    "nstvout = 5\n"
    "integrator = %s\n"
    "ld-seed = 234262\n"
    "tcoupl = %s\n"
    "ref-t = 298.0\n"
    "tau-t = 1.0\n"
    "tc-grps = System\n"
    "pcoupl = %s\n"
    "ref-p = 1\n"
    "compressibility = 5e-5\n"
};

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
class MdrunContinuationIsExact : public MdrunTestFixture
{
    public:
        //! Constructor
        MdrunContinuationIsExact() : testTolerance_(0, 0, 0, 0, true) {}
        /*! \brief Run a normal mdrun, the same mdrun stopping half
         * way, doing a continuation, and compare the results. */
        void runTest(bool isVelocityVerlet)
        {
            EXPECT_EQ(0, runner_.callGrompp());

            // prepare some names for files to use with the two mdrun calls
            std::string fullRunEdrFileName              = fileManager_.getTemporaryFilePath("full.edr");
            std::string firstHalfRunEdrFileName         = fileManager_.getTemporaryFilePath("firsthalf.edr");
            std::string firstHalfRunCheckpointFileName  = fileManager_.getTemporaryFilePath("firsthalf.cpt");
            std::string secondHalfRunEdrFileName        = fileManager_.getTemporaryFilePath("secondhalf");

            // do a normal mdrun
            {
                runner_.edrFileName_ = fullRunEdrFileName;
                CommandLine fullRunCaller;
                fullRunCaller.append("mdrun");
                ASSERT_EQ(0, runner_.callMdrun(fullRunCaller));
            }

            // do a repeat of the first half of the same mdrun
            {
                runner_.edrFileName_ = firstHalfRunEdrFileName;
                CommandLine firstHalfRunCaller;
                firstHalfRunCaller.append("mdrun");
                firstHalfRunCaller.addOption("-nsteps", 10);
                firstHalfRunCaller.addOption("-cpo", firstHalfRunCheckpointFileName);
                ASSERT_EQ(0, runner_.callMdrun(firstHalfRunCaller));
            }

            // do a continuation from the first half of that same mdrun
            {
                runner_.edrFileName_ = secondHalfRunEdrFileName;
                CommandLine secondHalfRunCaller;
                secondHalfRunCaller.append("mdrun");
                secondHalfRunCaller.append("-noappend");
                secondHalfRunCaller.addOption("-cpi", firstHalfRunCheckpointFileName);
                ASSERT_EQ(0, runner_.callMdrun(secondHalfRunCaller));
                // Cope with how -noappend works
                secondHalfRunEdrFileName += ".part0002.edr";
            }

            // Prepare objects to read the energy files
            EnergyFileReader  fullRunFile(fullRunEdrFileName);
            EnergyFileReader  firstHalfRunFile(firstHalfRunEdrFileName);
            EnergyFileReader  secondHalfRunFile(secondHalfRunEdrFileName);

            EnergyFrameReader fullRun(fullRunFile.openToReadFields(namesOfRequiredFields_));
            EnergyFrameReader firstHalfRun(firstHalfRunFile.openToReadFields(namesOfRequiredFields_));
            EnergyFrameReader secondHalfRun(secondHalfRunFile.openToReadFields(namesOfRequiredFields_));

            bool              isFirstHalf = true;
            EnergyFrameInfo   oldFullRunFrame;
            do
            {
                EnergyFrameInfo fullRunFrame;
                EnergyFrameInfo currentHalfRunFrame = isFirstHalf ? firstHalfRun.readNextFrame() : secondHalfRun.readNextFrame();
                if (isFirstHalf && !currentHalfRunFrame)
                {
                    fprintf(stderr, "run out of first half\n");
                    // Presumably we ran out of frames for the first half
                    isFirstHalf         = false;
                    currentHalfRunFrame = secondHalfRun.readNextFrame();
                    if (isVelocityVerlet)
                    {
                        // Get a full-run frame to compare with the new half-run frame.
                        fullRunFrame = fullRun.readNextFrame();
                    }
                    else
                    {
                        // The first frame in the second half is the
                        // same as the last of the first half, so we
                        // can compare it with the old frame from the
                        // full run.
                        fullRunFrame = oldFullRunFrame;
                    }
                }
                else
                {
                    // Get a full-run frame to compare with the new half-run frame.
                    fullRunFrame = fullRun.readNextFrame();
                }
                if (fullRunFrame && currentHalfRunFrame)
                {
                    for (auto const &name : namesOfRequiredFields_)
                    {
                        EXPECT_REAL_EQ_TOL(fullRunFrame.getValue(name), currentHalfRunFrame.getValue(name), testTolerance_)
                        << name << " didn't match between full run " << fullRunFrame.getFrameName() << " and restarted run " << currentHalfRunFrame.getFrameName();
                    }
                }
                else
                {
                    // At least one file is out of frames. Report if
                    // the either run had a frame when the other did
                    // not.
                    if (fullRunFrame != currentHalfRunFrame)
                    {
                        EXPECT_FALSE(currentHalfRunFrame) << "the energy files from the two-part run had at least one extra frame";
                        EXPECT_FALSE(fullRunFrame) << "the energy file from the full run had at least one extra frame";
                    }
                    break;
                }
                oldFullRunFrame = fullRunFrame;
            }
            while (true);
        }
        //! Help to prepare simulations (since we can't have a constructor with arguments in a test fixture)
        void prepareArgonSimulation(const char *integrator,
                                    const char *thermostat,
                                    const char *barostat)
        {
            runner_.useStringAsMdpFile(formatString(argonSimulationMdpFileFormatString, integrator, thermostat, barostat));

            // make the tpr file
            runner_.useTopGroAndNdxFromDatabase("argon");

            // Describe the energy-file fields that we wish to compare between
            // the two runs.
            namesOfRequiredFields_ = {{interaction_function[F_EPOT].longname,
                                       interaction_function[F_EKIN].longname,
                                       interaction_function[F_ETOT].longname,
                                       interaction_function[F_TEMP].longname,
                                       interaction_function[F_PRES].longname }};

            /* NOTE This is an arbitrary-but-sufficient choice for the
             * tolerance. */
            testTolerance_ = relativeToleranceAsPrecisionDependentUlp(1.0, 160, 1.2e5);
        }
        //! Names the energy fields required to compare equal in order to pass this test
        std::vector<std::string> namesOfRequiredFields_;
        //! Specifies the expected tolerance for this test
        FloatingPointTolerance   testTolerance_;
};

//! Helper function
bool isVelocityVerlet(const char *integrator)
{
    return (0 == std::strcmp(integrator, "md-vv"));
}

/* Listing all of these is tedious, but there's no other way to get a
 * usefully descriptive string into the test-case name, so that when
 * one breaks we can find out which one is broken without referring to
 * this source code file. */

typedef MdrunContinuationIsExact ArgonContinuationIsExactWith;

TEST_F(ArgonContinuationIsExactWith, Leapfrog_NVE)
{
    const char *integrator = "md";
    const char *thermostat = "no";
    const char *barostat   = "no";

    prepareArgonSimulation(integrator, thermostat, barostat);
    runTest(isVelocityVerlet(integrator));
}

TEST_F(ArgonContinuationIsExactWith, Leapfrog_NoThermostat_BerendsenBarostat)
{
    const char *integrator = "md";
    const char *thermostat = "no";
    const char *barostat   = "berendsen";

    prepareArgonSimulation(integrator, thermostat, barostat);
    runTest(isVelocityVerlet(integrator));
}

TEST_F(ArgonContinuationIsExactWith, Leapfrog_NoThermostat_ParrinelloRahmanBarostat)
{
    const char *integrator = "md";
    const char *thermostat = "no";
    const char *barostat   = "parrinello-rahman";

    prepareArgonSimulation(integrator, thermostat, barostat);
    runTest(isVelocityVerlet(integrator));
}

TEST_F(ArgonContinuationIsExactWith, Leapfrog_BerendsenThermostat_NoBarostat)
{
    const char *integrator = "md";
    const char *thermostat = "berendsen";
    const char *barostat   = "no";

    prepareArgonSimulation(integrator, thermostat, barostat);
    runTest(isVelocityVerlet(integrator));
}

TEST_F(ArgonContinuationIsExactWith, Leapfrog_BerendsenThermostat_BerendsenBarostat)
{
    const char *integrator = "md";
    const char *thermostat = "berendsen";
    const char *barostat   = "berendsen";

    prepareArgonSimulation(integrator, thermostat, barostat);
    runTest(isVelocityVerlet(integrator));
}

TEST_F(ArgonContinuationIsExactWith, Leapfrog_BerendsenThermostat_ParrinelloRahmanBarostat)
{
    const char *integrator = "md";
    const char *thermostat = "berendsen";
    const char *barostat   = "parrinello-rahman";

    prepareArgonSimulation(integrator, thermostat, barostat);
    runTest(isVelocityVerlet(integrator));
}

TEST_F(ArgonContinuationIsExactWith, Leapfrog_VrescaleThermostat_NoBarostat)
{
    const char *integrator = "md";
    const char *thermostat = "vrescale";
    const char *barostat   = "no";

    prepareArgonSimulation(integrator, thermostat, barostat);
    runTest(isVelocityVerlet(integrator));
}

TEST_F(ArgonContinuationIsExactWith, Leapfrog_VrescaleThermostat_BerendsenBarostat)
{
    const char *integrator = "md";
    const char *thermostat = "vrescale";
    const char *barostat   = "berendsen";

    prepareArgonSimulation(integrator, thermostat, barostat);
    runTest(isVelocityVerlet(integrator));
}

TEST_F(ArgonContinuationIsExactWith, Leapfrog_VrescaleThermostat_ParrinelloRahmanBarostat)
{
    const char *integrator = "md";
    const char *thermostat = "vrescale";
    const char *barostat   = "parrinello-rahman";

    prepareArgonSimulation(integrator, thermostat, barostat);
    runTest(isVelocityVerlet(integrator));
}

// md-vv tests

TEST_F(ArgonContinuationIsExactWith, VelocityVerlet_NVE)
{
    const char *integrator = "md-vv";
    const char *thermostat = "no";
    const char *barostat   = "no";

    prepareArgonSimulation(integrator, thermostat, barostat);
    runTest(isVelocityVerlet(integrator));
}

TEST_F(ArgonContinuationIsExactWith, VelocityVerlet_NoThermostat_BerendsenBarostat)
{
    const char *integrator = "md-vv";
    const char *thermostat = "no";
    const char *barostat   = "berendsen";

    prepareArgonSimulation(integrator, thermostat, barostat);
    runTest(isVelocityVerlet(integrator));
}

TEST_F(ArgonContinuationIsExactWith, VelocityVerlet_BerendsenThermostat_NoBarostat)
{
    const char *integrator = "md-vv";
    const char *thermostat = "berendsen";
    const char *barostat   = "no";

    prepareArgonSimulation(integrator, thermostat, barostat);
    runTest(isVelocityVerlet(integrator));
}

TEST_F(ArgonContinuationIsExactWith, VelocityVerlet_BerendsenThermostat_BerendsenBarostat)
{
    const char *integrator = "md-vv";
    const char *thermostat = "berendsen";
    const char *barostat   = "berendsen";

    prepareArgonSimulation(integrator, thermostat, barostat);
    runTest(isVelocityVerlet(integrator));
}

TEST_F(ArgonContinuationIsExactWith, VelocityVerlet_VrescaleThermostat_NoBarostat)
{
    const char *integrator = "md-vv";
    const char *thermostat = "v-rescale";
    const char *barostat   = "no";

    prepareArgonSimulation(integrator, thermostat, barostat);
    runTest(isVelocityVerlet(integrator));
}

TEST_F(ArgonContinuationIsExactWith, VelocityVerlet_VrescaleThermostat_BerendsenBarostat)
{
    const char *integrator = "md-vv";
    const char *thermostat = "v-rescale";
    const char *barostat   = "berendsen";

    prepareArgonSimulation(integrator, thermostat, barostat);
    runTest(isVelocityVerlet(integrator));
}

TEST_F(ArgonContinuationIsExactWith, VelocityVerlet_NoseHooverThermostat_NoBarostat)
{
    const char *integrator = "md-vv";
    const char *thermostat = "nose-hoover";
    const char *barostat   = "no";

    prepareArgonSimulation(integrator, thermostat, barostat);
    runTest(isVelocityVerlet(integrator));
}

TEST_F(ArgonContinuationIsExactWith, VelocityVerlet_NoseHooverThermostat_MTTKBarostat)
{
    const char *integrator = "md-vv";
    const char *thermostat = "nose-hoover";
    const char *barostat   = "mttk";

    prepareArgonSimulation(integrator, thermostat, barostat);
    runTest(isVelocityVerlet(integrator));
}

// stochastic dynamics tests

TEST_F(ArgonContinuationIsExactWith, StochasticDynamics_NoBarostat)
{
    const char *integrator = "sd";
    const char *thermostat = "no"; // required for sd
    const char *barostat   = "no";

    prepareArgonSimulation(integrator, thermostat, barostat);
    runTest(isVelocityVerlet(integrator));
}

TEST_F(ArgonContinuationIsExactWith, StochasticDynamics_BerendsenBarostat)
{
    const char *integrator = "sd";
    const char *thermostat = "no"; // required for sd
    const char *barostat   = "berendsen";

    prepareArgonSimulation(integrator, thermostat, barostat);
    runTest(isVelocityVerlet(integrator));
}

TEST_F(ArgonContinuationIsExactWith, BrownianDynamics_NoBarostat)
{
    const char *integrator = "bd";
    const char *thermostat = "no"; // required for bd
    const char *barostat   = "no";

    prepareArgonSimulation(integrator, thermostat, barostat);
    runTest(isVelocityVerlet(integrator));
}

TEST_F(ArgonContinuationIsExactWith, BrownianDynamics_BerendsenBarostat)
{
    const char *integrator = "bd";
    const char *thermostat = "no"; // required for bd
    const char *barostat   = "berendsen";

    prepareArgonSimulation(integrator, thermostat, barostat);
    runTest(isVelocityVerlet(integrator));
}

} // namespace
} // namespace
} // namespace
