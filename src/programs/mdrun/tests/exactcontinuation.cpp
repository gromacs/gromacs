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
class MdrunContinuationIsExact : public MdrunTestFixture,
                                 public ::testing::WithParamInterface<std::tuple<const char *, const char *, const char *> >
{
    public:
        /*! \brief Run a normal mdrun, the same mdrun stopping half
         * way, doing a continuation, and compare the results. */
        void runTest(std::vector<std::string> namesOfRequiredFields,
                     FloatingPointTolerance   testTolerance,
                     bool                     isVelocityVerlet)
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

            EnergyFrameReader fullRun(fullRunFile.openToReadFields(namesOfRequiredFields));
            EnergyFrameReader firstHalfRun(firstHalfRunFile.openToReadFields(namesOfRequiredFields));
            EnergyFrameReader secondHalfRun(secondHalfRunFile.openToReadFields(namesOfRequiredFields));

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
                    for (auto const &name : namesOfRequiredFields)
                    {
                        EXPECT_REAL_EQ_TOL(fullRunFrame.getValue(name), currentHalfRunFrame.getValue(name), testTolerance)
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
};

TEST_P(MdrunContinuationIsExact, ArgonSimulation)
{
    const char *integrator       = std::get<0>(GetParam());
    bool        isVelocityVerlet = (0 == std::strcmp(integrator, "md-vv"));
    const char *thermostat       = std::get<1>(GetParam());
    const char *barostat         = std::get<2>(GetParam());

    /* Do a simulation on argon
     * - writing frames from different kinds of steps: starting, ending, intermediate NS, intermediate non-NS
     * - with other steps between frame-writing steps
     * - with enough buffer that the rerun will compute the same potential energy even though it does NS every frame
     * - without constraints
     * whose restarted run will match with a tight tolerance.
     */
    runner_.useStringAsMdpFile(formatString("rcoulomb = 1.0\n"
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
                                            "compressibility = 5e-5\n",
                                            integrator,
                                            thermostat,
                                            barostat));

    // make the tpr file
    runner_.useTopGroAndNdxFromDatabase("argon");

    // Describe the energy-file fields that we wish to compare between
    // the two runs.
    std::vector<std::string> namesOfRequiredFields({{interaction_function[F_EPOT].longname,
                                                     interaction_function[F_EKIN].longname,
                                                     interaction_function[F_ETOT].longname,
                                                     interaction_function[F_TEMP].longname,
                                                     interaction_function[F_PRES].longname }});

    /* NOTE This is an arbitrary-but-sufficient choice for the
     * tolerance. */
    FloatingPointTolerance testTolerance = relativeToleranceAsPrecisionDependentUlp(1.0, 160, 1.2e5);

    runTest(namesOfRequiredFields, testTolerance, isVelocityVerlet);
}

//! Convenience typedef
typedef std::tuple<const char *, const char *, const char*> ContinuationTestType;

INSTANTIATE_TEST_CASE_P(WithVariousLeapFrogIntegrators,
                        MdrunContinuationIsExact,
                            ::testing::Values(ContinuationTestType("md", "no",        "no"),
                                              ContinuationTestType("md", "no",        "berendsen"),
                                              ContinuationTestType("md", "no",        "parrinello-rahman"),
                                              ContinuationTestType("md", "berendsen", "no"),
                                              ContinuationTestType("md", "berendsen", "berendsen"),
                                              ContinuationTestType("md", "berendsen", "parrinello-rahman"),
                                              ContinuationTestType("md", "v-rescale", "no"),
                                              ContinuationTestType("md", "v-rescale", "berendsen"),
                                              ContinuationTestType("md", "v-rescale", "parrinello-rahman")
                                              ));

INSTANTIATE_TEST_CASE_P(WithVariousVelocityVerletIntegrators,
                        MdrunContinuationIsExact,
                            ::testing::Values(ContinuationTestType("md-vv", "no",          "no"),
                                                                                                         // TODO As noted, some combinations are unsatisfactory. If we can't
                                                                                                         // find a fix, work out how to test that they're expected to be broken.
                                              ContinuationTestType("md-vv", "no",          "berendsen"), // slightly broken
                                              ContinuationTestType("md-vv", "berendsen",   "no"),        // bad
                                              ContinuationTestType("md-vv", "berendsen",   "berendsen"), // bad
                                              ContinuationTestType("md-vv", "v-rescale",   "no"),        // bad
                                              ContinuationTestType("md-vv", "v-rescale",   "berendsen"), // bad
                                              ContinuationTestType("md-vv", "nose-hoover", "no"),
                                              ContinuationTestType("md-vv", "nose-hoover", "mttk")       // slightly broken
                                              ));

INSTANTIATE_TEST_CASE_P(WithVariousStochasticIntegrators,
                        MdrunContinuationIsExact,
                            ::testing::Values(ContinuationTestType("sd", "no", "no"),
                                              ContinuationTestType("sd", "no", "berendsen"),
                                              ContinuationTestType("bd", "no", "no"),
                                              ContinuationTestType("bd", "no", "berendsen")
                                              ));

} // namespace
} // namespace
} // namespace
