/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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

/*! \brief Test fixture for mdrun -rerun
 *
 * This test ensures mdrun can run a simulation, writing a trajectory
 * and matching energies, and reproduce the same energies from a rerun
 * to within a tight tolerance. It says nothing about whether an mdrun
 * rerun can reproduce energies from a trajectory generated with older
 * code, since that is not a useful property. Whether mdrun produced
 * correct energies then and now needs different kinds of testing, but
 * if true, this test ensures the rerun has the expected property.
 *
 * Reproducing the same energies is currently only meaningful for
 * integration without thermostats or barostats, however the present
 * form of the test infrastructure has in-principle support for such,
 * if that is ever needed.
 *
 * \todo Add virtual sites peptide-in-water input case, so that
 * is also tested. */
class MdrunRerunCorrectlyReproduces : public MdrunTestFixture,
                                      public ::testing::WithParamInterface<std::tuple<const char *> >
{
    public:
        //! Run a normal mdrun, its rerun, and compare the results
        void runTest(std::vector<std::string> namesOfRequiredFields,
                     FloatingPointTolerance   testTolerance)
        {
            EXPECT_EQ(0, runner_.callGrompp());

            // prepare some names for files to use with the two mdrun calls
            std::string rerunTrajectoryFileName = fileManager_.getTemporaryFilePath("normal.trr");
            std::string normalRunEdrFileName    = fileManager_.getTemporaryFilePath("normal.edr");
            std::string rerunEdrFileName        = fileManager_.getTemporaryFilePath("rerun.edr");

            // do a normal mdrun
            {
                runner_.fullPrecisionTrajectoryFileName_ = rerunTrajectoryFileName;
                runner_.edrFileName_                     = normalRunEdrFileName;
                CommandLine normalRunCaller;
                normalRunCaller.append("mdrun");
                ASSERT_EQ(0, runner_.callMdrun(normalRunCaller));
            }

            // do a rerun on the .trr just produced
            {
                runner_.fullPrecisionTrajectoryFileName_.clear();
                runner_.edrFileName_ = rerunEdrFileName;
                CommandLine rerunCaller;
                rerunCaller.append("mdrun");
                rerunCaller.addOption("-rerun", rerunTrajectoryFileName);
                ASSERT_EQ(0, runner_.callMdrun(rerunCaller));
            }

            // Prepare objects to read the energy files
            EnergyFileReader  normalRunFile(normalRunEdrFileName);
            EnergyFileReader  rerunFile(rerunEdrFileName);

            EnergyFrameReader normalRun(normalRunFile.openToReadFields(namesOfRequiredFields));
            EnergyFrameReader rerun    (rerunFile.openToReadFields(namesOfRequiredFields));

            do
            {
                EnergyFrameInfo normalRunFrame = normalRun.readNextFrame();
                EnergyFrameInfo rerunFrame     = rerun.readNextFrame();
                if (normalRunFrame && rerunFrame)
                {
                    for (auto const &name : namesOfRequiredFields)
                    {
                        EXPECT_REAL_EQ_TOL(normalRunFrame.getValue(name), rerunFrame.getValue(name), testTolerance)
                        << name << " didn't match between normal run " << normalRunFrame.getFrameName() << " and rerun " << rerunFrame.getFrameName();
                    }
                }
                else
                {
                    // At least one file is out of frames. Report if the two
                    // files had different numbers of frames.
                    if (normalRunFrame != rerunFrame)
                    {
                        EXPECT_FALSE(rerunFrame) << "rerun energy file had at least one extra frame";
                        EXPECT_FALSE(normalRunFrame) << "normal run energy file had at least one extra frame";
                    }
                    break;
                }
            }
            while (true);
        }
};

TEST_P(MdrunRerunCorrectlyReproduces, ArgonSimulation)
{
    const char *integrator      = std::get<0>(GetParam());

    /* Do a simulation on argon
     * - writing frames from different kinds of steps: starting, ending, intermediate NS, intermediate non-NS
     * - with other steps between frame-writing steps
     * - with enough buffer that the rerun will compute the same potential energy even though it does NS every frame
     * - can generalize to using thermostats and barostats if that's ever useful
     * - without constraints
     * whose rerun will match with a tight tolerance.
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
                                            "tcoupl = no\n"
                                            "ref-t = 298.0\n"
                                            "tau-t = 1.0\n"
                                            "tc-grps = System\n"
                                            "pcoupl = no\n"
                                            "ref-p = 1\n"
                                            "compressibility = 5e-5\n",
                                            integrator));

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
     * tolerance. Rerun does pair search every frame, so it cannot in
     * general exactly reproduce quantities from a normal run. (Nor
     * does it reproduce pair-search frames exactly, either). */
    FloatingPointTolerance testTolerance = relativeToleranceAsPrecisionDependentUlp(1.0, 120, 1.2e5);

    runTest(namesOfRequiredFields, testTolerance);
}

TEST_P(MdrunRerunCorrectlyReproduces, WaterSimulation)
{
    const char *integrator      = std::get<0>(GetParam());

    /* Do a simulation
     * - writing frames from different kinds of steps: starting, ending, intermediate NS, intermediate non-NS
     * - with other steps between frame-writing steps
     * - with enough buffer that the rerun will compute the same potential energy even though it does NS every frame
     * - can generalize to using thermostats and barostats if that's ever useful
     * - with constraints (ie. SETTLE)
     * whose rerun will match with a tight tolerance.
     */
    runner_.useStringAsMdpFile(formatString("rcoulomb = 0.8\n"
                                            "rvdw = 0.8\n"
                                            "rlist = -1\n"
                                            "bd-fric = 1000\n"
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
                                            "tcoupl = no\n"
                                            "ref-t = 298.0\n"
                                            "tau-t = 1.0\n"
                                            "tc-grps = System\n"
                                            "pcoupl = no\n"
                                            "ref-p = 1\n"
                                            "compressibility = 5e-5\n",
                                            integrator));

    // make the tpr file
    runner_.useTopGroAndNdxFromDatabase("spc216");

    // Describe the energy-file fields that we wish to compare between
    // the two runs. We do not compare pressure, because with
    // constraints the non-search steps need a much larger tolerance,
    // and per Redmine 1858 we should stop computing pressure in
    // reruns anyway.
    std::vector<std::string> namesOfRequiredFields({{interaction_function[F_EPOT].longname,
                                                     interaction_function[F_EKIN].longname,
                                                     interaction_function[F_ETOT].longname,
                                                     interaction_function[F_TEMP].longname }});

    /* NOTE This is an arbitrary-but-sufficient choice for the
     * tolerance. Rerun does pair search every frame, so it cannot in
     * general exactly reproduce quantities from a normal run. (Nor
     * does it reproduce pair-search frames exactly, either). */
    FloatingPointTolerance testTolerance = relativeToleranceAsPrecisionDependentUlp(1.0, 30, 1.2e5);

    if (0 != std::strcmp("bd", integrator))
    {
        runTest(namesOfRequiredFields, testTolerance);
    }
    else
    {
        // TODO we can start testing bd integrator when we stop
        // expecting rerun to reproduce KE, per Redmine 1858
        fprintf(stderr, "bd + constraints rerun not tested because of known bug with first-frame KE\n");
    };
}

INSTANTIATE_TEST_CASE_P(WithVariousIntegrators,
                        MdrunRerunCorrectlyReproduces,
                            ::testing::Values("md", "md-vv", "bd", "sd"));

/*! \todo Add other tests for mdrun -rerun, e.g.
 *
 * - TpiExitsNormally (since it uses the -rerun machinery)
 */

} // namespace
} // namespace
} // namespace
