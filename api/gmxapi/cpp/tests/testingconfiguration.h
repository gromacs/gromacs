/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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

#ifndef GROMACS_TESTINGCONFIGURATION_H
#define GROMACS_TESTINGCONFIGURATION_H

#include "config.h"

#include <algorithm>
#include <string>
#include <vector>

#include "gromacs/gmxpreprocess/grompp.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/cmdlinetest.h"
#include "testutils/testfilemanager.h"

#include "programs/mdrun/tests/moduletest.h"
#if GMX_THREAD_MPI
#    include "testutils/mpitest.h"
#endif

namespace gmxapi
{

namespace testing
{

/*! \brief Helper function to get step size as floating point number.
 *
 * A step size of 0.002ps is suitable for the simulation.
 * We prefer to guarantee that the testing tools easily anticipate TPR time step size
 * and time value in trajectory outputs, so we explicitly choose a time step that is
 * exactly representable in binary. 2^-9 is just smaller than 0.002 and requires very
 * little floating point precision (mantissa == 0). This means that sums and multiples
 * of the timestep are also exactly representable, and thus particularly convenient for tests."
 *
 * For human readability we use the decimal representation, 1.0 x 2^-9 = 0.001953125.
 *
 * \returns Step size for tests.
 */
inline real getTestStepSize()
{
    return 0.001953125;
}

//! Provide command-line infrastructure for gmxapi tests.
class GmxApiTest : public gmx::test::MdrunTestFixture
{
public:
    GmxApiTest()
    {
        // For files that will always be generated, let's explicitly set the names so that
        // TestFileManager will clean them up properly.
        runner_.groOutputFileName_ = fileManager_.getTemporaryFilePath("confout.gro");
        runner_.cptOutputFileName_ = fileManager_.getTemporaryFilePath("state.cpt");
        runner_.reducedPrecisionTrajectoryFileName_ =
                fileManager_.getTemporaryFilePath("traj_comp.xtc");
    }

    /* \brief
     * Prepare a tpr to run the test with.
     *
     * Sets up the TPR to run a test of the GMXAPI with a set number of \p steps
     * defined in the test.
     *
     * \param[in] steps Number of steps for test to run.
     */
    void makeTprFile(int steps)
    {
        runner_.useTopGroAndNdxFromDatabase("spc_and_methane");
        runner_.useStringAsMdpFile(
                gmx::formatString("integrator = md\n"
                                  "cutoff-scheme = Verlet\n"
                                  "nsteps = %d\n"
                                  "dt = %11.9f\n"
                                  "nstxout = 2\n"
                                  "nstvout = 2\n"
                                  "nstfout = 4\n"
                                  "nstxout-compressed = 2\n"
                                  "tcoupl = v-rescale\n"
                                  "tc-grps = System\n"
                                  "tau-t = 1\n"
                                  "ref-t = 298\n"
                                  "compressed-x-grps = Sol\n",
                                  steps,
                                  getTestStepSize()));

        EXPECT_EQ(0, runner_.callGromppOnThisRank());
    }

    //! Make the md arguments to work with
    std::vector<std::string> makeMdArgs() const
    {
        std::vector<std::string> mdArgs;

        if (!runner_.fullPrecisionTrajectoryFileName_.empty())
        {
            EXPECT_TRUE(std::none_of(mdArgs.cbegin(), mdArgs.cend(), [](const std::string& arg) {
                return arg == "-o";
            }));
            mdArgs.emplace_back("-o");
            mdArgs.emplace_back(runner_.fullPrecisionTrajectoryFileName_);
        }
        if (!runner_.reducedPrecisionTrajectoryFileName_.empty())
        {
            EXPECT_TRUE(std::none_of(mdArgs.cbegin(), mdArgs.cend(), [](const std::string& arg) {
                return arg == "-x";
            }));
            mdArgs.emplace_back("-x");
            mdArgs.emplace_back(runner_.reducedPrecisionTrajectoryFileName_);
        }
        if (!runner_.groOutputFileName_.empty())
        {
            EXPECT_TRUE(std::none_of(mdArgs.cbegin(), mdArgs.cend(), [](const std::string& arg) {
                return arg == "-c";
            }));
            mdArgs.emplace_back("-c");
            mdArgs.emplace_back(runner_.groOutputFileName_);
        }
        if (!runner_.logFileName_.empty())
        {
            EXPECT_TRUE(std::none_of(mdArgs.cbegin(), mdArgs.cend(), [](const std::string& arg) {
                return arg == "-g";
            }));
            mdArgs.emplace_back("-g");
            mdArgs.emplace_back(runner_.logFileName_);
        }
        if (!runner_.edrFileName_.empty())
        {
            EXPECT_TRUE(std::none_of(mdArgs.cbegin(), mdArgs.cend(), [](const std::string& arg) {
                return arg == "-e";
            }));
            mdArgs.emplace_back("-e");
            mdArgs.emplace_back(runner_.edrFileName_);
        }
        if (!runner_.cptOutputFileName_.empty())
        {
            EXPECT_TRUE(std::none_of(mdArgs.cbegin(), mdArgs.cend(), [](const std::string& arg) {
                return arg == "-cpo";
            }));
            mdArgs.emplace_back("-cpo");
            mdArgs.emplace_back(runner_.cptOutputFileName_);
        }
#if GMX_THREAD_MPI
        mdArgs.emplace_back("-ntmpi");
        mdArgs.emplace_back(std::to_string(gmx::test::getNumberOfTestMpiRanks()));
#endif

#if GMX_OPENMP
        mdArgs.emplace_back("-ntomp");
        mdArgs.emplace_back(std::to_string(gmx::test::getNumberOfTestOpenMPThreads()));
#endif

        return mdArgs;
    }
};

} // namespace testing

} // end namespace gmxapi


#endif // GROMACS_TESTINGCONFIGURATION_H
