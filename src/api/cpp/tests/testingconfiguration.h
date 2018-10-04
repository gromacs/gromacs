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

#ifndef GROMACS_TESTINGCONFIGURATION_H
#define GROMACS_TESTINGCONFIGURATION_H

#include <gtest/gtest.h>

#include <string>
#include <vector>

#include "gromacs/gmxpreprocess/grompp.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/cmdlinetest.h"
#include "testutils/testfilemanager.h"
#include "programs/mdrun/tests/moduletest.h"

namespace gmxapi
{

namespace testing
{

/*! \brief
 * Helper function to get step size, currently hard coded to 2^-9.
 *
 * The step size was chosen to be exactly representable in binary
 * to avoid surprises when determining integer step numbers.
 *
 * \returns Step size for tests.
 */
inline real getTestStepSize() { return 0.001953125; }

//! Provide command-line infrastructure for gmxapi tests.
class GmxApiTest : public gmx::test::MdrunTestFixture
{
    public:
        GmxApiTest()
        {
            runner_.useTopGroAndNdxFromDatabase("spc-and-methanol");
            runner_.useStringAsMdpFile(gmx::formatString("integrator = md\n"
                                                         "cutoff-scheme = Verlet\n"
                                                         "nsteps = 10000\n"
                                                         "dt = %11.9f\n"
                                                         "nstxout = 2\n"
                                                         "nstvout = 2\n"
                                                         "nstfout = 4\n"
                                                         "nstxout-compressed = 5\n"
                                                         "tcoupl = v-rescale\n"
                                                         "tc-grps = System\n"
                                                         "tau-t = 1\n"
                                                         "ref-t = 298\n"
                                                         "compressed-x-grps = Sol\n", getTestStepSize()));
            EXPECT_EQ(0, runner_.callGromppOnThisRank());

            tprFilename_ = runner_.tprFileName_;
            // This is extremely ugly, but works when it comes to cleaning up the
            // generated temporary files.
        }

        //! Get the TPR file to work with
        std::string tprFileName() const
        {
            GMX_ASSERT(!tprFilename_.empty(), "Need to have declared filename for tpr file");
            return tprFilename_;
        }
        //! Make the md arguments to work with
        std::vector<std::string> makeMdArgs() const
        {
            std::vector<std::string> mdArgs;

            mdArgs.emplace_back("-o");
            mdArgs.emplace_back(runner_.fullPrecisionTrajectoryFileName_);
            mdArgs.emplace_back("-x");
            mdArgs.emplace_back(runner_.reducedPrecisionTrajectoryFileName_);
            mdArgs.emplace_back("-c");
            mdArgs.emplace_back(runner_.groOutputFileName_);
            mdArgs.emplace_back("-g");
            mdArgs.emplace_back(runner_.logFileName_);
            mdArgs.emplace_back("-e");
            mdArgs.emplace_back(runner_.edrFileName_);
            mdArgs.emplace_back("-cpo");
            mdArgs.emplace_back(runner_.cptFileName_);

            return mdArgs;
        }
    private:
        //! Name of input tpr file generated for tests.
        std::string                tprFilename_;
};

} // namespace testing

} // end namespace gmxapi


#endif //GROMACS_TESTINGCONFIGURATION_H
