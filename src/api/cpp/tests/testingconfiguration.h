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

//! Helper function to get step size, currently hard coded to 2^-9
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
            EXPECT_EQ(0, runner_.callGrompp());

            tprFilename_ = runner_.tprFileName_;
            // This is extremely ugly, but works when it comes to cleaning up the
            // generated temporary files.
            mdArgs_.emplace_back("-o");
            mdArgs_.emplace_back(runner_.fullPrecisionTrajectoryFileName_);
            mdArgs_.emplace_back("-x");
            mdArgs_.emplace_back(runner_.reducedPrecisionTrajectoryFileName_);
            mdArgs_.emplace_back("-c");
            mdArgs_.emplace_back(runner_.groOutputFileName_);
            mdArgs_.emplace_back("-g");
            mdArgs_.emplace_back(runner_.logFileName_);
            mdArgs_.emplace_back("-e");
            mdArgs_.emplace_back(runner_.edrFileName_);
            mdArgs_.emplace_back("-cpo");
            mdArgs_.emplace_back(runner_.cptFileName_);
        }

        //! Get the TPR file to work with
        std::string getTprFileName() const
        {
            GMX_ASSERT(!tprFilename_.empty(), "Need to have declared filename for tpr file");
            return tprFilename_;
        }
        //! Get the md arguments to work with
        std::vector<std::string> getMdArgs() const { return mdArgs_; }
        //! Get number of steps to run before stopping.
        int getStopStepNumber() const
        {
            GMX_ASSERT(haveSetStopStepNumber_, "Need to set stopStepNumber before using it");
            return stopStepNumber_;
        }
        //! Set number of steps to run before stopping.
        void setStopStepNumber(int number)
        {
            stopStepNumber_        = number;
            haveSetStopStepNumber_ = true;
        };
    private:
        //! Name of input tpr file generated for tests.
        std::string                tprFilename_;
        //! Arguments to run API driven simulation.
        std::vector<std::string>   mdArgs_;
        //! Time steps to run before stopping a run.
        int                        stopStepNumber_ = -1;
        //! True if we have set the number of steps before stopping.
        bool haveSetStopStepNumber_ = false;
};

} // namespace testing

} // end namespace gmxapi


#endif //GROMACS_TESTINGCONFIGURATION_H
