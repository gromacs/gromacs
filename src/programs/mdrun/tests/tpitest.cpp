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
 * This implements basic TPI sanity test.
 * It runs the input system with TPI for several steps, and checks the log output.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

// FIXME: we should be properly enabling (or, for TPI, disabling) FP exceptions in mdrun, not here.
// Remove this include and the exception toggling calls.
#include "gromacs/math/utilities.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textreader.h"

#include "testutils/refdata.h"

#include "moduletest.h"

namespace gmx
{
namespace test
{
namespace
{

/*! \brief A basic TPI runner */
class TpiTest : public MdrunTestFixture
{
    public:
        //! Runs the test with the given inputs
        void runTest();
};

void TpiTest::runTest()
{
    runner_.useTopGroAndNdxFromDatabase("spc216_with_methane");
    runner_.ndxFileName_ = "";
    ASSERT_EQ(0, runner_.callGrompp());

    gmx_fedisableexcept();

    auto        rerunFileName = fileManager_.getInputFilePath("spc216.gro");
    CommandLine commandLine;
    commandLine.append("-rerun");
    commandLine.append(rerunFileName);
    ASSERT_EQ(0, runner_.callMdrun(commandLine));

    gmx_feenableexcept();

    std::string logFileContexts = TextReader::readFileToString(runner_.logFileName_);
    std::string tpiOutputs      = logFileContexts.substr(logFileContexts.find("Started Test Particle Insertion"));
    auto        startIndex      = tpiOutputs.find("  <V>");
    auto        endIndex        = tpiOutputs.find("Finished mdrun");
    std::string expectedOutput  = tpiOutputs.substr(startIndex, endIndex - startIndex);


    TestReferenceData    refData;
    TestReferenceChecker checker(refData.rootChecker());
    checker.checkString(expectedOutput, "ExpectedOutput");
}

TEST_F(TpiTest, ReproducesOutput)
{
    const int         nsteps          = 200;
    const std::string mdpFileContents = formatString(R"(
        integrator               = tpi
        ld-seed                  = 1993
        rtpi                     = 0.05
        nstlog                   = 0
        nstenergy                = 0
        cutoff-scheme            = group
        nstlist                  = 10
        ns_type                  = grid
        rlist                    = 0.9
        coulombtype              = reaction-field
        rcoulomb                 = 0.9
        epsilon-r                = 1
        epsilon-rf               = 0
        vdw-type                 = cut-off
        rvdw                     = 0.9
        Tcoupl                   = no
        tc-grps                  = System
        tau_t                    = 0.5
        ref_t                    = 298
        nsteps                   = %d
    )", nsteps);

    runner_.useStringAsMdpFile(mdpFileContents);
    runTest();
}


}
}
}
