/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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

/*! \internal \file
 * \brief
 * Checks that expected output files are present
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_mdrun_integration_tests
 */

#include "gmxpre.h"

#include <filesystem>
#include <initializer_list>
#include <string>
#include <tuple>

#include <gtest/gtest.h>

#include "gromacs/utility/path.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/cmdlinetest.h"
#include "testutils/simulationdatabase.h"

#include "moduletest.h"

namespace gmx
{
namespace test
{
namespace
{

/*! \brief Test fixture base for presence of output files
 *
 * This test checks that some expected output files are actually written
 * in a very simple mdrun call. Currently, it tests for the presence of
 * full precision trajectory (.trr), log file (.log), energy trajectory (.edr),
 * final configuration (.gro) and checkpoint file (.cpt).
 *
 * The test only checks whether the files are existing, it does not check
 * their contents.
 */
using OutputFilesTestParams = std::tuple<std::string, std::string>;
class OutputFiles : public MdrunTestFixture, public ::testing::WithParamInterface<OutputFilesTestParams>
{
};

TEST_P(OutputFiles, FilesArePresent)
{
    auto params         = GetParam();
    auto simulationName = std::get<0>(params);
    auto integrator     = std::get<1>(params);

    SCOPED_TRACE(
            formatString("Checking for presence of expected output files using "
                         "simulation '%s' with integrator '%s'",
                         simulationName.c_str(),
                         integrator.c_str()));

    // Prepare the .tpr file
    {
        CommandLine caller;
        runner_.useTopGroAndNdxFromDatabase(simulationName);
        runner_.useStringAsMdpFile(prepareMdpFileContents(
                prepareMdpFieldValues(simulationName.c_str(), integrator.c_str(), "no", "no")));
        EXPECT_EQ(0, runner_.callGrompp(caller));
    }
    // Do mdrun
    {
        CommandLine mdrunCaller;
        ASSERT_EQ(0, runner_.callMdrun(mdrunCaller));
    }
    // Check if expected files are present
    {
        for (const auto& file : { runner_.fullPrecisionTrajectoryFileName_,
                                  runner_.logFileName_,
                                  runner_.edrFileName_,
                                  runner_.groOutputFileName_,
                                  runner_.cptOutputFileName_ })
        {
            EXPECT_TRUE(File::exists(file, File::returnFalseOnError))
                    << "File " << file << " was not found.";
        }
    }
}

INSTANTIATE_TEST_SUITE_P(Argon12,
                         OutputFiles,
                         ::testing::Combine(::testing::Values("argon12"),
                                            ::testing::Values("md", "md-vv")));

} // namespace
} // namespace test
} // namespace gmx
