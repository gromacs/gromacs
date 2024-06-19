/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2022- The GROMACS Authors
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
 * Tests for functionality of the "dssp" trajectory analysis module.
 *
 * \author Sergey Gorelov <gorelov_sv@pnpi.nrcki.ru>
 * \author Anatoly Titov <titov_ai@pnpi.nrcki.ru>
 * \author Alexey Shvetsov <alexxyum@gmail.com>
 * \ingroup module_trajectoryanalysis
 */

#include "gmxpre.h"

#include "gromacs/trajectoryanalysis/modules/dssp.h"

#include <filesystem>
#include <string>
#include <tuple>

#include <gtest/gtest-param-test.h>
#include <gtest/gtest.h>

#include "gromacs/utility/path.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/cmdlinetest.h"
#include "testutils/textblockmatchers.h"
#include "testutils/xvgtest.h"

#include "moduletest.h"

namespace gmx
{
namespace test
{
namespace
{

/********************************************************************
 * Tests for gmx::analysismodules::Dssp.
 */

using DsspTestParams =
        std::tuple<std::string, std::string, std::string, real, std::string, std::string, std::string>;

//! Test fixture for the dssp analysis module.
class DsspModuleTest :
    public TrajectoryAnalysisModuleTestFixture<gmx::analysismodules::DsspInfo>,
    public ::testing::WithParamInterface<DsspTestParams>
{
};

// See https://github.com/google/googletest/issues/2442 for the reason
// for this and following NOLINTNEXTLINE suppressions.

// NOLINTNEXTLINE(cppcoreguidelines-owning-memory)
TEST_P(DsspModuleTest, Works)
{
    DsspTestParams              params    = GetParam();
    const char* const           cmdline[] = { "dssp" };
    const std::string           inputFilename(std::get<0>(params));
    const std::filesystem::path inputBasename = stripExtension(inputFilename);
    CommandLine                 command(cmdline);
    setTopology(inputFilename.c_str());
    setTrajectory(inputFilename.c_str());
    setOutputFile("-o",
                  formatString("%s-%s-%s-%.1f-%s-%s.dat",
                               inputBasename.c_str(),
                               std::get<1>(params).c_str(),
                               std::get<2>(params).c_str(),
                               std::get<3>(params),
                               std::get<4>(params).c_str(),
                               std::get<5>(params).c_str())
                          .c_str(),
                  ExactTextMatch());
    command.addOption("-hmode", std::get<1>(params));
    command.addOption(std::string("-" + std::get<2>(params)).c_str());
    command.addOption("-cutoff", std::get<3>(params));
    command.addOption("-hbond", std::get<4>(params));
    command.addOption(std::string("-" + std::get<5>(params)).c_str());
    command.addOption(std::string("-" + std::get<6>(params)).c_str());
    setOutputFile("-num",
                  formatString("%s-%s-%s-%.1f-%s-%s.xvg",
                               inputBasename.c_str(),
                               std::get<1>(params).c_str(),
                               std::get<2>(params).c_str(),
                               std::get<3>(params),
                               std::get<4>(params).c_str(),
                               std::get<5>(params).c_str())
                          .c_str(),
                  test::XvgMatch());
    runTest(command);
}

INSTANTIATE_TEST_SUITE_P(
        MoleculeTests,
        DsspModuleTest,
        ::testing::Combine(::testing::Values("nonpdb.pdb", "RNAseA.pdb", "zyncfinger.pdb"),
                           ::testing::Values("dssp", "gromacs"),
                           ::testing::Values("nb", "nonb"),
                           ::testing::Values(0.9, 2.0),
                           ::testing::Values("energy", "geometry"),
                           ::testing::Values("clear", "noclear"),
                           ::testing::Values("polypro", "nopolypro")));
} // namespace
} // namespace test
} // namespace gmx
