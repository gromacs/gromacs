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

#include <string>

#include <gtest/gtest-param-test.h>
#include <gtest/gtest.h>

#include "gromacs/utility/path.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/cmdlinetest.h"
#include "testutils/textblockmatchers.h"

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

using DsspTestParamsDsspNB      = std::tuple<std::string, real>;
using DsspTestParamsGromacsNB   = std::tuple<std::string, real>;
using DsspTestParamsDsspNoNB    = std::string;
using DsspTestParamsGromacsNoNB = std::string;

//! Test fixture for the dssp analysis module.
class DsspModuleTestDsspNB :
    public TrajectoryAnalysisModuleTestFixture<gmx::analysismodules::DsspInfo>,
    public ::testing::WithParamInterface<DsspTestParamsDsspNB>
{
};

//! Test fixture for the dssp analysis module.
class DsspModuleTestGromacsNB :
    public TrajectoryAnalysisModuleTestFixture<gmx::analysismodules::DsspInfo>,
    public ::testing::WithParamInterface<DsspTestParamsGromacsNB>
{
};
//! Test fixture for the dssp analysis module.
class DsspModuleTestDsspNoNB :
    public TrajectoryAnalysisModuleTestFixture<gmx::analysismodules::DsspInfo>,
    public ::testing::WithParamInterface<DsspTestParamsDsspNoNB>
{
};

//! Test fixture for the dssp analysis module.
class DsspModuleTestGromacsNoNB :
    public TrajectoryAnalysisModuleTestFixture<gmx::analysismodules::DsspInfo>,
    public ::testing::WithParamInterface<DsspTestParamsGromacsNoNB>
{
};

// See https://github.com/google/googletest/issues/2442 for the reason
// for this and following NOLINTNEXTLINE suppressions.

// NOLINTNEXTLINE(cppcoreguidelines-owning-memory)
TEST_P(DsspModuleTestDsspNB, Works)
{
    std::tuple<std::string, real> params    = GetParam();
    const char* const             cmdline[] = { "dssp" };
    std::string                   inputFilename(std::get<0>(params));
    std::filesystem::path         inputBasename = stripExtension(inputFilename);
    CommandLine                   command(cmdline);
    setTopology(inputFilename.c_str());
    setTrajectory(inputFilename.c_str());
    setOutputFile("-o",
                  formatString("%s-dssp-nb-%.1f.dat", inputBasename.c_str(), std::get<1>(params)).c_str(),
                  ExactTextMatch());
    command.addOption("-hmode", "dssp");
    command.addOption("-nb");
    command.addOption("-cutoff", std::get<1>(params));
    runTest(command);
}

INSTANTIATE_TEST_SUITE_P(MoleculeTests,
                         DsspModuleTestDsspNB,
                         ::testing::Combine(::testing::Values("1cos.pdb",
                                                              "1hlc.pdb",
                                                              "1vzj.pdb",
                                                              "3byc.pdb",
                                                              "3kyy.pdb",
                                                              "4r80.pdb",
                                                              "4xjf.pdb",
                                                              "5u5p.pdb",
                                                              "7wgh.pdb",
                                                              "1gmc.pdb",
                                                              "1v3y.pdb",
                                                              "1yiw.pdb",
                                                              "2os3.pdb",
                                                              "3u04.pdb",
                                                              "4r6c.pdb",
                                                              "4wxl.pdb",
                                                              "5cvq.pdb",
                                                              "5i2b.pdb",
                                                              "5t8z.pdb",
                                                              "6jet.pdb"),
                                            ::testing::Values(0.9, 2.0)));


// NOLINTNEXTLINE(cppcoreguidelines-owning-memory)
TEST_P(DsspModuleTestGromacsNB, Works)
{
    std::tuple<std::string, real> params    = GetParam();
    const char* const             cmdline[] = { "dssp" };
    std::string                   inputFilename(std::get<0>(params));
    std::filesystem::path         inputBasename = stripExtension(inputFilename);
    CommandLine                   command(cmdline);
    setTopology(inputFilename.c_str());
    setTrajectory(inputFilename.c_str());
    setOutputFile("-o",
                  formatString("%s-gromacs-nb-%.1f.dat", inputBasename.c_str(), std::get<1>(params)).c_str(),
                  ExactTextMatch());
    command.addOption("-hmode", "gromacs");
    command.addOption("-nb");
    command.addOption("-cutoff", std::get<1>(params));
    runTest(command);
}

INSTANTIATE_TEST_SUITE_P(
        MoleculeTests,
        DsspModuleTestGromacsNB,
        ::testing::Combine(::testing::Values("hdac.pdb", "RNAseA.pdb", "zyncfinger.pdb"),
                           ::testing::Values(0.9, 2.0)));

// NOLINTNEXTLINE(cppcoreguidelines-owning-memory)
TEST_P(DsspModuleTestDsspNoNB, Works)
{
    const char* const     cmdline[] = { "dssp" };
    std::string           inputFilename(GetParam());
    std::filesystem::path inputBasename = stripExtension(inputFilename);
    CommandLine           command(cmdline);
    setTopology(inputFilename.c_str());
    setTrajectory(inputFilename.c_str());
    setOutputFile("-o", formatString("%s-dssp-nonb.dat", inputBasename.c_str()).c_str(), ExactTextMatch());
    command.addOption("-hmode", "dssp");
    command.addOption("-nonb");
    runTest(command);
}

INSTANTIATE_TEST_SUITE_P(MoleculeTests,
                         DsspModuleTestDsspNoNB,
                         ::testing::Values("1cos.pdb",
                                           "1hlc.pdb",
                                           "1vzj.pdb",
                                           "3byc.pdb",
                                           "3kyy.pdb",
                                           "4r80.pdb",
                                           "4xjf.pdb",
                                           "5u5p.pdb",
                                           "7wgh.pdb",
                                           "1gmc.pdb",
                                           "1v3y.pdb",
                                           "1yiw.pdb",
                                           "2os3.pdb",
                                           "3u04.pdb",
                                           "4r6c.pdb",
                                           "4wxl.pdb",
                                           "5cvq.pdb",
                                           "5i2b.pdb",
                                           "5t8z.pdb",
                                           "6jet.pdb"));


// NOLINTNEXTLINE(cppcoreguidelines-owning-memory)
TEST_P(DsspModuleTestGromacsNoNB, Works)
{
    const char* const     cmdline[] = { "dssp" };
    std::string           inputFilename(GetParam());
    std::filesystem::path inputBasename = stripExtension(inputFilename);
    CommandLine           command(cmdline);
    setTopology(inputFilename.c_str());
    setTrajectory(inputFilename.c_str());
    setOutputFile(
            "-o", formatString("%s-gromacs-nonb.dat", inputBasename.c_str()).c_str(), ExactTextMatch());
    command.addOption("-hmode", "gromacs");
    command.addOption("-nonb");
    runTest(command);
}

INSTANTIATE_TEST_SUITE_P(MoleculeTests,
                         DsspModuleTestGromacsNoNB,
                         ::testing::Values("hdac.pdb", "RNAseA.pdb", "zyncfinger.pdb"));

} // namespace
} // namespace test
} // namespace gmx
