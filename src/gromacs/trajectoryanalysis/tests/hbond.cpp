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
 * Tests for functionality of the "hbond" trajectory analysis module.
 *
 * \author Sergey Gorelov <gorelov_sv@pnpi.nrcki.ru>
 * \author Alexey Shvetsov <alexxyum@gmail.com>
 * \ingroup module_trajectoryanalysis
 */

#include "gmxpre.h"

#include "gromacs/trajectoryanalysis/modules/hbond.h"

#include <string>
#include <tuple>

#include <gtest/gtest-param-test.h>
#include <gtest/gtest.h>

#include "gromacs/utility/stringutil.h"

#include "testutils/cmdlinetest.h"
#include "testutils/testasserts.h"
#include "testutils/textblockmatchers.h"
#include "testutils/xvgtest.h"

#include "moduletest.h"

namespace gmx
{
namespace test
{

/********************************************************************
 * Tests for gmx::analysismodules::Hbond.
 */

using HbondTestParamsMerge = std::tuple<std::string, std::tuple<std::string, std::string>, std::string>;

//! Test fixture for the hbond analysis module.
class HbondModuleTest :
    public TrajectoryAnalysisModuleTestFixture<gmx::analysismodules::HbondInfo>,
    public ::testing::WithParamInterface<HbondTestParamsMerge>
{
};

// See https://github.com/google/googletest/issues/2442 for the reason
// for this and following NOLINTNEXTLINE suppressions.

// NOLINTNEXTLINE(cppcoreguidelines-owning-memory)
TEST_P(HbondModuleTest, Works)
{
    std::tuple<std::string, std::tuple<std::string, std::string>, std::string> params = GetParam();
    const char* const cmdline[]     = { "hbond2" };
    std::string       inputBasename = std::get<0>(params);
    CommandLine       command(cmdline);
    double            tolerance = 1e-2;
    test::XvgMatch    matcher;
    test::XvgMatch&   toleranceMatch =
            matcher.tolerance(gmx::test::relativeToleranceAsFloatingPoint(1, tolerance));
    setTrajectory((inputBasename + ".xtc").c_str());
    setTopology((inputBasename + ".tpr").c_str());
    std::string mergeCommand = "-" + std::get<2>(params);
    setOutputFile("-o",
                  formatString("%s-r%s-t%s%s.ndx",
                               inputBasename.c_str(),
                               std::get<0>(std::get<1>(params)).c_str(),
                               std::get<1>(std::get<1>(params)).c_str(),
                               mergeCommand.c_str())
                          .c_str(),
                  ExactTextMatch());
    command.addOption("-r", std::get<0>(std::get<1>(params)));
    command.addOption("-t", std::get<1>(std::get<1>(params)));
    setOutputFile("-num",
                  formatString("%s-r%s-t%s%s-num.xvg",
                               inputBasename.c_str(),
                               std::get<0>(std::get<1>(params)).c_str(),
                               std::get<1>(std::get<1>(params)).c_str(),
                               mergeCommand.c_str())
                          .c_str(),
                  test::XvgMatch());
    setOutputFile("-dist",
                  formatString("%s-r%s-t%s%s-dist.xvg",
                               inputBasename.c_str(),
                               std::get<0>(std::get<1>(params)).c_str(),
                               std::get<1>(std::get<1>(params)).c_str(),
                               mergeCommand.c_str())
                          .c_str(),
                  toleranceMatch);
    setOutputFile("-ang",
                  formatString("%s-r%s-t%s%s-ang.xvg",
                               inputBasename.c_str(),
                               std::get<0>(std::get<1>(params)).c_str(),
                               std::get<1>(std::get<1>(params)).c_str(),
                               mergeCommand.c_str())
                          .c_str(),
                  toleranceMatch);
    setOutputFile("-dan",
                  formatString("%s-r%s-t%s%s-dan.xvg",
                               inputBasename.c_str(),
                               std::get<0>(std::get<1>(params)).c_str(),
                               std::get<1>(std::get<1>(params)).c_str(),
                               mergeCommand.c_str())
                          .c_str(),
                  test::XvgMatch());
    command.addOption(mergeCommand.c_str());
    setDatasetTolerance("hbdist", gmx::test::relativeToleranceAsFloatingPoint(1, tolerance));
    setDatasetTolerance("histogram_dist", gmx::test::relativeToleranceAsFloatingPoint(1, tolerance));
    setDatasetTolerance("hbang", gmx::test::relativeToleranceAsFloatingPoint(1, tolerance));
    setDatasetTolerance("histogram_ang", gmx::test::relativeToleranceAsFloatingPoint(1, tolerance));
    runTest(command);
}

INSTANTIATE_TEST_SUITE_P(HBondTests,
                         HbondModuleTest,
                         ::testing::Combine(::testing::Values("trpcage"),
                                            ::testing::Values(std::make_tuple("Protein", "Protein"),
                                                              std::make_tuple("Protein", "Water"),
                                                              std::make_tuple("Water", "Water")),
                                            ::testing::Values("nom", "m")));


} // namespace test
} // namespace gmx
