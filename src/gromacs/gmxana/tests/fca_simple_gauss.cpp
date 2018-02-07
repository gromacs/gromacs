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
 * Tests special cases in domain decomposition
 *
 * \author Berenger Bramas <berenger.bramas@mpcdf.mpg.de>
 * \ingroup module_fca_integration_tests
 */
#include "gmxpre.h"

#include <gtest/gtest.h>

#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/utility/classhelpers.h"

#include "testutils/cmdlinetest.h"
#include "testutils/xvgtest.h"

namespace gmx
{

namespace
{

using test::CommandLine;
using test::XvgMatch;

//! Test fixture for fca gauss run in double
typedef test::CommandLineTestBase FcaModuleTestDouble;
//! Test fixture for fca gauss run in float
typedef test::CommandLineTestBase FcaModuleTest;
#if GMX_DOUBLE
//! The test group when using double
#define FcaModuleTest_DEF FcaModuleTestDouble
#else
//! The test group when using float
#define FcaModuleTest_DEF FcaModuleTest
#endif

TEST_F(FcaModuleTest_DEF, ComputesFca)
{
    setInputFile("-pos", "fca_gauss.xvg");

    const double tolerance = 1e-4;
    XvgMatch     xvg;
    XvgMatch    &toler     = xvg.tolerance(test::relativeToleranceAsFloatingPoint(1, tolerance));
    setOutputFile("-mi", ".xvg", toler);
    setOutputFile("-val", ".xvg", toler);
    setOutputFile("-ic", ".xvg", toler);

    auto &cmdline = commandLine();
    ASSERT_EQ(0, gmx_fca(cmdline.argc(), cmdline.argv()));

    checkOutputFiles();
}

TEST_F(FcaModuleTest_DEF, ComputesFcaBasic)
{
    setInputFile("-pos", "fca_gauss.xvg");

    const double tolerance = 1e-4;
    XvgMatch     xvg;
    XvgMatch    &toler     = xvg.tolerance(test::relativeToleranceAsFloatingPoint(1, tolerance));
    setOutputFile("-mi", ".xvg", toler);
    setOutputFile("-val", ".xvg", toler);
    setOutputFile("-ic", ".xvg", toler);

    auto &cmdline = commandLine();
    cmdline.addOption("-basic");
    ASSERT_EQ(0, gmx_fca(cmdline.argc(), cmdline.argv()));

    checkOutputFiles();
}

TEST_F(FcaModuleTest_DEF, ComputesFcaBasicKappa)
{
    setInputFile("-pos", "fca_gauss.xvg");

    const double tolerance = 1e-1;
    XvgMatch     xvg;
    XvgMatch    &toler     = xvg.tolerance(test::relativeToleranceAsFloatingPoint(1, tolerance));
    setOutputFile("-mi", ".xvg", toler);
    setOutputFile("-val", ".xvg", toler);
    setOutputFile("-ic", ".xvg", toler);

    auto &cmdline = commandLine();
    cmdline.addOption("-basic");
    cmdline.addOption("-kappa1", "1.");
    cmdline.addOption("-kappa2", "1.");
    ASSERT_EQ(0, gmx_fca(cmdline.argc(), cmdline.argv()));

    checkOutputFiles();
}

TEST_F(FcaModuleTest_DEF, ComputesFcaBasicKappaSortMin)
{
    setInputFile("-pos", "fca_gauss.xvg");

    const double tolerance = 1e-1;
    XvgMatch     xvg;
    XvgMatch    &toler     = xvg.tolerance(test::relativeToleranceAsFloatingPoint(1, tolerance));
    setOutputFile("-mi", ".xvg", toler);
    setOutputFile("-val", ".xvg", toler);
    setOutputFile("-ic", ".xvg", toler);

    auto &cmdline = commandLine();
    cmdline.addOption("-basic");
    cmdline.addOption("-kappa1", "1.");
    cmdline.addOption("-kappa2", "1.");
    cmdline.addOption("-sort");
    cmdline.addOption("-minimize");
    ASSERT_EQ(0, gmx_fca(cmdline.argc(), cmdline.argv()));

    checkOutputFiles();
}

TEST_F(FcaModuleTest_DEF, ComputesFcaKappaSortMin)
{
    setInputFile("-pos", "fca_gauss.xvg");

    const double tolerance = 1e-4;
    XvgMatch     xvg;
    XvgMatch    &toler     = xvg.tolerance(test::relativeToleranceAsFloatingPoint(1, tolerance));
    setOutputFile("-mi", ".xvg", toler);
    setOutputFile("-val", ".xvg", toler);
    setOutputFile("-ic", ".xvg", toler);

    auto &cmdline = commandLine();
    cmdline.addOption("-sort");
    cmdline.addOption("-minimize");
    ASSERT_EQ(0, gmx_fca(cmdline.argc(), cmdline.argv()));

    checkOutputFiles();
}



} // namespace

} // gmx namespace
