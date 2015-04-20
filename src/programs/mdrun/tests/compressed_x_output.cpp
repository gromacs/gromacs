/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
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
 * Tests for the mdrun -x functionality
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include <string>

#include <gtest/gtest.h>

#include "gromacs/options/filenameoption.h"
#include "gromacs/tools/check.h"

#include "testutils/cmdlinetest.h"

#include "moduletest.h"

namespace
{

//! Test fixture for mdrun -x
class CompressedXOutputTest : public gmx::test::MdrunTestFixture,
                              public testing::WithParamInterface<const char*>
{
};

/* Among other things, this test ensures mdrun can write a compressed trajectory. */
TEST_P(CompressedXOutputTest, ExitsNormally)
{
    std::string mdpFile("cutoff-scheme = Group\n"
                        "nsteps = 1\n"
                        "nstxout-compressed = 1\n");
    mdpFile += GetParam();
    runner_.useStringAsMdpFile(mdpFile.c_str());
    runner_.useTopGroAndNdxFromDatabase("spc2");
    ASSERT_EQ(0, runner_.callGrompp());

    runner_.reducedPrecisionTrajectoryFileName_ = fileManager_.getTemporaryFilePath(".xtc");
    ASSERT_EQ(0, runner_.callMdrun());

    ::gmx::test::CommandLine checkCaller;
    checkCaller.append("check");
    checkCaller.addOption("-f", runner_.reducedPrecisionTrajectoryFileName_);
    ASSERT_EQ(0, gmx_check(checkCaller.argc(), checkCaller.argv()));
}

#ifdef __INTEL_COMPILER
/* If we learn why this invocation triggers the "declared but not
   referenced" warning with ICC 12 on Windows and the case in
   gromacs/fft/tests/fft.cpp does not, then we can consider a more
   general solution. The pragma works also on MSVC, but so far is not
   required. */
#pragma warning( disable : 177 )
#endif

INSTANTIATE_TEST_CASE_P(WithDifferentMdpOptions, CompressedXOutputTest,
                            ::testing::Values
                            ( // Test writing the whole system via
                              // the default behaviour
                            "",

                            // Test writing the whole system
                            // explicitly
                            "compressed-x-grps = System\n",

                            // Test writing only part of the system.
                            // It would be nice to check that this test
                            // writes 3 atoms and the others write 6, but
                            // that's not yet easy.
                            "compressed-x-grps = SecondWaterMolecule\n"
                            ));

} // namespace
