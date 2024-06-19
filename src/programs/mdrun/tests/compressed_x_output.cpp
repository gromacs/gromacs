/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2014- The GROMACS Authors
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
 * Tests for the mdrun -x functionality
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include <filesystem>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/options/filenameoption.h"
#include "gromacs/tools/check.h"

#include "testutils/cmdlinetest.h"
#include "testutils/testfilemanager.h"

#include "moduletest.h"

namespace
{

//! Test fixture for mdrun -x
class MdrunCompressedXOutputTest :
    public gmx::test::MdrunTestFixture,
    public testing::WithParamInterface<const char*>
{
};

//! Helper typedef for naming test cases like sentences
typedef MdrunCompressedXOutputTest MdrunCompressedXOutput;

/* This test tests a few ways that mdrun can write a compressed trajectory. */
TEST_P(MdrunCompressedXOutput, ExitsNormally)
{
    std::string mdpFile(R"(cutoff-scheme = Verlet
                           verlet-buffer-tolerance = 0.005
                           nsteps = 1
                           nstxout-compressed = 1
                           )");
    std::string compressedXGrpsLine = GetParam();
    mdpFile += compressedXGrpsLine;
    runner_.useStringAsMdpFile(mdpFile.c_str());
    runner_.useTopGroAndNdxFromDatabase("spc2");
    ASSERT_EQ(0, runner_.callGrompp());

    runner_.reducedPrecisionTrajectoryFileName_ = fileManager_.getTemporaryFilePath(".xtc").string();
    ASSERT_EQ(0, runner_.callMdrun());

    ::gmx::test::CommandLine checkCaller;
    checkCaller.append("check");
    checkCaller.addOption("-f", runner_.reducedPrecisionTrajectoryFileName_);
    ASSERT_EQ(0, gmx_check(checkCaller.argc(), checkCaller.argv()));
}

INSTANTIATE_TEST_SUITE_P(WithDifferentOutputGroupSettings,
                         MdrunCompressedXOutput,
                         ::testing::Values( // Test writing the whole system via
                                            // the default behaviour
                                 "",

                                 // Test writing the whole system
                                 // explicitly
                                 "compressed-x-grps = System\n",

                                 // Test writing only part of the system.
                                 // It would be nice to check that this test
                                 // writes 3 atoms and the others write 6, but
                                 // that's not yet easy.
                                 "compressed-x-grps = SecondWaterMolecule\n"));

} // namespace
