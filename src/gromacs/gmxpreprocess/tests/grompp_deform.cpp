/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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
 * Tests for grompp with box deformation
 *
 * \ingroup module_gmxpreprocess
 */

#include "gmxpre.h"

#include <filesystem>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/fileio/tpxio.h"
#include "gromacs/gmxpreprocess/grompp.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textwriter.h"
#include "gromacs/utility/vectypes.h"

#include "testutils/cmdlinetest.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{
namespace
{

//! Test fixture for grompp with deform-init-flow, parametrized by the integrator
class GromppDeformInitFlowTest : public ::testing::TestWithParam<IntegrationAlgorithm>
{
};

//! Returns the name of a deform-init-flow test from its integrator
std::string nameOfDeformInitFlowTest(const ::testing::TestParamInfo<IntegrationAlgorithm>& info)
{
    // Test names may only contain alphanumeric characters and underscores
    return replaceAll(enumValueToString(info.param), "-", "");
}

TEST_P(GromppDeformInitFlowTest, AddsFlowProfileOnlyWithVelocities)
{
    const IntegrationAlgorithm integrator = GetParam();
    // Shear the box with a deformation velocity of b(x)
    const std::string mdpSettings =
            "nsteps = 10\n"
            "dt = 0.0005\n"
            "cutoff-scheme = Verlet\n"
            "coulombtype = Reaction-Field\n"
            "rcoulomb = 0.8\n"
            "rvdw = 0.8\n"
            "rlist = 0.8\n"
            "deform = 0 0 0 0.01 0 0\n"
            "deform-init-flow = yes\n";
    // Test-particle insertion needs a temperature and an insertion radius
    const std::string tpiSettings = EI_TPI(integrator) ? "rtpi = 0.05\n"
                                                         "tc-grps = System\n"
                                                         "tau-t = 0.1\n"
                                                         "ref-t = 300\n"
                                                       : "";
    const std::string mdpFileContents =
            formatString("integrator = %s\n", enumValueToString(integrator)) + mdpSettings + tpiSettings;

    TestFileManager fileManager;

    const std::filesystem::path mdpInputFileName = fileManager.getTemporaryFilePath("deform.mdp");
    TextWriter::writeFileFromString(mdpInputFileName, mdpFileContents);
    const std::filesystem::path tprFileName = fileManager.getTemporaryFilePath("deform.tpr");
    {
        SCOPED_TRACE("Calling grompp");
        CommandLine caller;
        caller.append("grompp");
        caller.addOption("-f", mdpInputFileName);
        caller.addOption("-p", TestFileManager::getInputFilePath("spc216-flexible.top").string());
        // This input file has no velocities
        caller.addOption("-c", TestFileManager::getInputFilePath("spc216.gro").string());
        caller.addOption("-o", tprFileName);
        // L-BFGS warns about the plain cut-off interactions used here
        const int maxWarningsTolerated = (integrator == IntegrationAlgorithm::LBFGS) ? 1 : 0;
        caller.addOption("-maxwarn", maxWarningsTolerated);
        ASSERT_EQ(0, gmx_grompp(caller.argc(), caller.argv()));
    }

    {
        SCOPED_TRACE("Checking output velocities");
        gmx_mtop_t top_after;
        t_inputrec ir_after;
        t_state    state_after;
        read_tpx_state(tprFileName, &ir_after, &state_after, &top_after);

        ASSERT_EQ(ir_after.eI, integrator);
        if (EI_DYNAMICS(integrator))
        {
            ASSERT_TRUE(state_after.hasEntry(StateEntry::V));
            ASSERT_EQ(state_after.v.size(), static_cast<size_t>(top_after.natoms));
            // The initial velocities are zero, so the velocities are the flow
            // profile of the shear in x along y
            ASSERT_NE(ir_after.deform[YY][XX], 0.0_real);
            const real coordToVelocity = ir_after.deform[YY][XX] / state_after.box[YY][YY];
            for (int i = 0; i < top_after.natoms; ++i)
            {
                EXPECT_REAL_EQ(coordToVelocity * state_after.x[i][YY], state_after.v[i][XX])
                        << "for atom " << i;
                EXPECT_EQ(state_after.v[i][YY], 0.0_real) << "for atom " << i;
                EXPECT_EQ(state_after.v[i][ZZ], 0.0_real) << "for atom " << i;
            }
        }
        else
        {
            EXPECT_FALSE(state_after.hasEntry(StateEntry::V))
                    << "the tpr has no velocities for a non-dynamical integrator";
            EXPECT_TRUE(state_after.v.empty());
        }
    }
}

INSTANTIATE_TEST_SUITE_P(Works,
                         GromppDeformInitFlowTest,
                         ::testing::Values(IntegrationAlgorithm::MD,
                                           IntegrationAlgorithm::Steep,
                                           IntegrationAlgorithm::CG,
                                           IntegrationAlgorithm::LBFGS,
                                           IntegrationAlgorithm::NM,
                                           IntegrationAlgorithm::TPI,
                                           IntegrationAlgorithm::TPIC),
                         nameOfDeformInitFlowTest);

} // namespace
} // namespace test
} // namespace gmx
