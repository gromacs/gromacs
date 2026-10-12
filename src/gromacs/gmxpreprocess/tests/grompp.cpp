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
 * Tests for end-to-end grompp behaviour
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 */

#include "gmxpre.h"

#include "gromacs/gmxpreprocess/grompp.h"

#include <algorithm>
#include <filesystem>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/fileio/tpxio.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/cmdlinetest.h"
#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{
namespace
{

//! Test parameters
using MdpParams = bool;

class GromppSpecialParticlesVelocitiesTest : public ::testing::TestWithParam<MdpParams>
{
};

TEST_P(GromppSpecialParticlesVelocitiesTest, ZeroWhenExpected)
{
    // Need nstcalcenergy = 1 for shells
    const bool        generateVelocities       = GetParam();
    const std::string generateVelocitiesString = generateVelocities ? "yes" : "no";
    const std::string mdpFileContents =
            "nstcalcenergy = 1\ngen-seed = 2623\ngen-vel = " + generateVelocitiesString;

    TestFileManager fileManager;

    // Use an input file containing special particles
    const int                   numSpecialParticles = 4;
    const std::string           name                = "sw-dimer";
    const std::filesystem::path mdpInputFileName = fileManager.getTemporaryFilePath(name + ".mdp");
    TextWriter::writeFileFromString(mdpInputFileName, mdpFileContents);
    const std::filesystem::path tprFileName = fileManager.getTemporaryFilePath(name + ".tpr");
    {
        SCOPED_TRACE("Calling grompp");
        CommandLine caller;
        caller.append("grompp");
        caller.addOption("-f", mdpInputFileName);
        caller.addOption("-p", TestFileManager::getInputFilePath(name + ".top").string());
        // This input file has non-zero velocities for shells and virtual sites
        caller.addOption("-c", TestFileManager::getInputFilePath(name + ".g96").string());
        caller.addOption("-o", tprFileName);
        EXPECT_EQ(0, gmx_grompp(caller.argc(), caller.argv()));
    }

    {
        SCOPED_TRACE("Checking output velocities");
        gmx_mtop_t top_after;
        t_inputrec ir_after;
        t_state    state_after;
        read_tpx_state(tprFileName, &ir_after, &state_after, &top_after);

        int numParticlesWithZeroedVelocityComponents = 0;
        for (const auto& v : state_after.v)
        {
            // Normally it is not advisable to compare floating-point
            // values for exact equality, but here we know the value
            // is the result of an assignment, not a computation.
            if (v[XX] == 0.0_real && v[YY] == 0.0_real && v[ZZ] == 0.0_real)
            {
                ++numParticlesWithZeroedVelocityComponents;
            }
        }
        EXPECT_EQ(numParticlesWithZeroedVelocityComponents, generateVelocities ? numSpecialParticles : 0)
                << "the tpr has zeroed velocity components of special particles only when "
                   "generating velocities";
    }
}

INSTANTIATE_TEST_SUITE_P(Works, GromppSpecialParticlesVelocitiesTest, testing::Bool());

//! Test fixture for grompp with frozen particles, parametrized by the integrator
class GromppFrozenParticlesTest : public ::testing::TestWithParam<IntegrationAlgorithm>
{
};

TEST_P(GromppFrozenParticlesTest, VelocitiesAreHandled)
{
    const IntegrationAlgorithm integrator = GetParam();
    // Freeze all dimensions of the methanol molecule, but not the water
    const std::string mdpSettings =
            "nsteps = 10\n"
            "cutoff-scheme = Verlet\n"
            "coulombtype = Reaction-Field\n"
            "rcoulomb = 0.9\n"
            "rvdw = 0.9\n"
            "freezegrps = Methanol\n"
            "freezedim = Y Y Y\n";
    // Test-particle insertion needs a temperature and an insertion radius
    const std::string tpiSettings = EI_TPI(integrator) ? "rtpi = 0.05\n"
                                                         "tc-grps = System\n"
                                                         "tau-t = 0.1\n"
                                                         "ref-t = 300\n"
                                                       : "";
    const std::string mdpFileContents =
            formatString("integrator = %s\n", enumValueToString(integrator)) + mdpSettings + tpiSettings;

    TestFileManager fileManager;

    // The first three atoms, the methanol molecule, will be frozen
    const int                   numFrozenParticles = 3;
    const std::string           name               = "spc-and-methanol";
    const std::filesystem::path mdpInputFileName = fileManager.getTemporaryFilePath(name + ".mdp");
    TextWriter::writeFileFromString(mdpInputFileName, mdpFileContents);
    const std::filesystem::path tprFileName = fileManager.getTemporaryFilePath(name + ".tpr");
    {
        SCOPED_TRACE("Calling grompp");
        CommandLine caller;
        caller.append("grompp");
        caller.addOption("-f", mdpInputFileName);
        caller.addOption("-p", TestFileManager::getInputFilePath(name + ".top").string());
        // This input file has non-zero velocities for all particles
        caller.addOption("-c", TestFileManager::getInputFilePath(name + ".gro").string());
        caller.addOption("-n", TestFileManager::getInputFilePath(name + ".ndx").string());
        caller.addOption("-o", tprFileName);
        // L-BFGS warns about the plain cut-off interactions used here
        const int maxWarningsTolerated = (integrator == IntegrationAlgorithm::LBFGS) ? 1 : 0;
        caller.addOption("-maxwarn", maxWarningsTolerated);
        ASSERT_EQ(0, gmx_grompp(caller.argc(), caller.argv()));
    }

    {
        SCOPED_TRACE("Checking frozen dimensions and output velocities");
        gmx_mtop_t top_after;
        t_inputrec ir_after;
        t_state    state_after;
        read_tpx_state(tprFileName, &ir_after, &state_after, &top_after);

        ASSERT_EQ(ir_after.eI, integrator);
        for (int i = 0; i < top_after.natoms; ++i)
        {
            const bool expectFrozen = (i < numFrozenParticles);
            const int freezeGroup = getGroupType(top_after.groups, SimulationAtomGroupType::Freeze, i);
            for (int d = 0; d < DIM; ++d)
            {
                EXPECT_EQ(ir_after.opts.nFreeze[freezeGroup][d] != 0, expectFrozen)
                        << "for atom " << i << " and dimension " << d;
            }
        }

        if (EI_DYNAMICS(integrator))
        {
            ASSERT_TRUE(state_after.hasEntry(StateEntry::V));
            ASSERT_EQ(state_after.v.size(), static_cast<size_t>(top_after.natoms));
            for (int i = 0; i < top_after.natoms; ++i)
            {
                const bool expectFrozen = (i < numFrozenParticles);
                // Normally it is not advisable to compare floating-point
                // values for exact equality, but here we know the value
                // is either the result of an assignment or read from
                // the input file.
                for (int d = 0; d < DIM; ++d)
                {
                    EXPECT_EQ(state_after.v[i][d] == 0.0_real, expectFrozen)
                            << "the tpr has zeroed velocity components only for frozen "
                               "dimensions, for atom "
                            << i << " and dimension " << d;
                }
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

//! Returns a test name for the integrator that is a valid GoogleTest identifier
std::string integratorTestName(const testing::TestParamInfo<IntegrationAlgorithm>& info)
{
    std::string name = enumValueToString(info.param);
    std::replace(name.begin(), name.end(), '-', '_');
    return name;
}

INSTANTIATE_TEST_SUITE_P(Works,
                         GromppFrozenParticlesTest,
                         testing::Values(IntegrationAlgorithm::MD,
                                         IntegrationAlgorithm::Steep,
                                         IntegrationAlgorithm::CG,
                                         IntegrationAlgorithm::LBFGS,
                                         IntegrationAlgorithm::NM,
                                         IntegrationAlgorithm::TPI,
                                         IntegrationAlgorithm::TPIC),
                         integratorTestName);

} // namespace
} // namespace test
} // namespace gmx
