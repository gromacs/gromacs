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
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trrio.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/mtop_atomloops.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
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

//! Settings for MD of spc-and-methanol with a frozen methanol molecule
const char c_frozenMethanolMdpSettings[] =
        "integrator = md\n"
        "nsteps = 10\n"
        "cutoff-scheme = Verlet\n"
        "coulombtype = Reaction-Field\n"
        "rcoulomb = 0.9\n"
        "rvdw = 0.9\n"
        "freezegrps = Methanol\n"
        "freezedim = Y Y Y\n";

/*! \brief Checks the velocities in \p tprFileName
 *
 * Checks that the velocities are zero for the fully frozen first
 * \p numFrozenParticles atoms and non-zero for all others. When
 * \p expectedVelocities is not empty, the velocities of the atoms that
 * are not frozen must be equal to it. */
void checkVelocitiesOfFrozenParticles(const std::filesystem::path& tprFileName,
                                      const int                    numFrozenParticles,
                                      ArrayRef<const RVec>         expectedVelocities)
{
    gmx_mtop_t top_after;
    t_inputrec ir_after;
    t_state    state_after;
    read_tpx_state(tprFileName, &ir_after, &state_after, &top_after);

    ASSERT_TRUE(state_after.hasEntry(StateEntry::V));
    ASSERT_EQ(state_after.v.size(), static_cast<size_t>(top_after.natoms));
    for (int i = 0; i < top_after.natoms; ++i)
    {
        const bool expectFrozen = (i < numFrozenParticles);
        const int  freezeGroup = getGroupType(top_after.groups, SimulationAtomGroupType::Freeze, i);
        for (int d = 0; d < DIM; ++d)
        {
            ASSERT_EQ(ir_after.opts.nFreeze[freezeGroup][d] != 0, expectFrozen)
                    << "for atom " << i << " and dimension " << d;
            // Normally it is not advisable to compare floating-point
            // values for exact equality, but here we know the value
            // is either the result of an assignment or read from a
            // file.
            EXPECT_EQ(state_after.v[i][d] == 0.0_real, expectFrozen)
                    << "the tpr has zeroed velocity components only for frozen "
                       "dimensions, for atom "
                    << i << " and dimension " << d;
            if (!expectFrozen && !expectedVelocities.empty())
            {
                EXPECT_EQ(state_after.v[i][d], expectedVelocities[i][d])
                        << "for atom " << i << " and dimension " << d;
            }
        }
    }
}

//! The velocities of the trajectory frame of writeSpcAndMethanolFrame(), non-zero for all particles
const std::vector<RVec> c_spcAndMethanolFrameVelocities = {
    { 0.1_real, -0.2_real, 0.3_real }, { -0.4_real, 0.5_real, -0.6_real },
    { 0.7_real, -0.8_real, 0.9_real }, { -1.0_real, 1.1_real, -1.2_real },
    { 1.3_real, -1.4_real, 1.5_real }, { -1.6_real, 1.7_real, -1.8_real }
};

/*! \brief Writes a trajectory frame of spc-and-methanol to \p trrFileName
 *
 * The frame has the coordinates and box of the configuration file, but
 * the velocities c_spcAndMethanolFrameVelocities for all particles,
 * including the ones that will be frozen. */
void writeSpcAndMethanolFrame(const std::filesystem::path& trrFileName)
{
    const std::vector<RVec> x = {
        { 1.970_real, 1.460_real, 1.209_real }, { 1.978_real, 1.415_real, 1.082_real },
        { 1.905_real, 1.460_real, 1.030_real }, { 1.555_real, 1.511_real, 0.703_real },
        { 1.498_real, 1.495_real, 0.784_real }, { 1.496_real, 1.521_real, 0.623_real }
    };
    const matrix box = { { 3.01_real, 0, 0 }, { 0, 3.01_real, 0 }, { 0, 0, 3.01_real } };
    gmx_trr_write_single_frame(trrFileName,
                               0,
                               0,
                               0,
                               box,
                               static_cast<int>(x.size()),
                               as_rvec_array(x.data()),
                               as_rvec_array(c_spcAndMethanolFrameVelocities.data()),
                               nullptr);
}

TEST(GromppFrozenParticlesChangedVelocitiesTest, VelocitiesFromTrajectoryFrameAreCleared)
{
    TestFileManager fileManager;

    // The first three atoms, the methanol molecule, will be frozen
    const int                   numFrozenParticles = 3;
    const std::string           name               = "spc-and-methanol";
    const std::filesystem::path mdpInputFileName = fileManager.getTemporaryFilePath(name + ".mdp");
    TextWriter::writeFileFromString(mdpInputFileName, c_frozenMethanolMdpSettings);

    const std::filesystem::path trrFileName = fileManager.getTemporaryFilePath(name + ".trr");
    writeSpcAndMethanolFrame(trrFileName);

    const std::filesystem::path tprFileName = fileManager.getTemporaryFilePath(name + ".tpr");
    {
        SCOPED_TRACE("Calling grompp");
        CommandLine caller;
        caller.append("grompp");
        caller.addOption("-f", mdpInputFileName);
        caller.addOption("-p", TestFileManager::getInputFilePath(name + ".top").string());
        caller.addOption("-c", TestFileManager::getInputFilePath(name + ".gro").string());
        caller.addOption("-n", TestFileManager::getInputFilePath(name + ".ndx").string());
        // This frame has non-zero velocities for all particles
        caller.addOption("-t", trrFileName);
        caller.addOption("-o", tprFileName);
        ASSERT_EQ(0, gmx_grompp(caller.argc(), caller.argv()));
    }

    {
        SCOPED_TRACE("Checking output velocities");
        checkVelocitiesOfFrozenParticles(tprFileName, numFrozenParticles, c_spcAndMethanolFrameVelocities);
    }
}

TEST(GromppFrozenParticlesChangedVelocitiesTest, VelocitiesOfDeformFlowProfileAreCleared)
{
    // Deform the box so that the flow profile adds velocity in x that
    // depends on the z coordinate, which is non-zero for all atoms
    const std::string mdpFileContents = std::string(c_frozenMethanolMdpSettings)
                                        + "deform = 0 0 0 0 0.01 0\n"
                                          "deform-init-flow = yes\n";

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
        ASSERT_EQ(0, gmx_grompp(caller.argc(), caller.argv()));
    }

    {
        SCOPED_TRACE("Checking output velocities");
        checkVelocitiesOfFrozenParticles(tprFileName, numFrozenParticles, {});
    }

    {
        SCOPED_TRACE("Checking that the flow profile was added");
        // Read the coordinates, velocities and box that grompp used
        // from the configuration file
        gmx_mtop_t topGro;
        bool       haveTopology = false;
        rvec*      xGro         = nullptr;
        rvec*      vGro         = nullptr;
        matrix     boxGro;
        readConfAndTopology(
                TestFileManager::getInputFilePath(name + ".gro"), &haveTopology, &topGro, nullptr, &xGro, &vGro, boxGro);

        gmx_mtop_t top_after;
        t_inputrec ir_after;
        t_state    state_after;
        read_tpx_state(tprFileName, &ir_after, &state_after, &top_after);
        ASSERT_EQ(state_after.v.size(), static_cast<size_t>(topGro.natoms));

        // The flow profile adds the deformation rate of the x
        // component of the z box vector, relative to the box length
        // in z, times the z coordinate to the x component of the
        // velocity
        const real flowPerZ = ir_after.deform[ZZ][XX] / boxGro[ZZ][ZZ];
        ASSERT_NE(flowPerZ, 0.0_real);
        for (int i = 0; i < topGro.natoms; ++i)
        {
            for (int d = 0; d < DIM; ++d)
            {
                if (i < numFrozenParticles)
                {
                    EXPECT_EQ(state_after.v[i][d], 0.0_real)
                            << "for frozen atom " << i << " and dimension " << d;
                }
                else if (d == XX)
                {
                    EXPECT_REAL_EQ_TOL(vGro[i][d] + flowPerZ * xGro[i][ZZ],
                                       state_after.v[i][d],
                                       relativeToleranceAsFloatingPoint(vGro[i][d], 10 * GMX_REAL_EPS))
                            << "for atom " << i << " and dimension " << d;
                }
                else
                {
                    EXPECT_EQ(state_after.v[i][d], vGro[i][d])
                            << "for atom " << i << " and dimension " << d;
                }
            }
        }
        sfree(xGro);
        sfree(vGro);
    }
}

TEST(GromppFrozenParticlesChangedVelocitiesTest, VelocitiesFromTrajectoryFrameWithDeformFlowAreCleared)
{
    // With -t, the flow profile is added to the generated velocities
    // before the coordinates of the trajectory frame replace those of
    // the configuration file (#5738), so only the frozen components
    // are checked here
    const std::string mdpFileContents = std::string(c_frozenMethanolMdpSettings)
                                        + "gen-vel = yes\n"
                                          "gen-temp = 300\n"
                                          "gen-seed = 2623\n"
                                          "deform = 0 0 0 0 0.01 0\n"
                                          "deform-init-flow = yes\n";

    TestFileManager fileManager;

    // The first three atoms, the methanol molecule, will be frozen
    const int                   numFrozenParticles = 3;
    const std::string           name               = "spc-and-methanol";
    const std::filesystem::path mdpInputFileName = fileManager.getTemporaryFilePath(name + ".mdp");
    TextWriter::writeFileFromString(mdpInputFileName, mdpFileContents);

    const std::filesystem::path trrFileName = fileManager.getTemporaryFilePath(name + ".trr");
    writeSpcAndMethanolFrame(trrFileName);

    const std::filesystem::path tprFileName = fileManager.getTemporaryFilePath(name + ".tpr");
    {
        SCOPED_TRACE("Calling grompp");
        CommandLine caller;
        caller.append("grompp");
        caller.addOption("-f", mdpInputFileName);
        caller.addOption("-p", TestFileManager::getInputFilePath(name + ".top").string());
        caller.addOption("-c", TestFileManager::getInputFilePath(name + ".gro").string());
        caller.addOption("-n", TestFileManager::getInputFilePath(name + ".ndx").string());
        caller.addOption("-t", trrFileName);
        caller.addOption("-o", tprFileName);
        ASSERT_EQ(0, gmx_grompp(caller.argc(), caller.argv()));
    }

    {
        SCOPED_TRACE("Checking output velocities");
        checkVelocitiesOfFrozenParticles(tprFileName, numFrozenParticles, {});
    }
}

//! Test fixture for the flow profile of special particles, parametrized by whether grompp generates the velocities
class GromppDeformFlowSpecialParticlesTest : public ::testing::TestWithParam<bool>
{
};

TEST_P(GromppDeformFlowSpecialParticlesTest, ShellsAndVirtualSitesGetTheFlowProfile)
{
    // Need nstcalcenergy = 1 for shells. Generated velocities at 0 K
    // are zero, so then all particles have only the flow profile.
    const bool        generateVelocities = GetParam();
    const std::string mdpFileContents =
            std::string(
                    "nstcalcenergy = 1\n"
                    "deform = 0 0 0 0 0.01 0\n"
                    "deform-init-flow = yes\n")
            + (generateVelocities ? "gen-vel = yes\ngen-temp = 0\ngen-seed = 2623\n" : "gen-vel = no\n");

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
        // grompp fails on any warning with -maxwarn 0
        caller.addOption("-maxwarn", 0);
        ASSERT_EQ(0, gmx_grompp(caller.argc(), caller.argv()));
    }

    // Read the coordinates, velocities and box that grompp used from
    // the configuration file
    gmx_mtop_t topConf;
    bool       haveTopology = false;
    rvec*      xConf        = nullptr;
    rvec*      vConf        = nullptr;
    matrix     boxConf;
    readConfAndTopology(
            TestFileManager::getInputFilePath(name + ".g96"), &haveTopology, &topConf, nullptr, &xConf, &vConf, boxConf);

    gmx_mtop_t top_after;
    t_inputrec ir_after;
    t_state    state_after;
    read_tpx_state(tprFileName, &ir_after, &state_after, &top_after);
    ASSERT_EQ(state_after.v.size(), static_cast<size_t>(topConf.natoms));

    int numSpecialParticlesFound = 0;
    for (const AtomProxy atomP : AtomRange(top_after))
    {
        const ParticleType ptype = atomP.atom().ptype;
        numSpecialParticlesFound +=
                (ptype == ParticleType::Shell || ptype == ParticleType::VSite) ? 1 : 0;
    }
    ASSERT_EQ(numSpecialParticlesFound, numSpecialParticles);

    // The flow profile adds the deformation rate of the x component of
    // the z box vector, relative to the box length in z, times the z
    // coordinate to the x component of the velocity of every particle,
    // shells and virtual sites included
    const real flowPerZ = ir_after.deform[ZZ][XX] / boxConf[ZZ][ZZ];
    ASSERT_NE(flowPerZ, 0.0_real);
    for (int i = 0; i < topConf.natoms; ++i)
    {
        for (int d = 0; d < DIM; ++d)
        {
            const real inputVelocity = generateVelocities ? 0.0_real : vConf[i][d];
            if (d == XX)
            {
                const real expected = inputVelocity + flowPerZ * xConf[i][ZZ];
                EXPECT_REAL_EQ_TOL(expected,
                                   state_after.v[i][d],
                                   relativeToleranceAsFloatingPoint(expected, 10 * GMX_REAL_EPS))
                        << "for atom " << i << " and dimension " << d;
            }
            else
            {
                EXPECT_EQ(state_after.v[i][d], inputVelocity)
                        << "for atom " << i << " and dimension " << d;
            }
        }
    }
    sfree(xConf);
    sfree(vConf);
}

INSTANTIATE_TEST_SUITE_P(WithAndWithoutGeneratedVelocities,
                         GromppDeformFlowSpecialParticlesTest,
                         testing::Bool());

} // namespace
} // namespace test
} // namespace gmx
