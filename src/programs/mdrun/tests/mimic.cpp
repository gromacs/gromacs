/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
 * Tests for MiMiC forces computation
 *
 * \author Viacheslav Bolnykh <v.bolnykh@hpc-leap.eu>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include <filesystem>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/topology/ifunc.h"
#include "gromacs/trajectory/energyframe.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/simulationdatabase.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

#include "energycomparison.h"
#include "energyreader.h"
#include "moduletest.h"

namespace gmx
{
namespace test
{

//! Test fixture for bonded interactions
class MimicTest : public gmx::test::MdrunTestFixture
{
public:
    //! Execute the trajectory writing test
    void setupGrompp(const char* index_file, const char* top_file, const char* gro_file)
    {
        runner_.topFileName_ = TestFileManager::getInputFilePath(top_file).string();
        runner_.groFileName_ = TestFileManager::getInputFilePath(gro_file).string();
        runner_.ndxFileName_ = TestFileManager::getInputFilePath(index_file).string();
        runner_.useStringAsMdpFile(
                "integrator                = mimic\n"
                "QMMM-grps                 = QMatoms");
    }
    //! Prepare an mdrun caller
    CommandLine setupRerun()
    {
        CommandLine rerunCaller;
        rerunCaller.append("mdrun");
        rerunCaller.addOption("-rerun", runner_.groFileName_);
        runner_.edrFileName_ = fileManager_.getTemporaryFilePath(".edr").string();
        return rerunCaller;
    }
    //! Check the output of mdrun
    void checkRerun()
    {
        EnergyTermsToCompare energyTermsToCompare{ {
                { interaction_function[F_EPOT].longname, relativeToleranceAsFloatingPoint(-20.1, 1e-4) },
        } };

        TestReferenceData refData;
        auto              checker = refData.rootChecker();
        checkEnergiesAgainstReferenceData(runner_.edrFileName_, energyTermsToCompare, &checker);
    }
};

// This test checks if the energies produced with one quantum molecule are reasonable
TEST_F(MimicTest, OneQuantumMol)
{
    setupGrompp("1quantum.ndx", "4water.top", "4water.gro");
    ASSERT_EQ(0, runner_.callGrompp());

    test::CommandLine rerunCaller = setupRerun();

    ASSERT_EQ(0, runner_.callMdrun(rerunCaller));
    if (gmx_node_rank() == 0)
    {
        checkRerun();
    }
}

// This test checks if the energies produced with all quantum molecules are reasonable (0)
TEST_F(MimicTest, AllQuantumMol)
{
    setupGrompp("allquantum.ndx", "4water.top", "4water.gro");
    ASSERT_EQ(0, runner_.callGrompp());

    test::CommandLine rerunCaller = setupRerun();
    ASSERT_EQ(0, runner_.callMdrun(rerunCaller));
    if (gmx_node_rank() == 0)
    {
        checkRerun();
    }
}

// This test checks if the energies produced with two quantum molecules are reasonable
// Needed to check the LJ intermolecular exclusions
TEST_F(MimicTest, TwoQuantumMol)
{
    setupGrompp("2quantum.ndx", "4water.top", "4water.gro");
    ASSERT_EQ(0, runner_.callGrompp());

    test::CommandLine rerunCaller = setupRerun();
    ASSERT_EQ(0, runner_.callMdrun(rerunCaller));
    if (gmx_node_rank() == 0)
    {
        checkRerun();
    }
}

// This test checks if the energies produced with QM/MM boundary cutting the bond are ok
TEST_F(MimicTest, BondCuts)
{
    setupGrompp("ala.ndx", "ala.top", "ala.gro");
    ASSERT_EQ(0, runner_.callGrompp());

    test::CommandLine rerunCaller = setupRerun();
    ASSERT_EQ(0, runner_.callMdrun(rerunCaller));
    if (gmx_node_rank() == 0)
    {
        checkRerun();
    }
}

} // namespace test

} // namespace gmx
