/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
 * Tests for functionality of the NNPot ForceProvider
 *
 * \author Lukas Müllender <lukas.muellender@gmail.com>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "gromacs/applied_forces/nnpot/nnpotforceprovider.h"

#include <gtest/gtest.h>

#include "gromacs/applied_forces/nnpot/nnpotoptions.h"
#include "gromacs/domdec/localatomset.h"
#include "gromacs/domdec/localatomsetmanager.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/gmxpreprocess/grompp.h"
#include "gromacs/hardware/device_information.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/setenv.h"
#include "testutils/test_hardware_environment.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace gmx
{

namespace test
{

class NNPotForceProviderTest : public ::testing::Test
{
public:
    /*! \brief Generates tpr file from *.top and *.gro existing in the simulation database directory
     * and loads gmx_mtop_t from it
     */
    void readDataFromFile(const std::string& fileName, const std::string& mdpContent)
    {
        const std::filesystem::path simData =
                gmx::test::TestFileManager::getTestSimulationDatabaseDirectory();

        // Generate empty mdp file
        const std::string mdpInputFileName =
                fileManager_.getTemporaryFilePath(fileName + ".mdp").string();
        gmx::TextWriter::writeFileFromString(mdpInputFileName, mdpContent);

        // Generate tpr file
        const std::string tprName = fileManager_.getTemporaryFilePath(fileName + ".tpr").string();
        {
            gmx::test::CommandLine caller;
            caller.append("grompp");
            caller.addOption("-f", mdpInputFileName);
            caller.addOption("-p", (simData / fileName).replace_extension(".top").string());
            caller.addOption("-c", (simData / fileName).replace_extension(".gro").string());
            caller.addOption("-o", tprName);
            ASSERT_EQ(0, gmx_grompp(caller.argc(), caller.argv()));
        }

        // Load topology
        bool fullTopology;

        readConfAndTopology(tprName.c_str(), &fullTopology, &mtop_, &pbcType_, &coords, nullptr, box_);

        x_ = gmx::constArrayRefFromArray(reinterpret_cast<gmx::RVec*>(coords), mtop_.natoms);
    }

    void setDefaultParameters()
    {
        params_.active_        = true;
        params_.modelFileName_ = gmx::test::TestFileManager::getInputFilePath("model.pt").string();
        params_.numAtoms_      = 6;
        params_.cr_            = &cr_;
    }

    void setParametersTopology()
    {
        params_.active_        = true;
        params_.modelFileName_ = gmx::test::TestFileManager::getInputFilePath("model.pt").string();
        params_.atoms_         = gmx_mtop_global_atoms(mtop_);
        params_.numAtoms_      = params_.atoms_.nr;
        params_.cr_            = &cr_;
        params_.modelInput_    = { "atom-positions", "atom-numbers", "box", "pbc" };

        std::vector<gmx::Index> inpIndices = { 0, 1, 2, 3, 4, 5 };
        LocalAtomSet            set = atomSetManager_.add(ArrayRef<const gmx::Index>(inpIndices));
        params_.inpAtoms_           = std::make_unique<LocalAtomSet>(set);
        params_.pbcType_            = std::make_unique<PbcType>(pbcType_);
    }

    void testCalculateForces(gmx::test::TestReferenceData& testData)
    {
        gmx::test::TestReferenceChecker checker(testData.rootChecker());

        std::unique_ptr<NNPotForceProvider> nnpotForceProvider;
        EXPECT_NO_THROW(nnpotForceProvider = std::make_unique<NNPotForceProvider>(params_, &logger_));
        EXPECT_NO_THROW(nnpotForceProvider->gatherAtomNumbersIndices());

        // Prepare input for force provider
        ForceProviderInput fInput(x_, params_.numAtoms_, {}, {}, 0.0, 0, box_, cr_);

        // Prepare output for force provider
        std::vector<RVec>   forces(params_.numAtoms_, RVec{ 0, 0, 0 });
        ForceWithVirial     forceWithVirial(forces, true);
        gmx_enerdata_t      enerdDummy(1, nullptr);
        ForceProviderOutput forceProviderOutput(&forceWithVirial, &enerdDummy);

        EXPECT_NO_THROW(nnpotForceProvider->calculateForces(fInput, &forceProviderOutput));

        checker.setDefaultTolerance(gmx::test::relativeToleranceAsFloatingPoint(100000.0, 5e-5));
        checker.checkReal(enerdDummy.term[F_ENNPOT], "Energy");
        checker.setDefaultTolerance(gmx::test::relativeToleranceAsFloatingPoint(100.0, 5e-5));
        checker.checkSequence(forces.begin(), forces.end(), "Forces");
    }

protected:
    NNPotParameters               params_;
    MDLogger                      logger_;
    LocalAtomSetManager           atomSetManager_;
    std::unique_ptr<LocalAtomSet> inpAtomSet_;
    t_commrec                     cr_;
    gmx::test::TestFileManager    fileManager_;

    rvec*                coords = nullptr;
    ArrayRef<const RVec> x_;
    gmx_mtop_t           mtop_;
    PbcType              pbcType_;
    matrix               box_;
};

TEST_F(NNPotForceProviderTest, CanConstruct)
{
    setDefaultParameters();

    // GMX_TORCH is defined by set_source_files_properties() in CMakeLists.txt
    if (GMX_TORCH)
    {
        {
            SCOPED_TRACE("Check construction on CPU");
            EXPECT_NO_THROW(NNPotForceProvider nnpotForceProvider(params_, &logger_));
        }

        {
            SCOPED_TRACE("Check construction with invalid model file name");
            params_.modelFileName_ = "model";
            EXPECT_THROW_GMX(NNPotForceProvider nnpotForceProvider(params_, &logger_), FileIOError);
            params_.modelFileName_ = gmx::test::TestFileManager::getInputFilePath("model.pt").string();
        }

        for (const auto& device : getTestHardwareEnvironment()->getTestDeviceList())
        {
            const DeviceInformation deviceInfo = device->deviceInfo();
            if (deviceInfo.deviceVendor == DeviceVendor::Nvidia)
            {
                SCOPED_TRACE("Check construction on default NVIDIA GPU");
                gmxSetenv("GMX_NN_DEVICE", "cuda", 1);
                EXPECT_NO_THROW(NNPotForceProvider nnpotForceProvider(params_, &logger_));
                gmxUnsetenv("GMX_NN_DEVICE");
                break; // Only test one GPU until we have a better way to hande device selection
            }
        }
    }
    else
    {
        EXPECT_ANY_THROW(NNPotForceProvider nnpotForceProvider(params_, &logger_));
    }
}

#if GMX_TORCH
/*! Test if the NNPotForceProvider can calculate forces.
 * This implicitly tests the NNPotModel as well as the input gathering functions.
 */
TEST_F(NNPotForceProviderTest, CanCalculateForces)
{
    gmx::test::TestReferenceData data;

    readDataFromFile("spc2", "");
    setParametersTopology();

    {
        SCOPED_TRACE("Check calculate forces on CPU");
        testCalculateForces(data);
    }
    if (!getTestHardwareEnvironment()->hasCompatibleDevices())
    {
        return;
    }
    for (const auto& device : getTestHardwareEnvironment()->getTestDeviceList())
    {
        const DeviceInformation deviceInfo = device->deviceInfo();
        if (deviceInfo.deviceVendor == DeviceVendor::Nvidia)
        {
            SCOPED_TRACE("Check calculate forces on default NVIDIA GPU");
            gmxSetenv("GMX_NN_DEVICE", "cuda", 1);
            testCalculateForces(data);
            gmxUnsetenv("GMX_NN_DEVICE");
            break; // Only test one GPU until we have a better way to hande device selection
        }
    }

    done_atom(&params_.atoms_);
}
#endif // GMX_TORCH

} // namespace test

} // namespace gmx
