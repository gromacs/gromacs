/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2023- The GROMACS Authors
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
 * Tests for ColvarsForceProvider class of Colvars MDModule.
 *
 * \author Hubert Santuz <hubert.santuz@gmail.com>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "gromacs/applied_forces/colvars/colvarsforceprovider.h"

#include <array>
#include <filesystem>
#include <map>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/applied_forces/colvars/colvarssimulationsparameters.h"
#include "gromacs/domdec/localatomsetmanager.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/gmxpreprocess/grompp.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/forcerec.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/textreader.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

enum class PbcType : int;

namespace gmx
{

class ColvarsForceProviderTest : public ::testing::Test
{
public:
    void PrepareInputForceProvider(const std::string& fileName)
    {

        gmx::test::TestFileManager  fileManager_;
        const std::filesystem::path simData =
                gmx::test::TestFileManager::getTestSimulationDatabaseDirectory();

        // Generate empty mdp file
        const std::string mdpInputFileName =
                fileManager_.getTemporaryFilePath(fileName + ".mdp").string();
        gmx::TextWriter::writeFileFromString(mdpInputFileName, "");

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


        bool fullTopology;

        // Load topology
        readConfAndTopology(tprName.c_str(), &fullTopology, &mtop, &pbcType_, &coords, nullptr, box_);

        x_ = gmx::constArrayRefFromArray(reinterpret_cast<gmx::RVec*>(coords), mtop.natoms);
    }

    void ColvarsConfigStringFromFile(const std::string& filename)
    {
        // Path to the sample colvars input file
        std::filesystem::path colvarsInputFile = gmx::test::TestFileManager::getInputFilePath(filename);

        colvarsConfigString_ = TextReader::readFileToString(colvarsInputFile);
    }

    void CorrectColvarsConfigString()
    {
        // atomNumbers start at 1 in colvars
        colvarsConfigString_ = R"(units gromacs
              colvar {
                  name d_atoms
                  distance {
                      group1 {
                          atomNumbers 1
                      }
                      group2 {
                          atomNumbers 5
                      }
                  }
              }
              harmonic {
                  colvars d_atoms
                  forceConstant 20000
                  centers 0.3
              })";
    }

    void IncorrectColvarsConfigString()
    {
        // Path to the sample colvars input file
        std::filesystem::path colvarsInputFile =
                gmx::test::TestFileManager::getInputFilePath("colvars_sample.dat");

        colvarsConfigString_ = TextReader::readFileToString(colvarsInputFile);
    }

    ~ColvarsForceProviderTest() override
    {
        if (coords)
        {
            sfree(coords);
        }
    }

protected:
    std::string                        colvarsConfigString_;
    std::vector<RVec>                  atomCoords_;
    LocalAtomSetManager                atomSetManager_;
    PbcType                            pbcType_;
    MDLogger                           logger_;
    t_atoms                            atoms_;
    t_commrec                          cr_;
    std::map<std::string, std::string> KVTInputs;
    ColvarsForceProviderState          colvarsState_;

    double      simulationTimeStep_ = 0.002;
    real        temperature_        = 300;
    int         seed_               = 123456;
    std::string prefixOutput_;

    rvec*                coords = nullptr;
    ArrayRef<const RVec> x_;
    gmx_mtop_t           mtop;
    matrix               box_;
};

TEST_F(ColvarsForceProviderTest, CanConstructOrNot)
{

    EXPECT_NO_THROW(ColvarsForceProvider forceProvider(colvarsConfigString_,
                                                       atoms_,
                                                       pbcType_,
                                                       &logger_,
                                                       KVTInputs,
                                                       temperature_,
                                                       seed_,
                                                       &atomSetManager_,
                                                       &cr_,
                                                       simulationTimeStep_,
                                                       atomCoords_,
                                                       prefixOutput_,
                                                       colvarsState_));
}

TEST_F(ColvarsForceProviderTest, SimpleInputs)
{
    PrepareInputForceProvider("4water");
    ColvarsConfigStringFromFile("colvars_sample.dat");
    atoms_ = gmx_mtop_global_atoms(mtop);

    // Indexes taken from the Colvars Config file.
    atomCoords_ = { RVec(x_[0]), RVec(x_[4]) };


    ColvarsForceProvider forceProvider(colvarsConfigString_,
                                       atoms_,
                                       pbcType_,
                                       &logger_,
                                       KVTInputs,
                                       temperature_,
                                       seed_,
                                       &atomSetManager_,
                                       &cr_,
                                       simulationTimeStep_,
                                       atomCoords_,
                                       prefixOutput_,
                                       colvarsState_);


    // Re-use the PreProcessorTest since the ForceProvider recalls colvars initilization and the input are identicals.
    gmx::test::TestReferenceData    data("ColvarsPreProcessorTest_CheckValuesFourWaters.xml");
    gmx::test::TestReferenceChecker checker(data.rootChecker());

    // Check colvars & atoms values are correctly read
    checker.setDefaultTolerance(gmx::test::absoluteTolerance(0.001));
    checker.checkVector(atomCoords_[1], "Coords Atom 4");

    const auto* const atomIds = forceProvider.get_atom_ids();
    checker.checkSequence(atomIds->begin(), atomIds->end(), "Index of colvars atoms");

    const auto* const atomMasses = forceProvider.get_atom_masses();
    checker.checkSequence(atomMasses->begin(), atomMasses->end(), "Masses of colvars atoms");

    const auto* const atomCharges = forceProvider.get_atom_charges();
    checker.checkSequence(atomCharges->begin(), atomCharges->end(), "Charges of colvars atoms");

    done_atom(&atoms_);
}

#if defined(__has_feature)
#    if !__has_feature(address_sanitizer)
TEST_F(ColvarsForceProviderTest, WrongColvarsInput)
{
    PrepareInputForceProvider("4water");
    atoms_ = gmx_mtop_global_atoms(mtop);

    // atom 100 does not exist
    colvarsConfigString_ = R"(units gromacs
              colvar {
                  name d_atoms
                  distance {
                      group1 {
                          atomNumbers 100
                      }
                      group2 {
                          atomNumbers 2
                      }
                  }
              }
              harmonic {
                  colvars d_atoms
                  forceConstant 2000
                  centers 0.6
              })";

    EXPECT_ANY_THROW(ColvarsForceProvider forceProvider(colvarsConfigString_,
                                                        atoms_,
                                                        pbcType_,
                                                        &logger_,
                                                        KVTInputs,
                                                        temperature_,
                                                        seed_,
                                                        &atomSetManager_,
                                                        &cr_,
                                                        simulationTimeStep_,
                                                        atomCoords_,
                                                        prefixOutput_,
                                                        colvarsState_));
    done_atom(&atoms_);
}
#    endif
#endif

TEST_F(ColvarsForceProviderTest, CalculateForces4water)
{
    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());

    PrepareInputForceProvider("4water");
    ColvarsConfigStringFromFile("colvars_sample.dat");
    atoms_ = gmx_mtop_global_atoms(mtop);

    // Indexes taken from the Colvars Config file.
    atomCoords_ = { RVec(x_[0]), RVec(x_[4]) };

    // Prepare a ForceProviderInput
    ForceProviderInput forceProviderInput(x_, atoms_.nr, {}, {}, 0.0, 0, box_, cr_);

    // Prepare a ForceProviderOutput
    std::vector<RVec>   forces(atoms_.nr, RVec{ 0, 0, 0 });
    ForceWithVirial     forceWithVirial(forces, true);
    gmx_enerdata_t      enerdDummy(1, nullptr);
    ForceProviderOutput forceProviderOutput(&forceWithVirial, &enerdDummy);

    ColvarsForceProvider forceProvider(colvarsConfigString_,
                                       atoms_,
                                       pbcType_,
                                       &logger_,
                                       KVTInputs,
                                       temperature_,
                                       seed_,
                                       &atomSetManager_,
                                       &cr_,
                                       simulationTimeStep_,
                                       atomCoords_,
                                       prefixOutput_,
                                       colvarsState_);

    forceProvider.calculateForces(forceProviderInput, &forceProviderOutput);

    checker.setDefaultTolerance(gmx::test::relativeToleranceAsFloatingPoint(100.0, 5e-5));
    checker.checkReal(enerdDummy.term[F_COM_PULL], "Bias Energy");
    checker.checkSequence(forces.begin(), forces.end(), "Forces");

    done_atom(&atoms_);
}

TEST_F(ColvarsForceProviderTest, CalculateForcesAlanine)
{
    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());

    PrepareInputForceProvider("ala");
    ColvarsConfigStringFromFile("colvars_sample_alanine.dat");
    atoms_ = gmx_mtop_global_atoms(mtop);

    // Indexes taken from the Colvars Config file.
    atomCoords_ = { RVec(x_[0]), RVec(x_[4]), RVec(x_[10]), RVec(x_[11]) };

    // Prepare a ForceProviderInput
    ForceProviderInput forceProviderInput(x_, atoms_.nr, {}, {}, 0.0, 0, box_, cr_);

    // Prepare a ForceProviderOutput
    std::vector<RVec>   forces(atoms_.nr, RVec{ 0, 0, 0 });
    ForceWithVirial     forceWithVirial(forces, true);
    gmx_enerdata_t      enerdDummy(1, nullptr);
    ForceProviderOutput forceProviderOutput(&forceWithVirial, &enerdDummy);

    ColvarsForceProvider forceProvider(colvarsConfigString_,
                                       atoms_,
                                       pbcType_,
                                       &logger_,
                                       KVTInputs,
                                       temperature_,
                                       seed_,
                                       &atomSetManager_,
                                       &cr_,
                                       simulationTimeStep_,
                                       atomCoords_,
                                       prefixOutput_,
                                       colvarsState_);

    forceProvider.calculateForces(forceProviderInput, &forceProviderOutput);

    checker.setDefaultTolerance(gmx::test::relativeToleranceAsFloatingPoint(10.0, 5e-5));
    checker.checkReal(enerdDummy.term[F_COM_PULL], "Bias Energy");
    checker.checkSequence(forces.begin(), forces.end(), "Forces");

    done_atom(&atoms_);
}

} // namespace gmx
