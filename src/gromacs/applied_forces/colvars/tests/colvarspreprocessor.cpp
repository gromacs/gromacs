/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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
 * Tests for ColvarsPreProcessor class for Colvars MDModule
 *
 * \author Hubert Santuz <hubert.santuz@gmail.com>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "gromacs/applied_forces/colvars/colvarspreprocessor.h"

#include <cstddef>

#include <filesystem>
#include <iostream>
#include <list>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/fileio/confio.h"
#include "gromacs/gmxpreprocess/grompp.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

enum class PbcType : int;

namespace gmx
{

class ColvarsPreProcessorTest : public ::testing::Test
{
public:
    /*! \brief Generates tpr file from *.top and *.gro existing in the simulation database directory
     * and loads gmx_mtop_t from it
     */
    void makeMtopFromFile(const std::string& fileName, const std::string& mdpContent)
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
        bool       fullTopology;
        gmx_mtop_t mtop_;
        readConfAndTopology(tprName.c_str(), &fullTopology, &mtop_, &pbcType_, &coords, nullptr, box_);
        atoms_ = gmx_mtop_global_atoms(mtop_);
        x_     = gmx::constArrayRefFromArray(reinterpret_cast<gmx::RVec*>(coords), atoms_.nr);
    }

    ~ColvarsPreProcessorTest() override
    {
        done_atom(&atoms_);
        sfree(coords);
    }

protected:
    rvec*                      coords;
    gmx::test::TestFileManager fileManager_;
    t_atoms                    atoms_;
    PbcType                    pbcType_;
    matrix                     box_;
    ArrayRef<const RVec>       x_;
};

TEST_F(ColvarsPreProcessorTest, CanConstructColvarsPreProcess)
{
    // Reference input 4x SPCE waters from database 4water.top
    makeMtopFromFile("4water", "");

    EXPECT_NO_THROW(ColvarsPreProcessor colvarsPreProcess("", atoms_, pbcType_, nullptr, 0, -1, box_, x_));
}

TEST_F(ColvarsPreProcessorTest, CheckValuesFourWaters)
{
    // Reference input 4x SPCE waters from database 4water.top
    makeMtopFromFile("4water", "");

    // atomNumbers start at 1 in colvars
    std::string colvarsInput = R"(units gromacs
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

    ColvarsPreProcessor colvarsPreProcess(colvarsInput, atoms_, pbcType_, nullptr, 0, -1, box_, x_);

    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());

    // Check colvars & atoms values are correctly read
    checker.setDefaultTolerance(gmx::test::absoluteTolerance(0.001));
    checker.checkVector(colvarsPreProcess.getColvarsCoords()[1], "Coords Atom 4");

    const auto* const atomIds = colvarsPreProcess.get_atom_ids();
    checker.checkSequence(atomIds->begin(), atomIds->end(), "Index of colvars atoms");

    const auto* const atomMasses = colvarsPreProcess.get_atom_masses();
    checker.checkSequence(atomMasses->begin(), atomMasses->end(), "Masses of colvars atoms");

    const auto* const atomCharges = colvarsPreProcess.get_atom_charges();
    checker.checkSequence(atomCharges->begin(), atomCharges->end(), "Charges of colvars atoms");
}


TEST_F(ColvarsPreProcessorTest, CheckNestedInputFiles)
{
    // Reference input 4x SPCE waters from database 4water.top
    makeMtopFromFile("4water", "");

    // TODO: add ref xyz file
    std::string colvarsInput = R"(units gromacs
              indexFile <template_ndx>
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
                  centers 0.6
              })";

    // Replace index.ndx by its absolute path so Colvars can parse it.
    std::string pathIndex = gmx::test::TestFileManager::getInputFilePath("index.ndx").string();
    size_t      index     = colvarsInput.find("<template_ndx>");
    colvarsInput.replace(index, std::string("<template_ndx>").length(), pathIndex);

    ColvarsPreProcessor colvarsPreProcess(colvarsInput, atoms_, pbcType_, nullptr, 0, -1, box_, x_);

    // Make sure the index file inside colvarsInput was correctly read
    auto listInputStreams = colvarsPreProcess.list_input_stream_names();
    auto inputStreamIndex = listInputStreams.begin();
    EXPECT_EQ(pathIndex, *inputStreamIndex);
}

#if defined(__has_feature)
#    if !__has_feature(address_sanitizer)
TEST_F(ColvarsPreProcessorTest, WrongColvarsInput)
{
    // Reference input 4x SPCE waters from database 4water.top
    makeMtopFromFile("4water", "");

    // atom 100 does not exist
    std::string colvarsInput = R"(units gromacs
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
                  forceConstant 20000
                  centers 0.6
              })";

    EXPECT_ANY_THROW(ColvarsPreProcessor colvarsPreProcess(
            colvarsInput, atoms_, pbcType_, nullptr, 0, -1, box_, x_));
}
#    endif
#endif

} // namespace gmx
